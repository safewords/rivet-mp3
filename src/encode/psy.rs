//! The psychoacoustic model: for each granule, the noise each scalefactor
//! band can carry without being heard (the masking threshold), computed on
//! the MDCT spectrum that is quantised, and the transient detector that
//! switches to short blocks.
//!
//! The model follows the outline of ISO/IEC 11172-3 Annex D (psychoacoustic
//! model 2) in a simplified form, partitioned by scalefactor band:
//!
//! - band energies and a tonality estimate from the spectral flatness of the
//!   band's lines (a flat band is noise-like, a peaked one tone-like);
//! - spreading across bands with the spreading function of Annex D
//!   (Schroeder's: 15.81 + 7.5 (dz + 0.474) - 17.5 sqrt(1 + (dz + 0.474)^2)
//!   dB, dz in Bark);
//! - a masking offset of 14.5 + z dB below a tonal masker and 5.5 dB below
//!   a noise masker, interpolated by tonality;
//! - the threshold in quiet (Terhardt's approximation of the absolute
//!   threshold of hearing), with a full-scale sine taken as 96 dB SPL;
//! - pre-echo control for long blocks: a band's threshold may not exceed
//!   twice its value in the previous granule.

use std::f64::consts::PI;

use super::analysis::Analysis;
use crate::layer3::imdct::mdct;

/// Bark scale (Zwicker and Terhardt).
fn bark(f: f64) -> f64 {
    13.0 * (0.00076 * f).atan() + 3.5 * (f / 7500.0).powi(2).atan()
}

/// The spreading function of 11172-3 Annex D, linear power, for a masker
/// `dz` Bark below the maskee.
fn spreading(dz: f64) -> f64 {
    let t = dz + 0.474;
    let db = 15.81 + 7.5 * t - 17.5 * (1.0 + t * t).sqrt();
    if db < -60.0 { 0.0 } else { 10f64.powf(db / 10.0) }
}

/// Absolute threshold of hearing, dB SPL (Terhardt).
fn ath_db(f: f64) -> f64 {
    let k = (f.max(20.0)) / 1000.0;
    (3.64 * k.powf(-0.8) - 6.5 * (-0.6 * (k - 3.3).powi(2)).exp() + 1e-3 * k.powi(4)).min(110.0)
}

/// Band layout with precomputed spreading and quiet thresholds.
#[derive(Clone, Debug)]
struct Bands {
    /// Line ranges.
    lo: Vec<usize>,
    hi: Vec<usize>,
    z: Vec<f64>,
    /// spread[b][j]: weight of masker band j on band b.
    spread: Vec<Vec<f64>>,
    /// Normalisation of the spread energy (see [`Bands::thresholds`]).
    norm: Vec<f64>,
    /// Threshold in quiet as band noise energy, in the spectrum's units.
    ath: Vec<f64>,
}

impl Bands {
    fn new(edges: &[u16], line_hz: f64, full_scale_tone: f64) -> Bands {
        let n = edges.len() - 1;
        let lo: Vec<usize> = edges[..n].iter().map(|&e| usize::from(e)).collect();
        let hi: Vec<usize> = edges[1..].iter().map(|&e| usize::from(e)).collect();
        let z: Vec<f64> = (0..n).map(|b| bark((lo[b] + hi[b]) as f64 * 0.5 * line_hz)).collect();
        let spread: Vec<Vec<f64>> =
            (0..n).map(|b| (0..n).map(|j| spreading(z[b] - z[j])).collect()).collect();
        let norm: Vec<f64> = (0..n)
            .map(|b| {
                let w = (hi[b] - lo[b]) as f64;
                (0..n).map(|j| spread[b][j] * (hi[j] - lo[j]) as f64).sum::<f64>() / w
            })
            .collect();
        let ath = (0..n)
            .map(|b| {
                let min_db = (lo[b]..hi[b]).map(|i| ath_db((i as f64 + 0.5) * line_hz)).fold(f64::MAX, f64::min);
                full_scale_tone * 10f64.powf((min_db - 96.0) / 10.0)
            })
            .collect();
        Bands { lo, hi, z, spread, norm, ath }
    }

    /// Thresholds for one spectrum (`x`, lines indexed by `lo..hi`),
    /// and the band energies.
    fn thresholds(&self, x: &[f32], out_thr: &mut [f64], out_energy: &mut [f64]) {
        let n = self.lo.len();
        let mut tonality = vec![0.0f64; n];
        for b in 0..n {
            let lines = &x[self.lo[b]..self.hi[b]];
            let e: f64 = lines.iter().map(|&v| f64::from(v) * f64::from(v)).sum();
            out_energy[b] = e;
            if e <= 1e-30 {
                continue;
            }
            // Spectral flatness: geometric over arithmetic mean of line
            // energies.
            let m = lines.len() as f64;
            let log_sum: f64 = lines.iter().map(|&v| (f64::from(v) * f64::from(v) + 1e-30).ln()).sum();
            let geo = (log_sum / m).exp();
            let sfm_db = 10.0 * (geo / (e / m)).log10();
            tonality[b] = (sfm_db / -25.0).clamp(0.0, 1.0);
        }
        for b in 0..n {
            let spread: f64 = (0..n).map(|j| out_energy[j] * self.spread[b][j]).sum();
            let offset = tonality[b] * (14.5 + self.z[b]) + (1.0 - tonality[b]) * 5.5;
            let thr = spread / self.norm[b] * 10f64.powf(-offset / 10.0);
            out_thr[b] = thr.max(self.ath[b]);
        }
    }
}

/// The allowed noise of one granule of one channel.
// The short variant is the larger (2 x 39 values against 2 x 22); masks
// live a frame at a time, so boxing it would buy nothing.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug)]
pub(crate) enum Mask {
    /// Long, start and stop blocks: 22 long bands.
    Long { thr: [f64; 22], energy: [f64; 22] },
    /// Short blocks: 13 short bands by 3 windows.
    Short { thr: [[f64; 3]; 13], energy: [[f64; 3]; 13] },
}

impl Mask {
    /// Perceptual entropy estimate (bits): sum over bands of
    /// width * log2(1 + sqrt(energy / threshold)).
    pub(crate) fn pe(&self, long_edges: &[u16; 23], short_edges: &[u16; 14]) -> f64 {
        let term = |w: f64, e: f64, t: f64| if e > t { w * (1.0 + (e / t).sqrt()).log2() } else { 0.0 };
        match self {
            Mask::Long { thr, energy } => (0..22)
                .map(|b| term(f64::from(long_edges[b + 1] - long_edges[b]), energy[b], thr[b]))
                .sum(),
            Mask::Short { thr, energy } => (0..13)
                .flat_map(|b| (0..3).map(move |w| (b, w)))
                .map(|(b, w)| term(f64::from(short_edges[b + 1] - short_edges[b]), energy[b][w], thr[b][w]))
                .sum(),
        }
    }

    /// The element-wise minimum of two masks of the same kind (the masks
    /// of M and S in mid/side coding are the smaller of L's and R's, so
    /// that noise kept below them in M and S stays below them in L and R).
    pub(crate) fn min(&self, other: &Mask) -> Mask {
        match (self, other) {
            (Mask::Long { thr: a, energy }, Mask::Long { thr: b, .. }) => {
                Mask::Long { thr: std::array::from_fn(|i| a[i].min(b[i])), energy: *energy }
            }
            (Mask::Short { thr: a, energy }, Mask::Short { thr: b, .. }) => Mask::Short {
                thr: std::array::from_fn(|i| std::array::from_fn(|w| a[i][w].min(b[i][w]))),
                energy: *energy,
            },
            _ => self.clone(),
        }
    }

    /// Replace the energies (keeping the thresholds) with those of another
    /// spectrum: used when the channel coded is M or S.
    pub(crate) fn with_energy_of(&self, x: &[f32; 576], long_edges: &[u16; 23], short_edges: &[u16; 14]) -> Mask {
        match self {
            Mask::Long { thr, .. } => {
                let energy = std::array::from_fn(|b| {
                    x[usize::from(long_edges[b])..usize::from(long_edges[b + 1])]
                        .iter()
                        .map(|&v| f64::from(v).powi(2))
                        .sum()
                });
                Mask::Long { thr: *thr, energy }
            }
            Mask::Short { thr, .. } => {
                let mut energy = [[0.0; 3]; 13];
                let mut i = 0;
                for (b, e) in energy.iter_mut().enumerate() {
                    let w = usize::from(short_edges[b + 1] - short_edges[b]);
                    for ew in e.iter_mut() {
                        *ew = x[i..i + w].iter().map(|&v| f64::from(v).powi(2)).sum();
                        i += w;
                    }
                }
                Mask::Short { thr: *thr, energy }
            }
        }
    }
}

/// The model for one sampling frequency.
#[derive(Clone, Debug)]
pub(crate) struct Psy {
    long: Bands,
    short: Bands,
    short_edges: [u16; 14],
}

impl Psy {
    pub(crate) fn new(sample_rate: u32, long_edges: &[u16; 23], short_edges: &[u16; 14]) -> Psy {
        let fs = f64::from(sample_rate);
        let (k_long, k_short) = full_scale_energy(sample_rate);
        Psy {
            long: Bands::new(long_edges, fs / 2.0 / 576.0, k_long),
            short: Bands::new(short_edges, fs / 2.0 / 192.0, k_short),
            short_edges: *short_edges,
        }
    }

    /// The mask of a long (or start / stop) block spectrum. `prev` is the
    /// previous granule's long thresholds, for pre-echo control.
    pub(crate) fn long(&self, x: &[f32; 576], prev: Option<&[f64; 22]>) -> Mask {
        let mut thr = [0.0; 22];
        let mut energy = [0.0; 22];
        self.long.thresholds(x, &mut thr, &mut energy);
        if let Some(p) = prev {
            for b in 0..22 {
                if p[b] > 0.0 {
                    thr[b] = thr[b].min(2.0 * p[b]).max(self.long.ath[b]);
                }
            }
        }
        Mask::Long { thr, energy }
    }

    /// The mask of a short block spectrum in bitstream order (band,
    /// window, line).
    pub(crate) fn short(&self, x: &[f32; 576]) -> Mask {
        // De-interleave into three window spectra of 192 lines.
        let mut win = [[0.0f32; 192]; 3];
        let mut i = 0;
        for b in 0..13 {
            let lo = usize::from(self.short_edges[b]);
            let w = usize::from(self.short_edges[b + 1]) - lo;
            for spec in win.iter_mut() {
                spec[lo..lo + w].copy_from_slice(&x[i..i + w]);
                i += w;
            }
        }
        let mut thr = [[0.0; 3]; 13];
        let mut energy = [[0.0; 3]; 13];
        for w in 0..3 {
            let mut t = [0.0; 13];
            let mut e = [0.0; 13];
            self.short.thresholds(&win[w], &mut t, &mut e);
            for b in 0..13 {
                thr[b][w] = t[b];
                energy[b][w] = e[b];
            }
        }
        Mask::Short { thr, energy }
    }
}

/// Energy, in the encoder's MDCT spectrum, of a full-scale sine (amplitude
/// 1.0) as seen by a long block and by one short window: the reference that
/// places the threshold in quiet (96 dB SPL). Measured by running a sine
/// through the analysis filterbank and the MDCT.
fn full_scale_energy(sample_rate: u32) -> (f64, f64) {
    let f = f64::from(sample_rate) / 64.0 * 2.5; // middle of subband 2
    let n = 576 * 8;
    let mut a = Analysis::default();
    let mut slots = Vec::new();
    for blk in 0..n / 32 {
        let input: Vec<f32> =
            (0..32).map(|i| ((2.0 * PI * f * (blk * 32 + i) as f64) / f64::from(sample_rate)).sin() as f32).collect();
        let mut s = [0.0; 32];
        a.run(&input, &mut s);
        slots.push(s);
    }
    let energy = |bt: u8| {
        let g = 5; // a granule well into the steady state
        let mut e = 0.0;
        for sb in 0..32 {
            let mut x = [0.0f64; 36];
            for (i, v) in x.iter_mut().enumerate() {
                let t = (g - 1) * 18 + i;
                *v = if sb % 2 == 1 && i % 2 == 1 { -slots[t][sb] } else { slots[t][sb] };
            }
            let mut c = [0.0f32; 18];
            mdct(&x, bt, &mut c);
            e += c.iter().map(|&v| f64::from(v).powi(2)).sum::<f64>();
        }
        e
    };
    let long = energy(0);
    // One short window holds a third of the long block's energy.
    let short = energy(2) / 3.0;
    (long, short)
}

/// Attack detection on the input: high-passed energy in short segments,
/// compared with the recent past.
#[derive(Clone, Debug, Default)]
pub(crate) struct Transient {
    /// Energy of recent segments.
    history: Vec<f64>,
    prev_in: f32,
    hp_state: f64,
}

/// Segment length for the attack detector (a third of a short block's
/// stride, so an attack is placed within a short window).
const SEG: usize = 64;

impl Transient {
    /// Feed the input samples of one granule's region; returns whether it
    /// holds an attack.
    pub(crate) fn granule(&mut self, x: &[f32]) -> bool {
        let mut attack = false;
        for seg in x.chunks(SEG) {
            let mut e = 0.0f64;
            for &v in seg {
                // First-order high-pass (differences, leaky), emphasising
                // the frequencies pre-echo is heard in.
                let d = f64::from(v - self.prev_in);
                self.prev_in = v;
                self.hp_state = 0.6 * self.hp_state + d;
                e += self.hp_state * self.hp_state;
            }
            let past = if self.history.is_empty() {
                e
            } else {
                self.history.iter().sum::<f64>() / self.history.len() as f64
            };
            let floor = 1e-6 * SEG as f64;
            if e > floor && e > 10.0 * past.max(floor / 10.0) {
                attack = true;
            }
            self.history.push(e);
            if self.history.len() > 9 {
                self.history.remove(0);
            }
        }
        attack
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spreading_peaks_at_zero_and_falls_off() {
        let at0 = spreading(0.0);
        assert!((10.0 * at0.log10()).abs() < 0.5, "{at0}");
        assert!(spreading(-3.0) < 0.01 * at0); // masker 3 Bark above
        assert!(spreading(3.0) < 0.01 * at0); // masker 3 Bark below: slower fall
        assert!(spreading(3.0) > spreading(-3.0));
    }

    #[test]
    fn quiet_threshold_shape() {
        assert!(ath_db(3300.0) < ath_db(1000.0));
        assert!(ath_db(100.0) > 20.0);
        assert!(ath_db(16000.0) > 40.0);
    }

    #[test]
    fn transient_detector() {
        let mut t = Transient::default();
        let quiet: Vec<f32> = (0..576).map(|i| 0.01 * (i as f32 * 0.05).sin()).collect();
        assert!(!t.granule(&quiet));
        assert!(!t.granule(&quiet));
        let mut hit = quiet.clone();
        for v in hit[300..].iter_mut() {
            *v += 0.5 * if (v.to_bits() & 1) == 0 { 1.0 } else { -1.0 };
        }
        assert!(t.granule(&hit));
    }
}
