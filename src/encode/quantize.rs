//! Quantisation of one granule of one channel against an allowed noise per
//! band, and the side-information fields that go with it.
//!
//! The requantiser of ISO/IEC 11172-3 2.4.3.4.7.1 gives every line of a
//! long band the step 2^(s/4) with s = global_gain - 210 -
//! 2 (1 + scalefac_scale)(scalefac + preflag pretab), and of a short band's
//! window s = global_gain - 210 - 8 subblock_gain - 2 (1 + scalefac_scale)
//! scalefac. The quantiser is its inverse,
//! ix = nint((|xr| / 2^(s/4))^(3/4) - 0.0946) (Annex C's rounding).
//!
//! For each band the coarsest step whose noise stays within the band's
//! allowance is found; global_gain is set by the coarsest band, and the
//! other bands are brought down to their steps with scalefactors (and, in
//! short blocks, subblock gains). A band whose scalefactor would exceed its
//! field keeps the largest one (more noise than allowed there).

use std::sync::OnceLock;

use super::huffman::{self, Coding, MAX_VALUE};
use super::psy::Mask;
use crate::bits::BitWriter;
use crate::tables::layer3::{PRETAB, SLEN};

/// |ix|^(4/3) for every codable magnitude.
fn pow43() -> &'static [f64] {
    static T: OnceLock<Vec<f64>> = OnceLock::new();
    T.get_or_init(|| (0..=MAX_VALUE as usize + 1).map(|i| (i as f64).powf(4.0 / 3.0)).collect())
}

/// One granule-channel's spectrum ready to quantise.
pub(crate) struct Spectrum<'a> {
    /// Lines in bitstream order (short blocks: band, window, line).
    pub(crate) xr: &'a [f32; 576],
    /// 0 normal, 1 start, 2 short, 3 stop.
    pub(crate) block_type: u8,
    pub(crate) mask: &'a Mask,
    pub(crate) long_edges: &'a [u16; 23],
    pub(crate) short_edges: &'a [u16; 14],
    pub(crate) lsf: bool,
}

/// The result: values and the side information that describes them.
#[derive(Clone, Debug)]
pub(crate) struct Quantised {
    pub(crate) ix: [i32; 576],
    pub(crate) global_gain: u8,
    pub(crate) scalefac_scale: bool,
    /// MPEG-1 long blocks: pretab adds to the high bands' scalefactors.
    pub(crate) preflag: bool,
    pub(crate) subblock_gain: [u8; 3],
    /// Long scalefactors (bands 0..21) or short ([band][window]).
    pub(crate) sf_l: [u8; 22],
    pub(crate) sf_s: [[u8; 3]; 13],
    pub(crate) scalefac_compress: u16,
    /// Scalefactor lengths for the writer (MPEG-1: slen1, slen2; LSF: four).
    pub(crate) slen: [u8; 4],
    pub(crate) coding: Coding,
    /// Bits of the scalefactors.
    pub(crate) part2_bits: u32,
}

impl Quantised {
    pub(crate) fn bits(&self) -> u32 {
        self.part2_bits + self.coding.bits
    }

    /// An all-zero granule.
    pub(crate) fn silent() -> Quantised {
        Quantised {
            ix: [0; 576],
            global_gain: 0,
            scalefac_scale: false,
            preflag: false,
            subblock_gain: [0; 3],
            sf_l: [0; 22],
            sf_s: [[0; 3]; 13],
            scalefac_compress: 0,
            slen: [0; 4],
            coding: Coding::default(),
            part2_bits: 0,
        }
    }
}

/// Bands' line ranges in bitstream order, with their allowed noise.
struct BandPlan {
    /// (start, end, allowed noise, window) per band; window 3 for long.
    bands: Vec<(usize, usize, f64, usize, usize)>,
}

fn plan(sp: &Spectrum, scale: f64) -> BandPlan {
    let mut bands = Vec::new();
    match sp.mask {
        Mask::Long { thr, .. } => {
            for b in 0..22 {
                let lo = usize::from(sp.long_edges[b]);
                let hi = usize::from(sp.long_edges[b + 1]);
                bands.push((lo, hi, thr[b] * scale, 3, b));
            }
        }
        Mask::Short { thr, .. } => {
            let mut i = 0;
            for b in 0..13 {
                let w = usize::from(sp.short_edges[b + 1] - sp.short_edges[b]);
                for win in 0..3 {
                    bands.push((i, i + w, thr[b][win] * scale, win, b));
                    i += w;
                }
            }
        }
    }
    BandPlan { bands }
}

/// Noise of quantising `xa` (magnitudes, with `x34` = xa^(3/4)) with step
/// 2^(s/4).
fn noise(xa: &[f64], x34: &[f64], s: i32) -> f64 {
    let p = pow43();
    let q = (-3.0 * f64::from(s) / 16.0).exp2();
    let g = (f64::from(s) / 4.0).exp2();
    let mut n = 0.0;
    for (&a, &c) in xa.iter().zip(x34) {
        let ix = ((c * q + 0.4054) as usize).min(MAX_VALUE as usize);
        let e = a - p[ix] * g;
        n += e * e;
    }
    n
}

/// The smallest s at which every value of the band is codable.
fn s_floor(xmax: f64) -> i32 {
    if xmax <= 0.0 {
        return i32::MIN / 4;
    }
    // (xmax / 2^(s/4))^(3/4) + 0.4054 < MAX_VALUE + 1
    (4.0 * xmax.log2() - (16.0 / 3.0) * (f64::from(MAX_VALUE) + 0.5).log2()).ceil() as i32
}

/// Quantise with every band's allowed noise multiplied by `scale`.
pub(crate) fn quantise(sp: &Spectrum, scale: f64) -> Quantised {
    let xa: Vec<f64> = sp.xr.iter().map(|&v| f64::from(v).abs()).collect();
    let x34: Vec<f64> = xa.iter().map(|&v| v.powf(0.75)).collect();
    let plan = plan(sp, scale);
    let short = sp.block_type == 2;
    // Required step per band: None for a band that may be quantised to
    // nothing (its energy is within its allowance) and for the highest band,
    // which has no scalefactor and takes the global step.
    let last_band = if short { 12 } else { 21 };
    let mut req: Vec<Option<i32>> = Vec::with_capacity(plan.bands.len());
    // The global step must keep every value of every band codable.
    let mut floor = i32::MIN / 4;
    for &(lo, hi, allowed, _, b) in &plan.bands {
        let xmax = xa[lo..hi].iter().copied().fold(0.0, f64::max);
        let energy: f64 = xa[lo..hi].iter().map(|v| v * v).sum();
        if xmax == 0.0 {
            req.push(None);
            continue;
        }
        let lo_s = s_floor(xmax);
        if b == last_band || energy <= allowed {
            floor = floor.max(lo_s);
            req.push(None);
            continue;
        }
        // Largest s in [lo_s, 60] with noise <= allowed (noise rises with s).
        let (mut a, mut z) = (lo_s, 60);
        if noise(&xa[lo..hi], &x34[lo..hi], a) > allowed {
            req.push(Some(a));
            continue;
        }
        while a < z {
            let m = (a + z + 1).div_euclid(2);
            if noise(&xa[lo..hi], &x34[lo..hi], m) <= allowed {
                a = m;
            } else {
                z = m - 1;
            }
        }
        req.push(Some(a));
    }
    if req.iter().all(Option::is_none) && floor == i32::MIN / 4 {
        return Quantised::silent();
    }
    // Band maxima of the scalefactor fields: 15 for the low bands (4 bits),
    // 7 for the high (3 bits) — MPEG-1 slen1 / slen2 and LSF slen0..3 alike.
    let sf_max = |b: usize| if (short && b < 6) || (!short && b < 11) { 15i32 } else { 7 };
    let mut best: Option<(bool, Quantised)> = None;
    // Every combination of scalefac_scale and (MPEG-1 long blocks) preflag:
    // each meets the allowances as far as its scalefactor range reaches; the
    // one that clips least, then costs least, wins.
    let variants: &[(bool, bool)] = if short || sp.lsf {
        &[(false, false), (true, false)]
    } else {
        &[(false, false), (true, false), (false, true), (true, true)]
    };
    for &(scale_bit, pre) in variants {
        let m = if scale_bit { 4 } else { 2 };
        let pretab = |b: usize| if pre && !short && b < 22 { i32::from(PRETAB[b]) } else { 0 };
        let mut q = Quantised::silent();
        q.scalefac_scale = scale_bit;
        q.preflag = pre;
        // Per window (one "window" for long blocks): the coarsest required
        // step, capped so that the finest band can still reach its step with
        // the largest scalefactor.
        let windows = if short { 3 } else { 1 };
        let mut base = [i32::MIN / 4; 3];
        for (w, bw) in base.iter_mut().enumerate().take(windows) {
            let mut hi = i32::MIN / 4;
            let mut cap = i32::MAX / 4;
            for (&(_, _, _, win, b), r) in plan.bands.iter().zip(&req) {
                if short && win != w {
                    continue;
                }
                if let Some(s) = *r {
                    hi = hi.max(s);
                    cap = cap.min(s + m * (sf_max(b) + pretab(b)));
                }
            }
            *bw = hi.min(cap);
        }
        // No band needs coding: as coarse as the values allow (bands that
        // need nothing must not pull the step down).
        let required = base[..windows].iter().copied().max().filter(|&b| b > i32::MIN / 4).unwrap_or(45);
        let top = required.max(floor);
        let mut gg = (top + 210).clamp(0, 255);
        let mut clipped = false;
        if short {
            for w in 0..3 {
                if base[w] > i32::MIN / 4 {
                    q.subblock_gain[w] = ((gg - 210 - base[w]) / 8).clamp(0, 7) as u8;
                }
            }
        }
        for (&(_, _, _, w, b), r) in plan.bands.iter().zip(&req) {
            let Some(s) = *r else { continue };
            let base = gg - 210 - if short { 8 * i32::from(q.subblock_gain[w]) } else { 0 };
            let need = ((base - s + m - 1).div_euclid(m) - pretab(b)).max(0);
            let sf = need.min(sf_max(b));
            clipped |= need > sf;
            if short {
                q.sf_s[b][w] = sf as u8;
            } else {
                q.sf_l[b] = sf as u8;
            }
        }
        // Values beyond the tables (a band left coarser or finer than its
        // floor allows): raise the global step until all are codable.
        loop {
            q.global_gain = gg as u8;
            if fill_values(&mut q, &plan, &x34, short) || gg >= 255 {
                break;
            }
            gg += 1;
        }
        choose_scalefac_compress(&mut q, short, sp.lsf);
        let region1_ws = if short { 3 * usize::from(sp.short_edges[3]) } else { usize::from(sp.long_edges[8]) };
        q.coding = huffman::choose(&q.ix, sp.long_edges, sp.block_type != 0, region1_ws);
        let better = match &best {
            None => true,
            Some((bc, bq)) => (clipped, q.bits()) < (*bc, bq.bits()),
        };
        if better {
            best = Some((clipped, q));
        }
    }
    best.expect("one variant tried").1

}

/// Quantise every line with the chosen parameters; false if a value
/// exceeds what the tables can code.
fn fill_values(q: &mut Quantised, plan: &BandPlan, x34: &[f64], short: bool) -> bool {
    let m = if q.scalefac_scale { 4 } else { 2 };
    let gg = i32::from(q.global_gain);
    for &(lo, hi, _, w, b) in &plan.bands {
        let s = if short {
            gg - 210 - 8 * i32::from(q.subblock_gain[w]) - m * i32::from(q.sf_s[b][w])
        } else {
            gg - 210 - m * (i32::from(q.sf_l[b]) + pre(q, b))
        };
        let qf = (-3.0 * f64::from(s) / 16.0).exp2();
        for i in lo..hi {
            let v = x34[i] * qf + 0.4054;
            if v >= f64::from(MAX_VALUE) + 1.0 {
                return false;
            }
            q.ix[i] = v as i32;
        }
    }
    true
}

/// The pretab amplification of long band `b` when preflag is set.
fn pre(q: &Quantised, b: usize) -> i32 {
    if q.preflag && b < 22 { i32::from(PRETAB[b]) } else { 0 }
}

/// Apply signs from the spectrum.
pub(crate) fn apply_signs(q: &mut Quantised, xr: &[f32; 576]) {
    for (v, &x) in q.ix.iter_mut().zip(xr) {
        if x < 0.0 {
            *v = -*v;
        }
    }
}

fn bits_for(max: u8) -> u8 {
    (8 - max.leading_zeros()) as u8
}

/// Pick the cheapest scalefac_compress that holds the scalefactors.
fn choose_scalefac_compress(q: &mut Quantised, short: bool, lsf: bool) {
    if lsf {
        // Table row 0 of 13818-3: four groups (6, 5, 5, 5 long bands; 3
        // short bands, i.e. 9 values, each); slen0, slen1 up to 4 bits,
        // slen2, slen3 up to 3.
        let groups: [(usize, usize); 4] = if short { [(0, 3), (3, 6), (6, 9), (9, 12)] } else { [(0, 6), (6, 11), (11, 16), (16, 21)] };
        let mut slen = [0u8; 4];
        for (g, &(a, b)) in groups.iter().enumerate() {
            let max = if short {
                (a..b).flat_map(|band| q.sf_s[band]).max().unwrap_or(0)
            } else {
                q.sf_l[a..b].iter().copied().max().unwrap_or(0)
            };
            slen[g] = bits_for(max);
        }
        let counts: [u32; 4] = if short { [9, 9, 9, 9] } else { [6, 5, 5, 5] };
        q.part2_bits = slen.iter().zip(counts).map(|(&s, c)| u32::from(s) * c).sum();
        q.scalefac_compress =
            ((u16::from(slen[0]) * 5 + u16::from(slen[1])) << 4) + (u16::from(slen[2]) << 2) + u16::from(slen[3]);
        q.slen = slen;
        return;
    }
    let (lo_max, hi_max) = if short {
        (
            (0..6).flat_map(|b| q.sf_s[b]).max().unwrap_or(0),
            (6..12).flat_map(|b| q.sf_s[b]).max().unwrap_or(0),
        )
    } else {
        (q.sf_l[..11].iter().copied().max().unwrap_or(0), q.sf_l[11..21].iter().copied().max().unwrap_or(0))
    };
    let (n1, n2) = if short { (18u32, 18u32) } else { (11, 10) };
    let (need1, need2) = (bits_for(lo_max), bits_for(hi_max));
    let mut best = (u32::MAX, 0usize);
    for (i, &(s1, s2)) in SLEN.iter().enumerate() {
        if s1 >= need1 && s2 >= need2 {
            let bits = n1 * u32::from(s1) + n2 * u32::from(s2);
            if bits < best.0 {
                best = (bits, i);
            }
        }
    }
    q.scalefac_compress = best.1 as u16;
    q.part2_bits = best.0;
    q.slen = [SLEN[best.1].0, SLEN[best.1].1, 0, 0];
}

/// Write the scalefactors (part 2) in the decoder's reading order.
pub(crate) fn write_scalefactors(w: &mut BitWriter, q: &Quantised, short: bool, lsf: bool) {
    if lsf {
        let counts: [usize; 4] = if short { [9, 9, 9, 9] } else { [6, 5, 5, 5] };
        let mut vals: Vec<u8> = Vec::with_capacity(36);
        if short {
            for b in 0..12 {
                vals.extend_from_slice(&q.sf_s[b]);
            }
        } else {
            vals.extend_from_slice(&q.sf_l[..21]);
        }
        let mut it = vals.into_iter();
        for (g, &n) in counts.iter().enumerate() {
            for _ in 0..n {
                let v = it.next().unwrap_or(0);
                w.put(u32::from(v), u32::from(q.slen[g]));
            }
        }
        return;
    }
    let (s1, s2) = (u32::from(q.slen[0]), u32::from(q.slen[1]));
    if short {
        for b in 0..12 {
            for win in 0..3 {
                w.put(u32::from(q.sf_s[b][win]), if b < 6 { s1 } else { s2 });
            }
        }
    } else {
        for b in 0..21 {
            w.put(u32::from(q.sf_l[b]), if b < 11 { s1 } else { s2 });
        }
    }
}

/// The decoder's reconstruction of `q` (for measuring noise).
pub(crate) fn dequantise(q: &Quantised, sp: &Spectrum) -> [f32; 576] {
    let plan = plan(sp, 1.0);
    let short = sp.block_type == 2;
    let m = if q.scalefac_scale { 4 } else { 2 };
    let gg = i32::from(q.global_gain);
    let p = pow43();
    let mut out = [0.0f32; 576];
    for &(lo, hi, _, w, b) in &plan.bands {
        let s = if short {
            gg - 210 - 8 * i32::from(q.subblock_gain[w]) - m * i32::from(q.sf_s[b][w])
        } else {
            gg - 210 - m * (i32::from(q.sf_l[b]) + pre(q, b))
        };
        let g = (f64::from(s) / 4.0).exp2();
        for i in lo..hi {
            let v = p[q.ix[i].unsigned_abs() as usize] * g;
            out[i] = if q.ix[i] < 0 { -v as f32 } else { v as f32 };
        }
    }
    out
}
