//! The polyphase synthesis subband filterbank of ISO/IEC 11172-3 2.4.3.2.2
//! / Annex A, Figure 3-A.2 ("synthesis subband filter flow chart"), shared
//! by all three layers.
//!
//! For each set of 32 subband samples S[k]:
//!
//! 1. shift the 1024-entry vector V by 64;
//! 2. V[i] = sum_k N[i][k] S[k], N[i][k] = cos((16 + i)(2k + 1) pi / 64),
//!    i = 0..63 (matrixing);
//! 3. build U from V: U[64i + j] = V[128i + j], U[64i + 32 + j] =
//!    V[128i + 96 + j];
//! 4. window: W[i] = U[i] D[i];
//! 5. output sample j = sum_{i=0}^{15} W[j + 32i].
//!
//! The matrixing uses the cosine matrix's symmetries: rows 0..16 and
//! 49..64 of N repeat rows 16..=48 (up to sign), so only 33 distinct dot
//! products are computed, directly, in double precision.
//!
//! Both sums run side by side over their outputs (the 33 rows; the 32
//! samples), each output's own terms added in the standard's order with
//! separate multiplies and adds, so the loops vectorise and every output
//! is the same to the bit as summing it alone, on any CPU.

use std::sync::OnceLock;

use crate::tables::window::synthesis_window;

struct Tables {
    /// cos((16 + i)(2k + 1) pi / 64) for i = 16..=48 (the 33 distinct
    /// rows), row-major [i - 16][k].
    #[cfg(test)]
    n: Vec<f64>,
    /// The same, column-major [k][i - 16], rows padded to 36 with zeros.
    n_t: [[f64; 36]; 32],
    /// D[0..512].
    d: [f64; 512],
}

fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| {
        let mut n = vec![0.0; 33 * 32];
        for i in 16..=48 {
            for k in 0..32 {
                n[(i - 16) * 32 + k] = (((16 + i) * (2 * k + 1)) as f64 * std::f64::consts::PI / 64.0).cos();
            }
        }
        let n_t = std::array::from_fn(|k| std::array::from_fn(|r| if r < 33 { n[r * 32 + k] } else { 0.0 }));
        Tables {
            #[cfg(test)]
            n,
            n_t,
            d: synthesis_window(),
        }
    })
}

/// One channel's synthesis filterbank state.
#[derive(Clone)]
pub(crate) struct Synth {
    /// V, 1024 entries, as a ring: entry `j` of the standard's V lives at
    /// `v[(off + j) % 1024]`.
    v: Box<[f64; 1024]>,
    off: usize,
}

impl Default for Synth {
    fn default() -> Self {
        Self { v: Box::new([0.0; 1024]), off: 0 }
    }
}

impl Synth {
    pub(crate) fn reset(&mut self) {
        self.v.fill(0.0);
        self.off = 0;
    }

    /// The 64 matrixed values for `s`, from the 33 distinct rows of N (see
    /// [`row_value`]; the unit test checks every entry of the reconstructed
    /// 64 x 32 matrix against its definition).
    #[cfg(test)]
    fn matrix(s: &[f32; 32], out: &mut [f64; 64]) {
        matrix_with(tables(), s, out);
    }

    /// Filter one set of 32 subband samples into 32 PCM samples.
    pub(crate) fn run(&mut self, s: &[f32; 32], pcm: &mut [f32]) {
        self.off = (self.off + 1024 - 64) % 1024;
        synthesise(&mut self.v, self.off, s, pcm);
    }
}

#[inline(always)]
fn matrix_with(t: &Tables, s: &[f32; 32], out: &mut [f64; 64]) {
    {
        // Each row's sum over k from -0.0, as `Iterator::sum` does.
        let mut acc = [-0.0f64; 36];
        for (col, &x) in t.n_t.iter().zip(s) {
            let x = f64::from(x);
            for (a, &c) in acc.iter_mut().zip(col) {
                *a += c * x;
            }
        }
        let mid: &[f64; 33] = acc[..33].try_into().expect("33 rows");
        for i in 0..64 {
            out[i] = row_value(i, mid);
        }
    }
}

crate::simd::multiversion! {
/// One step of the filterbank on the ring `v` whose start has just moved
/// to `off`: matrix `s` into it, window and sum into `pcm`.
fn synthesise(v: &mut [f64; 1024], off: usize, s: &[f32; 32], pcm: &mut [f32]) {
    let t = tables();
    let mut m = [0.0f64; 64];
    matrix_with(t, s, &mut m);
    for (i, &x) in m.iter().enumerate() {
        v[(off + i) % 1024] = x;
    }
        // `off` and the offsets are multiples of 32, so each run of 32
        // entries is contiguous in the ring.
        let mut sum = [0.0f64; 32];
        for i in 0..8 {
            let ra = (off + 128 * i) % 1024;
            let rb = (off + 128 * i + 96) % 1024;
            let a: &[f64; 32] = v[ra..ra + 32].try_into().expect("32 entries");
            let b: &[f64; 32] = v[rb..rb + 32].try_into().expect("32 entries");
            let da: &[f64; 32] = t.d[64 * i..64 * i + 32].try_into().expect("32 entries");
            let db: &[f64; 32] = t.d[64 * i + 32..64 * i + 64].try_into().expect("32 entries");
            for j in 0..32 {
                sum[j] += a[j] * da[j] + b[j] * db[j];
            }
        }
        for (out, &s) in pcm.iter_mut().zip(&sum) {
            *out = s as f32;
        }
}
}

/// N[i] · S for any i in 0..64 from the 33 computed rows (i = 16..=48).
/// With a = (2k+1)π/64, 32a is an odd multiple of π/2 and 64a an odd
/// multiple of π, so cos((32 - m)a) = -cos((32 + m)a) and
/// cos((64 - m)a) = cos((64 + m)a): row i < 16 is minus row 32 - i, and
/// row i > 48 equals row 96 - i.
fn row_value(i: usize, mid: &[f64; 33]) -> f64 {
    match i {
        0..=15 => -mid[16 - i],
        16..=48 => mid[i - 16],
        _ => mid[80 - i],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matrix_symmetries_reproduce_the_definition() {
        // Feed unit vectors and compare every V entry with N[i][k].
        for k in 0..32 {
            let mut s = [0.0f32; 32];
            s[k] = 1.0;
            let mut out = [0.0; 64];
            Synth::matrix(&s, &mut out);
            for (i, &o) in out.iter().enumerate() {
                let want = (((16 + i) * (2 * k + 1)) as f64 * std::f64::consts::PI / 64.0).cos();
                assert!((o - want).abs() < 1e-12, "N[{i}][{k}] = {o}, want {want}");
            }
        }
    }

    /// A float reference of the standard's flow chart, written out literally
    /// (shift, matrix with the full 64 x 32 N, build U, window, sum).
    struct Reference {
        v: Vec<f64>,
    }

    impl Reference {
        fn run(&mut self, s: &[f32; 32]) -> [f64; 32] {
            let d = synthesis_window();
            for i in (64..1024).rev() {
                self.v[i] = self.v[i - 64];
            }
            for i in 0..64 {
                self.v[i] = (0..32)
                    .map(|k| (((16 + i) * (2 * k + 1)) as f64 * std::f64::consts::PI / 64.0).cos() * f64::from(s[k]))
                    .sum();
            }
            let mut u = [0.0; 512];
            for i in 0..8 {
                for j in 0..32 {
                    u[64 * i + j] = self.v[128 * i + j];
                    u[64 * i + 32 + j] = self.v[128 * i + 96 + j];
                }
            }
            let mut out = [0.0; 32];
            for (j, o) in out.iter_mut().enumerate() {
                *o = (0..16).map(|i| u[j + 32 * i] * d[j + 32 * i]).sum();
            }
            out
        }
    }

    /// [`Synth::run`] as first written: one output's sum at a time, the
    /// ring indexed entry by entry; the side-by-side form must match it to
    /// the bit.
    fn run_by_outputs(syn: &mut Synth, s: &[f32; 32], pcm: &mut [f32; 32]) {
        let t = tables();
        syn.off = (syn.off + 1024 - 64) % 1024;
        let mut mid = [0.0f64; 33];
        for (r, m) in mid.iter_mut().enumerate() {
            let row = &t.n[r * 32..r * 32 + 32];
            *m = row.iter().zip(s.iter()).map(|(&a, &b)| a * f64::from(b)).sum();
        }
        for i in 0..64 {
            syn.v[(syn.off + i) % 1024] = row_value(i, &mid);
        }
        for (j, out) in pcm.iter_mut().enumerate() {
            let mut sum = 0.0f64;
            for i in 0..8 {
                let a = syn.v[(syn.off + 128 * i + j) % 1024];
                let b = syn.v[(syn.off + 128 * i + 96 + j) % 1024];
                sum += a * t.d[64 * i + j] + b * t.d[64 * i + 32 + j];
            }
            *out = sum as f32;
        }
    }

    #[test]
    fn side_by_side_sums_are_the_one_at_a_time_ones_to_the_bit() {
        let (mut fast, mut slow) = (Synth::default(), Synth::default());
        let mut seed = 4321u32;
        for n in 0..300 {
            let mut s = [0.0f32; 32];
            for (k, x) in s.iter_mut().enumerate() {
                seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12345);
                *x = match n % 4 {
                    0 if k > 12 => 0.0,
                    1 => -0.0,
                    _ => ((seed >> 8) as f32 / (1u32 << 24) as f32 - 0.5) * 0.7,
                };
            }
            let (mut a, mut b) = ([0.0f32; 32], [0.0f32; 32]);
            fast.run(&s, &mut a);
            run_by_outputs(&mut slow, &s, &mut b);
            assert_eq!(a.map(f32::to_bits), b.map(f32::to_bits), "set {n}");
        }
    }

    #[test]
    fn matches_the_flow_chart() {
        let mut fast = Synth::default();
        let mut slow = Reference { v: vec![0.0; 1024] };
        let mut seed = 12345u32;
        let mut worst = 0.0f64;
        for _ in 0..200 {
            let mut s = [0.0f32; 32];
            for x in s.iter_mut() {
                seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12345);
                *x = ((seed >> 8) as f32 / (1u32 << 24) as f32 - 0.5) * 0.5;
            }
            let mut pcm = [0.0f32; 32];
            fast.run(&s, &mut pcm);
            let want = slow.run(&s);
            for j in 0..32 {
                worst = worst.max((f64::from(pcm[j]) - want[j]).abs());
            }
        }
        assert!(worst < 1e-6, "worst difference {worst}");
    }

    #[test]
    fn a_single_subband_tone_comes_out_at_its_frequency() {
        // A constant in subband 0 is DC-ish (passband 0..fs/64); its output
        // reaches a steady level with gain ~1 for the right input sign.
        let mut syn = Synth::default();
        let mut s = [0.0f32; 32];
        s[0] = 1.0;
        let mut last = [0.0f32; 32];
        for _ in 0..40 {
            syn.run(&s, &mut last);
        }
        for &x in &last {
            assert!((x - last[0]).abs() < 0.01 * last[0].abs().max(1e-3), "{last:?}");
        }
        assert!(last[0].abs() > 0.5, "{}", last[0]);
    }
}
