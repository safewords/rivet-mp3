//! The Layer III IMDCT and its windows (ISO/IEC 11172-3 2.4.3.4.10.2,
//! "IMDCT"), and the forward MDCT the encoder uses with the same windows.
//!
//! x_i = sum_{k=0}^{n/2-1} X_k cos(pi / (2n) (2i + 1 + n/2)(2k + 1)),
//! i = 0..n-1, n = 36 for long blocks and 12 for each of the three short
//! windows. Windows by block_type:
//!
//! - 0 (normal): sin(pi/36 (i + 1/2)), i = 0..35
//! - 1 (start): the normal window's rising half, then 1 for i = 18..23,
//!   sin(pi/12 (i - 18 + 1/2)) for 24..29, 0 for 30..35
//! - 3 (stop): 0 for 0..5, sin(pi/12 (i - 6 + 1/2)) for 6..11, 1 for
//!   12..17, then the normal window's falling half
//! - 2 (short): sin(pi/12 (i + 1/2)), i = 0..11, for each short window; the
//!   three windowed outputs overlap at offsets 6, 12 and 18 of the 36.

use std::sync::OnceLock;

pub(crate) struct Tables {
    /// cos table for n = 36: [i * 18 + k].
    pub(crate) cos36: Vec<f64>,
    /// The same, transposed: [k][i], so the IMDCT runs its 36 outputs side
    /// by side.
    cos36_t: [[f64; 36]; 18],
    /// The n = 12 table transposed: [k][i].
    cos12_t: [[f64; 12]; 6],
    /// [`imdct`] of 18 zeros (+0.0) for each block type.
    zero: [[f64; 36]; 4],
    /// cos table for n = 12: [i * 6 + k].
    pub(crate) cos12: Vec<f64>,
    /// Long windows by block type (index 2 unused: short).
    pub(crate) win: [[f64; 36]; 4],
    pub(crate) win_short: [f64; 12],
}

pub(crate) fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| {
        use std::f64::consts::PI;
        let mut cos36 = vec![0.0; 36 * 18];
        for i in 0..36 {
            for k in 0..18 {
                cos36[i * 18 + k] = (PI / 72.0 * (2 * i + 1 + 18) as f64 * (2 * k + 1) as f64).cos();
            }
        }
        let mut cos12 = vec![0.0; 12 * 6];
        for i in 0..12 {
            for k in 0..6 {
                cos12[i * 6 + k] = (PI / 24.0 * (2 * i + 1 + 6) as f64 * (2 * k + 1) as f64).cos();
            }
        }
        let mut win = [[0.0; 36]; 4];
        for i in 0..36 {
            win[0][i] = (PI / 36.0 * (i as f64 + 0.5)).sin();
        }
        for i in 0..36 {
            win[1][i] = match i {
                0..=17 => (PI / 36.0 * (i as f64 + 0.5)).sin(),
                18..=23 => 1.0,
                24..=29 => (PI / 12.0 * ((i - 18) as f64 + 0.5)).sin(),
                _ => 0.0,
            };
            win[3][i] = match i {
                0..=5 => 0.0,
                6..=11 => (PI / 12.0 * ((i - 6) as f64 + 0.5)).sin(),
                12..=17 => 1.0,
                _ => (PI / 36.0 * (i as f64 + 0.5)).sin(),
            };
        }
        let mut win_short = [0.0; 12];
        for (i, w) in win_short.iter_mut().enumerate() {
            *w = (PI / 12.0 * (i as f64 + 0.5)).sin();
        }
        let cos36_t = std::array::from_fn(|k| std::array::from_fn(|i| cos36[i * 18 + k]));
        let cos12_t = std::array::from_fn(|k| std::array::from_fn(|i| cos12[i * 6 + k]));
        let mut t = Tables { cos36, cos12, cos36_t, cos12_t, zero: [[0.0; 36]; 4], win, win_short };
        for bt in 0..4u8 {
            let mut z = [0.0; 36];
            imdct_with(&t, &[0.0; 18], bt, &mut z);
            t.zero[usize::from(bt)] = z;
        }
        t
    })
}

/// [`imdct`] of 18 zeros for `block_type`: what a subband above the last
/// coded line transforms to (the sign of each zero included).
pub(crate) fn zero_output(block_type: u8) -> [f64; 36] {
    tables().zero[usize::from(block_type)]
}

/// Inverse transform and window one subband's 18 values for `block_type`
/// (for short blocks `x` holds the three windows' six values each,
/// window-major), giving 36 windowed samples.
///
/// Each output's sum runs over k in order from -0.0 (as `Iterator::sum`
/// does) with separate multiplies and adds; the outputs are computed side
/// by side, which vectorises without changing any of them.
pub(crate) fn imdct(x: &[f32], block_type: u8, out: &mut [f64; 36]) {
    imdct_dispatch(tables(), x, block_type, out);
}

crate::simd::multiversion! {
fn imdct_dispatch(t: &Tables, x: &[f32], block_type: u8, out: &mut [f64; 36]) {
    imdct_with(t, x, block_type, out)
}
}

#[inline(always)]
fn imdct_with(t: &Tables, x: &[f32], block_type: u8, out: &mut [f64; 36]) {
    if block_type == 2 {
        out.fill(0.0);
        for w in 0..3 {
            let xs = &x[w * 6..w * 6 + 6];
            let mut s = [-0.0f64; 12];
            for (col, &v) in t.cos12_t.iter().zip(xs) {
                let v = f64::from(v);
                for (s, &c) in s.iter_mut().zip(col) {
                    *s += c * v;
                }
            }
            for i in 0..12 {
                out[6 + 6 * w + i] += s[i] * t.win_short[i];
            }
        }
    } else {
        let win = &t.win[usize::from(block_type)];
        let mut s = [-0.0f64; 36];
        for (col, &v) in t.cos36_t.iter().zip(x) {
            let v = f64::from(v);
            for (s, &c) in s.iter_mut().zip(col) {
                *s += c * v;
            }
        }
        for ((o, &s), &w) in out.iter_mut().zip(&s).zip(win) {
            *o = s * w;
        }
    }
}

/// The forward MDCT matching [`imdct`]: 36 input samples (the previous and
/// the current granule's 18 of one subband) to 18 coefficients, for a long
/// block type; or, for short blocks, the three windows' six coefficients
/// each (window-major). Scaled so that `imdct` followed by overlap-add
/// reconstructs the input: X_k = (2 / (n/2))... see the unit test.
pub(crate) fn mdct(x: &[f64; 36], block_type: u8, out: &mut [f32; 18]) {
    let t = tables();
    if block_type == 2 {
        for w in 0..3 {
            for k in 0..6 {
                let mut s = 0.0;
                for i in 0..12 {
                    s += x[6 + 6 * w + i] * t.win_short[i] * t.cos12[i * 6 + k];
                }
                out[w * 6 + k] = (s * (2.0 / 6.0)) as f32;
            }
        }
    } else {
        let win = &t.win[usize::from(block_type)];
        for k in 0..18 {
            let mut s = 0.0;
            for i in 0..36 {
                s += x[i] * win[i] * t.cos36[i * 18 + k];
            }
            out[k] = (s * (2.0 / 18.0)) as f32;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The IMDCT as first written (one output's dot product at a time),
    /// which the side-by-side form must match to the bit.
    fn imdct_by_rows(x: &[f32], block_type: u8, out: &mut [f64; 36]) {
        let t = tables();
        if block_type == 2 {
            out.fill(0.0);
            for w in 0..3 {
                let xs = &x[w * 6..w * 6 + 6];
                for i in 0..12 {
                    let row = &t.cos12[i * 6..i * 6 + 6];
                    let s: f64 = row.iter().zip(xs).map(|(&c, &v)| c * f64::from(v)).sum();
                    out[6 + 6 * w + i] += s * t.win_short[i];
                }
            }
        } else {
            let win = &t.win[usize::from(block_type)];
            for i in 0..36 {
                let row = &t.cos36[i * 18..i * 18 + 18];
                let s: f64 = row.iter().zip(x).map(|(&c, &v)| c * f64::from(v)).sum();
                out[i] = s * win[i];
            }
        }
    }

    #[test]
    fn side_by_side_imdct_is_the_row_by_row_one_to_the_bit() {
        let mut seed = 9u32;
        for case in 0..400 {
            let x: Vec<f32> = (0..18)
                .map(|k| {
                    seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    match case % 5 {
                        // Zeros of both signs, and sparse spectra.
                        0 => 0.0,
                        1 => -0.0,
                        2 if k % 3 != 0 => 0.0,
                        _ => (seed >> 8) as f32 / (1u32 << 24) as f32 - 0.5,
                    }
                })
                .collect();
            for bt in 0..4u8 {
                let (mut a, mut b) = ([0.0; 36], [0.0; 36]);
                imdct(&x, bt, &mut a);
                imdct_by_rows(&x, bt, &mut b);
                assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits), "case {case} block type {bt}");
                if x.iter().all(|v| v.to_bits() == 0) {
                    assert_eq!(zero_output(bt).map(f64::to_bits), b.map(f64::to_bits));
                }
            }
        }
    }

    /// The IMDCT against its definition, evaluated literally.
    #[test]
    fn imdct_matches_the_formula() {
        use std::f64::consts::PI;
        let x: Vec<f32> = (0..18).map(|k| ((k * 7 % 11) as f32 - 5.0) / 7.0).collect();
        let mut out = [0.0; 36];
        imdct(&x, 0, &mut out);
        for i in 0..36 {
            let want: f64 = (0..18)
                .map(|k| f64::from(x[k]) * (PI / 72.0 * (2 * i + 1 + 18) as f64 * (2 * k + 1) as f64).cos())
                .sum::<f64>()
                * (PI / 36.0 * (i as f64 + 0.5)).sin();
            assert!((out[i] - want).abs() < 1e-12, "{i}");
        }
        imdct(&x, 2, &mut out);
        for i in 0..6 {
            assert_eq!(out[i], 0.0);
            assert_eq!(out[30 + i], 0.0);
        }
        // Window 0 alone contributes to 6..12.
        for i in 0..6 {
            let want: f64 = (0..6)
                .map(|k| f64::from(x[k]) * (PI / 24.0 * (2 * i + 1 + 6) as f64 * (2 * k + 1) as f64).cos())
                .sum::<f64>()
                * (PI / 12.0 * (i as f64 + 0.5)).sin();
            assert!((out[6 + i] - want).abs() < 1e-12, "{i}");
        }
    }

    #[test]
    fn windows_meet_princen_bradley() {
        let t = tables();
        // Long-long, start-short, short-stop and long transitions.
        for i in 0..18 {
            let a = t.win[0][i + 18].powi(2) + t.win[0][i].powi(2);
            assert!((a - 1.0).abs() < 1e-12);
        }
        for i in 0..6 {
            let a = t.win_short[i + 6].powi(2) + t.win_short[i].powi(2);
            assert!((a - 1.0).abs() < 1e-12);
            // Start block's tail against the first short window.
            assert!((t.win[1][24 + i].powi(2) + t.win_short[i].powi(2) - 1.0).abs() < 1e-12);
            // Stop block's head against the last short window.
            assert!((t.win[3][6 + i].powi(2) + t.win_short[6 + i].powi(2) - 1.0).abs() < 1e-12);
        }
    }

    /// MDCT then IMDCT with overlap-add reconstructs the middle granule, for
    /// every legal block-type sequence.
    #[test]
    fn perfect_reconstruction() {
        let signal: Vec<f64> = (0..18 * 6).map(|i| ((i * 37 % 23) as f64 - 11.0) / 13.0).collect();
        for seq in [[0, 0, 0, 0, 0], [0, 1, 2, 3, 0], [0, 1, 2, 2, 3], [3, 0, 1, 2, 3]] {
            let mut prev = [0.0f64; 18];
            let mut rec = vec![0.0; signal.len()];
            for (g, &bt) in seq.iter().enumerate() {
                let mut block = [0.0f64; 36];
                for i in 0..36 {
                    let n = g * 18 + i;
                    block[i] = if n >= 18 && n - 18 < signal.len() { signal[n - 18] } else { 0.0 };
                }
                let mut c = [0.0f32; 18];
                mdct(&block, bt, &mut c);
                let mut y = [0.0; 36];
                imdct(&c, bt, &mut y);
                for i in 0..18 {
                    if g * 18 + i < rec.len() {
                        rec[g * 18 + i] = y[i] + prev[i];
                    }
                    prev[i] = y[18 + i];
                }
            }
            // Granules 1..=3 of the output are fully overlapped: input
            // granules 0..=2.
            for n in 18..18 * 4 {
                assert!((rec[n] - signal[n - 18]).abs() < 1e-5, "{seq:?} {n}: {} vs {}", rec[n], signal[n - 18]);
            }
        }
    }
}
