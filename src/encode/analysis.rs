//! The polyphase analysis subband filterbank of ISO/IEC 11172-3 Annex C
//! (C.1.3, Figure 3-C.1 "analysis subband filter flow chart"):
//!
//! 1. shift 32 new samples into the 512-entry FIFO X (X[0] the newest);
//! 2. window: Z[i] = C[i] X[i], C[i] = D[i] / 32... with the analysis
//!    window C of Table 3-C.1, which is the synthesis window D scaled by
//!    1/32 (the standard prints both; their ratio is 32 to the printed
//!    precision);
//! 3. partial sums: Y[i] = sum_{j=0}^{7} Z[i + 64 j], i = 0..63;
//! 4. matrixing: S[k] = sum_{i=0}^{63} M[k][i] Y[i],
//!    M[k][i] = cos((2k + 1)(i - 16) pi / 64).

use std::sync::OnceLock;

use crate::tables::window::synthesis_window;

struct Tables {
    c: [f64; 512],
    /// M[k][i], row-major.
    m: Vec<f64>,
}

fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| {
        let d = synthesis_window();
        let mut c = [0.0; 512];
        for i in 0..512 {
            c[i] = d[i] / 32.0;
        }
        let mut m = vec![0.0; 32 * 64];
        for k in 0..32 {
            for i in 0..64 {
                m[k * 64 + i] = ((2 * k + 1) as f64 * (i as f64 - 16.0) * std::f64::consts::PI / 64.0).cos();
            }
        }
        Tables { c, m }
    })
}

/// One channel's analysis filterbank.
#[derive(Clone)]
pub(crate) struct Analysis {
    /// X as a ring: X[i] is `x[(head + i) % 512]`.
    x: Box<[f64; 512]>,
    head: usize,
}

impl Default for Analysis {
    fn default() -> Self {
        Self { x: Box::new([0.0; 512]), head: 0 }
    }
}

impl Analysis {
    /// Filter 32 input samples (in time order) into 32 subband samples.
    pub(crate) fn run(&mut self, input: &[f32], out: &mut [f64; 32]) {
        let t = tables();
        for &s in input.iter().take(32) {
            self.head = (self.head + 511) % 512;
            self.x[self.head] = f64::from(s);
        }
        let mut y = [0.0f64; 64];
        for (i, yi) in y.iter_mut().enumerate() {
            let mut sum = 0.0;
            for j in 0..8 {
                let n = i + 64 * j;
                sum += t.c[n] * self.x[(self.head + n) % 512];
            }
            *yi = sum;
        }
        for (k, o) in out.iter_mut().enumerate() {
            let row = &t.m[k * 64..k * 64 + 64];
            *o = row.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::Synth;

    /// Analysis followed by synthesis reproduces the input, delayed by
    /// 481 samples (512 - 31), to within the filterbank's small ripple and
    /// aliasing.
    #[test]
    fn analysis_then_synthesis_is_near_perfect() {
        let n = 32 * 400;
        let input: Vec<f32> = (0..n)
            .map(|i| {
                let t = i as f32;
                0.3 * (t * 0.0123).sin() + 0.2 * (t * 0.31).sin() + 0.1 * (t * 1.7).cos()
            })
            .collect();
        let mut a = Analysis::default();
        let mut s = Synth::default();
        let mut out = vec![0.0f32; n];
        for b in 0..n / 32 {
            let mut sb = [0.0f64; 32];
            a.run(&input[b * 32..b * 32 + 32], &mut sb);
            let sb32: [f32; 32] = std::array::from_fn(|k| sb[k] as f32);
            s.run(&sb32, &mut out[b * 32..b * 32 + 32]);
        }
        let delay = 481;
        let (mut err, mut sig) = (0.0f64, 0.0f64);
        for i in 2000..n - delay {
            let e = f64::from(out[i + delay]) - f64::from(input[i]);
            err += e * e;
            sig += f64::from(input[i]).powi(2);
        }
        let snr = 10.0 * (sig / err).log10();
        assert!(snr > 70.0, "analysis/synthesis SNR {snr:.1} dB");
    }
}
