//! Layer I and II tables: the scalefactors (ISO/IEC 11172-3 Table 3-B.1),
//! the classes of quantization (Table 3-B.4) and the possible quantization
//! per subband (Tables 3-B.2a–d, and ISO/IEC 13818-3 Table B.1 for the
//! lower sampling frequencies).

/// Scalefactor `index` (0..=62): 2^(1 - index/3). Table 3-B.1 prints these
/// values to 14 decimals; the formula gives them exactly.
pub(crate) fn scalefactor(index: usize) -> f32 {
    (2.0f64 * (-(index as f64) / 3.0).exp2()) as f32
}

/// One class of quantization (Table 3-B.4).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct QuantClass {
    /// Number of steps (3, 5, 7, 9, 15, …, 65535).
    pub(crate) steps: u32,
    /// Bits per sample, or per codeword of three grouped samples.
    pub(crate) bits: u32,
    /// Three samples share one codeword.
    pub(crate) grouped: bool,
}

const fn q(steps: u32) -> QuantClass {
    match steps {
        3 => QuantClass {
            steps,
            bits: 5,
            grouped: true,
        },
        5 => QuantClass {
            steps,
            bits: 7,
            grouped: true,
        },
        9 => QuantClass {
            steps,
            bits: 10,
            grouped: true,
        },
        _ => QuantClass {
            steps,
            bits: bits_for(steps),
            grouped: false,
        },
    }
}

const fn bits_for(steps: u32) -> u32 {
    // 2^bits - 1 = steps for the ungrouped classes.
    32 - (steps + 1).leading_zeros() - 1
}

impl QuantClass {
    /// Bits of each sample's code before grouping: ceil(log2(steps)).
    pub(crate) fn sample_bits(&self) -> u32 {
        32 - (self.steps - 1).leading_zeros()
    }

    /// C and D of Table 3-B.4: the dequantised value of code `v` is
    /// C * (v / 2^(nb-1) - 1 + D), nb = [`Self::sample_bits`]. C is
    /// 2^nb / steps; D is (2^nb - steps + 1) / 2^nb, which gives the
    /// printed D for every class (0.5 for 3, 5 and 9 steps).
    pub(crate) fn dequantise(&self, v: u32) -> f32 {
        let nb = self.sample_bits();
        let full = (1u64 << nb) as f64;
        let c = full / f64::from(self.steps);
        let d = (full - f64::from(self.steps) + 1.0) / full;
        let frac = f64::from(v) / (full / 2.0) - 1.0;
        (c * (frac + d)) as f32
    }
}

/// One subband's row of a Layer II allocation table: nbal bits, and the
/// class each allocation index (1..2^nbal) selects.
#[derive(Clone, Copy, Debug)]
pub(crate) struct AllocRow {
    pub(crate) nbal: u32,
    pub(crate) classes: &'static [QuantClass],
}

#[rustfmt::skip]
const R_A0: [QuantClass; 15] = [
    q(3), q(7), q(15), q(31), q(63), q(127), q(255), q(511), q(1023), q(2047), q(4095), q(8191), q(16383),
    q(32767), q(65535),
];
#[rustfmt::skip]
const R_A1: [QuantClass; 15] = [
    q(3), q(5), q(7), q(9), q(15), q(31), q(63), q(127), q(255), q(511), q(1023), q(2047), q(4095), q(8191),
    q(65535),
];
const R_A2: [QuantClass; 7] = [q(3), q(5), q(7), q(9), q(15), q(31), q(65535)];
const R_A3: [QuantClass; 3] = [q(3), q(5), q(65535)];
#[rustfmt::skip]
const R_C0: [QuantClass; 15] = [
    q(3), q(5), q(9), q(15), q(31), q(63), q(127), q(255), q(511), q(1023), q(2047), q(4095), q(8191), q(16383),
    q(32767),
];
const R_C1: [QuantClass; 7] = [q(3), q(5), q(9), q(15), q(31), q(63), q(127)];
#[rustfmt::skip]
const R_L0: [QuantClass; 15] = [
    q(3), q(5), q(7), q(9), q(15), q(31), q(63), q(127), q(255), q(511), q(1023), q(2047), q(4095), q(8191),
    q(16383),
];
const R_L1: [QuantClass; 7] = [q(3), q(5), q(9), q(15), q(31), q(63), q(127)];
const R_L2: [QuantClass; 3] = [q(3), q(5), q(9)];

const fn row(nbal: u32, classes: &'static [QuantClass]) -> AllocRow {
    AllocRow { nbal, classes }
}

/// Which allocation table a Layer II frame uses (11172-3 Annex B, B.2;
/// 13818-3 Annex B, B.1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum AllocTable {
    A,
    B,
    C,
    D,
    Lsf,
}

impl AllocTable {
    /// The table for an MPEG-1 frame: by bit rate per channel and sampling
    /// frequency; free format (`bitrate` 0) counts as a high rate.
    pub(crate) fn select(lsf: bool, sample_rate: u32, bitrate: u32, channels: usize) -> AllocTable {
        if lsf {
            return AllocTable::Lsf;
        }
        let per_channel = bitrate / channels as u32;
        if bitrate != 0 && per_channel <= 48_000 {
            if sample_rate == 32_000 {
                AllocTable::D
            } else {
                AllocTable::C
            }
        } else if (bitrate != 0 && per_channel <= 80_000) || sample_rate == 48_000 {
            AllocTable::A
        } else {
            AllocTable::B
        }
    }

    /// sblimit: the number of subbands coded.
    pub(crate) fn sblimit(self) -> usize {
        match self {
            AllocTable::A => 27,
            AllocTable::B => 30,
            AllocTable::C => 8,
            AllocTable::D => 12,
            AllocTable::Lsf => 30,
        }
    }

    pub(crate) fn row(self, sb: usize) -> AllocRow {
        match self {
            AllocTable::A | AllocTable::B => match sb {
                0..=2 => row(4, &R_A0),
                3..=10 => row(4, &R_A1),
                11..=22 => row(3, &R_A2),
                _ => row(2, &R_A3),
            },
            AllocTable::C | AllocTable::D => match sb {
                0..=1 => row(4, &R_C0),
                _ => row(3, &R_C1),
            },
            AllocTable::Lsf => match sb {
                0..=3 => row(4, &R_L0),
                4..=10 => row(3, &R_L1),
                _ => row(2, &R_L2),
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classes_match_table_3_b_4() {
        let c3 = q(3);
        assert_eq!(c3.sample_bits(), 2);
        // C = 1.33333333333, D = 0.5: codes 0, 1, 2 -> -2/3, 0, 2/3.
        assert!((c3.dequantise(0) + 2.0 / 3.0).abs() < 1e-6);
        assert!(c3.dequantise(1).abs() < 1e-6);
        assert!((c3.dequantise(2) - 2.0 / 3.0).abs() < 1e-6);
        // 7 steps: C = 1.14285714286, D = 0.25.
        let c7 = q(7);
        assert_eq!((c7.bits, c7.grouped), (3, false));
        assert!((c7.dequantise(0) - 1.142_857_1 * (-1.0 + 0.25)).abs() < 1e-6);
        assert_eq!(q(65535).bits, 16);
        assert_eq!(q(9).sample_bits(), 4);
        // Every class is symmetric about its middle code.
        for steps in [
            3, 5, 7, 9, 15, 31, 63, 127, 255, 511, 1023, 2047, 4095, 8191, 16383, 32767, 65535,
        ] {
            let c = q(steps);
            let a = c.dequantise(0);
            let b = c.dequantise(steps - 1);
            assert!((a + b).abs() < 1e-6, "{steps}: {a} {b}");
            assert!(c.dequantise(steps / 2).abs() < 1e-6, "{steps}");
            assert!(b < 1.0, "{steps}");
        }
    }
}
