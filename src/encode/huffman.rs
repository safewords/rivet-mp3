//! Huffman coding of quantised Layer III spectra: the division into
//! big_values pairs, count1 quadruples and zeros, the choice of regions and
//! tables, bit counting and writing (the inverse of ISO/IEC 11172-3
//! 2.4.3.4.4, with the tables of Annex B, Table 3-B.7).

use crate::bits::BitWriter;
use crate::tables::huffman::{QUAD_A_CODES, QUAD_A_LENS, QUAD_B_CODES, QUAD_B_LENS, TABLE_INFO, pair_table};

/// Largest magnitude a table can code: xlen - 1 without linbits, 15 +
/// 2^linbits - 1 with.
fn capacity(table: usize) -> u32 {
    match TABLE_INFO[table] {
        Some((0, _)) => 0,
        Some((base, 0)) => pair_table(base).expect("table").xlen as u32 - 1,
        Some((_, lin)) => 15 + (1u32 << lin) - 1,
        None => 0,
    }
}

/// Bits of one pair (absolute values `x`, `y`) in `table`, signs included.
fn pair_bits(table: usize, x: u32, y: u32) -> u32 {
    let Some((base, lin)) = TABLE_INFO[table] else { return u32::MAX / 4 };
    if base == 0 {
        return 0;
    }
    let t = pair_table(base).expect("table");
    let (cx, cy) = (x.min(15), y.min(15));
    let mut bits = u32::from(t.lens[cx as usize * t.xlen + cy as usize]);
    if lin > 0 {
        if x >= 15 {
            bits += lin;
        }
        if y >= 15 {
            bits += lin;
        }
    }
    bits + u32::from(x != 0) + u32::from(y != 0)
}

/// The tables worth considering: every code table except the unused 4
/// and 14 (and 0, which only zeros use).
const CANDIDATES: [usize; 29] = [
    1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
];

/// How a granule's spectrum is to be Huffman coded.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct Coding {
    pub(crate) big_values: u16,
    pub(crate) table_select: [u8; 3],
    pub(crate) region0_count: u8,
    pub(crate) region1_count: u8,
    pub(crate) count1_table_b: bool,
    /// Lines in the count1 region (a multiple of 4).
    pub(crate) count1_lines: u16,
    /// Bits of the big_values and count1 regions together.
    pub(crate) bits: u32,
}

/// Prefix sums of pair costs per table over a granule's pairs; a pair a
/// table cannot code costs [`INF`].
struct Costs {
    /// [candidate index][pair] cumulative bits, pairs + 1 entries each.
    prefix: Vec<Vec<u64>>,
}

const INF: u64 = 1 << 32;

impl Costs {
    fn new(abs: &[u32; 576], pairs: usize) -> Costs {
        let prefix = CANDIDATES
            .iter()
            .map(|&t| {
                let cap = capacity(t);
                let mut v = Vec::with_capacity(pairs + 1);
                v.push(0u64);
                // A table too small for one value in the granule still
                // serves regions without it.
                let mut acc = 0u64;
                for p in 0..pairs {
                    let (x, y) = (abs[2 * p], abs[2 * p + 1]);
                    acc += if x.max(y) > cap { INF } else { u64::from(pair_bits(t, x, y)) };
                    v.push(acc);
                }
                v
            })
            .collect();
        Costs { prefix }
    }

    /// The cheapest table for lines [a, b) (pairs a/2..b/2) and its bits;
    /// table 0 when every value there is zero.
    fn best(&self, a: usize, b: usize) -> (u8, u32) {
        let (pa, pb) = (a / 2, b / 2);
        if pa >= pb {
            return (0, 0);
        }
        let mut best = (0u8, u64::MAX);
        for (i, &t) in CANDIDATES.iter().enumerate() {
            let bits = self.prefix[i][pb] - self.prefix[i][pa];
            if bits < best.1 {
                best = (t as u8, bits);
            }
        }
        // All zero: table 0 codes it in no bits. (Every candidate counts
        // at least one bit per pair, so zero cost cannot come from them.)
        let zero = (pa..pb).all(|p| self.is_zero_pair(p));
        if zero {
            return (0, 0);
        }
        (best.0, best.1.min(u64::from(u32::MAX / 2)) as u32)
    }

    fn is_zero_pair(&self, p: usize) -> bool {
        // Table 1 codes (0, 0) in one bit and every other pair in more, with
        // no sign bits: a pair costs exactly 1 there only when it is zero.
        // (CANDIDATES[0] is table 1.)
        self.prefix[0][p + 1] - self.prefix[0][p] == 1
    }
}

/// Split a quantised granule (`ix`, signed) and choose regions and tables.
/// `sfb_long` are the long band boundaries; `short` selects the window-
/// switching region layout (region 1 from `short_region1`, no region 2).
pub(crate) fn choose(ix: &[i32; 576], sfb_long: &[u16; 23], window_switching: bool, region1_ws: usize) -> Coding {
    let mut abs = [0u32; 576];
    for (a, &v) in abs.iter_mut().zip(ix.iter()) {
        *a = v.unsigned_abs();
    }
    let mut n = 576;
    while n > 0 && abs[n - 1] == 0 {
        n -= 1;
    }
    let end = n.div_ceil(2) * 2;
    let mut c1 = end;
    while c1 >= 4 && abs[c1 - 4..c1].iter().all(|&v| v <= 1) {
        c1 -= 4;
    }
    // count1: try both tables.
    let (mut bits_a, mut bits_b) = (0u32, 0u32);
    for q in abs[c1..end].as_chunks::<4>().0 {
        let idx = (q[0] << 3 | q[1] << 2 | q[2] << 1 | q[3]) as usize;
        let signs = q.iter().filter(|&&v| v != 0).count() as u32;
        bits_a += u32::from(QUAD_A_LENS[idx]) + signs;
        bits_b += u32::from(QUAD_B_LENS[idx]) + signs;
    }
    let count1_table_b = bits_b < bits_a;
    let count1_bits = bits_a.min(bits_b);
    let pairs = c1 / 2;
    let costs = Costs::new(&abs, pairs);
    let mut coding = Coding {
        big_values: pairs as u16,
        count1_table_b,
        count1_lines: (end - c1) as u16,
        ..Default::default()
    };
    if window_switching {
        let r1 = region1_ws.min(c1);
        let (t0, b0) = costs.best(0, r1);
        let (t1, b1) = costs.best(r1, c1);
        coding.table_select = [t0, t1, 0];
        coding.bits = b0 + b1 + count1_bits;
        return coding;
    }
    let mut best = (u32::MAX, 0u8, 0u8, [0u8; 3]);
    for r0 in 0..16usize {
        let a = usize::from(sfb_long[r0 + 1]).min(c1);
        let (t0, b0) = costs.best(0, a);
        if b0 >= best.0 {
            continue;
        }
        for r1 in 0..8usize {
            let b = usize::from(sfb_long[(r0 + r1 + 2).min(22)]).min(c1);
            let (t1, b1) = costs.best(a, b);
            let (t2, b2) = costs.best(b, c1);
            let total = b0 + b1 + b2;
            if total < best.0 {
                best = (total, r0 as u8, r1 as u8, [t0, t1, t2]);
            }
        }
    }
    coding.region0_count = best.1;
    coding.region1_count = best.2;
    coding.table_select = best.3;
    coding.bits = best.0 + count1_bits;
    coding
}

/// Write the big_values and count1 regions of `ix` as `coding` says.
pub(crate) fn write(
    w: &mut BitWriter,
    ix: &[i32; 576],
    coding: &Coding,
    sfb_long: &[u16; 23],
    window_switching: bool,
    region1_ws: usize,
) {
    let big = usize::from(coding.big_values) * 2;
    let (r1, r2) = if window_switching {
        (region1_ws.min(big), big)
    } else {
        (
            usize::from(sfb_long[usize::from(coding.region0_count) + 1]).min(big),
            usize::from(sfb_long[(usize::from(coding.region0_count) + usize::from(coding.region1_count) + 2).min(22)])
                .min(big),
        )
    };
    let regions = [(0, r1, coding.table_select[0]), (r1, r2, coding.table_select[1]), (r2, big, coding.table_select[2])];
    for (a, b, table) in regions {
        let Some((base, lin)) = TABLE_INFO[usize::from(table)] else { continue };
        if base == 0 || a >= b {
            continue;
        }
        let t = pair_table(base).expect("table");
        for p in (a..b).step_by(2) {
            let (x, y) = (ix[p], ix[p + 1]);
            let (ax, ay) = (x.unsigned_abs(), y.unsigned_abs());
            let (cx, cy) = (ax.min(15) as usize, ay.min(15) as usize);
            let k = cx * t.xlen + cy;
            w.put(u32::from(t.codes[k]), u32::from(t.lens[k]));
            if lin > 0 && ax >= 15 {
                w.put(ax - 15, lin);
            }
            if ax != 0 {
                w.put(u32::from(x < 0), 1);
            }
            if lin > 0 && ay >= 15 {
                w.put(ay - 15, lin);
            }
            if ay != 0 {
                w.put(u32::from(y < 0), 1);
            }
        }
    }
    let (codes, lens) = if coding.count1_table_b { (&QUAD_B_CODES, &QUAD_B_LENS) } else { (&QUAD_A_CODES, &QUAD_A_LENS) };
    for q in ix[big..big + usize::from(coding.count1_lines)].as_chunks::<4>().0 {
        let idx = (q[0].unsigned_abs() << 3 | q[1].unsigned_abs() << 2 | q[2].unsigned_abs() << 1 | q[3].unsigned_abs())
            as usize;
        w.put(u32::from(codes[idx]), u32::from(lens[idx]));
        for &v in q {
            if v != 0 {
                w.put(u32::from(v < 0), 1);
            }
        }
    }
}

/// Largest magnitude any table can code (15 + 2^13 - 1).
pub(crate) const MAX_VALUE: u32 = 8206;
