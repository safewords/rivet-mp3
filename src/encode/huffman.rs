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

/// Bits of one pair (absolute values `x`, `y`) in `table`, signs included
/// (the definition [`Costs`] is tested against).
#[cfg(test)]
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
const CANDIDATES: [usize; 29] =
    [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31];

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

/// The code tables behind the candidates: tables 16–23 share the code of
/// 16 and 24–31 that of 24, differing only in linbits.
const BASES: [usize; 15] = [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 24];

/// Per candidate: its index in [`BASES`], its linbits and the largest
/// magnitude it codes.
fn candidate_info() -> &'static [(usize, u32, u32); 29] {
    static T: std::sync::OnceLock<[(usize, u32, u32); 29]> = std::sync::OnceLock::new();
    T.get_or_init(|| {
        std::array::from_fn(|i| {
            let t = CANDIDATES[i];
            let (base, lin) = TABLE_INFO[t].expect("a candidate is a code table");
            (BASES.iter().position(|&b| b == base).expect("a base table"), lin, capacity(t))
        })
    })
}

/// The bits of any range of a granule's pairs in any table, from prefix
/// sums: a pair costs its code word in the table's code (shared by the
/// tables of one base), linbits for each value of 15 or more, and a sign
/// bit for each nonzero value; a table that cannot code the range's
/// largest value is out. The same costs as pricing every pair in every
/// table ([`pair_bits`]), which the unit test holds it to, at a fraction
/// of the work.
struct Costs {
    /// Cumulative code-word bits, per pair boundary and base table (lane
    /// `i` is `BASES[i]`, lane 15 unused); pairs a table cannot code count
    /// 0 (the range maximum rules the table out there). Sixteen u16 lanes
    /// per pair, so a pair's step is one vector add.
    words: [[u16; 16]; 289],
    /// Cumulative values of 15 or more (each takes linbits).
    escapes: [u32; 289],
    /// Cumulative nonzero values (each takes a sign bit).
    signs: [u32; 289],
    /// Range maximum: `max[l][p]` is the largest value of pairs
    /// `p..p + 2^l`.
    max: Vec<[u32; 288]>,
}

/// Code-word lengths by `min(x, 15) * 16 + min(y, 15)`, one lane per base
/// table (0 where the table has no such pair).
fn word_lengths() -> &'static [[u16; 16]; 256] {
    static T: std::sync::OnceLock<[[u16; 16]; 256]> = std::sync::OnceLock::new();
    T.get_or_init(|| {
        let mut t = [[0u16; 16]; 256];
        for (i, &base) in BASES.iter().enumerate() {
            let pt = pair_table(base).expect("table");
            for x in 0..pt.xlen {
                for y in 0..pt.xlen {
                    t[x * 16 + y][i] = u16::from(pt.lens[x * pt.xlen + y]);
                }
            }
        }
        t
    })
}

impl Costs {
    fn new(abs: &[u32; 576], pairs: usize) -> Costs {
        let lens = word_lengths();
        let mut words = [[0u16; 16]; 289];
        let mut escapes = [0u32; 289];
        let mut signs = [0u32; 289];
        let mut level0 = [0u32; 288];
        for p in 0..pairs {
            let (x, y) = (abs[2 * p], abs[2 * p + 1]);
            level0[p] = x.max(y);
            escapes[p + 1] = escapes[p] + u32::from(x >= 15) + u32::from(y >= 15);
            signs[p + 1] = signs[p] + u32::from(x != 0) + u32::from(y != 0);
            let row = &lens[(x.min(15) * 16 + y.min(15)) as usize];
            let prev = words[p];
            for ((w, &a), &l) in words[p + 1].iter_mut().zip(&prev).zip(row) {
                // At most 19 bits a pair, 288 pairs: no overflow.
                *w = a + l;
            }
        }
        let mut max = vec![level0];
        let mut span = 1;
        while 2 * span <= pairs {
            let prev = max.last().expect("level 0");
            let mut next = [0u32; 288];
            for p in 0..=pairs - 2 * span {
                next[p] = prev[p].max(prev[p + span]);
            }
            max.push(next);
            span *= 2;
        }
        Costs { words, escapes, signs, max }
    }

    /// The largest value in pairs `pa..pb` (non-empty).
    fn range_max(&self, pa: usize, pb: usize) -> u32 {
        let l = (usize::BITS - 1 - (pb - pa).leading_zeros()) as usize;
        self.max[l][pa].max(self.max[l][pb - (1 << l)])
    }

    /// The cheapest table for lines [a, b) (pairs a/2..b/2) and its bits;
    /// table 0 when every value there is zero.
    fn best(&self, a: usize, b: usize) -> (u8, u32) {
        let (pa, pb) = (a / 2, b / 2);
        if pa >= pb {
            return (0, 0);
        }
        let top = self.range_max(pa, pb);
        // All zero: table 0 codes it in no bits.
        if top == 0 {
            return (0, 0);
        }
        let escapes = self.escapes[pb] - self.escapes[pa];
        let signs = self.signs[pb] - self.signs[pa];
        let mut best = (0u8, u32::MAX);
        // The candidates sharing a base table are adjacent and in rising
        // linbits (and capacity): within a base, the first that can code
        // the range costs least, and a later one could only tie, which the
        // first already wins. So one candidate per base is priced.
        let mut last_base = usize::MAX;
        for (&t, &(bi, lin, cap)) in CANDIDATES.iter().zip(candidate_info()) {
            if top > cap || bi == last_base {
                continue;
            }
            last_base = bi;
            let bits = u32::from(self.words[pb][bi] - self.words[pa][bi]) + lin * escapes + signs;
            if bits < best.1 {
                best = (t as u8, bits);
            }
        }
        (best.0, best.1.min(u32::MAX / 2))
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
    let mut coding =
        Coding { big_values: pairs as u16, count1_table_b, count1_lines: (end - c1) as u16, ..Default::default() };
    if window_switching {
        let r1 = region1_ws.min(c1);
        let (t0, b0) = costs.best(0, r1);
        let (t1, b1) = costs.best(r1, c1);
        coding.table_select = [t0, t1, 0];
        coding.bits = b0 + b1 + count1_bits;
        return coding;
    }
    let mut best = (u32::MAX, 0u8, 0u8, [0u8; 3]);
    // Region 2 runs from a band edge to the end: one price per edge.
    let mut tail: [Option<(u8, u32)>; 23] = [None; 23];
    for r0 in 0..16usize {
        let a = usize::from(sfb_long[r0 + 1]).min(c1);
        let (t0, b0) = costs.best(0, a);
        if b0 >= best.0 {
            continue;
        }
        for r1 in 0..8usize {
            let edge = (r0 + r1 + 2).min(22);
            let b = usize::from(sfb_long[edge]).min(c1);
            let (t1, b1) = costs.best(a, b);
            let (t2, b2) = *tail[edge].get_or_insert_with(|| costs.best(b, c1));
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
    let regions =
        [(0, r1, coding.table_select[0]), (r1, r2, coding.table_select[1]), (r2, big, coding.table_select[2])];
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
    let (codes, lens) =
        if coding.count1_table_b { (&QUAD_B_CODES, &QUAD_B_LENS) } else { (&QUAD_A_CODES, &QUAD_A_LENS) };
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

#[cfg(test)]
mod tests {
    use super::*;

    /// The definition: every pair priced in every table, a table that
    /// cannot code a pair of the range costing "infinity".
    fn best_by_pairs(abs: &[u32; 576], a: usize, b: usize) -> (u8, u32) {
        let (pa, pb) = (a / 2, b / 2);
        if pa >= pb || (pa..pb).all(|p| abs[2 * p] == 0 && abs[2 * p + 1] == 0) {
            return (0, 0);
        }
        let mut best = (0u8, u64::MAX);
        for &t in &CANDIDATES {
            let cap = capacity(t);
            let bits: u64 = (pa..pb)
                .map(|p| {
                    let (x, y) = (abs[2 * p], abs[2 * p + 1]);
                    if x.max(y) > cap { 1 << 32 } else { u64::from(pair_bits(t, x, y)) }
                })
                .sum();
            if bits < best.1 {
                best = (t as u8, bits);
            }
        }
        (best.0, best.1.min(u64::from(u32::MAX / 2)) as u32)
    }

    #[test]
    fn range_costs_match_pricing_every_pair() {
        let mut seed = 77u32;
        let mut next = move || {
            seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
            seed >> 8
        };
        for case in 0..60 {
            let mut abs = [0u32; 576];
            // Mostly small values with some large ones, like a spectrum
            // falling off with frequency.
            let pairs = 1 + (next() as usize % 288);
            for (i, v) in abs.iter_mut().enumerate().take(2 * pairs) {
                let r = next();
                *v = match (case % 4, r % 100) {
                    (0, _) => r % 2,
                    (_, 0..=2) => r % (MAX_VALUE + 1),
                    (_, 3..=10) => 15 + r % 40,
                    _ => (r % 16) >> (i * 4 / 576),
                };
            }
            let costs = Costs::new(&abs, pairs);
            for _ in 0..200 {
                let x = next() as usize % (2 * pairs + 1);
                let y = next() as usize % (2 * pairs + 1);
                let (a, b) = (x.min(y), x.max(y));
                assert_eq!(costs.best(a, b), best_by_pairs(&abs, a, b), "case {case} lines {a}..{b}");
            }
        }
    }
}
