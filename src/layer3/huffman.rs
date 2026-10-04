//! Huffman decoding of Layer III spectral values (ISO/IEC 11172-3
//! 2.4.2.7 "huffmancodebits", 2.4.3.4.4) from the code tables in
//! [`crate::tables::huffman`].

use std::sync::OnceLock;

use crate::bits::BitReader;
use crate::error::{Result, invalid};
use crate::tables::huffman::{
    QUAD_A_CODES, QUAD_A_LENS, QUAD_B_CODES, QUAD_B_LENS, TABLE_INFO, pair_table,
};

/// A binary decoding tree: `nodes[i]` holds the two children of node `i`;
/// a child >= 0x8000 is a leaf carrying symbol `child - 0x8000`. The first
/// [`LUT_BITS`] levels are also flattened into a table indexed by that
/// many bits, so a code word that short (nearly all of them) decodes in
/// one lookup; longer ones walk on from the node the table reaches.
pub(crate) struct Tree {
    nodes: Vec<[u16; 2]>,
    /// Per `LUT_BITS`-bit prefix: what the walk meets first. The low 16
    /// bits are the symbol (leaf) or node (internal), bits 16..24 the
    /// depth (the leaf's code length, or the 0-based bit at which an empty
    /// branch was met), bits 24.. the kind.
    lut: Vec<u32>,
}

const LEAF: u16 = 0x8000;
const EMPTY: u16 = 0;
/// Bits the lookup table resolves.
const LUT_BITS: u32 = 10;
const KIND_LEAF: u32 = 0;
const KIND_NODE: u32 = 1;
const KIND_EMPTY: u32 = 2;

impl Tree {
    fn build(codes: &[u32], lens: &[u8]) -> Tree {
        let mut nodes: Vec<[u16; 2]> = vec![[EMPTY; 2]];
        for (sym, (&code, &len)) in codes.iter().zip(lens).enumerate() {
            let mut node = 0usize;
            for b in (0..u32::from(len)).rev() {
                let bit = ((code >> b) & 1) as usize;
                if b == 0 {
                    nodes[node][bit] = LEAF + sym as u16;
                } else {
                    if nodes[node][bit] == EMPTY {
                        nodes.push([EMPTY; 2]);
                        nodes[node][bit] = (nodes.len() - 1) as u16;
                    }
                    node = usize::from(nodes[node][bit]);
                }
            }
        }
        let lut = (0..1u32 << LUT_BITS)
            .map(|prefix| {
                let mut node = 0usize;
                for i in 0..LUT_BITS {
                    let bit = ((prefix >> (LUT_BITS - 1 - i)) & 1) as usize;
                    let next = nodes[node][bit];
                    if next >= LEAF {
                        return (KIND_LEAF << 24) | ((i + 1) << 16) | u32::from(next - LEAF);
                    }
                    if next == EMPTY {
                        return (KIND_EMPTY << 24) | (i << 16);
                    }
                    node = usize::from(next);
                }
                (KIND_NODE << 24) | (LUT_BITS << 16) | node as u32
            })
            .collect();
        Tree { nodes, lut }
    }

    /// Decode one symbol.
    pub(crate) fn decode(&self, r: &mut BitReader) -> Result<usize> {
        // Peek 24 bits at once (no code is longer than 19) and walk.
        let avail = r.remaining().min(24) as u32;
        if avail == 0 {
            return Err(invalid("Huffman data ends inside a code word"));
        }
        let bits = r.peek_unchecked(24);
        // The walk's first LUT_BITS steps at once; the outcomes (and the
        // errors) are the bit-by-bit walk's.
        let e = self.lut[(bits >> (24 - LUT_BITS)) as usize];
        let depth = (e >> 16) & 0xFF;
        let (mut node, from) = match e >> 24 {
            KIND_LEAF if depth <= avail => {
                r.skip(depth as usize)?;
                return Ok((e & 0xFFFF) as usize);
            }
            KIND_EMPTY if depth < avail => {
                return Err(invalid("code word matches no Huffman table entry"));
            }
            KIND_NODE if depth < avail => ((e & 0xFFFF) as usize, LUT_BITS),
            _ => return Err(invalid("Huffman data ends inside a code word")),
        };
        for i in from..avail {
            let bit = ((bits >> (23 - i)) & 1) as usize;
            let next = self.nodes[node][bit];
            if next >= LEAF {
                r.skip(i as usize + 1)?;
                return Ok(usize::from(next - LEAF));
            }
            if next == EMPTY {
                return Err(invalid("code word matches no Huffman table entry"));
            }
            node = usize::from(next);
        }
        Err(invalid("Huffman data ends inside a code word"))
    }
}

pub(crate) struct Tables {
    /// Pair trees by code table number (index 0, 4, 14 and 17–23, 25–31
    /// unused: those reuse 16 and 24).
    pairs: Vec<Option<Tree>>,
    pub(crate) quad_a: Tree,
    pub(crate) quad_b: Tree,
}

pub(crate) fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(|| {
        let pairs = (0..=24)
            .map(|base| {
                pair_table(base).map(|t| {
                    let codes: Vec<u32> = t.codes.iter().map(|&c| u32::from(c)).collect();
                    Tree::build(&codes, t.lens)
                })
            })
            .collect();
        let a: Vec<u32> = QUAD_A_CODES.iter().map(|&c| u32::from(c)).collect();
        let b: Vec<u32> = QUAD_B_CODES.iter().map(|&c| u32::from(c)).collect();
        Tables {
            pairs,
            quad_a: Tree::build(&a, &QUAD_A_LENS),
            quad_b: Tree::build(&b, &QUAD_B_LENS),
        }
    })
}

/// Decode `count` pairs with table `table_select` into `out` (2 * count
/// values).
pub(crate) fn decode_pairs(r: &mut BitReader, table_select: u8, out: &mut [i32]) -> Result<()> {
    let Some((base, linbits)) = TABLE_INFO[usize::from(table_select)] else {
        return Err(invalid(format!("Huffman table {table_select} is not used")));
    };
    if base == 0 {
        out.fill(0);
        return Ok(());
    }
    let t = tables();
    let tree = t.pairs[base].as_ref().expect("every base table has a tree");
    let xlen = pair_table(base).expect("base table").xlen;
    for pair in out.as_chunks_mut::<2>().0 {
        let sym = tree.decode(r)?;
        let mut x = (sym / xlen) as i32;
        let mut y = (sym % xlen) as i32;
        if linbits > 0 && x == 15 {
            x += r.read(linbits)? as i32;
        }
        if x != 0 && r.bit()? {
            x = -x;
        }
        if linbits > 0 && y == 15 {
            y += r.read(linbits)? as i32;
        }
        if y != 0 && r.bit()? {
            y = -y;
        }
        pair[0] = x;
        pair[1] = y;
    }
    Ok(())
}

/// Decode one count1 quadruple (v, w, x, y) with table A or B.
pub(crate) fn decode_quad(r: &mut BitReader, table_b: bool) -> Result<[i32; 4]> {
    let t = tables();
    let sym = if table_b {
        t.quad_b.decode(r)?
    } else {
        t.quad_a.decode(r)?
    };
    let mut q = [
        ((sym >> 3) & 1) as i32,
        ((sym >> 2) & 1) as i32,
        ((sym >> 1) & 1) as i32,
        (sym & 1) as i32,
    ];
    for v in q.iter_mut() {
        if *v != 0 && r.bit()? {
            *v = -*v;
        }
    }
    Ok(q)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bits::BitWriter;

    /// The bit-by-bit walk alone (the decoder before the lookup table):
    /// the symbol and bits taken, or `None` for an error.
    fn walk(tree: &Tree, bits: u32, avail: u32) -> Option<(usize, u32)> {
        let mut node = 0usize;
        for i in 0..avail {
            let next = tree.nodes[node][((bits >> (23 - i)) & 1) as usize];
            if next >= LEAF {
                return Some((usize::from(next - LEAF), i + 1));
            }
            if next == EMPTY {
                return None;
            }
            node = usize::from(next);
        }
        None
    }

    #[test]
    fn the_lookup_table_decodes_as_the_walk_does() {
        let t = tables();
        let trees: Vec<&Tree> = t
            .pairs
            .iter()
            .flatten()
            .chain([&t.quad_a, &t.quad_b])
            .collect();
        let mut seed = 31u32;
        for tree in trees {
            for case in 0..3000 {
                seed = seed.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                // Long runs of zeros reach the long code words.
                let bits = if case % 3 == 0 { seed >> 12 } else { seed >> 8 };
                let avail = 1 + case % 24;
                let bytes = (bits << 8).to_be_bytes();
                let mut r = BitReader::new(&bytes);
                r.set_end(avail as usize);
                let got = tree.decode(&mut r).ok().map(|s| (s, r.pos() as u32));
                assert_eq!(
                    got,
                    walk(tree, bits & 0xFF_FFFF, avail),
                    "bits {bits:06x} avail {avail}"
                );
            }
        }
    }

    #[test]
    fn every_code_word_decodes_to_its_own_entry() {
        for base in [1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 24] {
            let t = pair_table(base).unwrap();
            let tree = tables().pairs[base].as_ref().unwrap();
            for (sym, (&c, &l)) in t.codes.iter().zip(t.lens).enumerate() {
                let mut w = BitWriter::new();
                w.put(u32::from(c), u32::from(l));
                let bytes = w.into_bytes();
                let mut r = BitReader::new(&bytes);
                r.set_end(usize::from(l));
                assert_eq!(tree.decode(&mut r).unwrap(), sym, "table {base}");
                assert_eq!(r.pos(), usize::from(l));
            }
        }
    }
}
