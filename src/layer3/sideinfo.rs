//! Layer III side information (ISO/IEC 11172-3 2.4.1.7; ISO/IEC 13818-3
//! 2.4.1.7 for the lower sampling frequencies): parsing, and writing for
//! the encoder.

use crate::bits::{BitReader, BitWriter};
use crate::error::{Result, invalid};
use crate::header::{FrameHeader, Mode};

/// One granule of one channel.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GranuleInfo {
    /// Bits of scalefactors and Huffman data.
    pub part2_3_length: u16,
    /// Pairs in the big_values region (at most 288).
    pub big_values: u16,
    /// The quantiser step of the granule.
    pub global_gain: u8,
    /// 4 bits (MPEG-1) or 9 bits (LSF).
    pub scalefac_compress: u16,
    /// A block type other than normal follows.
    pub window_switching: bool,
    /// 0 normal, 1 start, 2 short, 3 stop.
    pub block_type: u8,
    /// The two lowest subbands are long in a short block.
    pub mixed_block: bool,
    /// Huffman table per region.
    pub table_select: [u8; 3],
    /// Short windows: gain offset per window (8 quarter-steps each).
    pub subblock_gain: [u8; 3],
    /// Bands in region 0, less one.
    pub region0_count: u8,
    /// Bands in region 1, less one.
    pub region1_count: u8,
    /// MPEG-1: transmitted; LSF: scalefac_compress >= 500 (not for the
    /// intensity-coded right channel).
    pub preflag: bool,
    /// Scalefactors step by 2^1 rather than 2^0.5.
    pub scalefac_scale: bool,
    /// count1 quadruples use table B.
    pub count1_table_b: bool,
}

/// A frame's side information.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct SideInfo {
    /// Bytes before this frame's main data where it starts (the reservoir).
    pub main_data_begin: usize,
    /// Private bits.
    pub private_bits: u8,
    /// MPEG-1 only: `scfsi[ch][group]`.
    pub scfsi: [[bool; 4]; 2],
    /// `[granule][channel]`; LSF frames use granule 0 only.
    pub gr: [[GranuleInfo; 2]; 2],
}

impl SideInfo {
    /// Parse side information (`b` starts right after the header and CRC).
    pub fn parse(b: &[u8], h: &FrameHeader) -> Result<SideInfo> {
        let mut r = BitReader::new(b);
        let lsf = h.version.is_lsf();
        let nch = h.channels();
        let mut si = SideInfo::default();
        if lsf {
            si.main_data_begin = r.read(8)? as usize;
            si.private_bits = r.read(if nch == 1 { 1 } else { 2 })? as u8;
        } else {
            si.main_data_begin = r.read(9)? as usize;
            si.private_bits = r.read(if nch == 1 { 5 } else { 3 })? as u8;
            for ch in 0..nch {
                for g in 0..4 {
                    si.scfsi[ch][g] = r.bit()?;
                }
            }
        }
        let ngr = if lsf { 1 } else { 2 };
        for gr in 0..ngr {
            for ch in 0..nch {
                let mut g = GranuleInfo {
                    part2_3_length: r.read(12)? as u16,
                    big_values: r.read(9)? as u16,
                    global_gain: r.read(8)? as u8,
                    scalefac_compress: r.read(if lsf { 9 } else { 4 })? as u16,
                    window_switching: r.bit()?,
                    ..Default::default()
                };
                if g.big_values > 288 {
                    return Err(invalid("big_values above 288"));
                }
                if g.window_switching {
                    g.block_type = r.read(2)? as u8;
                    if g.block_type == 0 {
                        return Err(invalid("window switching with block_type 0 is forbidden"));
                    }
                    g.mixed_block = r.bit()?;
                    for t in 0..2 {
                        g.table_select[t] = r.read(5)? as u8;
                    }
                    for w in 0..3 {
                        g.subblock_gain[w] = r.read(3)? as u8;
                    }
                    // Implicit region counts (only used for reporting; the
                    // decoder derives the boundaries directly).
                    g.region0_count = if g.block_type == 2 && !g.mixed_block { 8 } else { 7 };
                    g.region1_count = 36;
                } else {
                    for t in 0..3 {
                        g.table_select[t] = r.read(5)? as u8;
                    }
                    g.region0_count = r.read(4)? as u8;
                    g.region1_count = r.read(3)? as u8;
                }
                if lsf {
                    let intensity_right = ch == 1 && h.mode == Mode::JointStereo && h.mode_extension & 1 == 1;
                    g.preflag = !intensity_right && g.scalefac_compress >= 500;
                } else {
                    g.preflag = r.bit()?;
                }
                g.scalefac_scale = r.bit()?;
                g.count1_table_b = r.bit()?;
                for &t in &g.table_select {
                    if t == 4 || t == 14 {
                        return Err(invalid(format!("Huffman table {t} is not used")));
                    }
                }
                si.gr[gr][ch] = g;
            }
        }
        Ok(si)
    }

    /// Write the side information for `h`.
    pub(crate) fn write(&self, h: &FrameHeader, w: &mut BitWriter) {
        let lsf = h.version.is_lsf();
        let nch = h.channels();
        if lsf {
            w.put(self.main_data_begin as u32, 8);
            w.put(u32::from(self.private_bits), if nch == 1 { 1 } else { 2 });
        } else {
            w.put(self.main_data_begin as u32, 9);
            w.put(u32::from(self.private_bits), if nch == 1 { 5 } else { 3 });
            for ch in 0..nch {
                for g in 0..4 {
                    w.put(u32::from(self.scfsi[ch][g]), 1);
                }
            }
        }
        let ngr = if lsf { 1 } else { 2 };
        for gr in 0..ngr {
            for ch in 0..nch {
                let g = &self.gr[gr][ch];
                w.put(u32::from(g.part2_3_length), 12);
                w.put(u32::from(g.big_values), 9);
                w.put(u32::from(g.global_gain), 8);
                w.put(u32::from(g.scalefac_compress), if lsf { 9 } else { 4 });
                w.put(u32::from(g.window_switching), 1);
                if g.window_switching {
                    w.put(u32::from(g.block_type), 2);
                    w.put(u32::from(g.mixed_block), 1);
                    for t in 0..2 {
                        w.put(u32::from(g.table_select[t]), 5);
                    }
                    for s in 0..3 {
                        w.put(u32::from(g.subblock_gain[s]), 3);
                    }
                } else {
                    for t in 0..3 {
                        w.put(u32::from(g.table_select[t]), 5);
                    }
                    w.put(u32::from(g.region0_count), 4);
                    w.put(u32::from(g.region1_count), 3);
                }
                if !lsf {
                    w.put(u32::from(g.preflag), 1);
                }
                w.put(u32::from(g.scalefac_scale), 1);
                w.put(u32::from(g.count1_table_b), 1);
            }
        }
    }
}
