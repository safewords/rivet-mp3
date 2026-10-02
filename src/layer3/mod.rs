//! Layer III (ISO/IEC 11172-3 2.4.1.7, 2.4.2.7, 2.4.3.4; ISO/IEC 13818-3
//! 2.4.3.2–2.4.3.4 for the lower sampling frequencies).
//!
//! A frame is the header, the side information and main data; the main
//! data of a frame can start in earlier frames (the bit reservoir,
//! `main_data_begin`), so the decoder keeps the tail of previous frames'
//! main data.

pub(crate) mod huffman;
pub(crate) mod imdct;
pub mod sideinfo;

use crate::bits::BitReader;
use crate::crc::frame_crc_bits;
use crate::error::{Result, invalid};
use crate::header::{FrameHeader, Mode};
use crate::tables::layer3::{NR_OF_SFB, PRETAB, SFB_LONG, SFB_SHORT, SLEN, rate_index};
use sideinfo::{GranuleInfo, SideInfo};

/// Bytes of previous main data kept: main_data_begin is at most 511.
const RESERVOIR_MAX: usize = 511;

/// Antialias butterfly coefficients c_i (Table 3-B.9).
const CI: [f64; 8] = [-0.6, -0.535, -0.33, -0.185, -0.095, -0.041, -0.0142, -0.0037];

/// Per-channel scalefactors of one granule.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Scalefactors {
    /// Long bands 0..=21 (21 has none and stays 0).
    pub(crate) l: [u8; 22],
    /// Short bands 0..=12 by window (12 has none).
    pub(crate) s: [[u8; 3]; 13],
    /// LSF intensity stereo: the value of each band's scalefactor that marks
    /// an illegal intensity position (2^slen - 1), long and short.
    pub(crate) max_l: [u8; 22],
    pub(crate) max_s: [u8; 13],
}

/// What the decoder found out about a frame's bit accounting (used by the
/// strict checks in the tests and by [`crate::FrameInfo`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Layer3Accounting {
    /// main_data_begin, bytes.
    pub main_data_begin: usize,
    /// Sum of part2_3_length over the frame's granules and channels, bits.
    pub part2_3_bits: usize,
    /// Bytes of main data this frame carries after its side information.
    pub main_data_bytes: usize,
    /// Every count1 region ended exactly at its part2_3_length (no code
    /// word ran past it), and no Huffman field ran out.
    pub exact: bool,
}

/// Layer III decoder state that persists across frames.
#[derive(Clone)]
pub(crate) struct Layer3 {
    reservoir: Vec<u8>,
    /// IMDCT overlap, [channel][subband][18].
    overlap: Box<[[[f64; 18]; 32]; 2]>,
    /// MPEG-1 granule 0 scalefactors, for scfsi in granule 1.
    prev_sf: [Scalefactors; 2],
}

impl Default for Layer3 {
    fn default() -> Self {
        Self { reservoir: Vec::new(), overlap: Box::new([[[0.0; 18]; 32]; 2]), prev_sf: Default::default() }
    }
}

/// One granule's decoded spectrum for both channels, before stereo
/// processing: requantised values in bitstream order.
struct GranuleData {
    xr: [[f32; 576]; 2],
    /// Index after the last nonzero value, per channel.
    nonzero: [usize; 2],
    sf: [Scalefactors; 2],
}

impl Layer3 {
    pub(crate) fn reset(&mut self) {
        self.reservoir.clear();
        self.reset_overlap();
    }

    /// Clear the IMDCT overlap only (the bit reservoir is kept).
    pub(crate) fn reset_overlap(&mut self) {
        for ch in self.overlap.iter_mut() {
            for sb in ch.iter_mut() {
                sb.fill(0.0);
            }
        }
    }

    /// Decode one frame into subband samples `out[slot][ch][sb]` (18 slots
    /// per granule). Returns the accounting; when the reservoir lacks the
    /// bytes main_data_begin points back to (decoding started mid-stream),
    /// the frame decodes to silence and `Ok(None)` comes back.
    pub(crate) fn decode(
        &mut self,
        frame: &[u8],
        h: &FrameHeader,
        out: &mut Vec<[[f32; 32]; 2]>,
        check_crc: bool,
    ) -> Result<Option<Layer3Accounting>> {
        let nch = h.channels();
        let lsf = h.version.is_lsf();
        let ngr = if lsf { 1 } else { 2 };
        let si_start = h.header_len();
        let si_len = h.side_info_len();
        if frame.len() < si_start + si_len {
            return Err(invalid("frame shorter than its side information"));
        }
        if check_crc && h.crc {
            let crc = frame_crc_bits(0xFFFF, frame, 16, 16);
            let crc = frame_crc_bits(crc, frame, si_start * 8, si_len * 8);
            let stored = u16::from_be_bytes([frame[4], frame[5]]);
            if crc != stored {
                return Err(invalid(format!("CRC mismatch: frame says {stored:04x}, data gives {crc:04x}")));
            }
        }
        let si = SideInfo::parse(&frame[si_start..si_start + si_len], h)?;
        let main = &frame[si_start + si_len..];
        let mut acct = Layer3Accounting {
            main_data_begin: si.main_data_begin,
            part2_3_bits: 0,
            main_data_bytes: main.len(),
            exact: true,
        };
        out.clear();
        out.resize(18 * ngr, [[0.0; 32]; 2]);

        let have = self.reservoir.len();
        let ok = si.main_data_begin <= have;
        let mut data = Vec::with_capacity(si.main_data_begin + main.len());
        if ok {
            data.extend_from_slice(&self.reservoir[have - si.main_data_begin..]);
        }
        data.extend_from_slice(main);
        // Keep the reservoir for the next frame.
        self.reservoir.extend_from_slice(main);
        if self.reservoir.len() > RESERVOIR_MAX {
            let cut = self.reservoir.len() - RESERVOIR_MAX;
            self.reservoir.drain(..cut);
        }
        if !ok {
            // Not enough history: output silence, but keep the overlap and
            // synthesis running so the next frames join cleanly.
            for gr in 0..ngr {
                for ch in 0..nch {
                    let xr = [0.0f32; 576];
                    let mut slots = [[0.0f32; 32]; 18];
                    self.hybrid(&xr, &si.gr[gr][ch], ch, h, &mut slots);
                    for (t, s) in slots.iter().enumerate() {
                        out[gr * 18 + t][ch] = *s;
                    }
                }
            }
            return Ok(None);
        }

        let rate = rate_index(h.version_index(), usize::from(h.sample_rate_index));
        let mut pos = 0usize;
        for gr in 0..ngr {
            let mut g = GranuleData { xr: [[0.0; 576]; 2], nonzero: [0; 2], sf: Default::default() };
            for ch in 0..nch {
                let gi = &si.gr[gr][ch];
                let end = pos + gi.part2_3_length as usize;
                acct.part2_3_bits += gi.part2_3_length as usize;
                if end > data.len() * 8 {
                    return Err(invalid("part2_3_length runs past the main data"));
                }
                let mut r = BitReader::at(&data, pos);
                r.set_end(end);
                let intensity_right = ch == 1 && h.mode == Mode::JointStereo && h.mode_extension & 1 == 1;
                let sf = if lsf {
                    read_scalefactors_lsf(&mut r, gi, intensity_right)?
                } else {
                    let prev = self.prev_sf[ch];
                    let sf = read_scalefactors_mpeg1(&mut r, gi, &si.scfsi[ch], gr, &prev)?;
                    self.prev_sf[ch] = sf;
                    sf
                };
                let mut is = [0i32; 576];
                let (nonzero, exact) = read_huffman(&mut r, gi, rate, &mut is)?;
                acct.exact &= exact;
                requantise(&is, gi, &sf, rate, nonzero, &mut g.xr[ch]);
                g.nonzero[ch] = nonzero;
                g.sf[ch] = sf;
                pos = end;
            }
            if h.mode == Mode::JointStereo && nch == 2 {
                stereo(&mut g, &si.gr[gr], h, rate);
            }
            for ch in 0..nch {
                let mut slots = [[0.0f32; 32]; 18];
                self.hybrid(&g.xr[ch], &si.gr[gr][ch], ch, h, &mut slots);
                for (t, s) in slots.iter().enumerate() {
                    out[gr * 18 + t][ch] = *s;
                }
            }
        }
        Ok(Some(acct))
    }

    /// Reorder, antialias, IMDCT with overlap-add and frequency inversion:
    /// one granule of one channel to 18 slots of 32 subband samples.
    fn hybrid(&mut self, xr: &[f32; 576], gi: &GranuleInfo, ch: usize, h: &FrameHeader, out: &mut [[f32; 32]; 18]) {
        let rate = rate_index(h.version_index(), usize::from(h.sample_rate_index));
        let mut x = *xr;
        let short = gi.block_type == 2;
        if short {
            reorder(&mut x, gi.mixed_block, rate);
        }
        // Antialias: all boundaries for long blocks, only the first for
        // mixed blocks, none for short.
        let bounds = if !short { 31 } else if gi.mixed_block { 1 } else { 0 };
        for sb in 1..=bounds {
            for (i, &c) in CI.iter().enumerate() {
                let cs = 1.0 / (1.0 + c * c).sqrt();
                let ca = c / (1.0 + c * c).sqrt();
                let bu = f64::from(x[18 * sb - 1 - i]);
                let bd = f64::from(x[18 * sb + i]);
                x[18 * sb - 1 - i] = (bu * cs - bd * ca) as f32;
                x[18 * sb + i] = (bd * cs + bu * ca) as f32;
            }
        }
        for sb in 0..32 {
            // The two lowest subbands of a block with mixed_block_flag set
            // use the normal window, whatever the block type (the flag can
            // accompany start and stop blocks too, and the reference
            // decoder's output for l3_10203 has it so).
            let bt = if gi.mixed_block && sb < 2 { 0 } else { gi.block_type };
            let mut y = [0.0f64; 36];
            imdct::imdct(&x[sb * 18..sb * 18 + 18], bt, &mut y);
            let ov = &mut self.overlap[ch][sb];
            for i in 0..18 {
                let mut v = y[i] + ov[i];
                ov[i] = y[18 + i];
                // Frequency inversion: odd time samples of odd subbands.
                if sb % 2 == 1 && i % 2 == 1 {
                    v = -v;
                }
                out[i][sb] = v as f32;
            }
        }
    }
}

/// MPEG-1 scalefactors (2.4.2.7, scale_factors()).
fn read_scalefactors_mpeg1(
    r: &mut BitReader,
    gi: &GranuleInfo,
    scfsi: &[bool; 4],
    gr: usize,
    prev: &Scalefactors,
) -> Result<Scalefactors> {
    let (slen1, slen2) = SLEN[usize::from(gi.scalefac_compress)];
    let (slen1, slen2) = (u32::from(slen1), u32::from(slen2));
    let mut sf = Scalefactors::default();
    if gi.block_type == 2 {
        if gi.mixed_block {
            for b in 0..8 {
                sf.l[b] = r.read(slen1)? as u8;
            }
            for b in 3..6 {
                for w in 0..3 {
                    sf.s[b][w] = r.read(slen1)? as u8;
                }
            }
        } else {
            for b in 0..6 {
                for w in 0..3 {
                    sf.s[b][w] = r.read(slen1)? as u8;
                }
            }
        }
        for b in 6..12 {
            for w in 0..3 {
                sf.s[b][w] = r.read(slen2)? as u8;
            }
        }
    } else {
        const GROUPS: [(usize, usize); 4] = [(0, 6), (6, 11), (11, 16), (16, 21)];
        for (g, &(a, b)) in GROUPS.iter().enumerate() {
            let len = if g < 2 { slen1 } else { slen2 };
            for band in a..b {
                sf.l[band] = if gr == 1 && scfsi[g] { prev.l[band] } else { r.read(len)? as u8 };
            }
        }
    }
    Ok(sf)
}

/// LSF scalefactors (13818-3 2.4.3.2): scalefac_compress (9 bits) gives
/// four lengths and a row of [`NR_OF_SFB`]; for the right channel of an
/// intensity-stereo frame the derivation differs.
fn read_scalefactors_lsf(r: &mut BitReader, gi: &GranuleInfo, intensity_right: bool) -> Result<Scalefactors> {
    let sfc = u32::from(gi.scalefac_compress);
    let (slen, row) = if !intensity_right {
        if sfc < 400 {
            ([(sfc >> 4) / 5, (sfc >> 4) % 5, (sfc % 16) >> 2, sfc % 4], 0)
        } else if sfc < 500 {
            let s = sfc - 400;
            ([(s >> 2) / 5, (s >> 2) % 5, s % 4, 0], 1)
        } else {
            let s = sfc - 500;
            ([s / 3, s % 3, 0, 0], 2)
        }
    } else {
        let s = sfc >> 1;
        if s < 180 {
            ([s / 36, (s % 36) / 6, (s % 36) % 6, 0], 3)
        } else if s < 244 {
            let s = s - 180;
            ([(s % 64) >> 4, (s % 16) >> 2, s % 4, 0], 4)
        } else {
            let s = s - 244;
            ([s / 3, s % 3, 0, 0], 5)
        }
    };
    let block_index = match (gi.block_type == 2, gi.mixed_block) {
        (false, _) => 0,
        (true, false) => 1,
        (true, true) => 2,
    };
    let counts = NR_OF_SFB[row][block_index];
    let mut sf = Scalefactors::default();
    // Values in transmission order, then distributed onto the bands.
    let mut vals: Vec<(u8, u8)> = Vec::with_capacity(39);
    for (g, &n) in counts.iter().enumerate() {
        let len = slen[g];
        let max = if len == 0 { 0 } else { ((1u32 << len) - 1) as u8 };
        for _ in 0..n {
            vals.push((r.read(len)? as u8, max));
        }
    }
    let mut it = vals.into_iter();
    let mut next = || it.next().unwrap_or((0, 0));
    match block_index {
        0 => {
            for b in 0..21 {
                (sf.l[b], sf.max_l[b]) = next();
            }
        }
        1 => {
            for b in 0..12 {
                for w in 0..3 {
                    let (v, m) = next();
                    sf.s[b][w] = v;
                    sf.max_s[b] = m;
                }
            }
        }
        _ => {
            for b in 0..6 {
                (sf.l[b], sf.max_l[b]) = next();
            }
            for b in 3..12 {
                for w in 0..3 {
                    let (v, m) = next();
                    sf.s[b][w] = v;
                    sf.max_s[b] = m;
                }
            }
        }
    }
    // The last band has no scalefactor; its intensity position is taken
    // from the band below.
    sf.max_l[21] = sf.max_l[20];
    sf.max_s[12] = sf.max_s[11];
    // preflag for LSF comes from scalefac_compress (the 500.. range).
    Ok(sf)
}

/// Decode the big_values and count1 regions into `is`. Returns the index
/// after the last decoded (possibly nonzero) value and whether the data
/// ended exactly at part2_3_length.
fn read_huffman(r: &mut BitReader, gi: &GranuleInfo, rate: usize, is: &mut [i32; 576]) -> Result<(usize, bool)> {
    let big = (gi.big_values as usize * 2).min(576);
    let (r1, r2) = region_bounds(gi, rate);
    let r1 = r1.min(big);
    let r2 = r2.min(big).max(r1);
    huffman::decode_pairs(r, gi.table_select[0], &mut is[..r1])?;
    huffman::decode_pairs(r, gi.table_select[1], &mut is[r1..r2])?;
    huffman::decode_pairs(r, gi.table_select[2], &mut is[r2..big])?;
    let mut i = big;
    let mut exact = true;
    while i + 4 <= 576 && r.remaining() > 0 {
        let save = r.pos();
        match huffman::decode_quad(r, gi.count1_table_b) {
            Ok(q) => {
                is[i..i + 4].copy_from_slice(&q);
                i += 4;
            }
            Err(_) => {
                // A code word that would run past part2_3_length: the
                // quadruple is discarded.
                exact = false;
                r.skip(r.end() - save).ok();
                break;
            }
        }
    }
    if r.remaining() > 0 {
        // Stuffing after 576 values (legal), but not "exact".
        exact = false;
    }
    let mut nz = i;
    while nz > 0 && is[nz - 1] == 0 {
        nz -= 1;
    }
    Ok((nz, exact))
}

/// region1 and region2 start lines (2.4.2.7 region0_count / region1_count;
/// with window switching: 36 lines for MPEG-1 short blocks — i.e. three
/// short bands — and the eighth long band otherwise; no region 2).
pub(crate) fn region_bounds(gi: &GranuleInfo, rate: usize) -> (usize, usize) {
    if gi.window_switching {
        let r1 = if gi.block_type == 2 && !gi.mixed_block {
            3 * usize::from(SFB_SHORT[rate][3])
        } else {
            usize::from(SFB_LONG[rate][8])
        };
        (r1, 576)
    } else {
        let a = (usize::from(gi.region0_count) + 1).min(22);
        let b = (usize::from(gi.region0_count) + usize::from(gi.region1_count) + 2).min(22);
        (usize::from(SFB_LONG[rate][a]), usize::from(SFB_LONG[rate][b]))
    }
}

/// The long bands of a mixed block: those ending at or below line 36.
pub(crate) fn mixed_long_bands(rate: usize) -> usize {
    (0..22).take_while(|&b| SFB_LONG[rate][b + 1] <= 36).count()
}

/// The first short band of a mixed block.
pub(crate) fn mixed_first_short(rate: usize) -> usize {
    (0..13).find(|&b| 3 * usize::from(SFB_SHORT[rate][b]) >= 36).unwrap_or(3)
}

/// 2^(x/4) for the integer quarter-steps of the requantiser.
fn pow2_quarter(q: i32) -> f64 {
    (f64::from(q) * 0.25).exp2()
}

/// |is|^(4/3) with the sign of `is`.
pub(crate) fn pow43(v: i32) -> f64 {
    let a = f64::from(v.unsigned_abs()).powf(4.0 / 3.0);
    if v < 0 { -a } else { a }
}

/// Requantisation (2.4.3.4.7.1): long bands
/// xr = sign(is)|is|^(4/3) 2^((global_gain - 210)/4) 2^-(sfm (sf + preflag pretab)),
/// short bands 2^((global_gain - 210 - 8 subblock_gain[w])/4) 2^-(sfm sf[w]),
/// sfm = 0.5 (1 + scalefac_scale).
fn requantise(is: &[i32; 576], gi: &GranuleInfo, sf: &Scalefactors, rate: usize, n: usize, xr: &mut [f32; 576]) {
    xr.fill(0.0);
    let shift = if gi.scalefac_scale { 2 } else { 1 }; // sfm in quarter-steps / 2
    let gg = i32::from(gi.global_gain) - 210;
    let long_end = if gi.block_type == 2 {
        if gi.mixed_block { mixed_long_bands(rate) } else { 0 }
    } else {
        22
    };
    // Long part.
    for b in 0..long_end {
        let lo = usize::from(SFB_LONG[rate][b]);
        let hi = usize::from(SFB_LONG[rate][b + 1]).min(n);
        if lo >= hi {
            continue;
        }
        let pre = if gi.preflag { i32::from(PRETAB[b]) } else { 0 };
        let q = gg - 2 * shift * (i32::from(sf.l[b]) + pre);
        let g = pow2_quarter(q);
        for i in lo..hi {
            xr[i] = (pow43(is[i]) * g) as f32;
        }
    }
    if gi.block_type != 2 {
        return;
    }
    // Short part, in bitstream order: band, window, line.
    let first = if gi.mixed_block { mixed_first_short(rate) } else { 0 };
    let mut i = if gi.mixed_block { usize::from(SFB_LONG[rate][long_end]) } else { 0 };
    for b in first..13 {
        let width = usize::from(SFB_SHORT[rate][b + 1] - SFB_SHORT[rate][b]);
        for w in 0..3 {
            let q = gg - 8 * i32::from(gi.subblock_gain[w]) - 2 * shift * i32::from(sf.s[b][w]);
            let g = pow2_quarter(q);
            for _ in 0..width {
                if i < n {
                    xr[i] = (pow43(is[i]) * g) as f32;
                }
                i += 1;
            }
        }
    }
}

/// Short blocks: from bitstream order (band, window, line) to the order
/// the IMDCT takes (subband, window, six lines).
fn reorder(x: &mut [f32; 576], mixed: bool, rate: usize) {
    let src = *x;
    let first = if mixed { mixed_first_short(rate) } else { 0 };
    let mut i = if mixed { usize::from(SFB_LONG[rate][mixed_long_bands(rate)]) } else { 0 };
    for b in first..13 {
        let lo = usize::from(SFB_SHORT[rate][b]);
        let width = usize::from(SFB_SHORT[rate][b + 1]) - lo;
        for w in 0..3 {
            for j in 0..width {
                let f = lo + j;
                let dest = (f / 6) * 18 + w * 6 + f % 6;
                if i < 576 && dest < 576 {
                    x[dest] = src[i];
                }
                i += 1;
            }
        }
    }
}

/// Joint stereo processing for one granule (2.4.3.4.9; 13818-3 2.4.3.4
/// for LSF intensity stereo).
fn stereo(g: &mut GranuleData, gi: &[GranuleInfo; 2], h: &FrameHeader, rate: usize) {
    let ms = h.mode_extension & 2 != 0;
    let intensity = h.mode_extension & 1 != 0;
    let lsf = h.version.is_lsf();
    let inv_sqrt2 = std::f64::consts::FRAC_1_SQRT_2;
    // Which lines are intensity-coded: per line, Some((kl, kr)).
    let mut is_lines: Vec<Option<(f64, f64)>> = vec![None; 576];
    if intensity {
        let right = &gi[1];
        let sf = &g.sf[1];
        let k = |pos: u8, max: u8| -> Option<(f64, f64)> {
            if lsf {
                // An illegal position, (2^slen) - 1, means the band is not
                // intensity-coded.
                if pos == max {
                    return None;
                }
                // 13818-3: io = 2^(-1/4) (intensity_scale 0) or 2^(-1/2).
                let io: f64 = if right.scalefac_compress & 1 == 1 { 0.5f64.sqrt() } else { 0.5f64.sqrt().sqrt() };
                let p = i32::from(pos);
                Some(if p == 0 {
                    (1.0, 1.0)
                } else if p % 2 == 1 {
                    (io.powi((p + 1) / 2), 1.0)
                } else {
                    (1.0, io.powi(p / 2))
                })
            } else {
                if pos >= 7 {
                    return None;
                }
                let ratio = (f64::from(pos) * std::f64::consts::PI / 12.0).tan();
                if pos == 6 {
                    // tan(pi/2): all to the left.
                    return Some((1.0, 0.0));
                }
                Some((ratio / (1.0 + ratio), 1.0 / (1.0 + ratio)))
            }
        };
        let xr_r = &g.xr[1];
        if right.block_type == 2 {
            let first = if right.mixed_block { mixed_first_short(rate) } else { 0 };
            // Bitstream positions of each short (band, window).
            let mut start = if right.mixed_block { usize::from(SFB_LONG[rate][mixed_long_bands(rate)]) } else { 0 };
            let mut pos = vec![[0usize; 3]; 13];
            for b in first..13 {
                let width = usize::from(SFB_SHORT[rate][b + 1] - SFB_SHORT[rate][b]);
                for w in 0..3 {
                    pos[b][w] = start;
                    start += width;
                }
            }
            let mut any_short_nonzero = false;
            for w in 0..3 {
                // The last band of this window with a nonzero right value.
                let mut last: Option<usize> = None;
                for b in first..13 {
                    let width = usize::from(SFB_SHORT[rate][b + 1] - SFB_SHORT[rate][b]);
                    let p = pos[b][w];
                    if (p..(p + width).min(576)).any(|i| xr_r[i] != 0.0) {
                        last = Some(b);
                    }
                }
                any_short_nonzero |= last.is_some();
                let from = last.map_or(first, |b| b + 1);
                for b in from..13 {
                    let width = usize::from(SFB_SHORT[rate][b + 1] - SFB_SHORT[rate][b]);
                    let kk = if b == 12 {
                        last_band(&k, from <= 11, sf.s[11][w], sf.max_s[11])
                    } else {
                        k(sf.s[b][w], sf.max_s[b])
                    };
                    let p = pos[b][w];
                    for line in is_lines.iter_mut().take((p + width).min(576)).skip(p) {
                        *line = kk;
                    }
                }
            }
            if right.mixed_block && !any_short_nonzero {
                let long_bands = mixed_long_bands(rate);
                let nz = g.nonzero[1].min(usize::from(SFB_LONG[rate][long_bands]));
                let from = (0..long_bands).find(|&b| usize::from(SFB_LONG[rate][b]) >= nz).unwrap_or(long_bands);
                for b in from..long_bands {
                    let kk = k(sf.l[b], sf.max_l[b]);
                    for line in
                        is_lines.iter_mut().take(usize::from(SFB_LONG[rate][b + 1])).skip(usize::from(SFB_LONG[rate][b]))
                    {
                        *line = kk;
                    }
                }
            }
        } else {
            let nz = g.nonzero[1];
            let from = (0..22).find(|&b| usize::from(SFB_LONG[rate][b]) >= nz).unwrap_or(22);
            for b in from..22 {
                let kk = if b == 21 {
                    last_band(&k, from <= 20, sf.l[20], sf.max_l[20])
                } else {
                    k(sf.l[b], sf.max_l[b])
                };
                for line in
                    is_lines.iter_mut().take(usize::from(SFB_LONG[rate][b + 1])).skip(usize::from(SFB_LONG[rate][b]))
                {
                    *line = kk;
                }
            }
        }
    }
    /// The highest band (long 21, short 12) has no scalefactor and so no
    /// transmitted intensity position. It continues the position of the
    /// band below when that band is itself intensity-coded; when intensity
    /// coding starts at the highest band, the position is 0. (The standards'
    /// text leaves this open; this is what the ISO reference decoder's
    /// output for the LSF test sequences l3_20245, l3_20246, l3_20248 and
    /// l3_25208 shows, which a decoder carrying the position from band to
    /// band, starting from 0, produces.)
    fn last_band(
        k: &dyn Fn(u8, u8) -> Option<(f64, f64)>,
        below_coded: bool,
        below: u8,
        below_max: u8,
    ) -> Option<(f64, f64)> {
        if below_coded { k(below, below_max) } else { k(0, u8::MAX) }
    }
    let (l, r) = g.xr.split_at_mut(1);
    let (l, r) = (&mut l[0], &mut r[0]);
    for i in 0..576 {
        match is_lines[i] {
            Some((kl, kr)) => {
                let x = f64::from(l[i]);
                l[i] = (x * kl) as f32;
                r[i] = (x * kr) as f32;
            }
            None if ms => {
                let m = f64::from(l[i]);
                let s = f64::from(r[i]);
                l[i] = ((m + s) * inv_sqrt2) as f32;
                r[i] = ((m - s) * inv_sqrt2) as f32;
            }
            None => {}
        }
    }
}
