//! The informational first frame many MP3 files carry: the Xing / Info
//! header with its LAME extension (encoder delay and padding, the figures
//! gapless playback needs), and Fraunhofer's VBRI header.
//!
//! None of these is part of ISO/IEC 11172-3; the layouts here follow their
//! public descriptions. The frame that carries one is a valid Layer III
//! frame whose main data decodes to silence; a player does not play it.
//!
//! Xing / Info, at the end of the side information of the first frame:
//!
//! | bytes | field |
//! |---|---|
//! | 4 | `Xing` (VBR) or `Info` (CBR) |
//! | 4 | flags: 1 frames, 2 bytes, 4 TOC, 8 quality |
//! | 4 | frames (if flagged): audio frames, the tag frame excluded |
//! | 4 | bytes (if flagged): the file's MPEG audio bytes, tag frame included |
//! | 100 | TOC (if flagged): seek points, 1/256 of the file per percent |
//! | 4 | quality (if flagged) |
//!
//! The LAME extension follows (36 bytes): encoder version (9), revision and
//! VBR method (1), lowpass / 100 Hz (1), peak amplitude (4), radio and
//! audiophile replay gain (2 + 2), encoding flags and ATH type (1),
//! bitrate (1), **encoder delay and padding (12 + 12 bits)**, misc (1),
//! MP3 gain (1), preset and surround (2), music length (4), music CRC (2)
//! and the CRC of the tag frame up to that field (2), both CRC-16/ARC.
//!
//! Gapless: the decoded stream holds `delay + 529` samples before the
//! first input sample (529 is a decoder's own delay, the synthesis
//! filterbank and the IMDCT overlap), and `padding - 529` after the last.
//!
//! VBRI, 32 bytes after the frame header: `VBRI`, version (2), delay (2),
//! quality (2), bytes (4), frames (4), TOC entries (2), scale (2), entry
//! size (2), frames per entry (2), TOC.

use crate::crc::crc16_arc;
use crate::header::{FrameHeader, Layer};

/// The decoder delay a LAME tag's figures assume: 528 samples of
/// filterbank delay plus one.
pub const DECODER_DELAY: u32 = 529;

/// A parsed Xing / Info header.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct XingHeader {
    /// `Info` (constant bit rate) rather than `Xing`.
    pub is_info: bool,
    /// Audio frames in the file, the tag frame excluded.
    pub frames: Option<u32>,
    /// MPEG audio bytes in the file, the tag frame included.
    pub bytes: Option<u32>,
    /// Seek table: entry i is the byte position of i% of the duration, in
    /// 1/256 of `bytes`.
    pub toc: Option<[u8; 100]>,
    /// Encoder quality indicator (0 best .. 100).
    pub quality: Option<u32>,
    /// The LAME extension, when present.
    pub lame: Option<LameTag>,
}

/// The LAME extension of a Xing / Info header.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LameTag {
    /// The nine-byte encoder string, e.g. `LAME3.100`.
    pub encoder: String,
    /// Info tag revision (high nibble of byte 9).
    pub revision: u8,
    /// VBR method (low nibble of byte 9): 1 CBR, 2 ABR, 3–6 VBR, …
    pub vbr_method: u8,
    /// Lowpass frequency in Hz (stored in units of 100 Hz).
    pub lowpass: u32,
    /// Samples added before the first input sample by the encoder.
    pub encoder_delay: u32,
    /// Samples added after the last input sample by the encoder (to fill
    /// the last frame).
    pub padding: u32,
    /// Bytes from the tag frame's first byte to the end of the last frame.
    pub music_length: u32,
    /// CRC-16/ARC of the audio frames after the tag frame.
    pub music_crc: u16,
    /// The stored CRC of the tag frame up to its CRC field.
    pub tag_crc: u16,
    /// Whether `tag_crc` matches the frame.
    pub tag_crc_ok: bool,
}

/// A parsed VBRI header.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct VbriHeader {
    /// Header version.
    pub version: u16,
    /// Encoder delay as the header states it.
    pub delay: u16,
    /// Quality indicator.
    pub quality: u16,
    /// MPEG audio bytes in the file.
    pub bytes: u32,
    /// Audio frames in the file.
    pub frames: u32,
}

/// Either informational header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum InfoHeader {
    /// Xing / Info (with or without a LAME extension).
    Xing(XingHeader),
    /// Fraunhofer VBRI.
    Vbri(VbriHeader),
}

impl InfoHeader {
    /// Audio frames the header says follow it.
    pub fn frames(&self) -> Option<u32> {
        match self {
            InfoHeader::Xing(x) => x.frames,
            InfoHeader::Vbri(v) => Some(v.frames),
        }
    }

    /// (encoder delay, padding) in samples, when the header gives them.
    /// VBRI states a delay and no padding.
    pub fn delay_padding(&self) -> Option<(u32, u32)> {
        match self {
            InfoHeader::Xing(x) => x.lame.as_ref().map(|l| (l.encoder_delay, l.padding)),
            InfoHeader::Vbri(v) => Some((u32::from(v.delay), 0)),
        }
    }
}

fn be32(b: &[u8]) -> u32 {
    u32::from_be_bytes([b[0], b[1], b[2], b[3]])
}

fn be16(b: &[u8]) -> u16 {
    u16::from_be_bytes([b[0], b[1]])
}

/// Byte offset of a Xing / Info header in a Layer III frame.
pub(crate) fn xing_offset(h: &FrameHeader) -> usize {
    h.header_len() + h.side_info_len()
}

/// Look for an informational header in one whole frame.
pub fn parse(frame: &[u8], h: &FrameHeader) -> Option<InfoHeader> {
    if h.layer != Layer::III {
        return None;
    }
    let at = xing_offset(h);
    if let Some(tag) = frame.get(at..at + 8)
        && (&tag[..4] == b"Xing" || &tag[..4] == b"Info")
    {
        return parse_xing(frame, at).map(InfoHeader::Xing);
    }
    let at = 4 + 32;
    if let Some(tag) = frame.get(at..at + 26)
        && &tag[..4] == b"VBRI"
    {
        return Some(InfoHeader::Vbri(VbriHeader {
            version: be16(&tag[4..]),
            delay: be16(&tag[6..]),
            quality: be16(&tag[8..]),
            bytes: be32(&tag[10..]),
            frames: be32(&tag[14..]),
        }));
    }
    None
}

fn parse_xing(frame: &[u8], at: usize) -> Option<XingHeader> {
    let mut x = XingHeader { is_info: &frame[at..at + 4] == b"Info", ..Default::default() };
    let flags = be32(frame.get(at + 4..at + 8)?);
    let mut p = at + 8;
    if flags & 1 != 0 {
        x.frames = Some(be32(frame.get(p..p + 4)?));
        p += 4;
    }
    if flags & 2 != 0 {
        x.bytes = Some(be32(frame.get(p..p + 4)?));
        p += 4;
    }
    if flags & 4 != 0 {
        let mut toc = [0u8; 100];
        toc.copy_from_slice(frame.get(p..p + 100)?);
        x.toc = Some(toc);
        p += 100;
    }
    if flags & 8 != 0 {
        x.quality = Some(be32(frame.get(p..p + 4)?));
        p += 4;
    }
    if let Some(t) = frame.get(p..p + 36) {
        // An encoder string of printable ASCII marks the extension (LAME
        // and encoders that copy its layout write one).
        let enc = &t[..9];
        if enc.iter().take(4).all(|&c| c.is_ascii_alphanumeric()) && enc.iter().all(|&c| c == 0 || (32..127).contains(&c))
        {
            let dp = &t[21..24];
            let stored = be16(&t[34..]);
            x.lame = Some(LameTag {
                encoder: String::from_utf8_lossy(enc).trim_end_matches('\0').to_string(),
                revision: t[9] >> 4,
                vbr_method: t[9] & 15,
                lowpass: u32::from(t[10]) * 100,
                encoder_delay: (u32::from(dp[0]) << 4) | (u32::from(dp[1]) >> 4),
                padding: ((u32::from(dp[1]) & 15) << 8) | u32::from(dp[2]),
                music_length: be32(&t[28..]),
                music_crc: be16(&t[32..]),
                tag_crc: stored,
                tag_crc_ok: crc16_arc(&frame[..p + 34]) == stored,
            });
        }
    }
    Some(x)
}

/// What the encoder puts in its tag frame.
pub(crate) struct TagContents<'a> {
    pub(crate) vbr: bool,
    pub(crate) frames: u32,
    pub(crate) bytes: u32,
    pub(crate) toc: [u8; 100],
    pub(crate) quality: u32,
    pub(crate) encoder: &'a [u8; 9],
    pub(crate) vbr_method: u8,
    pub(crate) lowpass_hz: u32,
    pub(crate) bitrate_kbps: u32,
    pub(crate) delay: u32,
    pub(crate) padding: u32,
    pub(crate) stereo_mode: u8,
    pub(crate) source_rate_code: u8,
    pub(crate) music_crc: u16,
}

/// Fill a tag frame (whose header and zeroed side information are already
/// in `frame`) with a Xing / Info header and LAME extension.
pub(crate) fn write_tag(frame: &mut [u8], h: &FrameHeader, c: &TagContents) {
    let at = xing_offset(h);
    frame[at..at + 4].copy_from_slice(if c.vbr { b"Xing" } else { b"Info" });
    frame[at + 4..at + 8].copy_from_slice(&15u32.to_be_bytes());
    frame[at + 8..at + 12].copy_from_slice(&c.frames.to_be_bytes());
    frame[at + 12..at + 16].copy_from_slice(&c.bytes.to_be_bytes());
    frame[at + 16..at + 116].copy_from_slice(&c.toc);
    frame[at + 116..at + 120].copy_from_slice(&c.quality.to_be_bytes());
    let t = at + 120;
    frame[t..t + 9].copy_from_slice(c.encoder);
    frame[t + 9] = c.vbr_method & 15; // revision 0
    frame[t + 10] = (c.lowpass_hz / 100).min(255) as u8;
    // Peak amplitude and replay gain: not measured (zero).
    frame[t + 19] = 0;
    frame[t + 20] = c.bitrate_kbps.min(255) as u8;
    let d = c.delay.min(4095);
    let p = c.padding.min(4095);
    frame[t + 21] = (d >> 4) as u8;
    frame[t + 22] = (((d & 15) << 4) | (p >> 8)) as u8;
    frame[t + 23] = (p & 255) as u8;
    frame[t + 24] = ((c.source_rate_code & 3) << 6) | ((c.stereo_mode & 7) << 2);
    frame[t + 28..t + 32].copy_from_slice(&c.bytes.to_be_bytes());
    frame[t + 32..t + 34].copy_from_slice(&c.music_crc.to_be_bytes());
    let crc = crc16_arc(&frame[..t + 34]);
    frame[t + 34..t + 36].copy_from_slice(&crc.to_be_bytes());
}
