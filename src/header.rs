//! The 32-bit frame header (ISO/IEC 11172-3 2.4.1.3 / 2.4.2.3, extended by
//! ISO/IEC 13818-3 2.4.2.3 for the lower sampling frequencies and, outside
//! either standard, by the "MPEG-2.5" extension for 8–12 kHz).

use crate::error::{Result, invalid};

/// The MPEG audio version, from the header's ID bit and the bit before it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Version {
    /// ISO/IEC 11172-3: 32, 44.1 and 48 kHz.
    Mpeg1,
    /// ISO/IEC 13818-3 lower sampling frequencies: 16, 22.05 and 24 kHz.
    Mpeg2,
    /// The MPEG-2.5 extension: 8, 11.025 and 12 kHz.
    Mpeg25,
}

impl Version {
    /// MPEG-2 and MPEG-2.5 share the lower-sampling-frequency syntax.
    pub fn is_lsf(self) -> bool {
        self != Version::Mpeg1
    }
}

/// The coding layer.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Layer {
    /// Layer I: 384 samples per frame.
    I,
    /// Layer II: 1152 samples per frame.
    II,
    /// Layer III ("MP3"): 1152 samples per frame (576 at the lower rates).
    III,
}

/// The header's mode field.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Mode {
    /// Two independent channels.
    Stereo,
    /// Joint stereo: intensity stereo (all layers) and, in Layer III,
    /// mid/side stereo, per the mode extension.
    JointStereo,
    /// Two channels carrying independent programmes.
    DualChannel,
    /// One channel.
    Mono,
}

/// The bit rates (kbit/s) by layer and bitrate_index; index 0 is free
/// format (15 is forbidden). ISO/IEC 11172-3 2.4.2.3.
pub const BITRATES_MPEG1: [[u32; 15]; 3] = [
    [
        0, 32, 64, 96, 128, 160, 192, 224, 256, 288, 320, 352, 384, 416, 448,
    ],
    [
        0, 32, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 320, 384,
    ],
    [
        0, 32, 40, 48, 56, 64, 80, 96, 112, 128, 160, 192, 224, 256, 320,
    ],
];
/// Bit rates for the lower sampling frequencies (MPEG-2, and MPEG-2.5),
/// ISO/IEC 13818-3 2.4.2.3.
pub const BITRATES_LSF: [[u32; 15]; 3] = [
    [
        0, 32, 48, 56, 64, 80, 96, 112, 128, 144, 160, 176, 192, 224, 256,
    ],
    [0, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160],
    [0, 8, 16, 24, 32, 40, 48, 56, 64, 80, 96, 112, 128, 144, 160],
];

/// Sampling frequencies by version (MPEG-1, MPEG-2, MPEG-2.5) and
/// sampling_frequency index.
pub const SAMPLE_RATES: [[u32; 3]; 3] = [
    [44_100, 48_000, 32_000],
    [22_050, 24_000, 16_000],
    [11_025, 12_000, 8_000],
];

/// A parsed frame header.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrameHeader {
    /// MPEG-1, MPEG-2 (LSF) or MPEG-2.5.
    pub version: Version,
    /// Layer I, II or III.
    pub layer: Layer,
    /// A 16-bit CRC follows the header (protection_bit 0).
    pub crc: bool,
    /// 1–14; 0 is free format.
    pub bitrate_index: u8,
    /// 0–2.
    pub sample_rate_index: u8,
    /// The frame carries one padding slot.
    pub padding: bool,
    /// The private bit.
    pub private: bool,
    /// The channel mode.
    pub mode: Mode,
    /// mode_extension (meaningful in joint stereo).
    pub mode_extension: u8,
    /// The copyright bit.
    pub copyright: bool,
    /// The original/copy bit.
    pub original: bool,
    /// 0 none, 1 50/15 µs, 3 CCITT J.17 (2 is reserved and refused).
    pub emphasis: u8,
}

impl FrameHeader {
    /// Parse the four header bytes. Refuses a missing syncword and every
    /// reserved value (version 01, layer 00, bitrate index 15, sampling
    /// index 3, emphasis 2).
    pub fn parse(b: &[u8]) -> Result<FrameHeader> {
        if b.len() < 4 {
            return Err(invalid("frame header is four bytes"));
        }
        let h = u32::from_be_bytes([b[0], b[1], b[2], b[3]]);
        if h >> 21 != 0x7FF {
            return Err(invalid("no syncword"));
        }
        let version = match (h >> 19) & 3 {
            0 => Version::Mpeg25,
            2 => Version::Mpeg2,
            3 => Version::Mpeg1,
            _ => return Err(invalid("reserved version")),
        };
        let layer = match (h >> 17) & 3 {
            1 => Layer::III,
            2 => Layer::II,
            3 => Layer::I,
            _ => return Err(invalid("reserved layer")),
        };
        let bitrate_index = ((h >> 12) & 15) as u8;
        if bitrate_index == 15 {
            return Err(invalid("forbidden bitrate index"));
        }
        let sample_rate_index = ((h >> 10) & 3) as u8;
        if sample_rate_index == 3 {
            return Err(invalid("reserved sampling frequency"));
        }
        let emphasis = (h & 3) as u8;
        if emphasis == 2 {
            return Err(invalid("reserved emphasis"));
        }
        let mode = match (h >> 6) & 3 {
            0 => Mode::Stereo,
            1 => Mode::JointStereo,
            2 => Mode::DualChannel,
            _ => Mode::Mono,
        };
        Ok(FrameHeader {
            version,
            layer,
            crc: (h >> 16) & 1 == 0,
            bitrate_index,
            sample_rate_index,
            padding: (h >> 9) & 1 == 1,
            private: (h >> 8) & 1 == 1,
            mode,
            mode_extension: ((h >> 4) & 3) as u8,
            copyright: (h >> 3) & 1 == 1,
            original: (h >> 2) & 1 == 1,
            emphasis,
        })
    }

    /// The four header bytes.
    pub fn to_bytes(&self) -> [u8; 4] {
        let version = match self.version {
            Version::Mpeg25 => 0,
            Version::Mpeg2 => 2,
            Version::Mpeg1 => 3,
        };
        let layer = match self.layer {
            Layer::III => 1,
            Layer::II => 2,
            Layer::I => 3,
        };
        let mode = match self.mode {
            Mode::Stereo => 0,
            Mode::JointStereo => 1,
            Mode::DualChannel => 2,
            Mode::Mono => 3,
        };
        let h: u32 = (0x7FF << 21)
            | (version << 19)
            | (layer << 17)
            | (u32::from(!self.crc) << 16)
            | (u32::from(self.bitrate_index) << 12)
            | (u32::from(self.sample_rate_index) << 10)
            | (u32::from(self.padding) << 9)
            | (u32::from(self.private) << 8)
            | (mode << 6)
            | (u32::from(self.mode_extension) << 4)
            | (u32::from(self.copyright) << 3)
            | (u32::from(self.original) << 2)
            | u32::from(self.emphasis);
        h.to_be_bytes()
    }

    /// 1 (mono) or 2.
    pub fn channels(&self) -> usize {
        if self.mode == Mode::Mono { 1 } else { 2 }
    }

    /// The sampling frequency in Hz.
    pub fn sample_rate(&self) -> u32 {
        SAMPLE_RATES[self.version_index()][usize::from(self.sample_rate_index)]
    }

    pub(crate) fn version_index(&self) -> usize {
        match self.version {
            Version::Mpeg1 => 0,
            Version::Mpeg2 => 1,
            Version::Mpeg25 => 2,
        }
    }

    fn layer_index(&self) -> usize {
        match self.layer {
            Layer::I => 0,
            Layer::II => 1,
            Layer::III => 2,
        }
    }

    /// The bit rate in bit/s; 0 for free format.
    pub fn bitrate(&self) -> u32 {
        let table = if self.version.is_lsf() {
            &BITRATES_LSF
        } else {
            &BITRATES_MPEG1
        };
        table[self.layer_index()][usize::from(self.bitrate_index)] * 1000
    }

    /// bitrate_index 0: the rate is fixed but not signalled.
    pub fn is_free_format(&self) -> bool {
        self.bitrate_index == 0
    }

    /// PCM samples per channel the frame decodes to.
    pub fn samples(&self) -> usize {
        match (self.layer, self.version.is_lsf()) {
            (Layer::I, _) => 384,
            (Layer::II, _) => 1152,
            (Layer::III, false) => 1152,
            (Layer::III, true) => 576,
        }
    }

    /// Bytes per slot: 4 in Layer I, 1 otherwise.
    pub fn slot_bytes(&self) -> usize {
        if self.layer == Layer::I { 4 } else { 1 }
    }

    /// Slots in a frame at `bitrate` bit/s, without padding (rounded down).
    pub fn slots_at(&self, bitrate: u32) -> usize {
        let fs = self.sample_rate() as u64;
        let br = u64::from(bitrate);
        (match self.layer {
            Layer::I => 12 * br / fs,
            _ => self.samples() as u64 / 8 * br / fs,
        }) as usize
    }

    /// The frame length in bytes at `bitrate` bit/s (pass [`Self::bitrate`],
    /// or the measured rate of a free-format stream), padding included.
    pub fn frame_len_at(&self, bitrate: u32) -> usize {
        (self.slots_at(bitrate) + usize::from(self.padding)) * self.slot_bytes()
    }

    /// The frame length in bytes; `None` for free format.
    pub fn frame_len(&self) -> Option<usize> {
        (!self.is_free_format()).then(|| self.frame_len_at(self.bitrate()))
    }

    /// Bytes of Layer III side information.
    pub fn side_info_len(&self) -> usize {
        match (self.version.is_lsf(), self.channels()) {
            (false, 1) => 17,
            (false, _) => 32,
            (true, 1) => 9,
            (true, _) => 17,
        }
    }

    /// Header (and CRC word) length in bytes.
    pub fn header_len(&self) -> usize {
        if self.crc { 6 } else { 4 }
    }

    /// Whether `other` can follow this header in one stream: the fields a
    /// stream keeps constant (version, layer, sampling frequency) agree,
    /// and both are free format or neither is.
    pub fn same_stream(&self, other: &FrameHeader) -> bool {
        self.version == other.version
            && self.layer == other.layer
            && self.sample_rate_index == other.sample_rate_index
            && self.is_free_format() == other.is_free_format()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lengths() {
        let h = FrameHeader::parse(&[0xFF, 0xFB, 0x90, 0x64]).unwrap();
        assert_eq!(h.version, Version::Mpeg1);
        assert_eq!(h.layer, Layer::III);
        assert_eq!(h.bitrate(), 128_000);
        assert_eq!(h.sample_rate(), 44_100);
        assert_eq!(h.frame_len(), Some(417));
        assert_eq!(h.mode, Mode::JointStereo);
        assert_eq!(h.to_bytes(), [0xFF, 0xFB, 0x90, 0x64]);
        let h = FrameHeader::parse(&[0xFF, 0xFB, 0xE4, 0x00]).unwrap();
        assert_eq!(h.frame_len(), Some(960));
        // MPEG-2 Layer III 24 kHz 64 kb/s: 72 * 64000 / 24000 = 192.
        let h = FrameHeader::parse(&[0xFF, 0xF3, 0x84, 0xC0]).unwrap();
        assert_eq!(h.version, Version::Mpeg2);
        assert_eq!(h.frame_len(), Some(192));
        assert_eq!(h.samples(), 576);
        // Layer I 48 kHz 384 kb/s: 12 * 384000 / 48000 * 4 = 384.
        let h = FrameHeader::parse(&[0xFF, 0xFF, 0xC4, 0x00]).unwrap();
        assert_eq!(h.layer, Layer::I);
        assert_eq!(h.frame_len(), Some(384));
        // Layer II 44.1 kHz 192 kb/s padded: 144 * 192000 / 44100 + 1 = 627.
        let h = FrameHeader::parse(&[0xFF, 0xFD, 0xA2, 0x00]).unwrap();
        assert_eq!(h.layer, Layer::II);
        assert_eq!(h.frame_len(), Some(627));
        assert!(FrameHeader::parse(&[0xFF, 0xFB, 0xF0, 0x00]).is_err());
        assert!(FrameHeader::parse(&[0xFF, 0xF9, 0x90, 0x00]).is_err());
        assert!(FrameHeader::parse(&[0xFF, 0xEB, 0x90, 0x00]).is_err());
    }
}
