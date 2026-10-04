//! Decoding: whole frames ([`FrameDecoder`]) or a byte stream in any
//! chunking ([`Decoder`]).
//!
//! Output is interleaved `f32` at full scale ±1.0 (a sample of the
//! standard's output range), one or two channels: mono streams decode to
//! one channel, stereo, joint stereo and dual channel streams to two.

use crate::error::{Result, invalid};
use crate::header::{FrameHeader, Layer};
use crate::layer3::{Layer3, Layer3Accounting};
use crate::layer12::{self, Subbands};
use crate::synth::Synth;
use crate::xing::{self, DECODER_DELAY, InfoHeader};

/// One decoded frame.
#[derive(Clone, Debug, PartialEq)]
pub struct Frame {
    /// Interleaved samples, `channels` per sample time, ±1.0 full scale.
    pub samples: Vec<f32>,
    /// 1 or 2.
    pub channels: u8,
    /// Hz.
    pub sample_rate: u32,
    /// The frame's header.
    pub header: FrameHeader,
    /// The frame's bit rate in bit/s (measured, for free format).
    pub bitrate: u32,
    /// Byte position of the frame in the stream (from the first byte given
    /// to the [`Decoder`]); 0 for [`FrameDecoder`].
    pub position: u64,
    /// Layer III bit accounting.
    pub layer3: Option<Layer3Accounting>,
    /// The frame could not be decoded and is silence standing in for it:
    /// a Layer III frame whose main data begins before the first byte the
    /// decoder received (decoding started mid-stream), or, when not in
    /// strict mode, a frame with corrupt data. It keeps the timeline; a
    /// caller that would rather drop it can.
    pub concealed: bool,
}

impl Frame {
    /// Samples per channel.
    pub fn len(&self) -> usize {
        self.samples.len() / usize::from(self.channels.max(1))
    }

    /// No samples (a frame trimmed away entirely by gapless trimming).
    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }
}

/// Decodes whole frames, one at a time, header included. Keeps the state
/// that crosses frames: the Layer III bit reservoir and IMDCT overlap, and
/// the synthesis filterbank.
#[derive(Clone, Default)]
pub struct FrameDecoder {
    synth: [Synth; 2],
    l3: Layer3,
    /// Check CRCs and refuse frames that fail.
    check_crc: bool,
    subbands: Subbands,
    /// Channels of the previous frame.
    channels: usize,
}

impl FrameDecoder {
    /// A decoder that ignores CRC words (as most players do).
    pub fn new() -> Self {
        Self::default()
    }

    /// Check each frame's CRC word, when it has one, and refuse frames whose
    /// CRC does not match ([`crate::Error::Invalid`]).
    pub fn check_crc(mut self, check: bool) -> Self {
        self.check_crc = check;
        self
    }

    /// Forget all state (after a seek).
    pub fn reset(&mut self) {
        for s in &mut self.synth {
            s.reset();
        }
        self.l3.reset();
    }

    /// Decode one frame. `frame` is exactly the frame's bytes (its length
    /// gives the bit rate of a free-format frame).
    pub fn decode_frame(&mut self, frame: &[u8]) -> Result<Frame> {
        let h = FrameHeader::parse(frame)?;
        let bitrate = if h.is_free_format() {
            let slots = frame.len() / h.slot_bytes() - usize::from(h.padding);
            let per = if h.layer == Layer::I { 12 } else { h.samples() / 8 };
            (slots as u64 * u64::from(h.sample_rate()) / per as u64) as u32
        } else {
            let want = h.frame_len().unwrap_or(0);
            if frame.len() < want {
                return Err(invalid(format!("frame is {} bytes, its header says {want}", frame.len())));
            }
            h.bitrate()
        };
        // A change between one and two channels starts the filterbanks
        // afresh (the reference decoder's output for l3_10201, whose mode
        // changes from mono to stereo and back, is that of a decoder started
        // anew at each change; the bit reservoir carries over).
        if self.channels != 0 && self.channels != h.channels() {
            for s in &mut self.synth {
                s.reset();
            }
            self.l3.reset_overlap();
        }
        self.channels = h.channels();
        let mut layer3 = None;
        match h.layer {
            Layer::I => layer12::decode_layer1(frame, &h, &mut self.subbands, self.check_crc)?,
            Layer::II => layer12::decode_layer2(frame, &h, bitrate, &mut self.subbands, self.check_crc)?,
            Layer::III => {
                let len = h.frame_len().unwrap_or(frame.len()).min(frame.len());
                layer3 = self.l3.decode(&frame[..len], &h, &mut self.subbands, self.check_crc)?;
            }
        }
        let nch = h.channels();
        let slots = self.subbands.len();
        let mut samples = vec![0.0f32; slots * 32 * nch];
        let mut pcm = [0.0f32; 32];
        for (t, slot) in self.subbands.iter().enumerate() {
            for ch in 0..nch {
                self.synth[ch].run(&slot[ch], &mut pcm);
                for (j, &v) in pcm.iter().enumerate() {
                    samples[(t * 32 + j) * nch + ch] = v;
                }
            }
        }
        Ok(Frame {
            samples,
            channels: nch as u8,
            sample_rate: h.sample_rate(),
            header: h,
            bitrate,
            position: 0,
            concealed: h.layer == Layer::III && layer3.is_none(),
            layer3,
        })
    }
}

/// Gapless figures for a stream: samples to drop at the start and the end
/// of the decoded output.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Gapless {
    /// Encoder delay as the tag states it.
    pub encoder_delay: u32,
    /// Encoder padding as the tag states it.
    pub padding: u32,
    /// Samples per channel to drop from the start of the decoded stream:
    /// `encoder_delay + 529`.
    pub skip_start: u32,
    /// Samples per channel of real audio, when the tag gives a frame count:
    /// `frames * samples_per_frame - encoder_delay - padding`.
    pub length: Option<u64>,
}

/// Options for [`Decoder`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DecoderOptions {
    /// Drop the encoder delay and padding a LAME tag (or VBRI header)
    /// declares, so the output is the encoder's input. Default on.
    pub trim_gapless: bool,
    /// Refuse frames whose CRC fails. Default off (decoded anyway).
    pub check_crc: bool,
    /// Treat every irregularity as an error instead of concealing it: bytes
    /// between frames, a CRC mismatch, a Layer III frame whose
    /// main_data_begin points before the data received or into the bits an
    /// earlier frame used, Huffman data that overruns part2_3_length.
    /// Default off. Meant for validating encoders.
    pub strict: bool,
}

impl Default for DecoderOptions {
    fn default() -> Self {
        Self { trim_gapless: true, check_crc: false, strict: false }
    }
}

/// Decodes an MPEG audio byte stream in any chunking: skips ID3v2 tags and
/// anything that is not a frame, confirms sync against the next frame's
/// header, handles free format, recognises (and does not play) a Xing /
/// Info / VBRI tag frame, and trims gapless delay and padding.
pub struct Decoder {
    opts: DecoderOptions,
    frames: FrameDecoder,
    buf: Vec<u8>,
    /// Bytes of `buf` already consumed: the data starts at `buf[head..]`
    /// (dropped from the front only now and then, so a large input is not
    /// moved down once per frame).
    head: usize,
    /// Stream position of `buf[0]`.
    base: u64,
    /// The header of the last frame accepted (sync is established).
    locked: Option<FrameHeader>,
    /// Free format: frame length without padding, in bytes.
    free_len: Option<usize>,
    /// The first frame has been examined for a tag.
    first_done: bool,
    info: Option<InfoHeader>,
    gapless: Option<Gapless>,
    /// Samples per channel output so far (before trimming).
    produced: u64,
    /// Unused main-data bytes after the last Layer III frame (strict mode).
    unused_main: Option<usize>,
    skipped_bytes: u64,
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}

/// The longest frame looked for: fixed-rate frames are under 2 KiB (Layer
/// II at 384 kb/s and 32 kHz is 1728 bytes); free-format frames can be
/// longer, at rates above the table's.
const MAX_FRAME: usize = 8192;

impl Decoder {
    /// A decoder with [`DecoderOptions::default`].
    pub fn new() -> Self {
        Self::with_options(DecoderOptions::default())
    }

    /// A decoder with the given options.
    pub fn with_options(opts: DecoderOptions) -> Self {
        Self {
            opts,
            frames: FrameDecoder::new().check_crc(opts.check_crc || opts.strict),
            buf: Vec::new(),
            head: 0,
            base: 0,
            locked: None,
            free_len: None,
            first_done: false,
            info: None,
            gapless: None,
            produced: 0,
            unused_main: None,
            skipped_bytes: 0,
        }
    }

    /// The Xing / Info / VBRI header of the stream, once the first frame has
    /// been seen.
    pub fn info(&self) -> Option<&InfoHeader> {
        self.info.as_ref()
    }

    /// The stream's gapless figures, once the first frame has been seen and
    /// if its tag gives them.
    pub fn gapless(&self) -> Option<Gapless> {
        self.gapless
    }

    /// Bytes skipped so far that were not frames (tags, garbage).
    pub fn skipped_bytes(&self) -> u64 {
        self.skipped_bytes
    }

    /// Feed bytes; returns the frames they complete. Frames are held back
    /// until the next frame's header confirms them (or [`Self::flush`]).
    pub fn decode(&mut self, bytes: &[u8]) -> Result<Vec<Frame>> {
        self.buf.extend_from_slice(bytes);
        self.drain(false)
    }

    /// No more bytes: decode what is buffered.
    pub fn flush(&mut self) -> Result<Vec<Frame>> {
        self.drain(true)
    }

    /// Decode a whole stream held in memory.
    pub fn decode_all(bytes: &[u8]) -> Result<Vec<Frame>> {
        let mut d = Decoder::new();
        let mut out = d.decode(bytes)?;
        out.extend(d.flush()?);
        Ok(out)
    }

    /// The buffered bytes not yet consumed.
    fn data(&self) -> &[u8] {
        &self.buf[self.head..]
    }

    /// Consume the first `n` buffered bytes.
    fn consume(&mut self, n: usize) {
        debug_assert!(n <= self.data().len());
        self.head += n;
        if self.head >= 1 << 16 && self.head * 2 >= self.buf.len() {
            self.buf.drain(..self.head);
            self.head = 0;
        }
    }

    fn clear_buf(&mut self) {
        self.buf.clear();
        self.head = 0;
    }

    fn skip(&mut self, n: usize) -> Result<()> {
        if self.opts.strict && n > 0 && self.locked.is_some() {
            return Err(invalid(format!("{n} bytes between frames at {}", self.base)));
        }
        self.consume(n);
        self.base += n as u64;
        self.skipped_bytes += n as u64;
        Ok(())
    }

    /// The length of the frame whose header `h` is at `buf[at..]`, or
    /// `None` if it cannot be told yet (free format before the next sync).
    fn frame_len(&mut self, h: &FrameHeader, at: usize, at_end: bool) -> Option<usize> {
        if let Some(len) = h.frame_len() {
            return Some(len);
        }
        if let Some(free) = self.free_len {
            return Some(free + usize::from(h.padding) * h.slot_bytes());
        }
        // Find the next header of the same free-format stream.
        let min = h.header_len() + 4;
        let mut i = at + min;
        while i + 4 <= self.data().len() && i < at + MAX_FRAME {
            if self.data()[i] == 0xFF
                && let Ok(n) = FrameHeader::parse(&self.data()[i..])
                && n.same_stream(h)
                && n.mode == h.mode
            {
                let len = i - at;
                let base = len - usize::from(h.padding) * h.slot_bytes();
                if base.is_multiple_of(h.slot_bytes()) {
                    self.free_len = Some(base);
                    return Some(len);
                }
            }
            i += 1;
        }
        if at_end && self.data().len() > at + min { Some(self.data().len() - at) } else { None }
    }

    fn drain(&mut self, at_end: bool) -> Result<Vec<Frame>> {
        let mut out = Vec::new();
        loop {
            // ID3v2 tag.
            if self.data().len() < 10 && self.data().starts_with(b"ID3") && !at_end {
                break;
            }
            if self.data().len() >= 10 && &self.data()[..3] == b"ID3" {
                let s = &self.data()[6..10];
                if s.iter().all(|&b| b < 0x80) {
                    let size = (usize::from(s[0]) << 21)
                        | (usize::from(s[1]) << 14)
                        | (usize::from(s[2]) << 7)
                        | usize::from(s[3]);
                    let footer = if self.data()[5] & 0x10 != 0 { 10 } else { 0 };
                    let total = 10 + size + footer;
                    if self.data().len() < total {
                        if at_end {
                            let n = self.data().len();
                            self.clear_buf();
                            self.base += n as u64;
                            self.skipped_bytes += n as u64;
                        }
                        break;
                    }
                    self.consume(total);
                    self.base += total as u64;
                    self.skipped_bytes += total as u64;
                    continue;
                }
            }
            if self.data().starts_with(b"ID3") {
                // Not a valid tag header: step over the three bytes.
                self.consume(3);
                self.base += 3;
                self.skipped_bytes += 3;
                continue;
            }
            // Sync search.
            let Some(start) = self.data().windows(2).position(|w| w[0] == 0xFF && w[1] & 0xE0 == 0xE0) else {
                let keep = usize::from(self.data().last() == Some(&0xFF));
                let n = self.data().len() - keep;
                if n > 0 && !(self.data().len() >= 3 && &self.data()[..3] == b"ID3") {
                    // Trailing ID3v1 / APE tags end up here too.
                    if self.opts.strict && self.locked.is_some() && !at_end {
                        return Err(invalid(format!("{n} bytes between frames at {}", self.base)));
                    }
                    self.consume(n);
                    self.base += n as u64;
                    self.skipped_bytes += n as u64;
                }
                break;
            };
            if start > 0 {
                // Allow an ID3 tag to be recognised.
                let id3 = self.data()[..start].windows(3).position(|w| w == b"ID3");
                let cut = id3.filter(|&p| p > 0).unwrap_or(start);
                if cut > 0 {
                    if self.opts.strict && self.locked.is_some() && !self.is_trailing_tag(at_end) {
                        return Err(invalid(format!("{cut} bytes between frames at {}", self.base)));
                    }
                    self.consume(cut);
                    self.base += cut as u64;
                    self.skipped_bytes += cut as u64;
                }
                continue;
            }
            if self.data().len() < 4 {
                if at_end {
                    let n = self.data().len();
                    self.clear_buf();
                    self.base += n as u64;
                    self.skipped_bytes += n as u64;
                }
                break;
            }
            let h = match FrameHeader::parse(self.data()) {
                Ok(h) if self.locked.is_none_or(|l| l.same_stream(&h)) => h,
                _ => {
                    self.skip(1)?;
                    continue;
                }
            };
            let Some(len) = self.frame_len(&h, 0, at_end) else {
                if self.data().len() > MAX_FRAME {
                    self.skip(1)?;
                    continue;
                }
                break;
            };
            if len < h.header_len() + 2 {
                self.skip(1)?;
                continue;
            }
            if self.data().len() < len {
                if at_end {
                    // A truncated last frame: decode what there is, padded.
                    if self.locked.is_some() && self.data().len() > h.header_len() + 8 {
                        let mut f = self.data().to_vec();
                        f.resize(len, 0);
                        let n = self.data().len();
                        self.clear_buf();
                        self.base += n as u64;
                        if let Ok(frame) = self.frames.decode_frame(&f) {
                            self.emit(frame, &mut out)?;
                        }
                    } else {
                        let n = self.data().len();
                        self.clear_buf();
                        self.base += n as u64;
                        self.skipped_bytes += n as u64;
                    }
                }
                break;
            }
            // Confirm against the next header unless already locked.
            if self.locked.is_none() {
                if self.data().len() < len + 4 {
                    if !at_end {
                        break;
                    }
                } else {
                    match FrameHeader::parse(&self.data()[len..]) {
                        Ok(n) if n.same_stream(&h) => {}
                        _ => {
                            self.skip(1)?;
                            continue;
                        }
                    }
                }
            }
            let frame_bytes: Vec<u8> = self.data()[..len].to_vec();
            let position = self.base;
            self.consume(len);
            self.base += len as u64;
            self.locked = Some(h);
            if !self.first_done {
                self.first_done = true;
                if let Some(info) = xing::parse(&frame_bytes, &h) {
                    self.set_info(info, &h);
                    // The tag frame is not audio; for Layer III keep its
                    // (empty) main data out of the reservoir by not decoding.
                    continue;
                }
            }
            match self.frames.decode_frame(&frame_bytes) {
                Ok(mut frame) => {
                    frame.position = position;
                    if self.opts.strict {
                        self.check_accounting(&frame)?;
                    }
                    self.emit(frame, &mut out)?;
                }
                Err(e) => {
                    if self.opts.strict {
                        return Err(e);
                    }
                    // Conceal: a frame of silence keeps the timeline.
                    let nch = h.channels();
                    let frame = Frame {
                        samples: vec![0.0; h.samples() * nch],
                        channels: nch as u8,
                        sample_rate: h.sample_rate(),
                        header: h,
                        bitrate: h.bitrate(),
                        position,
                        layer3: None,
                        concealed: true,
                    };
                    self.emit(frame, &mut out)?;
                }
            }
        }
        Ok(out)
    }

    fn is_trailing_tag(&self, _at_end: bool) -> bool {
        self.data().starts_with(b"TAG") || self.data().starts_with(b"APETAGEX") || self.data().starts_with(b"LYRICS")
    }

    fn set_info(&mut self, info: InfoHeader, h: &FrameHeader) {
        if let Some((delay, padding)) = info.delay_padding() {
            let spf = h.samples() as u64;
            let length =
                info.frames().map(|f| (u64::from(f) * spf).saturating_sub(u64::from(delay) + u64::from(padding)));
            self.gapless = Some(Gapless { encoder_delay: delay, padding, skip_start: delay + DECODER_DELAY, length });
        }
        self.info = Some(info);
    }

    fn check_accounting(&mut self, frame: &Frame) -> Result<()> {
        let Some(a) = frame.layer3 else {
            if frame.header.layer == Layer::III {
                return Err(invalid("main_data_begin points before the start of the stream"));
            }
            return Ok(());
        };
        if !a.exact {
            return Err(invalid(format!("Huffman data does not end at part2_3_length (frame at {})", frame.position)));
        }
        if let Some(unused) = self.unused_main
            && a.main_data_begin > unused
        {
            return Err(invalid(format!(
                "main_data_begin {} reaches into the previous frame's data ({} bytes unused) at {}",
                a.main_data_begin, unused, frame.position
            )));
        }
        if self.unused_main.is_none() && a.main_data_begin > 0 {
            return Err(invalid("the first frame's main_data_begin is not 0"));
        }
        let avail = (a.main_data_begin + a.main_data_bytes) * 8;
        if a.part2_3_bits > avail {
            return Err(invalid("frame uses more main data than it has"));
        }
        // Unused bytes at the end of this frame's main data, capped by the
        // reservoir limit the next frame can reach.
        self.unused_main = Some((avail - a.part2_3_bits) / 8);
        Ok(())
    }

    /// Apply gapless trimming and push the frame.
    fn emit(&mut self, mut frame: Frame, out: &mut Vec<Frame>) -> Result<()> {
        let n = frame.len() as u64;
        let start = self.produced;
        self.produced += n;
        if self.opts.trim_gapless
            && let Some(g) = self.gapless
        {
            let lo = u64::from(g.skip_start);
            let hi = g.length.map_or(u64::MAX, |len| lo + len);
            let a = lo.clamp(start, start + n) - start;
            let b = hi.clamp(start, start + n) - start;
            let ch = usize::from(frame.channels);
            if a > 0 || b < n {
                let b = b.max(a);
                frame.samples = frame.samples[a as usize * ch..b as usize * ch].to_vec();
            }
        }
        out.push(frame);
        Ok(())
    }
}
