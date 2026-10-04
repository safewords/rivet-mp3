//! A Layer III (MP3) encoder.
//!
//! - **Rates**: every Layer III sampling frequency — 32, 44.1, 48 kHz
//!   (MPEG-1), 16, 22.05, 24 kHz (MPEG-2 LSF) and 8, 11.025, 12 kHz
//!   (MPEG-2.5); [`coding_rate`] names the nearest one for other input.
//! - **Bit rate**: constant ([`BitrateMode::Cbr`], any bit rate of the
//!   version's table) or variable ([`BitrateMode::Vbr`], quality 0 best –
//!   9 smallest), with the bit reservoir.
//! - **Channels**: one or two; two are coded as joint stereo (mid/side,
//!   chosen frame by frame) unless [`EncoderConfig::joint_stereo`] is off.
//! - **Tools**: the polyphase analysis filterbank and MDCT of 11172-3
//!   Annex C, long, start, short and stop blocks switched by a transient
//!   detector, a psychoacoustic model, quantisation against the
//!   masking threshold with one noise-to-mask offset per frame chosen by
//!   bisection against the bit budget, Huffman table and region selection.
//! - **Gapless**: [`Encoder::tag_frame`] is a Xing (VBR) or Info (CBR)
//!   frame with a LAME-style extension carrying the encoder delay
//!   ([`Encoder::delay`]) and padding, which a gapless decoder (this
//!   crate's among them) uses to give back exactly the input.
//!
//! Input is interleaved `f32` at ±1.0 full scale.

pub(crate) mod analysis;
mod huffman;
pub(crate) mod psy;
mod quantize;
#[cfg(test)]
mod tests;

use std::collections::VecDeque;
use std::sync::{Arc, Mutex, mpsc};

use crate::bits::BitWriter;
use crate::crc::frame_crc_bits;
use crate::error::{Result, config};
use crate::header::{BITRATES_LSF, BITRATES_MPEG1, FrameHeader, Layer, Mode, SAMPLE_RATES, Version};
use crate::layer3::imdct::mdct;
use crate::layer3::sideinfo::{GranuleInfo, SideInfo};
use crate::tables::layer3::{SFB_LONG, SFB_SHORT, rate_index};
use crate::xing::{self, DECODER_DELAY};
use analysis::Analysis;
use psy::{Mask, Psy, Transient};
use quantize::{Quantised, Spectrum};

/// How the bit rate is set.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BitrateMode {
    /// Constant bit rate, bit/s: one of the version's Layer III rates
    /// (MPEG-1: 32 000 – 320 000; MPEG-2 / 2.5: 8 000 – 160 000).
    Cbr(u32),
    /// Variable bit rate by quality, 0 (best) to 9 (smallest).
    Vbr(u8),
}

/// What to encode.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EncoderConfig {
    /// One of the nine Layer III sampling frequencies.
    pub sample_rate: u32,
    /// 1 or 2.
    pub channels: u8,
    /// Constant or variable bit rate.
    pub bitrate: BitrateMode,
    /// Code two channels as joint stereo (mid/side where it pays, frame by
    /// frame). Off: plain stereo.
    pub joint_stereo: bool,
    /// Protect each frame with a CRC word.
    pub crc: bool,
    /// Low-pass cutoff in Hz; `None` picks one from the bit rate.
    pub lowpass: Option<u32>,
}

impl Default for EncoderConfig {
    fn default() -> Self {
        Self {
            sample_rate: 44_100,
            channels: 2,
            bitrate: BitrateMode::Cbr(128_000),
            joint_stereo: true,
            crc: false,
            lowpass: None,
        }
    }
}

/// Encoder delay written in the tag: samples the decoded stream holds,
/// beyond the decoder's own 529, before the first input sample. The
/// analysis filterbank delays by 481 samples and the hybrid filterbank by a
/// granule (576); 481 + 576 = 529 + 528.
pub const ENCODER_DELAY: u32 = 528;

/// The Layer III sampling frequency to resample other input to: the input
/// rate if it is one, otherwise the nearest.
pub fn coding_rate(input: u32) -> u32 {
    let all = SAMPLE_RATES.iter().flatten().copied();
    if all.clone().any(|r| r == input) {
        return input;
    }
    all.min_by_key(|&r| (i64::from(r) - i64::from(input)).abs()).unwrap_or(44_100)
}

/// The encoder string of the tag (nine bytes).
const ENCODER_NAME: &[u8; 9] = b"rivetmp3 ";

/// Default low-pass cutoff (Hz) for a bit rate per channel (kbit/s).
fn default_lowpass(kbps_per_channel: f64) -> f64 {
    const POINTS: [(f64, f64); 12] = [
        (8.0, 3000.0),
        (16.0, 5500.0),
        (24.0, 7500.0),
        (32.0, 10_000.0),
        (40.0, 12_000.0),
        (48.0, 14_000.0),
        (56.0, 15_000.0),
        (64.0, 16_000.0),
        (80.0, 17_500.0),
        (96.0, 18_500.0),
        (112.0, 19_500.0),
        (128.0, 20_000.0),
    ];
    if kbps_per_channel <= POINTS[0].0 {
        return POINTS[0].1;
    }
    for w in POINTS.windows(2) {
        if kbps_per_channel <= w[1].0 {
            let t = (kbps_per_channel - w[0].0) / (w[1].0 - w[0].0);
            return w[0].1 + t * (w[1].1 - w[0].1);
        }
    }
    20_500.0
}

/// VBR quality: noise-to-mask offset (dB) and low-pass cutoff.
fn vbr_setting(quality: u8) -> (f64, f64) {
    let q = f64::from(quality.min(9));
    (-4.0 + 1.6 * q, 20_000.0 - 700.0 * q)
}

/// A frame's quantised granules, [granule][channel].
type FrameQuant = Vec<Vec<Quantised>>;

/// One channel's analysis state.
struct Channel {
    analysis: Analysis,
    /// The previous granule's 18 slots of 32 subband samples.
    prev_slots: [[f64; 32]; 18],
    transient: Transient,
    /// Pre-echo control: the last long granule's thresholds.
    prev_thr: Option<[f64; 22]>,
}

/// One granule, analysed: spectra and masks per channel.
struct Granule {
    block_type: u8,
    xr: Vec<[f32; 576]>,
    mask: Vec<Mask>,
}

/// A frame waiting for later frames' main data to fill its slot.
struct Pending {
    bytes: Vec<u8>,
    /// Offset of the main-data slot in `bytes`.
    slot_start: usize,
    filled: usize,
}

/// The Layer III encoder.
pub struct Encoder {
    cfg: EncoderConfig,
    version: Version,
    sample_rate_index: u8,
    lsf: bool,
    ngr: usize,
    spf: usize,
    nch: usize,
    long_edges: [u16; 23],
    short_edges: [u16; 14],
    psy: Psy,
    /// Lines (long) at and above which the spectrum is cut.
    cutoff_long: usize,
    cutoff_short: usize,
    vbr_offset_db: f64,
    channels: Vec<Channel>,
    /// Input not yet analysed, per channel; `input_base` is the absolute
    /// index of its first sample.
    input: Vec<Vec<f32>>,
    input_base: u64,
    /// Samples per channel given so far.
    samples_in: u64,
    /// Next granule to analyse.
    next_granule: u64,
    /// Attack flags computed ahead (granule index -> attack).
    attacks: VecDeque<(u64, bool)>,
    prev_block: u8,
    frame_granules: Vec<Granule>,
    /// Reservoir: unused main-data bytes at the end of the pending frames.
    free_tail: usize,
    pending: VecDeque<Pending>,
    /// CBR padding accumulator (in 1/sample_rate of a slot).
    pad_acc: u64,
    frames_out: u64,
    bytes_out: u64,
    frame_offsets: Vec<u64>,
    music_crc: u16,
    flushed: bool,
    padding: u32,
    bitrate_sum: u64,
    /// Threads for a batch of frames; 0 is the machine's count.
    threads: usize,
    /// The quantiser threads, once started.
    pool: Option<QuantPool>,
}

/// A frame's work that does not depend on the frames before it: the
/// stereo decision, the spectra and the first quantisation the rate
/// control starts from. Frames of one batch prepare in parallel.
struct Prepared {
    grs: Vec<Granule>,
    ms: bool,
    work: Arc<FrameWork>,
    /// The quantisation at the mask (CBR) or at the quality's offset (VBR).
    first: FrameQuant,
}

/// What quantising a frame needs, shared with the quantiser threads: the
/// spectra to code (granule-major, then channel), their masks and lines.
struct FrameWork {
    items: Vec<WorkItem>,
    /// Channels coded per granule.
    nch: usize,
    long_edges: [u16; 23],
    short_edges: [u16; 14],
    lsf: bool,
}

struct WorkItem {
    xr: [f32; 576],
    mask: Mask,
    lines: Box<quantize::Lines>,
    block_type: u8,
}

impl FrameWork {
    /// Item `i` quantised with every band's allowance times `scale`.
    fn quantise(&self, i: usize, scale: f64) -> Quantised {
        let it = &self.items[i];
        let sp = Spectrum {
            lines: &it.lines,
            block_type: it.block_type,
            mask: &it.mask,
            long_edges: &self.long_edges,
            short_edges: &self.short_edges,
            lsf: self.lsf,
        };
        let mut q = quantize::quantise(&sp, scale);
        quantize::apply_signs(&mut q, &it.xr);
        q
    }
}

/// Threads that quantise a frame's granules and channels side by side
/// during the rate control's search (each trial quantises them all, and
/// the trials themselves follow one from the other). Started on first
/// use; they end when the encoder is dropped.
struct QuantPool {
    jobs: mpsc::Sender<(Arc<FrameWork>, usize, f64)>,
    /// Behind a lock only so the encoder stays `Sync`; one thread reads it.
    done: Mutex<mpsc::Receiver<(usize, Quantised)>>,
}

impl QuantPool {
    fn new(workers: usize) -> QuantPool {
        let (jobs, job_rx) = mpsc::channel::<(Arc<FrameWork>, usize, f64)>();
        let (done_tx, done) = mpsc::channel();
        let job_rx = Arc::new(Mutex::new(job_rx));
        for _ in 0..workers {
            let job_rx = Arc::clone(&job_rx);
            let done_tx = done_tx.clone();
            std::thread::spawn(move || {
                loop {
                    // The lock is held only to take a job, not to run it.
                    let job = job_rx.lock().map_or(Err(mpsc::RecvError), |rx| rx.recv());
                    let Ok((work, i, scale)) = job else { return };
                    if done_tx.send((i, work.quantise(i, scale))).is_err() {
                        return;
                    }
                }
            });
        }
        QuantPool { jobs, done: Mutex::new(done) }
    }
}

fn crc16_arc_update(mut crc: u16, data: &[u8]) -> u16 {
    for &b in data {
        crc ^= u16::from(b);
        for _ in 0..8 {
            crc = if crc & 1 == 1 { (crc >> 1) ^ 0xA001 } else { crc >> 1 };
        }
    }
    crc
}

impl Encoder {
    /// An encoder for `cfg`; refuses a sampling frequency Layer III does not
    /// code, a channel count other than 1 or 2, or a bit rate off the
    /// version's table.
    pub fn new(cfg: EncoderConfig) -> Result<Encoder> {
        let (vi, si) = SAMPLE_RATES
            .iter()
            .enumerate()
            .find_map(|(v, row)| row.iter().position(|&r| r == cfg.sample_rate).map(|s| (v, s)))
            .ok_or_else(|| {
                config(format!(
                    "{} Hz is not a Layer III sampling frequency (8, 11.025, 12, 16, 22.05, 24, 32, 44.1, 48 kHz; \
                     see encode::coding_rate)",
                    cfg.sample_rate
                ))
            })?;
        if !(1..=2).contains(&cfg.channels) {
            return Err(config(format!("MP3 carries one or two channels, not {}", cfg.channels)));
        }
        let version = [Version::Mpeg1, Version::Mpeg2, Version::Mpeg25][vi];
        let lsf = vi != 0;
        let table = if lsf { &BITRATES_LSF[2] } else { &BITRATES_MPEG1[2] };
        let vbr_offset_db;
        let lowpass;
        let nch = usize::from(cfg.channels);
        match cfg.bitrate {
            BitrateMode::Cbr(b) => {
                if b % 1000 != 0 || !table[1..].contains(&(b / 1000)) {
                    return Err(config(format!(
                        "{b} bit/s is not a {} Layer III bit rate ({} kbit/s)",
                        if lsf { "MPEG-2/2.5" } else { "MPEG-1" },
                        table[1..].iter().map(u32::to_string).collect::<Vec<_>>().join(", ")
                    )));
                }
                vbr_offset_db = 0.0;
                let per = f64::from(b) / 1000.0 / nch as f64 * if nch == 2 && cfg.joint_stereo { 1.15 } else { 1.0 };
                lowpass = default_lowpass(per);
            }
            BitrateMode::Vbr(q) => {
                if q > 9 {
                    return Err(config(format!("VBR quality {q} is outside 0..=9")));
                }
                let (off, lp) = vbr_setting(q);
                vbr_offset_db = off;
                lowpass = lp;
            }
        }
        let lowpass = cfg.lowpass.map_or(lowpass, f64::from).min(f64::from(cfg.sample_rate) * 0.5);
        let rate = rate_index(vi, si);
        let long_edges = SFB_LONG[rate];
        let short_edges = SFB_SHORT[rate];
        let nyquist = f64::from(cfg.sample_rate) / 2.0;
        let cutoff_long = ((lowpass / nyquist) * 576.0).round().clamp(1.0, 576.0) as usize;
        let cutoff_short = ((lowpass / nyquist) * 192.0).round().clamp(1.0, 192.0) as usize;
        Ok(Encoder {
            cfg,
            version,
            sample_rate_index: si as u8,
            lsf,
            ngr: if lsf { 1 } else { 2 },
            spf: if lsf { 576 } else { 1152 },
            nch,
            long_edges,
            short_edges,
            psy: Psy::new(cfg.sample_rate, &long_edges, &short_edges),
            cutoff_long,
            cutoff_short,
            vbr_offset_db,
            channels: (0..nch)
                .map(|_| Channel {
                    analysis: Analysis::default(),
                    prev_slots: [[0.0; 32]; 18],
                    transient: Transient::default(),
                    prev_thr: None,
                })
                .collect(),
            input: vec![Vec::new(); nch],
            input_base: 0,
            samples_in: 0,
            next_granule: 0,
            attacks: VecDeque::new(),
            prev_block: 0,
            frame_granules: Vec::new(),
            free_tail: 0,
            pending: VecDeque::new(),
            pad_acc: 0,
            frames_out: 0,
            bytes_out: 0,
            frame_offsets: Vec::new(),
            music_crc: 0,
            flushed: false,
            padding: 0,
            bitrate_sum: 0,
            threads: 0,
            pool: None,
        })
    }

    /// How many threads prepare a batch of frames (the frames one
    /// [`encode`](Self::encode) call completes): 0, the default, is one per
    /// CPU; 1 does everything on the caller's thread. The stream is the
    /// same byte for byte whatever the count.
    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    /// Samples per channel per frame: 1152 (MPEG-1) or 576.
    pub fn frame_samples(&self) -> usize {
        self.spf
    }

    /// The sampling frequency coded.
    pub fn sample_rate(&self) -> u32 {
        self.cfg.sample_rate
    }

    /// Encoder delay (samples per channel) as the tag states it: the
    /// decoded stream starts `delay() + 529` samples before the input.
    pub fn delay(&self) -> u32 {
        ENCODER_DELAY
    }

    /// Samples added after the input to fill the last frame (known after
    /// [`Self::flush`]).
    pub fn padding(&self) -> u32 {
        self.padding
    }

    /// Frames emitted so far (the tag frame not included).
    pub fn frames(&self) -> u64 {
        self.frames_out
    }

    /// Average bit rate of the frames emitted so far, bit/s.
    pub fn average_bitrate(&self) -> f64 {
        if self.frames_out == 0 {
            return 0.0;
        }
        self.bitrate_sum as f64 / self.frames_out as f64
    }

    /// Encode interleaved samples; returns the frames completed (one
    /// packet per frame). Frames can be held back while later frames may
    /// still put main data into them (the bit reservoir).
    pub fn encode(&mut self, pcm: &[f32]) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        if self.flushed {
            return out;
        }
        for (i, &v) in pcm.iter().enumerate() {
            self.input[i % self.nch].push(v);
        }
        self.samples_in += (pcm.len() / self.nch) as u64;
        self.run(false, &mut out);
        out
    }

    /// No more input: encode the rest, padded with silence, and emit every
    /// remaining frame.
    pub fn flush(&mut self) -> Vec<Vec<u8>> {
        let mut out = Vec::new();
        if self.flushed {
            return out;
        }
        self.run(true, &mut out);
        // Frames still waiting for main data: their slots stay unused.
        while let Some(p) = self.pending.pop_front() {
            self.emit(p.bytes, &mut out);
        }
        self.flushed = true;
        out
    }

    /// The input sample with absolute index `n` of channel `ch` (zero
    /// before the start and past the end).
    fn sample(&self, ch: usize, n: i64) -> f32 {
        if n < self.input_base as i64 {
            return 0.0;
        }
        let i = (n - self.input_base as i64) as usize;
        self.input[ch].get(i).copied().unwrap_or(0.0)
    }

    fn run(&mut self, at_end: bool, out: &mut Vec<Vec<u8>>) {
        let total_frames = if at_end {
            let need = self.samples_in + u64::from(ENCODER_DELAY + DECODER_DELAY);
            need.div_ceil(self.spf as u64)
        } else {
            u64::MAX
        };
        let mut batch: Vec<Vec<Granule>> = Vec::new();
        loop {
            let g = self.next_granule;
            if at_end && g >= total_frames * self.ngr as u64 {
                break;
            }
            // Input needed: the granule's own 576 samples and the look-ahead
            // of the next granule's attack region.
            let need = 576 * (g + 1) + 351;
            if !at_end && self.samples_in < need {
                break;
            }
            self.analyse_granule(g);
            self.next_granule += 1;
            if self.frame_granules.len() == self.ngr {
                batch.push(std::mem::take(&mut self.frame_granules));
            }
            // Drop input no longer needed: the next granule's analysis starts
            // at 576 (g + 1), and every later attack region after it.
            let keep_from = 576 * (g + 1);
            if keep_from > self.input_base {
                let cut = ((keep_from - self.input_base) as usize).min(self.input[0].len());
                for ch in &mut self.input {
                    ch.drain(..cut.min(ch.len()));
                }
                self.input_base += cut as u64;
            }
        }
        // The frames' independent work in parallel, then the rate control
        // and the reservoir in order.
        let threads =
            if self.threads == 0 { std::thread::available_parallelism().map_or(1, usize::from) } else { self.threads };
        let prepared = parallel_map(batch, threads, |grs| self.prepare(grs));
        if threads > 1 && self.pool.is_none() && !prepared.is_empty() {
            // One thread per granule-channel beyond the caller's.
            self.pool = Some(QuantPool::new((self.ngr * self.nch - 1).min(threads - 1)));
        }
        for p in prepared {
            self.encode_frame(p, out);
        }
        if at_end {
            self.padding = (total_frames * self.spf as u64 - u64::from(ENCODER_DELAY) - self.samples_in) as u32;
        }
    }

    /// Whether granule `g` has an attack in any channel.
    fn attack(&mut self, g: u64) -> bool {
        if let Some(&(_, a)) = self.attacks.iter().find(|(i, _)| *i == g) {
            return a;
        }
        let start = 576 * g as i64 - 225;
        let mut any = false;
        for ch in 0..self.nch {
            let x: Vec<f32> = (0..576).map(|i| self.sample(ch, start + i)).collect();
            any |= self.channels[ch].transient.granule(&x);
        }
        self.attacks.push_back((g, any));
        while self.attacks.len() > 4 {
            self.attacks.pop_front();
        }
        any
    }

    fn analyse_granule(&mut self, g: u64) {
        let now = self.attack(g);
        let next = self.attack(g + 1);
        let bt = match (self.prev_block, now, next) {
            (_, true, _) if matches!(self.prev_block, 1 | 2) => 2,
            (1, _, _) => 2, // a start block must be followed by a short one
            (2, false, true) => 2,
            (2, false, false) => 3,
            (_, true, _) => 1, // attack without a start block before it: start now, short next
            (_, false, true) => 1,
            _ => 0,
        };
        self.prev_block = bt;
        let mut xr_all = Vec::with_capacity(self.nch);
        let mut masks = Vec::with_capacity(self.nch);
        for ch in 0..self.nch {
            let base = 576 * g as i64;
            let mut slots = [[0.0f64; 32]; 18];
            for (t, slot) in slots.iter_mut().enumerate() {
                let input: Vec<f32> = (0..32).map(|i| self.sample(ch, base + 32 * t as i64 + i)).collect();
                self.channels[ch].analysis.run(&input, slot);
            }
            let mut xr = [0.0f32; 576];
            let prev = self.channels[ch].prev_slots;
            let mut coeffs = [[0.0f32; 18]; 32];
            for (sb, c) in coeffs.iter_mut().enumerate() {
                let mut x = [0.0f64; 36];
                for i in 0..36 {
                    let v = if i < 18 { prev[i][sb] } else { slots[i - 18][sb] };
                    x[i] = if sb % 2 == 1 && i % 2 == 1 { -v } else { v };
                }
                mdct(&x, bt, c);
            }
            self.channels[ch].prev_slots = slots;
            if bt == 2 {
                // Bitstream order: band, window, line.
                let mut i = 0;
                for b in 0..13 {
                    let lo = usize::from(self.short_edges[b]);
                    let hi = usize::from(self.short_edges[b + 1]);
                    for w in 0..3 {
                        for f in lo..hi {
                            xr[i] = if f < self.cutoff_short { coeffs[f / 6][w * 6 + f % 6] } else { 0.0 };
                            i += 1;
                        }
                    }
                }
                masks.push(self.psy.short(&xr));
                self.channels[ch].prev_thr = None;
            } else {
                for (sb, c) in coeffs.iter().enumerate() {
                    xr[sb * 18..sb * 18 + 18].copy_from_slice(c);
                }
                // Forward antialias butterflies (the inverse of the
                // decoder's).
                const CI: [f64; 8] = [-0.6, -0.535, -0.33, -0.185, -0.095, -0.041, -0.0142, -0.0037];
                for sb in 1..32 {
                    for (i, &c) in CI.iter().enumerate() {
                        let cs = 1.0 / (1.0 + c * c).sqrt();
                        let ca = c / (1.0 + c * c).sqrt();
                        let bu = f64::from(xr[18 * sb - 1 - i]);
                        let bd = f64::from(xr[18 * sb + i]);
                        xr[18 * sb - 1 - i] = (bu * cs + bd * ca) as f32;
                        xr[18 * sb + i] = (bd * cs - bu * ca) as f32;
                    }
                }
                for v in xr[self.cutoff_long..].iter_mut() {
                    *v = 0.0;
                }
                let mask = self.psy.long(&xr, self.channels[ch].prev_thr.as_ref());
                if let Mask::Long { thr, .. } = &mask {
                    self.channels[ch].prev_thr = Some(*thr);
                }
                masks.push(mask);
            }
            xr_all.push(xr);
        }
        self.frame_granules.push(Granule { block_type: bt, xr: xr_all, mask: masks });
    }

    fn header(&self, bitrate_index: u8, padding: bool, ms: bool) -> FrameHeader {
        let mode = match (self.nch, self.cfg.joint_stereo) {
            (1, _) => Mode::Mono,
            (_, true) => Mode::JointStereo,
            _ => Mode::Stereo,
        };
        FrameHeader {
            version: self.version,
            layer: Layer::III,
            crc: self.cfg.crc,
            bitrate_index,
            sample_rate_index: self.sample_rate_index,
            padding,
            private: false,
            mode,
            mode_extension: if ms { 2 } else { 0 },
            copyright: false,
            original: true,
            emphasis: 0,
        }
    }

    fn bitrate_table(&self) -> &'static [u32; 15] {
        if self.lsf { &BITRATES_LSF[2] } else { &BITRATES_MPEG1[2] }
    }

    /// Main-data bytes of a frame at `index` with `padding`.
    fn slot_bytes(&self, index: u8, padding: bool) -> usize {
        let h = self.header(index, padding, false);
        h.frame_len().unwrap_or(0) - h.header_len() - h.side_info_len()
    }

    /// Largest main_data_begin allowed: the 9- (8-) bit field, and the
    /// decoder's 7680-bit input buffer less the frame.
    fn reservoir_max(&self, frame_bytes: usize) -> usize {
        let field = if self.lsf { 255 } else { 511 };
        field.min((7680 / 8usize).saturating_sub(frame_bytes))
    }

    fn prepare(&self, grs: Vec<Granule>) -> Prepared {
        let nch = self.nch;
        // Mid/side decision for the frame.
        let ms = nch == 2 && self.cfg.joint_stereo && self.choose_ms(&grs);
        let mut spectra: Vec<Vec<([f32; 576], Mask)>> = Vec::new();
        for g in &grs {
            let mut v = Vec::new();
            if ms {
                let s2 = std::f32::consts::FRAC_1_SQRT_2;
                let m: [f32; 576] = std::array::from_fn(|i| (g.xr[0][i] + g.xr[1][i]) * s2);
                let s: [f32; 576] = std::array::from_fn(|i| (g.xr[0][i] - g.xr[1][i]) * s2);
                let mask = g.mask[0].min(&g.mask[1]);
                let mm = mask.with_energy_of(&m, &self.long_edges, &self.short_edges);
                let ms_ = mask.with_energy_of(&s, &self.long_edges, &self.short_edges);
                v.push((m, mm));
                v.push((s, ms_));
            } else {
                for ch in 0..nch {
                    v.push((g.xr[ch], g.mask[ch].clone()));
                }
            }
            spectra.push(v);
        }
        let items = spectra
            .into_iter()
            .zip(&grs)
            .flat_map(|(v, g)| {
                v.into_iter().map(|(xr, mask)| WorkItem {
                    lines: quantize::Lines::new(&xr),
                    xr,
                    mask,
                    block_type: g.block_type,
                })
            })
            .collect();
        let work = Arc::new(FrameWork {
            items,
            nch,
            long_edges: self.long_edges,
            short_edges: self.short_edges,
            lsf: self.lsf,
        });
        let offset = match self.cfg.bitrate {
            BitrateMode::Cbr(_) => 0.0,
            BitrateMode::Vbr(_) => self.vbr_offset_db,
        };
        let first = quantise_serial(&work, offset);
        Prepared { grs, ms, work, first }
    }

    /// Quantise every granule and channel of a frame with every band's
    /// allowed noise raised by `offset_db`: on the quantiser threads too
    /// when there are any, the result the same either way.
    fn quantise_frame(&self, work: &Arc<FrameWork>, offset_db: f64) -> FrameQuant {
        let Some(pool) = &self.pool else { return quantise_serial(work, offset_db) };
        let scale = 10f64.powf(offset_db / 10.0);
        let n = work.items.len();
        // Items 1.. go to the threads; this thread takes item 0.
        for i in 1..n {
            pool.jobs.send((Arc::clone(work), i, scale)).expect("quantiser threads run while the encoder lives");
        }
        let mut out: Vec<Option<Quantised>> = vec![None; n];
        out[0] = Some(work.quantise(0, scale));
        let done = pool.done.lock().expect("one reader");
        for _ in 1..n {
            let (i, q) = done.recv().expect("quantiser threads run while the encoder lives");
            out[i] = Some(q);
        }
        let mut items = out.into_iter().map(|q| q.expect("every item quantised"));
        (0..n / work.nch).map(|_| items.by_ref().take(work.nch).collect()).collect()
    }

    fn encode_frame(&mut self, p: Prepared, out: &mut Vec<Vec<u8>>) {
        let Prepared { grs, ms, work, first } = p;
        let quantise_all = |enc: &Encoder, offset_db: f64| -> FrameQuant { enc.quantise_frame(&work, offset_db) };
        let total = |q: &Vec<Vec<Quantised>>| -> (u32, bool) {
            let mut sum = 0;
            let mut ok = true;
            for gq in q {
                for c in gq {
                    let b = c.bits();
                    ok &= b <= 4095;
                    sum += b;
                }
            }
            (sum, ok)
        };
        // Choose the bit rate (VBR) and the budget.
        let (index, padding) = match self.cfg.bitrate {
            BitrateMode::Cbr(b) => {
                let index = self.bitrate_table().iter().position(|&r| r == b / 1000).unwrap_or(1) as u8;
                // Padding keeps the average frame length exact.
                let num = (self.spf as u64 / 8) * u64::from(b);
                let fs = u64::from(self.cfg.sample_rate);
                self.pad_acc += num % fs;
                let padding = self.pad_acc >= fs;
                if padding {
                    self.pad_acc -= fs;
                }
                (index, padding)
            }
            BitrateMode::Vbr(_) => (0, false),
        };
        let (q, index) = match self.cfg.bitrate {
            BitrateMode::Cbr(_) => {
                let slot = self.slot_bytes(index, padding);
                let frame_bytes = self.header(index, padding, ms).frame_len().unwrap_or(0);
                let resv_max = self.reservoir_max(frame_bytes);
                let avail = 8 * (self.free_tail + slot) as u32;
                let mean = 8 * slot as u32;
                // Demand: the bits that just meet the mask.
                let at_mask = first;
                let demand = total(&at_mask).0;
                let room = 8 * resv_max.saturating_sub(self.free_tail) as u32; // reservoir space left
                let target = if demand > mean {
                    demand.min(mean + (8 * self.free_tail as u32) * 6 / 10)
                } else {
                    // Save into the reservoir while it has room; spend the
                    // surplus on quality once it is full.
                    demand.max(mean.saturating_sub(room))
                }
                .min(avail);
                let q = self.fit(&quantise_all, &total, at_mask, target);
                (q, index)
            }
            BitrateMode::Vbr(_) => {
                let q0 = first;
                let need = total(&q0).0;
                // Smallest bit rate whose frame, with half the reservoir,
                // holds the frame's bits.
                let mut pick = 14u8;
                for i in 1..15u8 {
                    let slot = self.slot_bytes(i, false);
                    if 8 * slot as u32 + 8 * self.free_tail as u32 / 2 >= need {
                        pick = i;
                        break;
                    }
                }
                let slot = self.slot_bytes(pick, false);
                let avail = 8 * (self.free_tail + slot) as u32;
                let q = if need <= avail && total(&q0).1 { q0 } else { self.fit(&quantise_all, &total, q0, avail) };
                (q, pick)
            }
        };
        self.write_frame(&q, &grs, index, padding, ms, out);
    }

    /// Bisect the noise-to-mask offset for the finest quantisation within
    /// `target` bits (and every granule within 4095).
    fn fit(
        &self,
        quantise_all: &dyn Fn(&Encoder, f64) -> FrameQuant,
        total: &dyn Fn(&FrameQuant) -> (u32, bool),
        at_zero: FrameQuant,
        target: u32,
    ) -> FrameQuant {
        let fits = |q: &Vec<Vec<Quantised>>| {
            let (b, ok) = total(q);
            b <= target && ok
        };
        let (mut lo, mut hi);
        let mut best;
        if fits(&at_zero) {
            // Room to spare: look for a finer quantisation.
            best = at_zero;
            lo = -30.0;
            hi = 0.0;
            let q = quantise_all(self, lo);
            if fits(&q) {
                return q;
            }
        } else {
            lo = 0.0;
            hi = 60.0;
            best = quantise_all(self, hi);
            if !fits(&best) {
                // Even far above the mask the frame does not fit: coarser
                // still, until it does.
                let mut off = hi;
                while !fits(&best) && off < 200.0 {
                    off += 20.0;
                    best = quantise_all(self, off);
                }
                return best;
            }
        }
        for _ in 0..9 {
            let mid = 0.5 * (lo + hi);
            let q = quantise_all(self, mid);
            if fits(&q) {
                best = q;
                hi = mid;
            } else {
                lo = mid;
            }
        }
        best
    }

    /// Mid/side when it saves perceptual entropy.
    fn choose_ms(&self, grs: &[Granule]) -> bool {
        let mut lr = 0.0;
        let mut ms = 0.0;
        let s2 = std::f32::consts::FRAC_1_SQRT_2;
        for g in grs {
            lr += g.mask[0].pe(&self.long_edges, &self.short_edges) + g.mask[1].pe(&self.long_edges, &self.short_edges);
            let mask = g.mask[0].min(&g.mask[1]);
            let m: [f32; 576] = std::array::from_fn(|i| (g.xr[0][i] + g.xr[1][i]) * s2);
            let s: [f32; 576] = std::array::from_fn(|i| (g.xr[0][i] - g.xr[1][i]) * s2);
            ms += mask.with_energy_of(&m, &self.long_edges, &self.short_edges).pe(&self.long_edges, &self.short_edges)
                + mask.with_energy_of(&s, &self.long_edges, &self.short_edges).pe(&self.long_edges, &self.short_edges);
        }
        ms < lr
    }

    fn write_frame(
        &mut self,
        q: &[Vec<Quantised>],
        grs: &[Granule],
        index: u8,
        padding: bool,
        ms: bool,
        out: &mut Vec<Vec<u8>>,
    ) {
        let h = self.header(index, padding, ms);
        let frame_len = h.frame_len().unwrap_or(0);
        let slot = frame_len - h.header_len() - h.side_info_len();
        let rate = rate_index(h.version_index(), usize::from(h.sample_rate_index));
        let region1_short = 3 * usize::from(SFB_SHORT[rate][3]);
        let region1_long = usize::from(SFB_LONG[rate][8]);
        // Main data.
        let mut main = BitWriter::new();
        let mut si = SideInfo { main_data_begin: self.free_tail, ..Default::default() };
        for (gr, gq) in q.iter().enumerate() {
            let bt = grs[gr].block_type;
            for (ch, qc) in gq.iter().enumerate() {
                let start = main.len();
                quantize::write_scalefactors(&mut main, qc, bt == 2, self.lsf);
                let r1 = if bt == 2 { region1_short } else { region1_long };
                huffman::write(&mut main, &qc.ix, &qc.coding, &self.long_edges, bt != 0, r1);
                let len = main.len() - start;
                debug_assert_eq!(len as u32, qc.bits(), "bit count");
                si.gr[gr][ch] = GranuleInfo {
                    part2_3_length: len as u16,
                    big_values: qc.coding.big_values,
                    global_gain: qc.global_gain,
                    scalefac_compress: qc.scalefac_compress,
                    window_switching: bt != 0,
                    block_type: bt,
                    mixed_block: false,
                    table_select: qc.coding.table_select,
                    subblock_gain: qc.subblock_gain,
                    region0_count: qc.coding.region0_count,
                    region1_count: qc.coding.region1_count,
                    preflag: qc.preflag,
                    scalefac_scale: qc.scalefac_scale,
                    count1_table_b: qc.coding.count1_table_b,
                };
            }
        }
        let mut main = main.into_bytes();
        let resv_max = self.reservoir_max(frame_len);
        let mut free = self.free_tail + slot - main.len().min(self.free_tail + slot);
        if free > resv_max {
            // Stuffing (ancillary bytes) keeps main_data_begin in range.
            main.resize(main.len() + free - resv_max, 0);
            free = resv_max;
        }
        // Header, CRC and side information.
        let mut w = BitWriter::new();
        for b in h.to_bytes() {
            w.put(u32::from(b), 8);
        }
        if h.crc {
            w.put(0, 16);
        }
        si.write(&h, &mut w);
        let mut bytes = w.into_bytes();
        if h.crc {
            let crc = frame_crc_bits(0xFFFF, &bytes, 16, 16);
            let crc = frame_crc_bits(crc, &bytes, 48, h.side_info_len() * 8);
            bytes[4..6].copy_from_slice(&crc.to_be_bytes());
        }
        let slot_start = bytes.len();
        bytes.resize(frame_len, 0);
        self.pending.push_back(Pending { bytes, slot_start, filled: 0 });
        // Place the main data: first the free tails of earlier frames, then
        // this frame's slot.
        let mut data = &main[..];
        let skip_tail = {
            // Bytes of pending slots before the free region.
            let total_slots: usize = self.pending.iter().map(|p| p.bytes.len() - p.slot_start).sum();
            let used: usize = self.pending.iter().map(|p| p.filled).sum();
            total_slots - used - (self.free_tail + slot)
        };
        debug_assert_eq!(skip_tail, 0);
        for p in self.pending.iter_mut() {
            let cap = p.bytes.len() - p.slot_start - p.filled;
            if cap == 0 || data.is_empty() {
                continue;
            }
            let n = cap.min(data.len());
            let at = p.slot_start + p.filled;
            p.bytes[at..at + n].copy_from_slice(&data[..n]);
            p.filled += n;
            data = &data[n..];
        }
        debug_assert!(data.is_empty(), "main data overflows the reservoir");
        self.free_tail = free;
        self.bitrate_sum += u64::from(h.bitrate());
        // Emit every frame whose slot is full (free bytes are only ever at the
        // end, within the reservoir's reach, so older frames fill up).
        while let Some(front) = self.pending.front() {
            let unfilled = front.bytes.len() - front.slot_start - front.filled;
            if unfilled == 0 {
                let p = self.pending.pop_front().expect("front");
                self.emit(p.bytes, out);
            } else {
                break;
            }
        }
    }

    fn emit(&mut self, bytes: Vec<u8>, out: &mut Vec<Vec<u8>>) {
        self.frame_offsets.push(self.bytes_out);
        self.bytes_out += bytes.len() as u64;
        self.frames_out += 1;
        self.music_crc = crc16_arc_update(self.music_crc, &bytes);
        out.push(bytes);
    }

    /// The Xing (VBR) or Info (CBR) tag frame for the start of an `.mp3`
    /// file: frame count, byte count, seek table and the LAME-style
    /// extension with encoder delay and padding. Complete after
    /// [`Self::flush`]; it belongs before the first audio frame and is not
    /// audio (decoders that know it skip it; others decode it as one frame
    /// of silence).
    pub fn tag_frame(&self) -> Vec<u8> {
        let vbr = matches!(self.cfg.bitrate, BitrateMode::Vbr(_));
        let table = self.bitrate_table();
        let want = match self.cfg.bitrate {
            BitrateMode::Cbr(b) => table.iter().position(|&r| r == b / 1000).unwrap_or(1),
            BitrateMode::Vbr(_) => 1,
        };
        let probe = self.header(1, false, false);
        let need = xing::xing_offset(&probe) + 120 + 36;
        let index = (want..15)
            .chain(1..want)
            .find(|&i| self.header(i as u8, false, false).frame_len().unwrap_or(0) >= need)
            .unwrap_or(14) as u8;
        let mut h = self.header(index, false, false);
        h.crc = false;
        if h.mode == Mode::JointStereo {
            h.mode_extension = 0;
        }
        let len = h.frame_len().unwrap_or(0);
        let mut frame = vec![0u8; len];
        frame[..4].copy_from_slice(&h.to_bytes());
        let total_bytes = self.bytes_out + len as u64;
        let mut toc = [0u8; 100];
        if self.frames_out > 0 {
            for (i, t) in toc.iter_mut().enumerate() {
                let f = ((i as u64 * self.frames_out) / 100) as usize;
                let off = self.frame_offsets.get(f).copied().unwrap_or(0) + len as u64;
                *t = ((off * 256) / total_bytes.max(1)).min(255) as u8;
            }
        }
        let quality = match self.cfg.bitrate {
            BitrateMode::Vbr(q) => u32::from(q) * 10,
            BitrateMode::Cbr(_) => 50,
        };
        let stereo_mode = match (self.nch, self.cfg.joint_stereo) {
            (1, _) => 0,
            (_, true) => 3,
            _ => 1,
        };
        let source_rate_code = match self.cfg.sample_rate {
            r if r <= 32_000 => 0,
            44_100 => 1,
            48_000 => 2,
            _ => 3,
        };
        let bitrate_kbps = match self.cfg.bitrate {
            BitrateMode::Cbr(b) => b / 1000,
            BitrateMode::Vbr(_) => table[1],
        };
        xing::write_tag(
            &mut frame,
            &h,
            &xing::TagContents {
                vbr,
                frames: self.frames_out as u32,
                bytes: total_bytes.min(u64::from(u32::MAX)) as u32,
                toc,
                quality,
                encoder: ENCODER_NAME,
                vbr_method: if vbr { 4 } else { 1 },
                lowpass_hz: (self.cutoff_long as f64 / 576.0 * f64::from(self.cfg.sample_rate) / 2.0) as u32,
                bitrate_kbps,
                delay: ENCODER_DELAY,
                padding: self.padding,
                stereo_mode,
                source_rate_code,
                music_crc: self.music_crc,
            },
        );
        frame
    }
}

/// `f` over `items`, in order, on up to `threads` threads (the calling
/// thread among them), each taking the next item as it finishes one.
fn parallel_map<T: Send, R: Send>(items: Vec<T>, threads: usize, f: impl Fn(T) -> R + Sync) -> Vec<R> {
    let n = items.len();
    if threads <= 1 || n <= 1 {
        return items.into_iter().map(f).collect();
    }
    let queue = std::sync::Mutex::new(items.into_iter().enumerate());
    let work = || {
        let mut done = Vec::new();
        loop {
            let next = queue.lock().expect("work queue").next();
            let Some((i, item)) = next else { return done };
            done.push((i, f(item)));
        }
    };
    let mut all: Vec<(usize, R)> = std::thread::scope(|s| {
        let helpers: Vec<_> = (1..threads.min(n)).map(|_| s.spawn(work)).collect();
        let mut all = work();
        for h in helpers {
            all.extend(h.join().expect("an encoder thread panicked"));
        }
        all
    });
    all.sort_unstable_by_key(|(i, _)| *i);
    all.into_iter().map(|(_, r)| r).collect()
}

/// [`Encoder::quantise_frame`] on the calling thread alone.
fn quantise_serial(work: &FrameWork, offset_db: f64) -> FrameQuant {
    let scale = 10f64.powf(offset_db / 10.0);
    let mut items = (0..work.items.len()).map(|i| work.quantise(i, scale));
    (0..work.items.len() / work.nch).map(|_| items.by_ref().take(work.nch).collect()).collect()
}
