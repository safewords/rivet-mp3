//! Round trips through this crate's decoder, run strictly (every CRC
//! checked, the bit reservoir accounted exactly, every Huffman region ending
//! at its part2_3_length), measured by SNR and by noise-to-mask ratio per
//! scalefactor band under the encoder's own psychoacoustic model.

use super::*;
use crate::decode::{Decoder, DecoderOptions};
use crate::header::FrameHeader;
use crate::xing::InfoHeader;

/// A deterministic generator (64-bit LCG).
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }
}

/// A test programme: plucked harmonic notes with sharp attacks, a sustained
/// chord, coloured noise, castanet-like clicks; partly correlated stereo.
fn programme(seconds: f64, rate: u32, nch: usize) -> Vec<f32> {
    use std::f64::consts::PI;
    let n = (seconds * f64::from(rate)) as usize;
    let fs = f64::from(rate);
    let mut rng = Rng(0x5eed);
    let notes = [220.0, 277.18, 329.63, 440.0, 392.0, 293.66, 246.94, 523.25];
    let clicks = [0.9, 1.45, 1.8, 2.6];
    let mut lp = [0.0f64; 2];
    let mut out = Vec::with_capacity(n * nch);
    for i in 0..n {
        let t = i as f64 / fs;
        // Notes every 0.35 s.
        let k = (t / 0.35) as usize;
        let tn = t - k as f64 * 0.35;
        let f0 = notes[k % notes.len()];
        let mut note = 0.0;
        for h in 1..12 {
            let fh = f0 * h as f64;
            if fh > fs * 0.45 {
                break;
            }
            note += (2.0 * PI * fh * tn).sin() * (-tn * (3.0 + h as f64)).exp() / h as f64;
        }
        let pad = 0.06 * ((2.0 * PI * 130.81 * t).sin() + (2.0 * PI * 196.0 * t).sin() + (2.0 * PI * 311.13 * t).sin());
        let mut click = 0.0;
        for &c in &clicks {
            let d = t - c;
            if (0.0..0.03).contains(&d) {
                click += 0.6 * rng.next() * (-d * 200.0).exp();
            }
        }
        for c in 0..nch {
            let white = rng.next();
            lp[c] = 0.9 * lp[c] + 0.1 * white;
            let noise = 0.05 * lp[c] + 0.004 * white;
            let pan = if c == 0 { 0.8 } else { 0.55 };
            let x = 0.25 * note * pan + pad + noise + click * if c == 0 { 1.0 } else { 0.7 };
            out.push(x.clamp(-1.0, 1.0) as f32);
        }
    }
    out
}

/// Encode `pcm` and return the complete stream (tag frame first) and the
/// encoder.
fn encode(cfg: EncoderConfig, pcm: &[f32]) -> (Vec<u8>, Encoder, Vec<Vec<u8>>) {
    let mut enc = Encoder::new(cfg).unwrap();
    let mut frames = Vec::new();
    // Feed in uneven chunks.
    let nch = usize::from(cfg.channels);
    let mut at = 0;
    let mut step = 1000;
    while at < pcm.len() / nch {
        let end = (at + step).min(pcm.len() / nch);
        frames.extend(enc.encode(&pcm[at * nch..end * nch]));
        at = end;
        step = step * 7 % 3001 + 100;
    }
    frames.extend(enc.flush());
    let mut stream = enc.tag_frame();
    for f in &frames {
        stream.extend_from_slice(f);
    }
    (stream, enc, frames)
}

/// Decode strictly, gapless-trimmed.
fn decode_strict(stream: &[u8]) -> (Vec<f32>, Decoder) {
    let mut d = Decoder::with_options(DecoderOptions { strict: true, check_crc: true, trim_gapless: true });
    let mut out = Vec::new();
    for chunk in stream.chunks(777) {
        out.extend(d.decode(chunk).unwrap_or_else(|e| panic!("strict decode: {e}")));
    }
    out.extend(d.flush().unwrap_or_else(|e| panic!("strict decode: {e}")));
    (out.iter().flat_map(|f| f.samples.iter().copied()).collect(), d)
}

/// Long-block spectra of a channel, granule by granule, as the encoder
/// computes them (analysis filterbank, MDCT, antialias butterflies).
fn spectra(x: &[f32]) -> Vec<[f32; 576]> {
    let mut a = Analysis::default();
    let mut prev = [[0.0f64; 32]; 18];
    let mut out = Vec::new();
    for g in 0..x.len() / 576 {
        let mut slots = [[0.0f64; 32]; 18];
        for (t, s) in slots.iter_mut().enumerate() {
            a.run(&x[g * 576 + 32 * t..g * 576 + 32 * t + 32], s);
        }
        let mut xr = [0.0f32; 576];
        for sb in 0..32 {
            let mut v = [0.0f64; 36];
            for i in 0..36 {
                let s = if i < 18 { prev[i][sb] } else { slots[i - 18][sb] };
                v[i] = if sb % 2 == 1 && i % 2 == 1 { -s } else { s };
            }
            let mut c = [0.0f32; 18];
            mdct(&v, 0, &mut c);
            xr[sb * 18..sb * 18 + 18].copy_from_slice(&c);
        }
        prev = slots;
        out.push(xr);
    }
    out
}

/// What a round trip measured.
#[derive(Debug, Default)]
struct Measure {
    snr_db: f64,
    /// Mean NMR over (band, granule) pairs below the cutoff, dB.
    nmr_mean_db: f64,
    /// 95th percentile NMR, dB.
    nmr_p95_db: f64,
    /// Share of (band, granule) pairs whose noise exceeds the mask.
    over_mask: f64,
    /// Mean NMR per long band, dB.
    per_band: [f64; 22],
}

fn measure(orig: &[f32], dec: &[f32], nch: usize, rate: u32, cutoff_hz: f64) -> Measure {
    let (vi, si) = SAMPLE_RATES
        .iter()
        .enumerate()
        .find_map(|(v, row)| row.iter().position(|&r| r == rate).map(|s| (v, s)))
        .unwrap();
    let r = rate_index(vi, si);
    let psy = Psy::new(rate, &SFB_LONG[r], &SFB_SHORT[r]);
    let mut snr = f64::MAX;
    let mut nmrs = Vec::new();
    let mut band_sum = [0.0f64; 22];
    let mut band_n = [0usize; 22];
    let line_hz = f64::from(rate) / 2.0 / 576.0;
    for c in 0..nch {
        let x: Vec<f32> = orig.iter().skip(c).step_by(nch).copied().collect();
        let y: Vec<f32> = dec.iter().skip(c).step_by(nch).copied().collect();
        let e: Vec<f32> = x.iter().zip(&y).map(|(a, b)| b - a).collect();
        let sig: f64 = x.iter().map(|&v| f64::from(v).powi(2)).sum();
        let err: f64 = e.iter().map(|&v| f64::from(v).powi(2)).sum();
        snr = snr.min(10.0 * (sig / err.max(1e-30)).log10());
        let sx = spectra(&x);
        let se = spectra(&e);
        let mut prev: Option<[f64; 22]> = None;
        for (g, (a, b)) in sx.iter().zip(&se).enumerate() {
            let mask = psy.long(a, prev.as_ref());
            let Mask::Long { thr, .. } = mask else { unreachable!() };
            prev = Some(thr);
            if g < 2 || g + 2 >= sx.len() {
                continue;
            }
            for band in 0..22 {
                let lo = usize::from(SFB_LONG[r][band]);
                let hi = usize::from(SFB_LONG[r][band + 1]);
                if hi as f64 * line_hz > cutoff_hz {
                    continue;
                }
                let n: f64 = b[lo..hi].iter().map(|&v| f64::from(v).powi(2)).sum();
                let nmr = 10.0 * (n.max(1e-30) / thr[band]).log10();
                nmrs.push(nmr);
                band_sum[band] += nmr;
                band_n[band] += 1;
            }
        }
    }
    nmrs.sort_by(f64::total_cmp);
    let mean = nmrs.iter().sum::<f64>() / nmrs.len().max(1) as f64;
    let p95 = nmrs.get(nmrs.len() * 95 / 100).copied().unwrap_or(0.0);
    let over = nmrs.iter().filter(|&&v| v > 0.0).count() as f64 / nmrs.len().max(1) as f64;
    Measure {
        snr_db: snr,
        nmr_mean_db: mean,
        nmr_p95_db: p95,
        over_mask: over,
        per_band: std::array::from_fn(|b| if band_n[b] > 0 { band_sum[b] / band_n[b] as f64 } else { f64::NAN }),
    }
}

/// Every frame of `frames` parses, and what its headers say.
fn headers(frames: &[Vec<u8>]) -> Vec<FrameHeader> {
    frames
        .iter()
        .map(|f| {
            let h = FrameHeader::parse(f).unwrap();
            assert_eq!(h.frame_len(), Some(f.len()), "frame length matches its header");
            h
        })
        .collect()
}

struct Case {
    rate: u32,
    nch: u8,
    mode: BitrateMode,
    joint: bool,
    crc: bool,
}

fn case(rate: u32, nch: u8, mode: BitrateMode) -> Case {
    Case { rate, nch, mode, joint: true, crc: false }
}

/// Run one configuration: strict decode, exact length, tag fields, and the
/// figures.
fn run_case(c: &Case, seconds: f64) -> (Measure, f64, usize, usize) {
    let pcm = programme(seconds, c.rate, usize::from(c.nch));
    let cfg = EncoderConfig {
        sample_rate: c.rate,
        channels: c.nch,
        bitrate: c.mode,
        joint_stereo: c.joint,
        crc: c.crc,
        lowpass: None,
    };
    let (stream, enc, frames) = encode(cfg, &pcm);
    let hs = headers(&frames);
    let (dec, d) = decode_strict(&stream);
    assert_eq!(dec.len(), pcm.len(), "gapless: decoded length equals the input");
    // The tag.
    let Some(InfoHeader::Xing(x)) = d.info() else { panic!("no Xing/Info tag") };
    assert_eq!(x.is_info, matches!(c.mode, BitrateMode::Cbr(_)));
    assert_eq!(x.frames, Some(frames.len() as u32));
    assert_eq!(x.bytes, Some(stream.len() as u32));
    let lame = x.lame.as_ref().unwrap();
    assert!(lame.tag_crc_ok, "tag CRC");
    assert_eq!(lame.encoder_delay, ENCODER_DELAY);
    assert_eq!(lame.padding, enc.padding());
    assert!(enc.padding() >= DECODER_DELAY);
    let music: Vec<u8> = frames.concat();
    assert_eq!(lame.music_crc, crc16_arc_update(0, &music), "music CRC");
    let g = d.gapless().unwrap();
    assert_eq!(g.length, Some(pcm.len() as u64 / u64::from(c.nch)));
    if let BitrateMode::Cbr(b) = c.mode {
        assert!(hs.iter().all(|h| h.bitrate() == b), "CBR: every frame at {b}");
    }
    assert!(hs.iter().all(|h| h.crc == c.crc));
    let ms_frames = hs.iter().filter(|h| h.mode_extension & 2 != 0).count();
    let cutoff = enc.cutoff_long as f64 / 576.0 * f64::from(c.rate) / 2.0;
    let m = measure(&pcm, &dec, usize::from(c.nch), c.rate, cutoff);
    let short = count_short_granules(&frames);
    let _ = ms_frames;
    (m, enc.average_bitrate(), short, ms_frames)
}

/// Granules coded with short blocks, from the side information.
fn count_short_granules(frames: &[Vec<u8>]) -> usize {
    let mut n = 0;
    for f in frames {
        let h = FrameHeader::parse(f).unwrap();
        let si = SideInfo::parse(&f[h.header_len()..], &h).unwrap();
        let ngr = if h.version.is_lsf() { 1 } else { 2 };
        for gr in 0..ngr {
            for ch in 0..h.channels() {
                if si.gr[gr][ch].block_type == 2 {
                    n += 1;
                }
            }
        }
    }
    n
}

fn label(c: &Case) -> String {
    let mode = match c.mode {
        BitrateMode::Cbr(b) => format!("CBR {:>3}k", b / 1000),
        BitrateMode::Vbr(q) => format!("VBR q{q}  "),
    };
    format!(
        "{:>5} Hz {} {} {}{}",
        c.rate,
        if c.nch == 1 { "mono  " } else if c.joint { "joint " } else { "stereo" },
        mode,
        if c.crc { "CRC" } else { "   " },
        ""
    )
}

/// The configuration matrix. Every case: strict decode (CRC, reservoir
/// accounting, part2_3_length), exact gapless length, tag fields; and
/// quality floors. `MP3_ENCODER_REPORT=1` prints the table.
#[test]
fn round_trips() {
    use BitrateMode::{Cbr, Vbr};
    let mut cases = vec![
        case(44_100, 2, Cbr(32_000)),
        case(44_100, 2, Cbr(64_000)),
        case(44_100, 2, Cbr(96_000)),
        case(44_100, 2, Cbr(128_000)),
        case(44_100, 2, Cbr(160_000)),
        case(44_100, 2, Cbr(192_000)),
        case(44_100, 2, Cbr(256_000)),
        case(44_100, 2, Cbr(320_000)),
        Case { joint: false, ..case(44_100, 2, Cbr(128_000)) },
        Case { crc: true, ..case(44_100, 2, Cbr(128_000)) },
        case(44_100, 1, Cbr(64_000)),
        case(44_100, 1, Cbr(128_000)),
        case(48_000, 2, Cbr(128_000)),
        case(48_000, 2, Cbr(320_000)),
        case(32_000, 2, Cbr(64_000)),
        case(32_000, 2, Cbr(128_000)),
        case(44_100, 2, Vbr(0)),
        case(44_100, 2, Vbr(2)),
        case(44_100, 2, Vbr(4)),
        case(44_100, 2, Vbr(6)),
        case(44_100, 2, Vbr(9)),
        case(48_000, 2, Vbr(4)),
        case(32_000, 1, Vbr(4)),
        case(22_050, 2, Cbr(32_000)),
        case(22_050, 2, Cbr(64_000)),
        case(22_050, 2, Cbr(160_000)),
        case(24_000, 2, Cbr(96_000)),
        Case { crc: true, ..case(16_000, 1, Cbr(32_000)) },
        case(22_050, 2, Vbr(4)),
        case(11_025, 2, Cbr(32_000)),
        case(12_000, 1, Cbr(24_000)),
        case(8_000, 1, Cbr(8_000)),
        case(8_000, 1, Cbr(16_000)),
        case(11_025, 2, Vbr(4)),
    ];
    if cfg!(debug_assertions) && std::env::var_os("MP3_ENCODER_REPORT").is_none() {
        // Unoptimised builds: a representative subset.
        let keep = [3, 7, 9, 10, 18, 25, 27, 31];
        cases = cases.into_iter().enumerate().filter(|(i, _)| keep.contains(i)).map(|(_, c)| c).collect();
    }
    let report = std::env::var_os("MP3_ENCODER_REPORT").is_some();
    if report {
        println!(
            "{:<36} {:>8} {:>8} {:>9} {:>9} {:>7} {:>6} {:>5}",
            "configuration", "kbit/s", "SNR dB", "NMR mean", "NMR p95", ">mask", "short", "M/S"
        );
    }
    for c in &cases {
        let (m, br, short, ms) = run_case(c, 3.0);
        if report {
            println!(
                "{:<36} {:>8.1} {:>8.1} {:>9.1} {:>9.1} {:>6.1}% {:>6} {:>5}",
                label(c),
                br / 1000.0,
                m.snr_db,
                m.nmr_mean_db,
                m.nmr_p95_db,
                100.0 * m.over_mask,
                short,
                ms
            );
        }
        // Quality floors: every configuration keeps its signal; the
        // transparent-range rates keep the mean noise below the mask.
        assert!(m.snr_db > 2.0, "{}: SNR {:.1}", label(c), m.snr_db);
        let rich = match c.mode {
            BitrateMode::Cbr(b) => b / u32::from(c.nch) >= 96_000 && c.rate >= 32_000,
            BitrateMode::Vbr(q) => q <= 4 && c.rate >= 32_000,
        };
        if rich {
            assert!(m.nmr_mean_db < 0.0, "{}: mean NMR {:.1} dB", label(c), m.nmr_mean_db);
        }
        if report && std::env::var_os("MP3_ENCODER_BANDS").is_some() {
            let v: Vec<String> = m.per_band.iter().map(|x| format!("{x:.0}")).collect();
            println!("    NMR by band: {}", v.join(" "));
        }
    }
}

/// The codec's delay is what the tag says: an impulse comes back
/// `ENCODER_DELAY + 529` samples late in the untrimmed output.
#[test]
fn delay_matches_the_tag() {
    let rate = 44_100;
    let mut pcm = vec![0.0f32; rate as usize / 2];
    pcm[5000] = 0.9;
    let (stream, _, _) = encode(
        EncoderConfig { channels: 1, bitrate: BitrateMode::Cbr(320_000), ..Default::default() },
        &pcm,
    );
    let mut d = Decoder::with_options(DecoderOptions { trim_gapless: false, ..Default::default() });
    let mut out = d.decode(&stream).unwrap();
    out.extend(d.flush().unwrap());
    let dec: Vec<f32> = out.iter().flat_map(|f| f.samples.iter().copied()).collect();
    let peak = dec.iter().enumerate().max_by(|a, b| a.1.abs().total_cmp(&b.1.abs())).unwrap().0;
    assert_eq!(peak, 5000 + (ENCODER_DELAY + DECODER_DELAY) as usize);
}

/// Transients switch to short blocks, and the noise before an attack
/// stays far below the attack (no audible pre-echo smear across a long
/// window).
#[test]
fn attacks_use_short_blocks() {
    let rate = 44_100u32;
    let n = rate as usize;
    let mut rng = Rng(7);
    let mut pcm = vec![0.0f32; n];
    let at = 20_000;
    for (i, v) in pcm.iter_mut().enumerate() {
        let quiet = 0.001 * (i as f32 * 0.03).sin();
        *v = if i >= at { quiet + 0.7 * rng.next() as f32 * (-((i - at) as f32) / 3000.0).exp() } else { quiet };
    }
    let (stream, _, frames) =
        encode(EncoderConfig { channels: 1, bitrate: BitrateMode::Cbr(128_000), ..Default::default() }, &pcm);
    assert!(count_short_granules(&frames) >= 1, "no short blocks for an attack");
    let (dec, _) = decode_strict(&stream);
    let pre: f64 = (at - 576..at - 64).map(|i| f64::from(dec[i] - pcm[i]).powi(2)).sum::<f64>() / 512.0;
    let post: f64 = (at..at + 512).map(|i| f64::from(pcm[i]).powi(2)).sum::<f64>() / 512.0;
    let ratio = 10.0 * (post / pre.max(1e-30)).log10();
    assert!(ratio > 30.0, "pre-echo only {ratio:.1} dB below the attack");
}

/// Correlated stereo is coded mid/side; independent channels mostly not.
#[test]
fn mid_side_follows_correlation() {
    let rate = 44_100u32;
    let n = rate as usize;
    let mut rng = Rng(3);
    let mut same = Vec::with_capacity(2 * n);
    let mut apart = Vec::with_capacity(2 * n);
    for i in 0..n {
        let t = i as f32 / rate as f32;
        let a = 0.3 * (2.0 * std::f32::consts::PI * 440.0 * t).sin() + 0.02 * rng.next() as f32;
        let b = 0.3 * (2.0 * std::f32::consts::PI * 1234.0 * t).sin() + 0.02 * rng.next() as f32;
        same.extend([a, a * 0.98]);
        apart.extend([a, b]);
    }
    let cfg = EncoderConfig::default();
    let (_, _, f1) = encode(cfg, &same);
    let (_, _, f2) = encode(cfg, &apart);
    let ms = |f: &[Vec<u8>]| headers(f).iter().filter(|h| h.mode_extension & 2 != 0).count() as f64 / f.len() as f64;
    assert!(ms(&f1) > 0.9, "correlated: {}", ms(&f1));
    assert!(ms(&f2) < 0.5, "independent: {}", ms(&f2));
}

/// The bit reservoir is accounted exactly: walking the frames, each
/// frame's main_data_begin never exceeds the bytes the frames before it
/// left unused, nor the field's range.
#[test]
fn reservoir_accounting() {
    for mode in [BitrateMode::Cbr(64_000), BitrateMode::Cbr(256_000), BitrateMode::Vbr(3)] {
        let pcm = programme(2.0, 44_100, 2);
        let (_, _, frames) = encode(EncoderConfig { bitrate: mode, ..Default::default() }, &pcm);
        let mut unused = 0usize;
        let mut used_reservoir = false;
        for f in &frames {
            let h = FrameHeader::parse(f).unwrap();
            let si = SideInfo::parse(&f[h.header_len()..], &h).unwrap();
            assert!(si.main_data_begin <= unused, "main_data_begin {} > {unused} unused", si.main_data_begin);
            assert!(si.main_data_begin <= 511);
            used_reservoir |= si.main_data_begin > 0;
            let bits: usize = si.gr.iter().flatten().map(|g| usize::from(g.part2_3_length)).sum();
            let have = (si.main_data_begin + f.len() - h.header_len() - h.side_info_len()) * 8;
            assert!(bits <= have);
            unused = (have - bits) / 8;
        }
        assert!(used_reservoir, "{mode:?}: the reservoir was never used");
    }
}

#[test]
fn configuration_errors() {
    let bad = |c: EncoderConfig| Encoder::new(c).err().unwrap().to_string();
    assert!(bad(EncoderConfig { sample_rate: 44_000, ..Default::default() }).contains("sampling frequency"));
    assert!(bad(EncoderConfig { channels: 3, ..Default::default() }).contains("one or two"));
    assert!(bad(EncoderConfig { bitrate: BitrateMode::Cbr(100_000), ..Default::default() }).contains("bit rate"));
    assert!(
        bad(EncoderConfig { sample_rate: 22_050, bitrate: BitrateMode::Cbr(192_000), ..Default::default() })
            .contains("bit rate")
    );
    assert!(bad(EncoderConfig { bitrate: BitrateMode::Vbr(10), ..Default::default() }).contains("quality"));
    assert_eq!(coding_rate(96_000), 48_000);
    assert_eq!(coding_rate(22_050), 22_050);
    assert_eq!(coding_rate(44_000), 44_100);
}

/// The quantiser keeps each band that needs coding within its allowed
/// noise. Noise does not fall monotonically with the step for bands of a
/// few lines, and the step a scalefactor lands on is finer than the one the
/// search found, so a few bands come out slightly above; the bound is on
/// how many.
#[test]
fn quantiser_meets_the_allowance() {
    let pcm = programme(1.0, 44_100, 1);
    let mut enc = Encoder::new(EncoderConfig { channels: 1, ..Default::default() }).unwrap();
    enc.input[0].extend_from_slice(&pcm);
    enc.samples_in = pcm.len() as u64;
    for g in 0..40 {
        enc.analyse_granule(g);
    }
    let (mut over, mut total) = (0, 0);
    for gr in enc.frame_granules.iter().filter(|g| g.block_type == 0) {
        let lines = quantize::Lines::new(&gr.xr[0]);
        let sp = quantize::Spectrum {
            lines: &lines,
            block_type: 0,
            mask: &gr.mask[0],
            long_edges: &SFB_LONG[0],
            short_edges: &SFB_SHORT[0],
            lsf: false,
        };
        let mut q = quantize::quantise(&sp, 1.0);
        quantize::apply_signs(&mut q, &gr.xr[0]);
        let r = quantize::dequantise(&q, &sp);
        let Mask::Long { thr, energy } = &gr.mask[0] else { unreachable!() };
        for b in 0..21 {
            let (lo, hi) = (usize::from(SFB_LONG[0][b]), usize::from(SFB_LONG[0][b + 1]));
            let n: f64 = (lo..hi).map(|i| f64::from(r[i] - gr.xr[0][i]).powi(2)).sum();
            if energy[b] > thr[b] {
                total += 1;
                over += usize::from(n > thr[b] * 1.26); // more than 1 dB above
            }
        }
    }
    assert!(total > 100);
    assert!(over * 20 < total, "{over} of {total} bands more than 1 dB above their allowance");
}
