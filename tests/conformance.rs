//! ISO conformance: the MPEG-1/2 audio test sequences (Layer I, II and III)
//! and their reference waveforms, from ISO/IEC 14496-26's public
//! conformance package (the sequences descend from ISO/IEC 11172-4 and
//! 13818-4), judged by ISO/IEC 11172-4's accuracy criterion:
//!
//! - **full accuracy**: RMS of (decoded - reference) below 2^-15 / sqrt(12)
//!   of full scale, and no sample differing by more than 2^-14;
//! - **limited accuracy**: RMS below 2^-11 / sqrt(12).
//!
//! The data is not in the repository (ISO's copyright, ~45 MB):
//! `python tools/fetch_conformance.py DIR` fetches it, and
//! `MP3_CONFORMANCE_DIR=DIR cargo test --release --test conformance -- --nocapture`
//! runs this. Without the variable the test reports that it skipped.

use std::path::{Path, PathBuf};

use mp3::{Decoder, DecoderOptions};

/// Interleaved samples at ±1.0 full scale, channels, rate.
fn read_wav(path: &Path) -> (Vec<f64>, usize, u32) {
    let b = std::fs::read(path).unwrap();
    assert_eq!(&b[..4], b"RIFF");
    assert_eq!(&b[8..12], b"WAVE");
    let mut p = 12;
    let (mut ch, mut rate, mut bits, mut float) = (0usize, 0u32, 0u16, false);
    let mut data: &[u8] = &[];
    while p + 8 <= b.len() {
        let id = &b[p..p + 4];
        let len = u32::from_le_bytes(b[p + 4..p + 8].try_into().unwrap()) as usize;
        let body = &b[p + 8..(p + 8 + len).min(b.len())];
        if id == b"fmt " {
            let tag = u16::from_le_bytes([body[0], body[1]]);
            ch = usize::from(u16::from_le_bytes([body[2], body[3]]));
            rate = u32::from_le_bytes(body[4..8].try_into().unwrap());
            bits = u16::from_le_bytes([body[14], body[15]]);
            let tag = if tag == 0xFFFE { u16::from_le_bytes([body[24], body[25]]) } else { tag };
            float = tag == 3;
        } else if id == b"data" {
            data = body;
        }
        p += 8 + len + (len & 1);
    }
    let bps = usize::from(bits / 8);
    let samples = data
        .chunks_exact(bps)
        .map(|s| match (bps, float) {
            (2, _) => f64::from(i16::from_le_bytes([s[0], s[1]])) / 32768.0,
            (3, _) => f64::from(i32::from_le_bytes([0, s[0], s[1], s[2]]) >> 8) / 8_388_608.0,
            (4, false) => f64::from(i32::from_le_bytes(s.try_into().unwrap())) / 2_147_483_648.0,
            (4, true) => f64::from(f32::from_le_bytes(s.try_into().unwrap())),
            _ => panic!("{bits}-bit WAV"),
        })
        .collect();
    (samples, ch, rate)
}

struct Outcome {
    name: String,
    detail: String,
    rms: f64,
    max: f64,
    verdict: &'static str,
}

/// Decode a stream into segments of constant channel count (a stream may
/// change mode, and ISO gives one reference per segment), leaving out
/// frames the decoder had to conceal: a Layer III frame whose main data
/// starts before the stream does produces no output in the reference
/// decoder, where this decoder emits silence to keep the timeline.
fn decode_segments(bytes: &[u8]) -> Result<Vec<Vec<mp3::Frame>>, mp3::Error> {
    let mut d = Decoder::with_options(DecoderOptions { trim_gapless: false, ..Default::default() });
    let mut frames = d.decode(bytes)?;
    frames.extend(d.flush()?);
    let mut segments: Vec<Vec<mp3::Frame>> = Vec::new();
    for f in frames.into_iter().filter(|f| !f.concealed) {
        match segments.last_mut() {
            Some(s) if s[0].channels == f.channels => s.push(f),
            _ => segments.push(vec![f]),
        }
    }
    Ok(segments)
}

fn judge(name: &str, frames: &[mp3::Frame], reference: &Path) -> Outcome {
    let name = name.to_string();
    let (want, wch, wrate) = read_wav(reference);
    let first = &frames[0];
    let ch = usize::from(first.channels);
    let rate = first.sample_rate;
    let got: Vec<f32> = frames.iter().flat_map(|f| f.samples.iter().copied()).collect();
    let h = first.header;
    let mut detail = format!(
        "{:?} L{} {} Hz {:?} {} kb/s{}{}",
        h.version,
        match h.layer {
            mp3::Layer::I => 1,
            mp3::Layer::II => 2,
            mp3::Layer::III => 3,
        },
        rate,
        h.mode,
        first.bitrate / 1000,
        if h.crc { " CRC" } else { "" },
        if h.is_free_format() { " free" } else { "" },
    );
    if wrate != rate {
        detail += &format!(" (reference rate {wrate})");
    }
    // Compare channel by channel where the reference has the same count; a
    // mono reference of a two-channel decode is compared with channel 0.
    let n_got = got.len() / ch;
    let n_want = want.len() / wch;
    let n = n_got.min(n_want);
    if n_got != n_want {
        detail += &format!(" (lengths {n_got} vs {n_want})");
    }
    let cmp_ch = wch.min(ch);
    let (mut sum, mut max) = (0.0f64, 0.0f64);
    for i in 0..n {
        for c in 0..cmp_ch {
            // The reference is PCM, clipped to [-1, 1); the decoder's float
            // output is not clipped, so clip it the same way first.
            let g = f64::from(got[i * ch + c]).clamp(-1.0, 1.0 - 2f64.powi(-23));
            let e = g - want[i * wch + c];
            sum += e * e;
            max = max.max(e.abs());
        }
    }
    let rms = (sum / (n * cmp_ch) as f64).sqrt();
    let full = 2f64.powi(-15) / 12f64.sqrt();
    let limited = 2f64.powi(-11) / 12f64.sqrt();
    let verdict = if rms < full && max <= 2f64.powi(-14) {
        "full"
    } else if rms < limited {
        "limited"
    } else {
        "FAIL"
    };
    Outcome { name, detail, rms, max, verdict }
}

/// The reference waveforms for a stream: `<name>.wav`, or for streams with
/// several, `<name>_<segment>[_<variant>].wav`, grouped by segment.
fn references(dir: &Path, name: &str) -> Vec<(usize, Vec<PathBuf>)> {
    let exact = dir.join(format!("{name}.wav"));
    if exact.exists() {
        return vec![(0, vec![exact])];
    }
    let mut by_seg: std::collections::BTreeMap<usize, Vec<PathBuf>> = Default::default();
    for e in std::fs::read_dir(dir).unwrap().filter_map(|e| e.ok()) {
        let p = e.path();
        let file = p.file_name().unwrap().to_string_lossy().to_string();
        let Some(rest) = file.strip_prefix(&format!("{name}_")).and_then(|r| r.strip_suffix(".wav")) else {
            continue;
        };
        if let Some(seg) = rest.split('_').next().and_then(|s| s.parse::<usize>().ok()) {
            by_seg.entry(seg).or_default().push(p);
        }
    }
    by_seg.into_iter().map(|(k, mut v)| {
        v.sort();
        (k, v)
    }).collect()
}

#[test]
fn iso_conformance_sequences() {
    let Some(dir) = std::env::var_os("MP3_CONFORMANCE_DIR") else {
        eprintln!("skipping: MP3_CONFORMANCE_DIR is not set (see tools/fetch_conformance.py)");
        return;
    };
    let dir = PathBuf::from(dir);
    let mut streams: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| p.extension().is_some_and(|e| e == "mpa"))
        .collect();
    streams.sort();
    let mut failures = Vec::new();
    println!("{:<22} {:<56} {:>10} {:>10}  verdict", "reference", "stream", "RMS", "max |e|");
    for s in &streams {
        let name = s.file_stem().unwrap().to_string_lossy().to_string();
        let refs = references(&dir, &name);
        if refs.is_empty() {
            continue;
        }
        let segments = match decode_segments(&std::fs::read(s).unwrap()) {
            Ok(seg) => seg,
            Err(e) => {
                println!("{name:<22} decode error: {e}");
                failures.push(name);
                continue;
            }
        };
        for (seg, variants) in refs {
            let Some(frames) = segments.get(seg) else {
                println!("{name:<22} segment {seg} missing ({} decoded)", segments.len());
                failures.push(format!("{name}_{seg}"));
                continue;
            };
            // Several variants of one segment's reference: the best counts.
            let best = variants
                .iter()
                .map(|r| judge(&r.file_stem().unwrap().to_string_lossy(), frames, r))
                .min_by(|a, b| a.rms.total_cmp(&b.rms))
                .unwrap();
            println!("{:<22} {:<56} {:>10.3e} {:>10.3e}  {}", best.name, best.detail, best.rms, best.max, best.verdict);
            if best.verdict == "FAIL" {
                failures.push(best.name);
            }
        }
    }
    assert!(failures.is_empty(), "failed: {failures:?}");
}
