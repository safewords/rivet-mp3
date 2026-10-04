//! Encoder and decoder throughput:
//! `cargo run --release --example bench -- <pcm.raw> [seconds] [runs] [file.mp3 …]`.
//!
//! `pcm.raw` is 16-bit little-endian stereo PCM at 44.1 kHz. Its first
//! `seconds` (default 60) are encoded at 128 kb/s CBR, 320 kb/s CBR and
//! VBR quality 2, and each stream decoded back, each the best of `runs`
//! (default 3). Each extra MP3 file named is decoded too. Every encoded
//! stream's FNV-1a hash is printed, so a change to the encoder's output
//! shows. Encoding runs on one thread, then on all of them (which must give
//! the same stream).

use std::time::Instant;

use mp3::{BitrateMode, Decoder, Encoder, EncoderConfig};

fn best<T>(runs: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut out = None;
    let mut t = f64::MAX;
    for _ in 0..runs {
        let s = Instant::now();
        let v = f();
        t = t.min(s.elapsed().as_secs_f64());
        out = Some(v);
    }
    (t, out.expect("at least one run"))
}

fn decode(bytes: &[u8]) -> (Vec<f32>, u32, u8) {
    let frames = Decoder::decode_all(bytes).expect("decode");
    let rate = frames.first().map_or(44_100, |f| f.sample_rate);
    let ch = frames.first().map_or(2, |f| f.channels);
    (frames.into_iter().flat_map(|f| f.samples).collect(), rate, ch)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = args.get(1).expect("usage: bench <pcm.raw> [seconds] [runs] [file.mp3 ...]");
    let seconds: f64 = args.get(2).map_or(60.0, |s| s.parse().expect("seconds"));
    let runs: usize = args.get(3).map_or(3, |s| s.parse().expect("runs"));
    let raw = std::fs::read(path).expect("read");
    let take = ((seconds * 44_100.0) as usize * 2).min(raw.len() / 2);
    let pcm: Vec<f32> =
        raw.as_chunks::<2>().0.iter().take(take).map(|b| f32::from(i16::from_le_bytes(*b)) / 32768.0).collect();
    let dur = pcm.len() as f64 / 2.0 / 44_100.0;
    println!("{dur:.1} s of stereo 44.1 kHz");
    for (name, bitrate) in [
        ("CBR 128", BitrateMode::Cbr(128_000)),
        ("CBR 320", BitrateMode::Cbr(320_000)),
        ("VBR q2", BitrateMode::Vbr(2)),
    ] {
        let encode = |threads: usize| {
            let mut e = Encoder::new(EncoderConfig { bitrate, ..EncoderConfig::default() }).expect("encoder");
            e.set_threads(threads);
            let mut out: Vec<u8> = e.encode(&pcm).into_iter().flatten().collect();
            out.extend(e.flush().into_iter().flatten());
            out
        };
        let (te, stream) = best(runs, || encode(1));
        let (tm, threaded) = best(runs, || encode(0));
        assert!(threaded == stream, "the threaded encoder's stream differs");
        let mut h = 0xcbf2_9ce4_8422_2325u64;
        for &b in &stream {
            h = (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3);
        }
        let (td, (out, _, _)) = best(runs, || decode(&stream));
        println!(
            "{name:8} encode {:7.1} x realtime, {:7.1} x threaded ({} kB, hash {h:016x}); decode {:7.1} x realtime ({} samples)",
            dur / te,
            dur / tm,
            stream.len() / 1000,
            dur / td,
            out.len()
        );
    }
    for file in args.iter().skip(4) {
        let bytes = std::fs::read(file).expect("read mp3");
        let (td, (out, rate, ch)) = best(runs, || decode(&bytes));
        let d = out.len() as f64 / f64::from(ch) / f64::from(rate);
        println!("decode {file}: {:7.1} x realtime ({d:.1} s)", d / td);
    }
}
