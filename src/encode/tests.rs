//! Round trips through this crate's decoder.
use super::*;
use crate::decode::{Decoder, DecoderOptions};

fn tone(n: usize, nch: usize, rate: u32) -> Vec<f32> {
    let mut v = Vec::with_capacity(n * nch);
    for i in 0..n {
        let t = i as f64 / f64::from(rate);
        for c in 0..nch {
            let x = 0.3 * (2.0 * std::f64::consts::PI * 440.0 * t * (1.0 + 0.5 * c as f64)).sin()
                + 0.1 * (2.0 * std::f64::consts::PI * 3150.0 * t).sin();
            v.push(x as f32);
        }
    }
    v
}

#[test]
fn smoke() {
    let cfg = EncoderConfig::default();
    let mut enc = Encoder::new(cfg).unwrap();
    let pcm = tone(44_100, 2, 44_100);
    let mut frames = enc.encode(&pcm);
    frames.extend(enc.flush());
    let mut stream = enc.tag_frame();
    for f in &frames {
        stream.extend_from_slice(f);
    }
    let mut d = Decoder::with_options(DecoderOptions { strict: true, check_crc: true, ..Default::default() });
    let mut out = d.decode(&stream).unwrap();
    out.extend(d.flush().unwrap());
    let dec: Vec<f32> = out.iter().flat_map(|f| f.samples.iter().copied()).collect();
    assert_eq!(dec.len(), pcm.len());
    let (mut e, mut s) = (0.0f64, 0.0f64);
    for (a, b) in dec.iter().zip(&pcm) {
        e += f64::from(a - b).powi(2);
        s += f64::from(*b).powi(2);
    }
    eprintln!("SNR {:.1} dB, {} frames, avg {:.0}", 10.0 * (s / e).log10(), frames.len(), enc.average_bitrate());
}

#[test]
fn find_lag() {
    let cfg = EncoderConfig { channels: 1, bitrate: BitrateMode::Cbr(320_000), ..Default::default() };
    let mut enc = Encoder::new(cfg).unwrap();
    let pcm = tone(44_100, 1, 44_100);
    let mut frames = enc.encode(&pcm);
    frames.extend(enc.flush());
    let stream: Vec<u8> = frames.concat();
    let mut d = Decoder::with_options(DecoderOptions { trim_gapless: false, ..Default::default() });
    let mut out = d.decode(&stream).unwrap();
    out.extend(d.flush().unwrap());
    let dec: Vec<f32> = out.iter().flat_map(|f| f.samples.iter().copied()).collect();
    let mut best = (0usize, f64::MAX);
    for lag in 900..1300 {
        let e: f64 = (5000..30000).map(|i| f64::from(dec[i + lag] - pcm[i]).powi(2)).sum();
        if e < best.1 {
            best = (lag, e);
        }
    }
    let s: f64 = (5000..30000).map(|i| f64::from(pcm[i]).powi(2)).sum();
    eprintln!("best lag {} SNR {:.1}", best.0, 10.0 * (s / best.1).log10());
}

