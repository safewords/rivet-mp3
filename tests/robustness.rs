//! The decoder on input that is not a clean stream: garbage, damage,
//! truncation, tags around the frames, free format, any chunking. Errors
//! are allowed (in strict mode); panics are not.

use mp3::{
    BitrateMode, Decoder, DecoderOptions, Encoder, EncoderConfig, FrameDecoder, FrameHeader,
};

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }
}

fn signal(n: usize, nch: usize) -> Vec<f32> {
    (0..n * nch)
        .map(|i| {
            let t = (i / nch) as f32;
            0.3 * (t * 0.031).sin() + 0.2 * (t * 0.17 * (1.0 + (i % nch) as f32)).sin()
        })
        .collect()
}

/// Encoded frames (no tag frame).
fn frames(cfg: EncoderConfig, seconds: f32) -> Vec<Vec<u8>> {
    let mut enc = Encoder::new(cfg).unwrap();
    let pcm = signal(
        (seconds * cfg.sample_rate as f32) as usize,
        usize::from(cfg.channels),
    );
    let mut f = enc.encode(&pcm);
    f.extend(enc.flush());
    f
}

fn decode_all(bytes: &[u8], opts: DecoderOptions, chunk: usize) -> mp3::Result<Vec<mp3::Frame>> {
    let mut d = Decoder::with_options(opts);
    let mut out = Vec::new();
    for c in bytes.chunks(chunk.max(1)) {
        out.extend(d.decode(c)?);
    }
    out.extend(d.flush()?);
    Ok(out)
}

fn pcm(frames: &[mp3::Frame]) -> Vec<f32> {
    frames
        .iter()
        .flat_map(|f| f.samples.iter().copied())
        .collect()
}

#[test]
fn garbage_never_panics() {
    let mut rng = Rng(1);
    for len in [0usize, 1, 3, 4, 100, 4096, 20_000] {
        let bytes: Vec<u8> = (0..len).map(|_| rng.next() as u8).collect();
        for strict in [false, true] {
            let _ = decode_all(
                &bytes,
                DecoderOptions {
                    strict,
                    ..Default::default()
                },
                333,
            );
        }
        // Garbage that is all sync words.
        let syncs: Vec<u8> = (0..len)
            .map(|i| {
                if i % 2 == 0 {
                    0xFF
                } else {
                    0xE0 | rng.next() as u8
                }
            })
            .collect();
        let _ = decode_all(&syncs, DecoderOptions::default(), 50);
        let _ = FrameDecoder::new().decode_frame(&syncs);
    }
}

#[test]
fn damaged_streams_never_panic() {
    let mut rng = Rng(2);
    for cfg in [
        EncoderConfig::default(),
        EncoderConfig {
            sample_rate: 22_050,
            bitrate: BitrateMode::Cbr(64_000),
            ..Default::default()
        },
        EncoderConfig {
            sample_rate: 8_000,
            channels: 1,
            bitrate: BitrateMode::Vbr(5),
            ..Default::default()
        },
    ] {
        let clean: Vec<u8> = frames(cfg, 1.0).concat();
        for round in 0..60 {
            let mut b = clean.clone();
            match round % 4 {
                0 => {
                    for _ in 0..20 {
                        let i = rng.next() as usize % b.len();
                        b[i] ^= 1 << (rng.next() % 8);
                    }
                }
                1 => b.truncate(rng.next() as usize % b.len()),
                2 => {
                    let at = rng.next() as usize % b.len();
                    let junk: Vec<u8> = (0..rng.next() % 500).map(|_| rng.next() as u8).collect();
                    b.splice(at..at, junk);
                }
                _ => {
                    for _ in 0..200 {
                        let i = rng.next() as usize % b.len();
                        b[i] = rng.next() as u8;
                    }
                }
            }
            let out = decode_all(
                &b,
                DecoderOptions::default(),
                1 + rng.next() as usize % 2000,
            );
            assert!(
                out.is_ok(),
                "the lenient decoder conceals damage: {:?}",
                out.err()
            );
            let _ = decode_all(
                &b,
                DecoderOptions {
                    strict: true,
                    check_crc: true,
                    ..Default::default()
                },
                700,
            );
        }
    }
}

#[test]
fn chunking_does_not_change_the_output() {
    let bytes: Vec<u8> = frames(EncoderConfig::default(), 1.0).concat();
    let whole = pcm(&decode_all(&bytes, DecoderOptions::default(), bytes.len()).unwrap());
    for chunk in [1, 7, 418, 1000] {
        assert_eq!(
            pcm(&decode_all(&bytes, DecoderOptions::default(), chunk).unwrap()),
            whole,
            "chunk {chunk}"
        );
    }
}

#[test]
fn tags_around_the_frames_are_skipped() {
    let body: Vec<u8> = frames(EncoderConfig::default(), 0.5).concat();
    let plain = pcm(&decode_all(&body, DecoderOptions::default(), 4096).unwrap());
    let mut tagged = b"ID3\x04\x00\x00\x00\x00\x02\x01".to_vec(); // ID3v2, 257 bytes of payload
    tagged.extend(std::iter::repeat_n(0xFFu8, 257)); // a payload full of false syncs
    tagged.extend_from_slice(&body);
    let mut v1 = b"TAG".to_vec();
    v1.resize(128, b' ');
    tagged.extend_from_slice(&v1);
    let mut d = Decoder::new();
    let mut out = d.decode(&tagged).unwrap();
    out.extend(d.flush().unwrap());
    assert_eq!(pcm(&out), plain);
    assert_eq!(d.skipped_bytes(), 10 + 257 + 128);
}

/// Free format: the same frames with bitrate_index 0 decode identically,
/// the frame length found from the distance between syncs.
#[test]
fn free_format() {
    let cfg = EncoderConfig {
        sample_rate: 48_000,
        bitrate: BitrateMode::Cbr(160_000),
        ..Default::default()
    };
    let f = frames(cfg, 1.0);
    let fixed: Vec<u8> = f.concat();
    let free: Vec<u8> = f
        .iter()
        .flat_map(|fr| {
            let mut fr = fr.clone();
            fr[2] &= 0x0F; // bitrate_index 0
            fr
        })
        .collect();
    let h = FrameHeader::parse(&free).unwrap();
    assert!(h.is_free_format());
    let a = decode_all(&fixed, DecoderOptions::default(), 999).unwrap();
    let b = decode_all(
        &free,
        DecoderOptions {
            strict: true,
            ..Default::default()
        },
        999,
    )
    .unwrap();
    assert_eq!(pcm(&a), pcm(&b));
    assert!(
        b.iter().all(|fr| fr.bitrate == 160_000),
        "measured free-format bit rate"
    );
}

/// A VBRI header (Fraunhofer's) is recognised, and its frame not played.
#[test]
fn vbri_header() {
    let f = frames(
        EncoderConfig {
            bitrate: BitrateMode::Cbr(128_000),
            ..Default::default()
        },
        0.5,
    );
    let mut tag = vec![0u8; f[0].len()];
    tag[..4].copy_from_slice(&f[0][..4]);
    tag[2] &= !0x02; // no padding
    tag.truncate(FrameHeader::parse(&tag).unwrap().frame_len().unwrap());
    let v = &mut tag[36..];
    v[..4].copy_from_slice(b"VBRI");
    v[4..6].copy_from_slice(&1u16.to_be_bytes());
    v[6..8].copy_from_slice(&1105u16.to_be_bytes());
    v[8..10].copy_from_slice(&75u16.to_be_bytes());
    v[14..18].copy_from_slice(&(f.len() as u32).to_be_bytes());
    let mut stream = tag;
    stream.extend(f.concat());
    let mut d = Decoder::with_options(DecoderOptions {
        trim_gapless: false,
        ..Default::default()
    });
    let mut out = d.decode(&stream).unwrap();
    out.extend(d.flush().unwrap());
    match d.info() {
        Some(mp3::xing::InfoHeader::Vbri(v)) => {
            assert_eq!(v.delay, 1105);
            assert_eq!(v.frames, f.len() as u32);
        }
        other => panic!("{other:?}"),
    }
    assert_eq!(out.len(), f.len(), "the VBRI frame is not audio");
}

/// The gapless figures and the trimming they drive.
#[test]
fn gapless_trimming() {
    let cfg = EncoderConfig {
        channels: 1,
        bitrate: BitrateMode::Vbr(4),
        ..Default::default()
    };
    let n = 12_345;
    let input = signal(n, 1);
    let mut enc = Encoder::new(cfg).unwrap();
    let mut f = enc.encode(&input);
    f.extend(enc.flush());
    let mut stream = enc.tag_frame();
    stream.extend(f.concat());
    let trimmed = pcm(&decode_all(&stream, DecoderOptions::default(), 512).unwrap());
    assert_eq!(trimmed.len(), n);
    let untrimmed = pcm(&decode_all(
        &stream,
        DecoderOptions {
            trim_gapless: false,
            ..Default::default()
        },
        512,
    )
    .unwrap());
    assert_eq!(untrimmed.len(), f.len() * 1152);
    let skip = (enc.delay() + mp3::xing::DECODER_DELAY) as usize;
    assert_eq!(&untrimmed[skip..skip + n], &trimmed[..]);
    let mut d = Decoder::new();
    d.decode(&stream).unwrap();
    let g = d.gapless().unwrap();
    assert_eq!(
        (g.encoder_delay, g.padding, g.length),
        (enc.delay(), enc.padding(), Some(n as u64))
    );
}
