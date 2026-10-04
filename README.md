# rivet-mp3

[![CI](https://github.com/safewords/rivet-mp3/actions/workflows/ci.yml/badge.svg)](https://github.com/safewords/rivet-mp3/actions/workflows/ci.yml)

An **MPEG audio decoder** (Layers I, II and III; MPEG-1, MPEG-2 lower
sampling frequencies and MPEG-2.5) and an **MP3 (Layer III) encoder** in
Rust: no C, no system libraries, no build script, nothing to install on a
build host. Written from ISO/IEC 11172-3 and ISO/IEC 13818-3, not translated
from any other implementation. The decoder meets ISO's *full accuracy*
criterion on every Layer I, II and III conformance sequence it was checked
on — 64 of 64 (the figures are [below](#how-it-is-checked)).

Written for the **[rivet](https://github.com/safewords/rivet)**
transcoder, where it replaces minimp3 (decoding) and LAME (encoding). Usable
on its own by anything that has MPEG audio bytes and wants PCM back, or PCM
and wants MP3.

Published as `rivet-mp3`; **imported as `mp3`** (`use mp3::…`). One
dependency (`thiserror`), no features, no build script.

```toml
[dependencies]
mp3 = { package = "rivet-mp3", git = "https://github.com/safewords/rivet-mp3", branch = "develop" }
```

## What it decodes

| | supported |
|---|---|
| **Layers** | I, II, III |
| **Versions and rates** | MPEG-1 (32, 44.1, 48 kHz), MPEG-2 LSF (16, 22.05, 24 kHz), MPEG-2.5 (8, 11.025, 12 kHz) |
| **Bit rates** | every table rate; free format (the frame length found from the next sync) |
| **Channels** | mono, stereo, dual channel, joint stereo: intensity stereo (all layers; MPEG-1 and LSF rules) and mid/side |
| **Layer III** | long, start, short and stop blocks, mixed blocks, scfsi, preflag, both count1 tables, every Huffman table, the bit reservoir |
| **Integrity** | CRC words checked on request (`check_crc`); a strict mode that refuses any irregularity, for validating encoders |
| **Tags** | ID3v2 (skipped), ID3v1 / APE at the end (skipped), Xing / Info with the LAME extension and VBRI (parsed; the tag frame is not played) |
| **Gapless** | the LAME tag's encoder delay and padding (VBRI's delay) trim the output to the encoder's input; on by default |
| **Not decoded** | the MPEG-2 multichannel extension (its MPEG-1 compatible stereo is decoded); emphasis is reported in the header, not undone |

Output is interleaved `f32` at ±1.0 full scale, one or two channels, not
clipped. Bytes may arrive in any chunking; sync is confirmed against the
next frame's header before it is trusted, garbage between frames is
skipped, and damage is concealed with silence (`Frame::concealed`) so the
timeline holds — or reported, in strict mode.

## What it encodes

Layer III at all nine sampling frequencies (`encode::coding_rate` names the
one to resample other input to), mono or stereo:

- **constant bit rate** at any rate of the version's table, or **variable
  bit rate** by quality 0 (best) – 9;
- joint stereo with **mid/side** chosen frame by frame (or plain stereo);
- the polyphase analysis filterbank and MDCT of 11172-3 Annex C with
  **long, start, short and stop blocks** switched by a transient detector;
- a **psychoacoustic model** on the MDCT spectrum (band energies,
  spectral-flatness tonality, the Annex D spreading function, tonal and
  noise masking offsets, the threshold in quiet, pre-echo control);
- quantisation against the mask with one noise-to-mask offset per frame
  found by bisection against the bit budget; scalefactors,
  scalefac_scale, preflag and subblock gains chosen per granule; Huffman
  regions and tables chosen for the fewest bits;
- the **bit reservoir**, main_data_begin kept within its field and the
  decoder's 7680-bit buffer;
- optional CRC words;
- a **Xing (VBR) or Info (CBR) tag frame** with a LAME-style extension:
  frame and byte counts, seek table, encoder delay (528) and padding, music
  and tag CRCs — what a gapless decoder needs to give back exactly the
  input.

Left out, all optional for an encoder: intensity stereo, mixed blocks,
scfsi, free format, Layer I and II encoding.

## Speed

On a Ryzen 9 9950X (Windows, a shared machine, best of five), in multiples
of real time, for 60 s of a 16-bit stereo 44.1 kHz album track;
`cargo run --release --example bench -- <pcm.raw>` measures it on any
16-bit stereo 44.1 kHz raw PCM.

| | before | now |
|---|---|---|
| decode 128 kb/s CBR | 264 | 670 |
| decode 320 kb/s CBR | 198 | 647 |
| encode 128 kb/s CBR (one thread / threaded) | 6.7 | 13.3 / 20.8 |
| encode 320 kb/s CBR | 4.9 | 11.5 / 19.2 |
| encode VBR quality 2 | 22 | 29 / 39 |

Decoding: the Huffman codes resolve their first 10 bits in one table
lookup; |is|^(4/3) comes from a table; subbands above the last coded line
skip the IMDCT (their output, zeros' transform, is kept); the IMDCT and the
synthesis filterbank's matrixing and windowing compute their outputs side
by side, which vectorises (x86-64 builds them for the baseline and for
AVX2, picked at run time; aarch64 uses NEON); and `Decoder` no longer moves
its whole input down once per frame, which made decoding a large buffer
at once quadratic. Encoding: the Huffman table and region search prices
ranges from per-table prefix sums instead of pricing every pair in every
table; the 3/4 powers and step factors are computed once; and the frames
an `encode` call completes are analysed and first quantised in parallel,
with each rate-control trial quantising its granules and channels on
separate threads (`Encoder::set_threads`; 1 keeps everything on the
caller's thread).

**Same output everywhere.** No fused multiply-add, and every vectorised
sum adds its terms in the order of the plain loop, so the decoder's output
and the encoder's stream are the same to the bit on every CPU, on every
code path (the `force-scalar` feature compiles the run-time selection
out; CI tests both ways, and on arm64), at every thread count — and the
same as before this work: the encoded streams' hashes and the decoded
output of every ISO conformance stream are unchanged, and unit tests hold
each rewritten kernel to the loop it replaced.

## How it is checked

- **ISO conformance** (`tests/conformance.rs`): the MPEG-1/2 audio
  sequences of ISO/IEC 14496-26's public conformance package (descended
  from 11172-4 and 13818-4) against ISO's reference waveforms, by 11172-4's
  criterion. The data is ISO's and not in the repository:
  `python tools/fetch_conformance.py DIR` fetches it (about 45 MB, by HTTP
  range requests from the 10 GB package), then
  `MP3_CONFORMANCE_DIR=DIR cargo test --release --test conformance -- --nocapture`.
- **Encoder round trips** (`src/encode/tests.rs`): 34 configurations —
  32–320 kbit/s CBR, VBR quality 0–9, joint and plain stereo, mono, CRC,
  32 / 44.1 / 48 kHz, the LSF and MPEG-2.5 rates — each decoded strictly
  (CRCs, exact reservoir accounting, every Huffman region ending at its
  part2_3_length) with gapless trimming to exactly the input's length, the
  tag's fields checked, SNR and per-band noise-to-mask ratio measured
  (`MP3_ENCODER_REPORT=1` prints the table).
- **Robustness** (`tests/robustness.rs`): garbage, flipped bits, cut and
  spliced streams in any chunking — concealed, never a panic — ID3 tags,
  free format, VBRI, gapless trimming.
- **Tables and transforms**: every Huffman table a complete prefix code
  whose words decode to their own entries; the synthesis window's printed
  values and symmetry; the synthesis filterbank against a literal float
  rendering of the standard's flow chart; the IMDCT against its formula;
  MDCT/IMDCT perfect reconstruction across every block-type transition;
  analysis followed by synthesis.

| | result |
|---|---|
| ISO conformance | 64 of 64 references at full accuracy (9 Layer I, 16 Layer II, 39 Layer III; MPEG-1, LSF and MPEG-2.5); worst RMS error 1.3e-7 of full scale against a bound of 8.8e-6, worst sample 1.3e-6 against 6.1e-5 |
| Encoder, 44.1 kHz joint stereo | mean noise-to-mask +1.4 dB at 128 kbit/s, −1.4 dB at 160, −4.0 dB at 192, −12.0 dB at 320, on a demanding test programme; VBR quality 0–9 spans 246–74 kbit/s |

[docs/VALIDATION.md](docs/VALIDATION.md) lists every sequence and every
configuration; [docs/PROVENANCE.md](docs/PROVENANCE.md) where each part came
from, and the points where the standards leave a decoder a choice that
ISO's reference output settled.

## Provenance and licensing

Written from the standards' text; **no MPEG audio implementation's source
was read** — not minimp3, LAME, mpg123, libmad, FFmpeg's, Helix, dist10 or
any other — and none was run. The normative tables are the standards'
data, verified as above.

**Patents.** The principal MP3 patents are generally regarded as expired.
Nothing here is a licence to any patent, and the authors make no claim
about whether anyone needs one.

## Using it

```rust
// Decoding: bytes in any chunking (an .mp3 file, or an MP4 / Matroska /
// AVI stream's packets).
let mut dec = mp3::Decoder::new();
for packet in packets {
    for frame in dec.decode(packet)? {
        // frame.samples: interleaved f32; frame.channels; frame.sample_rate
    }
}
let tail = dec.flush()?;
// dec.gapless(): encoder delay, padding and length, when a tag gave them.

// Encoding: interleaved f32 at a Layer III rate, one packet per frame.
let mut enc = mp3::Encoder::new(mp3::EncoderConfig {
    sample_rate: 44_100,
    channels: 2,
    bitrate: mp3::BitrateMode::Cbr(192_000),
    ..Default::default()
})?;
let mut packets = enc.encode(&pcm);
packets.extend(enc.flush());
// A bare .mp3 file: the tag frame first, then the packets.
let tag = enc.tag_frame();
// In a container: skip enc.delay() + mp3::xing::DECODER_DELAY samples.
```

## License

Open Encoding Attribution License v1.0 — a source-available (not OSI open-source)
license, royalty-free, with a commercial-attribution requirement. See
[LICENSE.md](LICENSE.md) and [NOTICE](NOTICE).
