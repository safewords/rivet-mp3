# Validation

Figures measured 2026-10-02 on Windows x86-64, Rust 1.99, release build.

## Decoder: ISO conformance sequences

The MPEG-1/2 audio sequences of ISO/IEC 14496-26 (2nd edition) — the
Layer I, II and III streams that descend from ISO/IEC 11172-4 and 13818-4 —
and ISO's reference waveforms (24-bit WAV), fetched by
`tools/fetch_conformance.py` and compared by `tests/conformance.rs`.

Criterion (ISO/IEC 11172-4): *full accuracy* is an RMS difference below
2^-15 / sqrt(12) = 8.81e-6 of full scale with no sample off by more than
2^-14 = 6.10e-5; *limited accuracy* an RMS below 2^-11 / sqrt(12) = 1.41e-4.

Result: **64 of 64 references at full accuracy**; the largest RMS difference
is 1.30e-7 (68 times below the bound), the largest single-sample difference
1.31e-6 (47 times below).

How the comparison treats three things the standard leaves to the decoder:

- The reference is PCM clipped to [-1, 1); the decoder's float output is not
  clipped, so it is clipped the same way before comparing (l3_10204 peaks at
  1.03).
- A Layer III frame whose main data begins before the stream does (l3_10204
  starts mid-stream) produces no output from the reference decoder; this
  decoder emits silence for it, flagged `Frame::concealed`, to keep the
  timeline. Concealed frames are left out of the comparison.
- l3_10201 changes between mono and stereo; ISO gives one reference per
  segment (`_0`, `_1`, `_2`), and each is compared with its segment. Where a
  segment has two reference variants (`_1_AA`, `_1_noAA`), the closer counts.

l2_test27 has no reference waveform in the package; l3_10501 is a
multichannel (MPEG-2 BC) stream whose MPEG-1 compatible stereo base is
compared (`l3_10501_f00`); the multichannel extension is not decoded.

| reference | stream | RMS error | max \|error\| | verdict |
|---|---|---|---|---|
| l1_fl1 | Mpeg1 L1 32000 Hz Stereo 384 kb/s CRC | 3.445e-8 | 6.706e-8 | full |
| l1_fl2 | Mpeg1 L1 44100 Hz JointStereo 384 kb/s CRC | 3.457e-8 | 6.706e-8 | full |
| l1_fl3 | Mpeg1 L1 48000 Hz JointStereo 384 kb/s CRC | 3.453e-8 | 6.706e-8 | full |
| l1_fl4 | Mpeg1 L1 32000 Hz Mono 32 kb/s | 3.440e-8 | 6.706e-8 | full |
| l1_fl5 | Mpeg1 L1 48000 Hz DualChannel 448 kb/s CRC | 3.955e-8 | 1.192e-7 | full |
| l1_fl6 | Mpeg1 L1 44100 Hz Stereo 384 kb/s CRC | 3.714e-8 | 1.192e-7 | full |
| l1_fl7 | Mpeg1 L1 44100 Hz Stereo 384 kb/s CRC | 3.454e-8 | 8.941e-8 | full |
| l1_fl8 | Mpeg1 L1 44100 Hz Stereo 384 kb/s | 3.447e-8 | 6.706e-8 | full |
| l2_fl10 | Mpeg1 L2 32000 Hz Stereo 192 kb/s CRC | 3.444e-8 | 6.706e-8 | full |
| l2_fl11 | Mpeg1 L2 44100 Hz Stereo 192 kb/s CRC | 3.452e-8 | 6.706e-8 | full |
| l2_fl12 | Mpeg1 L2 48000 Hz Stereo 192 kb/s CRC | 3.443e-8 | 6.706e-8 | full |
| l2_fl13 | Mpeg1 L2 32000 Hz Mono 32 kb/s | 3.445e-8 | 6.706e-8 | full |
| l2_fl14 | Mpeg1 L2 48000 Hz DualChannel 384 kb/s CRC | 3.965e-8 | 7.451e-8 | full |
| l2_fl15 | Mpeg1 L2 48000 Hz Stereo 384 kb/s CRC | 3.563e-8 | 1.192e-7 | full |
| l2_fl16 | Mpeg1 L2 48000 Hz Stereo 256 kb/s CRC | 3.483e-8 | 1.192e-7 | full |
| l2_test24 | Mpeg2 L2 16000 Hz JointStereo 96 kb/s CRC | 3.715e-8 | 1.192e-7 | full |
| l2_test25 | Mpeg2 L2 22050 Hz JointStereo 96 kb/s CRC | 3.743e-8 | 1.192e-7 | full |
| l2_test26 | Mpeg2 L2 24000 Hz JointStereo 96 kb/s CRC | 3.716e-8 | 1.192e-7 | full |
| l2_test28 | Mpeg2 L2 24000 Hz JointStereo 96 kb/s CRC | 3.716e-8 | 1.192e-7 | full |
| l2_test29 | Mpeg2 L2 24000 Hz JointStereo 160 kb/s CRC | 3.739e-8 | 1.192e-7 | full |
| l2_test30 | Mpeg2 L2 24000 Hz Mono 96 kb/s CRC | 3.746e-8 | 1.192e-7 | full |
| l2_test31 | Mpeg2 L2 24000 Hz DualChannel 96 kb/s CRC | 3.755e-8 | 1.192e-7 | full |
| l2_test32 | Mpeg2 L2 24000 Hz Stereo 128 kb/s | 3.958e-8 | 1.192e-7 | full |
| l2_test33 | Mpeg2 L1 22050 Hz JointStereo 128 kb/s CRC | 3.661e-8 | 8.941e-8 | full |
| l2_test34 | Mpeg2 L2 24000 Hz Stereo 160 kb/s CRC | 3.442e-8 | 6.706e-8 | full |
| l3_10101 | Mpeg1 L3 32000 Hz Mono 32 kb/s (lengths 172800 vs 171648) | 3.554e-8 | 8.941e-8 | full |
| l3_10102 | Mpeg1 L3 44100 Hz Mono 32 kb/s (lengths 472320 vs 471168) | 3.554e-8 | 8.941e-8 | full |
| l3_10103 | Mpeg1 L3 48000 Hz Mono 32 kb/s (lengths 172800 vs 171648) | 3.554e-8 | 8.941e-8 | full |
| l3_10104 | Mpeg1 L3 44100 Hz Mono 64 kb/s (lengths 135936 vs 134784) | 3.261e-8 | 1.043e-7 | full |
| l3_10105 | Mpeg1 L3 44100 Hz Mono 64 kb/s (lengths 73728 vs 72576) | 3.456e-8 | 6.519e-8 | full |
| l3_10106 | Mpeg1 L3 44100 Hz Mono 64 kb/s | 3.397e-8 | 8.941e-8 | full |
| l3_10107 | Mpeg1 L3 48000 Hz Mono 160 kb/s | 3.394e-8 | 6.752e-8 | full |
| l3_10201_0 | Mpeg1 L3 44100 Hz Mono 128 kb/s | 3.403e-8 | 6.519e-8 | full |
| l3_10201_1_AA | Mpeg1 L3 44100 Hz DualChannel 128 kb/s | 3.444e-8 | 7.451e-8 | full |
| l3_10201_2 | Mpeg1 L3 44100 Hz Mono 128 kb/s | 3.331e-8 | 6.414e-8 | full |
| l3_10202 | Mpeg1 L3 44100 Hz Stereo 128 kb/s | 3.936e-8 | 1.490e-7 | full |
| l3_10203 | Mpeg1 L3 44100 Hz JointStereo 128 kb/s | 3.262e-8 | 7.451e-8 | full |
| l3_10204 | Mpeg1 L3 44100 Hz JointStereo 128 kb/s | 4.895e-8 | 1.788e-7 | full |
| l3_10205 | Mpeg1 L3 44100 Hz JointStereo 320 kb/s | 2.475e-8 | 6.007e-8 | full |
| l3_10206 | Mpeg1 L3 48000 Hz Mono 64 kb/s (lengths 249984 vs 248832) | 3.454e-8 | 6.706e-8 | full |
| l3_10207 | Mpeg1 L3 48000 Hz JointStereo 160 kb/s | 3.397e-8 | 8.941e-8 | full |
| l3_10501_f00 | Mpeg1 L3 48000 Hz Mono 160 kb/s | 3.394e-8 | 6.752e-8 | full |
| l3_20140 | Mpeg2 L3 16000 Hz Mono 8 kb/s | 3.549e-8 | 8.941e-8 | full |
| l3_20141 | Mpeg2 L3 22050 Hz Mono 8 kb/s | 3.560e-8 | 8.941e-8 | full |
| l3_20142 | Mpeg2 L3 24000 Hz Mono 8 kb/s | 3.531e-8 | 8.941e-8 | full |
| l3_20143 | Mpeg2 L3 24000 Hz Mono 128 kb/s | 3.456e-8 | 7.078e-8 | full |
| l3_20244 | Mpeg2 L3 22050 Hz JointStereo 96 kb/s | 3.496e-8 | 8.941e-8 | full |
| l3_20245 | Mpeg2 L3 22050 Hz JointStereo 160 kb/s | 3.418e-8 | 8.941e-8 | full |
| l3_20246 | Mpeg2 L3 22050 Hz JointStereo 160 kb/s | 3.361e-8 | 9.686e-8 | full |
| l3_20247 | Mpeg2 L3 22050 Hz JointStereo 96 kb/s | 2.871e-8 | 7.451e-8 | full |
| l3_20248 | Mpeg2 L3 22050 Hz JointStereo 96 kb/s | 3.445e-8 | 6.892e-8 | full |
| l3_25101 | Mpeg25 L3 12000 Hz Mono 64 kb/s | 5.559e-8 | 2.682e-7 | full |
| l3_25102 | Mpeg25 L3 8000 Hz Mono 8 kb/s | 6.208e-8 | 2.682e-7 | full |
| l3_25103 | Mpeg25 L3 11025 Hz Mono 8 kb/s | 7.095e-8 | 4.023e-7 | full |
| l3_25104 | Mpeg25 L3 12000 Hz Mono 8 kb/s | 6.817e-8 | 4.768e-7 | full |
| l3_25201 | Mpeg25 L3 11025 Hz JointStereo 64 kb/s | 6.021e-8 | 3.153e-7 | full |
| l3_25202 | Mpeg25 L3 11025 Hz JointStereo 64 kb/s | 3.946e-8 | 2.235e-7 | full |
| l3_25203 | Mpeg25 L3 11025 Hz JointStereo 64 kb/s | 5.502e-8 | 2.356e-7 | full |
| l3_25204 | Mpeg25 L3 8000 Hz JointStereo 64 kb/s | 6.175e-8 | 4.359e-7 | full |
| l3_25205 | Mpeg25 L3 11025 Hz JointStereo 64 kb/s | 6.453e-8 | 4.321e-7 | full |
| l3_25206 | Mpeg25 L3 12000 Hz JointStereo 64 kb/s | 6.581e-8 | 3.390e-7 | full |
| l3_25207 | Mpeg25 L3 8000 Hz JointStereo 32 kb/s | 1.196e-7 | 1.267e-6 | full |
| l3_25208 | Mpeg25 L3 11025 Hz JointStereo 32 kb/s | 1.297e-7 | 1.311e-6 | full |
| l3_25209 | Mpeg25 L3 12000 Hz JointStereo 32 kb/s | 1.250e-7 | 1.192e-6 | full |

## Encoder: round trips

`src/encode/tests.rs` encodes a three-second programme (plucked harmonic
notes with sharp attacks, a sustained chord, coloured noise, castanet-like
clicks, partly correlated stereo) in each configuration below, decodes it
with this crate's decoder in strict mode — every CRC checked, the bit
reservoir accounted exactly (no frame's main_data_begin reaching into bits
an earlier frame used), every Huffman region ending exactly at its
part2_3_length — with gapless trimming, and checks that the output has
exactly the input's length, that the Info/Xing tag's frame and byte counts,
delay, padding, tag CRC and music CRC are right, and that CBR frames all
carry the asked-for rate.

SNR is the worst channel's, over the full band (so it counts the low-pass
too). NMR is noise-to-mask ratio per long scalefactor band and granule,
noise and mask both computed by the encoder's own psychoacoustic model on
the encoder's long-block spectrum of the original and of the error, over
the bands below the configuration's low-pass; ">mask" is the share of
(band, granule) pairs whose noise exceeds the mask. "short" counts granules
coded with short blocks, "M/S" frames coded mid/side (of 116 frames at
44.1 kHz).

| configuration | kbit/s | SNR dB | NMR mean dB | NMR p95 dB | bands over mask | short granules | M/S frames |
|---|---|---|---|---|---|---|---|
| 44100 Hz joint  CBR  32k | 32.0 | 2.5 | 4.1 | 18.5 | 70.3% | 28 | 87 |
| 44100 Hz joint  CBR  64k | 64.0 | 7.1 | 4.6 | 12.7 | 77.2% | 28 | 43 |
| 44100 Hz joint  CBR  96k | 96.0 | 10.0 | 3.3 | 11.4 | 74.5% | 28 | 33 |
| 44100 Hz joint  CBR 128k | 128.0 | 10.9 | 1.4 | 8.5 | 71.1% | 28 | 32 |
| 44100 Hz joint  CBR 160k | 160.0 | 13.6 | -1.4 | 5.7 | 50.2% | 28 | 32 |
| 44100 Hz joint  CBR 192k | 192.0 | 15.5 | -4.0 | 3.7 | 16.4% | 28 | 32 |
| 44100 Hz joint  CBR 256k | 256.0 | 19.0 | -8.4 | -1.1 | 4.3% | 28 | 31 |
| 44100 Hz joint  CBR 320k | 320.0 | 21.0 | -12.0 | -5.3 | 2.4% | 28 | 31 |
| 44100 Hz stereo CBR 128k | 128.0 | 9.7 | 1.9 | 9.0 | 74.0% | 28 | 0 |
| 44100 Hz joint  CBR 128k CRC | 128.0 | 10.6 | 1.5 | 8.5 | 71.8% | 28 | 32 |
| 44100 Hz mono   CBR  64k | 64.0 | 10.5 | 1.8 | 9.8 | 73.6% | 14 | 0 |
| 44100 Hz mono   CBR 128k | 128.0 | 17.8 | -7.8 | 0.7 | 5.8% | 14 | 0 |
| 48000 Hz joint  CBR 128k | 128.0 | 11.0 | 0.8 | 7.7 | 67.9% | 28 | 30 |
| 48000 Hz joint  CBR 320k | 320.0 | 20.7 | -12.9 | -6.4 | 1.7% | 28 | 30 |
| 32000 Hz joint  CBR  64k | 64.0 | 7.1 | 2.9 | 16.6 | 68.5% | 28 | 38 |
| 32000 Hz joint  CBR 128k | 128.0 | 13.7 | -1.3 | 7.8 | 60.4% | 28 | 22 |
| 44100 Hz joint  VBR q0 | 245.8 | 18.8 | -7.9 | -0.4 | 4.7% | 28 | 32 |
| 44100 Hz joint  VBR q2 | 201.7 | 16.3 | -4.8 | 2.5 | 11.0% | 28 | 32 |
| 44100 Hz joint  VBR q4 | 162.8 | 13.8 | -1.6 | 5.5 | 47.8% | 28 | 32 |
| 44100 Hz joint  VBR q6 | 129.0 | 11.3 | 1.2 | 8.4 | 69.8% | 28 | 32 |
| 44100 Hz joint  VBR q9 | 73.8 | 9.8 | 4.7 | 12.2 | 77.1% | 28 | 36 |
| 48000 Hz joint  VBR q4 | 155.2 | 13.4 | -1.5 | 5.2 | 47.6% | 28 | 30 |
| 32000 Hz mono   VBR q4 | 74.8 | 13.6 | -2.4 | 7.4 | 49.5% | 14 | 0 |
| 22050 Hz joint  CBR  32k | 32.0 | 4.6 | 2.1 | 19.4 | 63.3% | 26 | 88 |
| 22050 Hz joint  CBR  64k | 64.0 | 9.8 | 3.1 | 13.0 | 76.4% | 26 | 54 |
| 22050 Hz joint  CBR 160k | 160.0 | 20.4 | -7.2 | 9.1 | 8.7% | 26 | 54 |
| 24000 Hz joint  CBR  96k | 96.0 | 12.6 | 0.9 | 11.9 | 65.4% | 28 | 56 |
| 16000 Hz mono   CBR  32k CRC | 32.0 | 12.3 | 0.3 | 13.2 | 58.8% | 11 | 0 |
| 22050 Hz joint  VBR q4 | 96.5 | 13.4 | -0.6 | 11.9 | 47.9% | 26 | 54 |
| 11025 Hz joint  CBR  32k | 32.0 | 11.8 | 0.1 | 12.3 | 59.0% | 10 | 54 |
| 12000 Hz mono   CBR  24k | 24.0 | 12.0 | 0.3 | 12.2 | 53.5% | 10 | 0 |
| 8000 Hz mono   CBR   8k | 8.0 | 8.1 | 2.8 | 17.0 | 62.9% | 6 | 0 |
| 8000 Hz mono   CBR  16k | 16.0 | 11.2 | 0.7 | 12.9 | 58.5% | 6 | 0 |
| 11025 Hz joint  VBR q4 | 46.3 | 14.8 | -2.9 | 9.1 | 35.1% | 10 | 54 |

The rate control spends a constant-rate frame's bits by moving one
noise-to-mask offset per frame, so the NMR is near-uniform across bands: at
128 kbit/s this demanding programme (broadband noise in every band) sits
about 1.4 dB above its mask on average; from 160 kbit/s joint stereo, and
from VBR quality 4, the mean noise is below the mask. VBR quality 0–9 spans
246–74 kbit/s on this programme.

Further encoder tests: the codec delay is the tag's (an impulse returns
528 + 529 samples late); an attack after near-silence switches to short
blocks and the noise before it stays more than 30 dB below the attack;
correlated stereo is coded mid/side and independent channels mostly not;
the reservoir is used and its accounting holds frame by frame at 64 and
256 kbit/s and VBR; 95 % or more of the bands that need coding meet their
allowance to within 1 dB.
