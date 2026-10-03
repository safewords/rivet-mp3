# Provenance

Where every part of this crate came from. The short version: the code is
this repository's own, written from ISO/IEC 11172-3 and ISO/IEC 13818-3 and
the public descriptions of the MPEG-2.5 extension and of the Xing, LAME and
VBRI tags; the normative tables are the standards' data, reproduced from the
author's knowledge of them and verified mechanically and against ISO's
conformance data; no MPEG audio implementation was read or run.

## Clean-room rules

- **No MPEG audio implementation's source was opened or read**, nor
  searched for: not minimp3, LAME, mpg123, libmad, FFmpeg's mpegaudio*,
  Helix, the ISO dist10 reference code, symphonia or any other. None was
  run either: there is no oracle decoder in the tests, and no C.
- The only external data is ISO's own conformance package (below), used as
  data: compressed streams in, reference waveforms to compare against.

## The standards, by part

**ISO/IEC 11172-3** (MPEG-1 audio):
- 2.4.1–2.4.2: frame header, error check, audio data syntax of Layers I, II
  and III, side information, scalefactors, Huffman code bits, ancillary data.
- 2.4.3.1: the CRC (generator x^16 + x^15 + x^2 + 1, preset all ones).
- 2.4.3.2–2.4.3.3: Layer I and II requantisation, scalefactors, grouping,
  the joint-stereo bound.
- 2.4.3.4: Layer III — main_data_begin and the bit reservoir, scalefactor
  selection information, Huffman decoding, requantisation, reordering,
  M/S and intensity stereo, alias reduction, IMDCT, windowing and
  overlap-add, frequency inversion.
- Annex A: the synthesis filterbank flow chart.
- Annex B tables: 3-B.1 (Layer I/II scalefactors; generated as
  2^(1 - i/3)), 3-B.2a–d (Layer II allocation), 3-B.3 (synthesis window,
  257 values, the rest by the table's symmetry), 3-B.4 (quantisation
  classes; C and D generated from the number of steps and checked against
  the printed values), 3-B.6 (pretab), 3-B.7 (Huffman codes), 3-B.8
  (scalefactor bands), 3-B.9 (alias-reduction coefficients).
- Annex C (informative, encoder): the analysis filterbank flow chart and
  window (D / 32), the quantiser's rounding (nint(x^(3/4) - 0.0946)).
- Annex D (informative): the psychoacoustic model's spreading function and
  the tonal / noise masking offsets, used in simplified form.

**ISO/IEC 13818-3** (MPEG-2 audio, lower sampling frequencies): the 16,
22.05 and 24 kHz rates and bit-rate tables, Layer II allocation table B.1,
Layer III LSF side information (one granule, 9-bit scalefac_compress, no
scfsi), the scalefactor length derivation and nr_of_sfb table, LSF
intensity stereo (intensity_scale and the is_pos pairing), the LSF
scalefactor bands.

**MPEG-2.5** (8, 11.025, 12 kHz; not part of either standard): the header's
version bits 00 and the band tables as commonly documented — 11.025 and
12 kHz use the 16 kHz tables; 8 kHz its own.

**Xing / Info, LAME extension, VBRI**: the publicly documented byte layouts
(field order and widths, the 12 + 12-bit delay and padding, CRC-16/ARC for
the music and tag checksums, the 529-sample decoder delay convention).

**Psychoacoustics**: Zwicker and Terhardt's Bark formula and Terhardt's
threshold-in-quiet approximation (published literature).

## The tables

Reproduced as data from the standards. They are checked by tests that do
not depend on any other implementation:

- every Huffman table is a complete prefix code (Kraft sum exactly 1, no
  word a prefix of another) and every word decodes to its own entry;
- the synthesis window's printed values are multiples of 2^-16, its
  prototype is symmetric and smooth (third differences bounded);
- the scalefactor band tables rise to 576 and 192;
- the quantisation classes are symmetric and match the printed C and D.

And by ISO's conformance data: every Layer I, II and III sequence of the
package decodes at full accuracy (docs/VALIDATION.md), at all nine
sampling frequencies, which a wrong code word, band edge or window
coefficient in anything those streams use would prevent.

## Where the standards leave the decoder a choice

Settled by ISO's reference waveforms, which are the reference decoder's
output; in each case the crate's comment names the sequence:

- **The highest scalefactor band's intensity position** (long band 21,
  short band 12, which have no scalefactor): it continues the band below's
  position when that band is itself intensity-coded, and is 0 when
  intensity coding starts at the highest band (l3_20245, l3_20246,
  l3_20248, l3_25208). The rule is applied to MPEG-1 as well; no MPEG-1
  sequence exercises it.
- **mixed_block_flag on start and stop blocks**: the two lowest subbands use
  the normal window whenever the flag is set, whatever the block type
  (l3_10203, l3_10105).
- **A change between mono and stereo**: the filterbanks start afresh, the
  bit reservoir carries over (l3_10201).
- **LSF illegal intensity positions**: a position equal to 2^slen - 1 marks
  the band as not intensity-coded, including slen 0.

## The conformance data

ISO/IEC 14496-26 (2nd edition) conformance package, freely available
from ISO at https://standards.iso.org/iso-iec/14496/-26/ed-2/en/:
`compressedMp4.zip` (Layer I and II streams, whose MP4 samples are plain
frames), `compressedMpeg12.zip` (Layer III streams in native format) and
`referencesWav.zip` (reference waveforms). `tools/fetch_conformance.py`
reads only the members needed, by HTTP range requests, and unwraps the MP4
files with a box parser of its own. The data is ISO's and is not in the
repository.
