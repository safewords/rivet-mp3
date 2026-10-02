//! Layer III tables other than the Huffman codes: scalefactor bands
//! (ISO/IEC 11172-3 Table 3-B.8; ISO/IEC 13818-3 Table B.2 for 16, 22.05
//! and 24 kHz; the MPEG-2.5 rates as publicly documented), pretab
//! (Table 3-B.6), the scalefactor lengths of scalefac_compress (2.4.2.7)
//! and the LSF nr_of_sfb table (13818-3 2.4.3.2).

/// Long-block scalefactor band boundaries (23 entries, 0..=576) by
/// sampling frequency, in [`crate::FrameHeader::sample_rate`] order:
/// 44.1, 48, 32 kHz; 22.05, 24, 16 kHz; 11.025, 12, 8 kHz.
pub(crate) const SFB_LONG: [[u16; 23]; 9] = [
    [0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 52, 62, 74, 90, 110, 134, 162, 196, 238, 288, 342, 418, 576],
    [0, 4, 8, 12, 16, 20, 24, 30, 36, 42, 50, 60, 72, 88, 106, 128, 156, 190, 230, 276, 330, 384, 576],
    [0, 4, 8, 12, 16, 20, 24, 30, 36, 44, 54, 66, 82, 102, 126, 156, 194, 240, 296, 364, 448, 550, 576],
    [0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396, 464, 522, 576],
    [0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 114, 136, 162, 194, 232, 278, 332, 394, 464, 540, 576],
    [0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396, 464, 522, 576],
    [0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396, 464, 522, 576],
    [0, 6, 12, 18, 24, 30, 36, 44, 54, 66, 80, 96, 116, 140, 168, 200, 238, 284, 336, 396, 464, 522, 576],
    [0, 12, 24, 36, 48, 60, 72, 88, 108, 132, 160, 192, 232, 280, 336, 400, 476, 566, 568, 570, 572, 574, 576],
];

/// Short-block scalefactor band boundaries (14 entries, 0..=192, per
/// window), same order as [`SFB_LONG`].
pub(crate) const SFB_SHORT: [[u16; 14]; 9] = [
    [0, 4, 8, 12, 16, 22, 30, 40, 52, 66, 84, 106, 136, 192],
    [0, 4, 8, 12, 16, 22, 28, 38, 50, 64, 80, 100, 126, 192],
    [0, 4, 8, 12, 16, 22, 30, 42, 58, 78, 104, 138, 180, 192],
    [0, 4, 8, 12, 18, 24, 32, 42, 56, 74, 100, 132, 174, 192],
    [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 136, 180, 192],
    [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 134, 174, 192],
    [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 134, 174, 192],
    [0, 4, 8, 12, 18, 26, 36, 48, 62, 80, 104, 134, 174, 192],
    [0, 8, 16, 24, 36, 52, 72, 96, 124, 160, 162, 164, 166, 192],
];

/// pretab (Table 3-B.6), by long scalefactor band 0..=21.
pub(crate) const PRETAB: [u8; 22] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 2, 0];

/// MPEG-1 scalefac_compress -> (slen1, slen2), 2.4.2.7.
pub(crate) const SLEN: [(u8, u8); 16] = [
    (0, 0),
    (0, 1),
    (0, 2),
    (0, 3),
    (3, 0),
    (1, 1),
    (1, 2),
    (1, 3),
    (2, 1),
    (2, 2),
    (2, 3),
    (3, 1),
    (3, 2),
    (3, 3),
    (4, 2),
    (4, 3),
];

/// nr_of_sfb_block[table][block_index][group] (13818-3 2.4.3.2): rows are
/// the six slen derivations of scalefac_compress (three without intensity
/// stereo on the channel, three with), block_index 0 for long blocks, 1
/// short, 2 mixed.
pub(crate) const NR_OF_SFB: [[[u8; 4]; 3]; 6] = [
    [[6, 5, 5, 5], [9, 9, 9, 9], [6, 9, 9, 9]],
    [[6, 5, 7, 3], [9, 9, 12, 6], [6, 9, 12, 6]],
    [[11, 10, 0, 0], [18, 18, 0, 0], [15, 18, 0, 0]],
    [[7, 7, 7, 0], [12, 12, 12, 0], [6, 15, 12, 0]],
    [[6, 6, 6, 3], [12, 9, 9, 6], [6, 12, 9, 6]],
    [[8, 8, 5, 0], [15, 12, 9, 0], [6, 18, 9, 0]],
];

/// Index into the 9-row band tables for a frame.
pub(crate) fn rate_index(version_index: usize, sample_rate_index: usize) -> usize {
    version_index * 3 + sample_rate_index
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bands_rise_to_the_granule() {
        for t in SFB_LONG {
            assert_eq!(t[0], 0);
            assert_eq!(t[22], 576);
            assert!(t.windows(2).all(|w| w[0] < w[1] && (w[1] - w[0]) % 2 == 0));
        }
        for t in SFB_SHORT {
            assert_eq!(t[0], 0);
            assert_eq!(t[13], 192);
            assert!(t.windows(2).all(|w| w[0] < w[1] && (w[1] - w[0]) % 2 == 0));
        }
    }

    #[test]
    fn lsf_groups_cover_the_bands() {
        for row in NR_OF_SFB {
            assert!(row[0].iter().map(|&n| n as usize).sum::<usize>() <= 21);
            assert!(row[1].iter().map(|&n| n as usize).sum::<usize>() <= 36);
        }
    }
}
