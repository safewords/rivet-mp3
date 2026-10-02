//! The two CRC-16s an MP3 file carries.
//!
//! - The frame CRC of ISO/IEC 11172-3 2.4.3.1 ("CRC check"): generator
//!   x^16 + x^15 + x^2 + 1, register preset to all ones, the bits fed most
//!   significant first.
//! - The checksums of the LAME tag (music CRC and tag CRC), which use the
//!   same polynomial bit-reflected (0xA001), register preset to zero and
//!   bytes fed least significant bit first — the "CRC-16/ARC" form the
//!   publicly documented LAME tag layout names.

/// The frame CRC over `bits` bits read MSB-first from `data`, starting at
/// bit `from`, continuing a register `crc`.
pub(crate) fn frame_crc_bits(mut crc: u16, data: &[u8], from: usize, bits: usize) -> u16 {
    for i in from..from + bits {
        let bit = (data[i / 8] >> (7 - i % 8)) & 1;
        let top = (crc >> 15) as u8 & 1;
        crc <<= 1;
        if top ^ bit == 1 {
            crc ^= 0x8005;
        }
    }
    crc
}

/// CRC-16/ARC over whole bytes.
pub(crate) fn crc16_arc(data: &[u8]) -> u16 {
    let mut crc: u16 = 0;
    for &b in data {
        crc ^= u16::from(b);
        for _ in 0..8 {
            crc = if crc & 1 == 1 { (crc >> 1) ^ 0xA001 } else { crc >> 1 };
        }
    }
    crc
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arc_check_value() {
        // The standard check input "123456789".
        assert_eq!(crc16_arc(b"123456789"), 0xBB3D);
    }

    #[test]
    fn frame_crc_is_crc16_with_all_ones_preset() {
        // MSB-first CRC-16 (poly 0x8005, init 0xFFFF, no reflection, no final
        // xor) over "123456789" is 0xAEE7 (the "CRC-16/CMS" check value).
        let d = b"123456789";
        assert_eq!(frame_crc_bits(0xFFFF, d, 0, 72), 0xAEE7);
    }
}
