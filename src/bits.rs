//! MSB-first bit reading and writing, the order of every field in an MPEG
//! audio frame (ISO/IEC 11172-3 2.1, "bslbf" / "uimsbf").

use crate::error::{Result, invalid};

/// Reads bits, most significant first, from a byte slice. Reads past the end
/// fail with [`crate::Error::Invalid`] rather than panic.
#[derive(Clone, Debug)]
pub(crate) struct BitReader<'a> {
    data: &'a [u8],
    /// Position in bits from the start of `data`.
    pos: usize,
    /// Reading stops here (in bits); at most `data.len() * 8`.
    end: usize,
}

impl<'a> BitReader<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0, end: data.len() * 8 }
    }

    /// A reader over `data` starting at bit `pos`.
    pub(crate) fn at(data: &'a [u8], pos: usize) -> Self {
        Self { data, pos, end: data.len() * 8 }
    }

    pub(crate) fn pos(&self) -> usize {
        self.pos
    }

    /// Limit reading to `end` bits from the start of the data.
    pub(crate) fn set_end(&mut self, end: usize) {
        self.end = end.min(self.data.len() * 8);
    }

    pub(crate) fn end(&self) -> usize {
        self.end
    }

    pub(crate) fn remaining(&self) -> usize {
        self.end.saturating_sub(self.pos)
    }

    /// Read `n` (0..=32) bits.
    pub(crate) fn read(&mut self, n: u32) -> Result<u32> {
        if n == 0 {
            return Ok(0);
        }
        if self.pos + n as usize > self.end {
            return Err(invalid("bitstream ends inside a field"));
        }
        let v = self.peek_unchecked(n);
        self.pos += n as usize;
        Ok(v)
    }

    pub(crate) fn bit(&mut self) -> Result<bool> {
        Ok(self.read(1)? == 1)
    }

    /// The next `n` (1..=32) bits without consuming them; bits past the end
    /// of the data read as zero.
    pub(crate) fn peek_unchecked(&self, n: u32) -> u32 {
        let mut v: u64 = 0;
        let byte = self.pos / 8;
        for i in 0..5 {
            v = (v << 8) | u64::from(*self.data.get(byte + i).unwrap_or(&0));
        }
        let shift = 40 - (self.pos % 8) as u32 - n;
        ((v >> shift) & ((1u64 << n) - 1)) as u32
    }

    pub(crate) fn skip(&mut self, n: usize) -> Result<()> {
        if self.pos + n > self.end {
            return Err(invalid("bitstream ends inside a field"));
        }
        self.pos += n;
        Ok(())
    }
}

/// Writes bits, most significant first.
#[derive(Clone, Debug, Default)]
pub(crate) struct BitWriter {
    pub(crate) bytes: Vec<u8>,
    /// Bits used in the last byte (0 when byte-aligned).
    used: u32,
}

impl BitWriter {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    /// Bits written so far.
    pub(crate) fn len(&self) -> usize {
        if self.used == 0 { self.bytes.len() * 8 } else { (self.bytes.len() - 1) * 8 + self.used as usize }
    }

    /// Write the low `n` (0..=32) bits of `v`.
    pub(crate) fn put(&mut self, v: u32, n: u32) {
        debug_assert!(n == 32 || v >> n == 0, "{v} does not fit {n} bits");
        let mut n = n;
        while n > 0 {
            if self.used == 0 {
                self.bytes.push(0);
            }
            let room = 8 - self.used;
            let take = room.min(n);
            let chunk = ((v >> (n - take)) & ((1 << take) - 1)) as u8;
            let last = self.bytes.len() - 1;
            self.bytes[last] |= chunk << (room - take);
            self.used = (self.used + take) % 8;
            n -= take;
        }
    }

    pub(crate) fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip() {
        let mut w = BitWriter::new();
        w.put(0b101, 3);
        w.put(0x1234, 16);
        w.put(0, 0);
        w.put(0x7fff_ffff, 31);
        w.put(1, 1);
        let n = w.len();
        assert_eq!(n, 51);
        let bytes = w.into_bytes();
        let mut r = BitReader::new(&bytes);
        assert_eq!(r.read(3).unwrap(), 0b101);
        assert_eq!(r.read(16).unwrap(), 0x1234);
        assert_eq!(r.read(31).unwrap(), 0x7fff_ffff);
        assert!(r.bit().unwrap());
        assert_eq!(r.read(5).unwrap(), 0);
        assert!(r.read(1).is_err());
    }
}
