//! Layer I and Layer II audio data (ISO/IEC 11172-3 2.4.1.5–2.4.1.6,
//! 2.4.2.5–2.4.2.6, 2.4.3.2–2.4.3.3; ISO/IEC 13818-3 for the lower
//! sampling frequencies) decoded to subband samples.

use crate::bits::BitReader;
use crate::crc::frame_crc_bits;
use crate::error::{Result, invalid};
use crate::header::{FrameHeader, Mode};
use crate::tables::layer12::{AllocTable, QuantClass, scalefactor};

/// Subband samples of one frame: `[channel][time slot][subband]`.
pub(crate) type Subbands = Vec<[[f32; 32]; 2]>;

/// The intensity-stereo bound (first subband coded jointly).
fn bound(h: &FrameHeader, sblimit: usize) -> usize {
    if h.mode == Mode::JointStereo {
        (4 + 4 * usize::from(h.mode_extension)).min(sblimit)
    } else {
        sblimit
    }
}

/// Check the frame CRC when present: it covers header bits 16..32 and the
/// `bits` audio-data bits that follow the CRC word.
fn check_crc(frame: &[u8], h: &FrameHeader, bits: usize) -> Result<()> {
    if !h.crc {
        return Ok(());
    }
    if frame.len() * 8 < 48 + bits {
        return Err(invalid("frame too short for its CRC-protected fields"));
    }
    let crc = frame_crc_bits(0xFFFF, frame, 16, 16);
    let crc = frame_crc_bits(crc, frame, 48, bits);
    let stored = u16::from_be_bytes([frame[4], frame[5]]);
    if crc != stored {
        return Err(invalid(format!(
            "CRC mismatch: frame says {stored:04x}, data gives {crc:04x}"
        )));
    }
    Ok(())
}

/// Layer I: 12 time slots of 32 subbands per channel.
pub(crate) fn decode_layer1(
    frame: &[u8],
    h: &FrameHeader,
    out: &mut Subbands,
    check: bool,
) -> Result<()> {
    let nch = h.channels();
    let bound = bound(h, 32);
    let mut r = BitReader::at(frame, h.header_len() * 8);
    let start = r.pos();
    let mut alloc = [[0u8; 32]; 2];
    for sb in 0..32 {
        let chans = if sb < bound { nch } else { 1 };
        for ch in 0..chans {
            let a = r.read(4)? as u8;
            if a == 15 {
                return Err(invalid("Layer I allocation 15 is forbidden"));
            }
            alloc[ch][sb] = a;
        }
        if sb >= bound {
            alloc[1][sb] = alloc[0][sb];
        }
    }
    if check {
        check_crc(frame, h, r.pos() - start)?;
    }
    let mut scf = [[0.0f32; 32]; 2];
    for sb in 0..32 {
        for ch in 0..nch {
            if alloc[ch][sb] != 0 {
                let i = r.read(6)? as usize;
                if i == 63 {
                    return Err(invalid("scalefactor index 63 is not defined"));
                }
                scf[ch][sb] = scalefactor(i);
            }
        }
    }
    out.clear();
    out.resize(12, [[0.0; 32]; 2]);
    for slot in out.iter_mut() {
        for sb in 0..32 {
            if sb < bound {
                for ch in 0..nch {
                    let a = alloc[ch][sb];
                    if a != 0 {
                        let nb = u32::from(a) + 1;
                        let v = r.read(nb)?;
                        slot[ch][sb] = layer1_class(nb).dequantise(v) * scf[ch][sb];
                    }
                }
            } else {
                let a = alloc[0][sb];
                if a != 0 {
                    let nb = u32::from(a) + 1;
                    let v = layer1_class(nb).dequantise(r.read(nb)?);
                    for ch in 0..nch {
                        slot[ch][sb] = v * scf[ch][sb];
                    }
                }
            }
        }
    }
    Ok(())
}

fn layer1_class(nb: u32) -> QuantClass {
    QuantClass {
        steps: (1 << nb) - 1,
        bits: nb,
        grouped: false,
    }
}

/// Layer II: 36 time slots of 32 subbands per channel.
pub(crate) fn decode_layer2(
    frame: &[u8],
    h: &FrameHeader,
    free_bitrate: u32,
    out: &mut Subbands,
    check: bool,
) -> Result<()> {
    let nch = h.channels();
    let bitrate = if h.is_free_format() {
        free_bitrate
    } else {
        h.bitrate()
    };
    let table = AllocTable::select(h.version.is_lsf(), h.sample_rate(), bitrate, nch);
    let sblimit = table.sblimit();
    let bound = bound(h, sblimit);
    let mut r = BitReader::at(frame, h.header_len() * 8);
    let start = r.pos();
    let mut alloc: [[Option<QuantClass>; 32]; 2] = [[None; 32]; 2];
    for sb in 0..sblimit {
        let row = table.row(sb);
        let chans = if sb < bound { nch } else { 1 };
        for ch in 0..chans {
            let a = r.read(row.nbal)? as usize;
            alloc[ch][sb] = if a == 0 {
                None
            } else {
                Some(row.classes[a - 1])
            };
        }
        if sb >= bound {
            alloc[1][sb] = alloc[0][sb];
        }
    }
    let mut scfsi = [[0u8; 32]; 2];
    for sb in 0..sblimit {
        for ch in 0..nch {
            if alloc[ch][sb].is_some() {
                scfsi[ch][sb] = r.read(2)? as u8;
            }
        }
    }
    if check {
        check_crc(frame, h, r.pos() - start)?;
    }
    let mut scf = [[[0.0f32; 3]; 32]; 2];
    for sb in 0..sblimit {
        for ch in 0..nch {
            if alloc[ch][sb].is_none() {
                continue;
            }
            let mut read = || -> Result<f32> {
                let i = r.read(6)? as usize;
                if i == 63 {
                    return Err(invalid("scalefactor index 63 is not defined"));
                }
                Ok(scalefactor(i))
            };
            scf[ch][sb] = match scfsi[ch][sb] {
                0 => [read()?, read()?, read()?],
                1 => {
                    let a = read()?;
                    let b = read()?;
                    [a, a, b]
                }
                2 => {
                    let a = read()?;
                    [a, a, a]
                }
                _ => {
                    let a = read()?;
                    let b = read()?;
                    [a, b, b]
                }
            };
        }
    }
    out.clear();
    out.resize(36, [[0.0; 32]; 2]);
    for gr in 0..12 {
        let part = gr / 4;
        for sb in 0..sblimit {
            let chans = if sb < bound { nch } else { 1 };
            for ch in 0..chans {
                let Some(class) = alloc[ch][sb] else { continue };
                let v = read_triple(&mut r, &class)?;
                for (s, &code) in v.iter().enumerate() {
                    let x = class.dequantise(code);
                    if sb < bound {
                        out[gr * 3 + s][ch][sb] = x * scf[ch][sb][part];
                    } else {
                        for c in 0..nch {
                            out[gr * 3 + s][c][sb] = x * scf[c][sb][part];
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

/// Three consecutive sample codes of one subband: one grouped codeword or
/// three separate codes.
fn read_triple(r: &mut BitReader, class: &QuantClass) -> Result<[u32; 3]> {
    if class.grouped {
        let mut c = r.read(class.bits)?;
        let n = class.steps;
        if c >= n * n * n {
            return Err(invalid("grouped codeword out of range"));
        }
        let mut v = [0; 3];
        for x in v.iter_mut() {
            *x = c % n;
            c /= n;
        }
        Ok(v)
    } else {
        let a = r.read(class.bits)?;
        let b = r.read(class.bits)?;
        let c = r.read(class.bits)?;
        Ok([a, b, c])
    }
}
