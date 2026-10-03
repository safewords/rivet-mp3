//! MPEG audio, both ways, written from ISO/IEC 11172-3 and ISO/IEC 13818-3.
//!
//! - [`decode`]: a decoder for Layers I, II and III at every MPEG-1, MPEG-2
//!   (lower sampling frequencies) and MPEG-2.5 rate: free format, intensity
//!   and mid/side stereo, every block type including mixed blocks, the bit
//!   reservoir, CRC checking; Xing / Info / LAME and VBRI tags, with
//!   gapless trimming.
//! - [`encode`]: a Layer III (MP3) encoder: constant and variable bit rate,
//!   a psychoacoustic model, mid/side stereo, the bit reservoir, short
//!   blocks for transients, and a Xing / Info tag frame with LAME-style
//!   gapless delay and padding.
//! - [`header`]: the frame header; [`xing`]: the tag frames.
//!
//! PCM on both sides is interleaved `f32` at full scale ±1.0, one or two
//! channels.

// Filterbanks, transforms and band loops index several arrays with one
// counter; written as index loops they read like the standard's formulas.
#![allow(clippy::needless_range_loop)]
#![warn(missing_docs)]

pub mod decode;
pub mod encode;
mod bits;
mod crc;
mod error;
pub mod header;
mod layer12;
pub mod layer3;
mod synth;
mod tables;
pub mod xing;

pub use decode::{Decoder, DecoderOptions, Frame, FrameDecoder, Gapless};
pub use encode::{BitrateMode, Encoder, EncoderConfig};
pub use error::{Error, Result};
pub use header::{FrameHeader, Layer, Mode, Version};
