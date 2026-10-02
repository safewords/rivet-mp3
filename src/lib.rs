//! MPEG audio.
mod bits;
mod crc;
mod error;
pub mod header;
mod tables;
pub use error::{Error, Result};
