//! The crate's one error type.

/// What went wrong. Every malformed input comes back as one of these; the
/// decoder never panics on bytes it is given.
#[derive(thiserror::Error, Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// The bitstream breaks the syntax or the semantics of the standard: a
    /// reserved field value, a codeword that matches nothing, a CRC that
    /// does not match, data that runs past the end of the frame.
    #[error("invalid MPEG audio data: {0}")]
    Invalid(String),
    /// Valid MPEG audio this crate does not implement, named.
    #[error("unsupported MPEG audio feature: {0}")]
    Unsupported(String),
    /// A configuration the caller asked for that cannot be coded: a channel
    /// count, sample rate or bit rate outside what the encoder supports.
    #[error("invalid MP3 encoder configuration: {0}")]
    Config(String),
}

pub(crate) fn invalid(msg: impl Into<String>) -> Error {
    Error::Invalid(msg.into())
}

pub(crate) fn config(msg: impl Into<String>) -> Error {
    Error::Config(msg.into())
}

/// `std::result::Result` with this crate's [`Error`].
pub type Result<T> = std::result::Result<T, Error>;
