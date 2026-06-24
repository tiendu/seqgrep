use std::error::Error as StdError;
use std::fmt;
use std::io;

/// Result type used throughout the crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Errors produced by parsing, validation, and search operations.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    /// An underlying I/O operation failed.
    Io(io::Error),
    /// A sequence contained a symbol outside the selected alphabet.
    InvalidSymbol {
        /// Human-readable alphabet name.
        alphabet: &'static str,
        /// Invalid symbol.
        symbol: char,
        /// One-based symbol position.
        position: usize,
    },
    /// A search or execution option was invalid.
    InvalidOption(String),
    /// A FASTA or FASTQ record was malformed.
    Parse {
        /// Input format name.
        format: &'static str,
        /// Approximate one-based input line.
        line: usize,
        /// Parse failure description.
        message: String,
    },
    /// Input format could not be inferred from the path.
    UnknownFormat(String),
    /// A worker thread terminated unexpectedly.
    WorkerPanicked,
}

impl fmt::Display for Error {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => write!(formatter, "{error}"),
            Self::InvalidSymbol {
                alphabet,
                symbol,
                position,
            } => write!(
                formatter,
                "Unsupported {alphabet} symbol {symbol:?} at position {position}"
            ),
            Self::InvalidOption(message) => formatter.write_str(message),
            Self::Parse {
                format,
                line,
                message,
            } => write!(formatter, "Invalid {format} near line {line}: {message}"),
            Self::UnknownFormat(message) => formatter.write_str(message),
            Self::WorkerPanicked => formatter.write_str("a search worker thread panicked"),
        }
    }
}

impl StdError for Error {
    fn source(&self) -> Option<&(dyn StdError + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            _ => None,
        }
    }
}

impl From<io::Error> for Error {
    fn from(error: io::Error) -> Self {
        Self::Io(error)
    }
}
