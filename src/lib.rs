//! Sequence search primitives and FASTA/FASTQ readers.
//!
//! The crate exposes a small backend API suitable for command-line tools,
//! services, and larger bioinformatics applications. Implementation details
//! such as packed target representations and parallel chunk scheduling remain
//! private so callers are not coupled to storage internals.
//!
//! # Example
//!
//! ```
//! use seqgrep::{AmbigMode, SearchEngine, SearchQuery, SequenceRecord};
//!
//! # fn main() -> seqgrep::Result<()> {
//! let engine = SearchEngine::default();
//! let prepared = engine.prepare_query(
//!     SearchQuery::new("ATGN").with_ambig_mode(AmbigMode::Query),
//! )?;
//! let record = SequenceRecord::new("example", "ATGNATGA");
//! let matches = engine.search_prepared(&record, &prepared)?;
//! assert_eq!(matches[0].start, 5);
//! # Ok(())
//! # }
//! ```

#![forbid(unsafe_code)]
#![deny(missing_docs)]
#![deny(rustdoc::broken_intra_doc_links)]

mod alphabet;
mod codec;
mod error;
pub mod fastx;
mod packed;
mod search;

/// Crate version reported by the bundled command-line tool.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

pub use error::{Error, Result};
pub use fastx::{FastxReader, InputFormat};
pub use search::{
    AmbigMode, ExecutionConfig, Match, PreparedQuery, SearchEngine, SearchQuery, SequenceRecord,
    SequenceType, Strand,
};

#[cfg(test)]
mod trait_tests {
    use super::*;

    fn assert_send_sync<T: Send + Sync>() {}
    fn assert_send<T: Send>() {}

    #[test]
    fn backend_types_are_thread_safe() {
        assert_send_sync::<SearchEngine>();
        assert_send_sync::<PreparedQuery>();
        assert_send_sync::<SequenceRecord>();
        assert_send::<FastxReader>();
    }
}
