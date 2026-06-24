# seqgrep

`seqgrep` is a Rust library and command-line tool for searching nucleotide and
amino-acid records in FASTA and FASTQ files.

The Rust rewrite keeps the existing CLI behavior while providing a small,
intentionally narrow backend API for other tools. It uses no unsafe code.

## Install

From a source checkout:

```bash
cargo install --path .
```

For development:

```bash
make check
make release
```

## CLI

Exact nucleotide search is the default:

```bash
seqgrep ATG genome.fa
seqgrep ATGN genome.fa.gz --with-header
```

Choose how IUPAC ambiguity is interpreted:

```bash
# Query-side ambiguity; recommended for assembled reference genomes
seqgrep ATGNNRY genome.fa --ambig-mode query

# Symmetric compatibility between query and target symbols
seqgrep AAAAA consensus.fa --ambig-mode both

# Backward-compatible alias for query mode
seqgrep ATGNNRY genome.fa --ambig
```

Search proteins with the amino-acid alphabet:

```bash
seqgrep MTEYK proteins.faa --sequence-type amino-acid
seqgrep LVVVG proteins.fa.gz -t amino-acid --jobs 4
```

Complete usage:

```text
seqgrep [OPTIONS] PATTERN INPUT

Options:
  -t, --sequence-type TYPE    nucleotide or amino-acid
      --ambig-mode MODE       none, query, or both
      --ambig                 alias for --ambig-mode query
      --revcomp               also search the nucleotide reverse complement
      --circular              allow matches to cross the sequence boundary
      --with-header           print a TSV header
      --format FORMAT         auto, fasta, or fastq
  -j, --jobs JOBS             worker threads per record
      --chunk-size N          candidate starts per work chunk
```

Each result is one TSV row:

```text
record  strand  start  end  matched  circular
```

Coordinates are one-based and inclusive. Circular matches report a wrapped end
coordinate and `circular=true`.

## Ambiguity modes

### `none`

Exact nucleotide matching. IUPAC letters remain literal, `T` and `U` are
equivalent, and `-` and `.` remain distinct.

### `query`

Ambiguity expands only in the query:

```text
query N  matches target A, C, G, or T
query R  matches target A or G
query A  does not match target N
query N  does not match target N
```

This avoids false wildcard matches inside assembly `N` blocks.

### `both`

Ambiguity expands in both query and target. Two positions match when their
IUPAC base sets overlap:

```text
query AAAAA  matches target NNANN
query NNNNN  matches target AAAAA
```

Use this for consensus or deliberately ambiguous references. Whole-genome
searches may produce very large outputs inside unknown regions.

In `query` and `both` modes, `-` and `.` are equivalent gaps. Gaps do not match
nucleotide bases.

## Amino-acid mode

Protein search is exact. Supported symbols are:

```text
A C D E F G H I K L M N P Q R S T V W Y
B J Z X U O * - .
```

Extended symbols are literal. Protein ambiguity expansion and reverse
complement are intentionally unsupported.

## Backend API

The public API is intentionally small. One-off searches can use the collecting
convenience method:

```rust
use seqgrep::{AmbigMode, SearchEngine, SearchQuery, SequenceRecord};

fn main() -> seqgrep::Result<()> {
    let engine = SearchEngine::default();
    let record = SequenceRecord::new("example", "ATGNATGA");
    let query = SearchQuery::new("ATGN").with_ambig_mode(AmbigMode::Query);

    for matched in engine.search_record(&record, &query)? {
        println!("{}:{}-{}", matched.record, matched.start, matched.end);
    }

    Ok(())
}
```

For FASTA/FASTQ streams, prepare the query once and visit matches without
collecting every result in memory:

```rust
use seqgrep::{FastxReader, InputFormat, SearchEngine, SearchQuery};

fn main() -> seqgrep::Result<()> {
    let engine = SearchEngine::default();
    let prepared = engine.prepare_query(SearchQuery::new("ATG"))?;

    for record in FastxReader::from_path("genome.fa.gz", InputFormat::Auto)? {
        let record = record?;
        engine.visit_prepared_matches(&record, &prepared, |matched| {
            println!("{}\t{}", matched.record, matched.start);
            Ok(())
        })?;
    }

    Ok(())
}
```

`SequenceRecord`, `SearchQuery`, `PreparedQuery`, and `SearchEngine` are
`Send + Sync`; `FastxReader` is `Send`. Their internal invariants and packed
storage are private so downstream tools cannot accidentally construct invalid
backend state. Public enums and errors are non-exhaustive so compatible variants
can be added without forcing a major-version release.

See `examples/backend.rs` and `docs/architecture.md`.

## Reliability and execution model

- safe Rust only (`#![forbid(unsafe_code)]`)
- deterministic serial and parallel result ordering
- prepared queries reusable across many records
- streaming match visitor with bounded parallel buffering
- scoped worker threads with no global thread pool
- no multiprocessing, serialization, or shared-memory cleanup
- streaming FASTA and multiline FASTQ parsing
- transparent plain-text and gzip input
- FASTQ quality-length validation
- overlapping matches
- circular patterns longer than the target
- adaptive packed target storage
- explicit error types suitable for backend callers

Search path selection:

```text
exact nucleotide, jobs=1     linear KMP search
exact amino acid, jobs=1     linear KMP search
IUPAC nucleotide, jobs=1     packed compatibility scan
any mode, jobs>1             deterministic packed threaded scan
```

Packed targets remain internal:

```text
canonical nucleotide          2 bits/base
exact mixed nucleotide        5 bits/symbol
query-mode mixed target       2-bit bases + validity/gap bitmaps
both-mode mixed target        5-bit IUPAC masks
protein                       5 bits/residue
```

## Development

```bash
make format
make lint
make test
make doc
make check
```

Real chromosome integration test:

```bash
make test-chr21
```

The test downloads and caches GRCh38 chromosome 21 under `tests/data/`, then
compares serial and threaded output and verifies query-only versus symmetric
ambiguity semantics.
