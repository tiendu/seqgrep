# Changelog

## 0.5.0

- Reimplemented seqgrep as a Rust library and CLI.
- Preserved exact nucleotide, IUPAC query-only, IUPAC symmetric, and exact
  amino-acid semantics.
- Replaced multiprocessing and shared memory with deterministic scoped threads.
- Preserved adaptive 2-bit and 5-bit target representations.
- Added a streaming FASTA/FASTQ reader with transparent gzip support.
- Added a narrow backend API with private invariants for embedding in other Rust tools.
- Added reusable prepared queries and bounded streaming match visitation.
- Added deterministic reference-parity tests and cross-platform CI with an MSRV check.
