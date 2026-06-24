# Architecture

The public API is deliberately smaller than the implementation:

```text
SequenceRecord
SearchQuery
PreparedQuery
ExecutionConfig
SearchEngine
Match
FastxReader
```

Internal modules own alphabet validation, packed storage, codec semantics,
exact matching, and parallel chunk scheduling.

## Search paths

```text
exact + jobs=1        KMP over normalized bytes
IUPAC + jobs=1        packed compatibility scan
any mode + jobs>1     packed compatibility scan with scoped threads
```

Exact serial search uses KMP rather than a regular-expression engine. It is
linear, supports overlapping matches, and searches circular virtual sequences
without materializing repeated chromosomes.

Parallel work is split into deterministic coordinate chunks. Threads borrow a
single immutable packed target and send chunk-local start positions through
bounded channels. The caller consumes contiguous worker partitions in coordinate
order, so output is deterministic without retaining all chromosome hits in
memory. No global thread pool, shared-memory segment, serialization, or unsafe
code is used.

`PreparedQuery` validates symbols and computes forward/reverse query encodings
once. This matters for FASTQ and multi-record FASTA inputs, where re-encoding the
same pattern for every record would be unnecessary work. `search_record` remains
a convenience wrapper; backend integrations should prepare once and use
`search_prepared` or `visit_prepared_matches`.

## Packed target representations

```text
canonical nucleotide          2 bits/base
exact mixed nucleotide        5 bits/symbol
query-mode mixed target       2-bit bases + validity/gap bitmaps
both-mode mixed target        5-bit IUPAC masks
protein                       5 bits/residue
```

Packing is private. Downstream tools depend on biological semantics, not the
physical representation, which allows storage strategies to evolve without a
breaking API change.
