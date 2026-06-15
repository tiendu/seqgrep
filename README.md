# seqgrep

DNA/protein-aware grep for FASTA/FASTQ files.

Default mode is literal matching, so protein FASTA works safely:

```bash
seqgrep MTEYK proteins.fa
seqgrep ATGN genome.fa
```

Use `--ambig` to enable IUPAC nucleotide ambiguity:

```bash
seqgrep ATGNNRY genome.fa --ambig
seqgrep GGGAAA plasmid.fa --ambig --revcomp --circular
```

Supports:

- FASTA and FASTQ
- plain text and `.gz`
- literal search by default
- IUPAC nucleotide ambiguity with `--ambig`
- reverse-complement search with `--revcomp`
- circular DNA search with `--circular`
- multiprocessing over a single long sequence with `--jobs`

## Install for development

```bash
make install
make check
```

## Examples

```bash
seqgrep ATG tests/example.fa --with-header
seqgrep ATGN tests/example.fa --ambig --with-header
seqgrep GGGAAA tests/example.fa --ambig --revcomp --circular --with-header
seqgrep ATG huge.fa.gz --ambig --jobs 8 --chunk-size 1000000
```
