# seqgrep

Compact sequence search for nucleotide and amino-acid FASTA/FASTQ files.

`seqgrep` defaults to **exact nucleotide mode**:

```bash
seqgrep ATG genome.fa
seqgrep ATGN genome.fa
```

Use IUPAC nucleotide ambiguity explicitly:

```bash
seqgrep ATGNNRY genome.fa --ambig
seqgrep GGGAAA plasmid.fa --ambig --revcomp --circular
```

Select amino-acid mode for protein sequences:

```bash
seqgrep MTEYK proteins.fa --sequence-type amino-acid
seqgrep LVVVG proteins.fa.gz -t amino-acid --jobs 4
```

## Sequence modes

### Nucleotide — default

```bash
seqgrep PATTERN INPUT
seqgrep PATTERN INPUT --sequence-type nucleotide
```

Without `--ambig`, symbols are matched exactly. For example, `N` matches a
literal `N`, not every base. `T` and `U` are treated as equivalent nucleotide
symbols.

Canonical targets containing only `A`, `C`, `G`, `T`, or `U` are stored using
two bits per base. Targets containing IUPAC symbols fall back to one compact
exact code per byte.

### Nucleotide ambiguity

```bash
seqgrep PATTERN INPUT --ambig
```

`--ambig` enables IUPAC compatibility for both pattern and target:

- `R` = A or G
- `Y` = C or T
- `N` = any canonical nucleotide
- and the remaining standard IUPAC nucleotide symbols

Canonical targets remain packed at two bits per base. Ambiguous targets use
one four-bit mask per byte.

`--ambig` and `--revcomp` are valid only in nucleotide mode.

### Amino acid

```bash
seqgrep PATTERN INPUT --sequence-type amino-acid
```

Protein matching is exact and uses five-bit packed target storage. Supported
symbols are:

```text
A C D E F G H I K L M N P Q R S T V W Y
B J Z X U O * - .
```

`B`, `J`, `Z`, and `X` are literal symbols in amino-acid mode; they are not
expanded as ambiguity codes.

## Features

- nucleotide and amino-acid sequence types
- two-bit canonical nucleotide targets
- five-bit protein targets
- optional IUPAC nucleotide ambiguity
- reverse-complement nucleotide search
- circular sequence search
- FASTA and FASTQ input
- plain text and `.gz` input
- serial and shared-memory multiprocessing modes
- TSV output with 1-based inclusive coordinates
- no runtime dependencies outside the Python standard library

## Usage

```text
seqgrep [-h]
        [-t {nucleotide,amino-acid}]
        [--ambig]
        [--revcomp]
        [--circular]
        [--with-header]
        [--format {auto,fasta,fastq}]
        [-j JOBS]
        [--chunk-size CHUNK_SIZE]
        pattern input
```

### Exact nucleotide search

```bash
seqgrep ATGN genome.fa --with-header
```

Here `N` matches only a literal `N`.

### IUPAC nucleotide search

```bash
seqgrep ATGN genome.fa --ambig --with-header
```

Here `N` matches any compatible nucleotide.

### Reverse-complement and circular search

```bash
seqgrep GGGAAA plasmid.fa \
  --ambig \
  --revcomp \
  --circular \
  --with-header
```

### Protein search

```bash
seqgrep MTEYK proteins.fa \
  --sequence-type amino-acid \
  --with-header
```

### Multiprocessing for a long sequence

```bash
seqgrep ATGNNRY genome.fa.gz \
  --ambig \
  --jobs 8 \
  --chunk-size 1000000
```

Multiprocessing splits candidate start positions while sharing one encoded
target buffer between workers.

## Output

Without `--with-header`, each match is one TSV row:

```text
record  strand  start  end  matched  circular
```

Coordinates are 1-based and inclusive. For a circular match that crosses the
boundary, `end` wraps to the beginning of the record and `circular` is `true`.

## Install for development

```bash
make install
make check
```

Or install directly:

```bash
python -m pip install -e ".[dev]"
```

## Validation

```bash
make test
make lint
make typecheck
make check
```

