# seqgrep

`seqgrep` searches nucleotide and amino-acid records in FASTA or FASTQ files.
It has no runtime dependencies outside the Python standard library.

The default is **exact nucleotide search**:

```bash
seqgrep ATG genome.fa
seqgrep ATGN genome.fa.gz --with-header
```

Enable nucleotide ambiguity explicitly:

```bash
seqgrep ATGNNRY genome.fa --ambig
seqgrep GGGAAA plasmid.fa --ambig --revcomp --circular
```

Select amino-acid mode for proteins:

```bash
seqgrep MTEYK proteins.fa --sequence-type amino-acid
seqgrep LVVVG proteins.fa.gz -t amino-acid --jobs 4
```

## Matching modes

### Exact nucleotide — default

```bash
seqgrep PATTERN INPUT
```

IUPAC letters are literal unless `--ambig` is present. For example, `N`
matches only `N`. `T` and `U` are treated as equivalent nucleotide symbols.
The gap characters `-` and `.` remain distinct in exact mode.

Serial exact searches use Python's native string-search implementation. With
`--jobs > 1`, canonical `A/C/G/T/U` targets are packed at two bits per base in
shared memory; targets containing ambiguity symbols or gaps use one exact code
per byte.

### IUPAC nucleotide ambiguity

```bash
seqgrep PATTERN INPUT --ambig
```

Compatibility is enabled for both pattern and target:

```text
A  A                 B  C/G/T
C  C                 D  A/G/T
G  G                 H  A/C/T
T  T                 V  A/C/G
U  T                 N  A/C/G/T
R  A/G               -  gap
Y  C/T               .  gap
S  G/C
W  A/T
K  G/T
M  A/C
```

In ambiguity mode, `-` and `.` are equivalent gap symbols. Canonical targets
remain two-bit packed; ambiguous targets use one mask byte per symbol.

### Amino acid

```bash
seqgrep PATTERN INPUT --sequence-type amino-acid
```

Protein search is exact. Supported symbols are:

```text
A C D E F G H I K L M N P Q R S T V W Y
B J Z X U O * - .
```

`B`, `J`, `Z`, and `X` are accepted as literal symbols. Protein ambiguity
expansion is intentionally not enabled, so `--ambig` is nucleotide-only.
Serial exact protein search uses native string search. Multiprocessing stores
targets at five bits per residue in shared memory.

## Execution strategy

`seqgrep` chooses the backend automatically:

```text
exact nucleotide, jobs=1     native string search
exact amino acid, jobs=1     native string search
IUPAC nucleotide, jobs=1     encoded compatibility scanner
any mode, jobs>1             packed shared-memory scanner
```

Packing therefore reduces memory where it matters without imposing bit
unpacking on ordinary serial exact searches.

## Features

- exact nucleotide search by default
- optional IUPAC nucleotide compatibility
- exact amino-acid search
- reverse-complement nucleotide search
- circular sequence search, including patterns longer than the record
- FASTA and multiline FASTQ input
- plain text and gzip input
- overlapping matches
- shared-memory multiprocessing for long records
- TSV output with 1-based inclusive coordinates
- compact two-bit nucleotide and five-bit protein targets

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
        [--version]
        pattern input
```

### Output

Each match is a TSV row:

```text
record  strand  start  end  matched  circular
```

Coordinates are 1-based and inclusive. For a circular match crossing the
record boundary, `end` wraps and `circular` is `true`.

## Development

```bash
make install
make check
```

Individual checks:

```bash
make test
make lint
make format-check
make typecheck
```

## Internal organization

```text
alphabets.py  symbol tables, validation, masks, packing, unpacking
codecs.py     matching semantics and target representation selection
exact.py      fast native exact serial search
window.py     serial IUPAC compatibility search
chunked.py    shared-memory multiprocessing search
planner.py    backend selection
```
