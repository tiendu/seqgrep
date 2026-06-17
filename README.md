# seqgrep

`seqgrep` searches nucleotide and amino-acid records in FASTA or FASTQ files.
It has no runtime dependencies outside the Python standard library.

The default is exact nucleotide search:

```bash
seqgrep ATG genome.fa
seqgrep ATGN genome.fa.gz --with-header
```

Choose how IUPAC nucleotide ambiguity is interpreted:

```bash
# Safe default for reference genomes: ambiguity only in the query
seqgrep ATGNNRY genome.fa --ambig-mode query

# Symmetric compatibility: ambiguity in both query and target
seqgrep AAAAA consensus.fa --ambig-mode both

# Backward-compatible alias for --ambig-mode query
seqgrep ATGNNRY genome.fa --ambig
```

Select amino-acid mode for proteins:

```bash
seqgrep MTEYK proteins.faa --sequence-type amino-acid
seqgrep LVVVG proteins.fa.gz -t amino-acid --jobs 4
```

## Ambiguity modes

### `none` — exact nucleotide matching

This is the default:

```bash
seqgrep PATTERN INPUT
seqgrep PATTERN INPUT --ambig-mode none
```

IUPAC letters remain literal. For example, `N` matches only `N`. `T` and
`U` are treated as equivalent nucleotide spellings. The gap characters `-`
and `.` remain distinct.

### `query` — query-side IUPAC ambiguity

```bash
seqgrep PATTERN INPUT --ambig-mode query
seqgrep PATTERN INPUT --ambig
```

Ambiguity expands only in the query:

```text
query N  matches target A, C, G, or T
query R  matches target A or G
query A  does not match target N
query N  does not match target N
```

This is the recommended mode for assembled reference genomes. Target `N`
usually means that the reference base is unknown, not that every possible
query should be reported as a match. Long assembly gaps therefore produce no
false wildcard hits.

### `both` — symmetric IUPAC compatibility

```bash
seqgrep PATTERN INPUT --ambig-mode both
```

Ambiguity expands in both query and target. Two symbols match when their IUPAC
base sets overlap:

```text
query AAAAA  matches target NNANN
query NNNNN  matches target AAAAA
query R      matches target A or G
query A      matches target R or N
```

Use this mode for consensus sequences, degenerate references, uncertain base
calls, and compatibility questions. On whole reference genomes it may report
very large numbers of possible matches inside `N` blocks.

In both ambiguity modes, `-` and `.` are equivalent gap symbols. Gaps match
gaps but do not match nucleotide bases.

Supported nucleotide symbols:

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

## Amino-acid mode

```bash
seqgrep PATTERN INPUT --sequence-type amino-acid
```

Protein search is exact. Supported symbols are:

```text
A C D E F G H I K L M N P Q R S T V W Y
B J Z X U O * - .
```

`B`, `J`, `Z`, and `X` are literal symbols. Protein ambiguity expansion is not
enabled, so `--ambig` and `--ambig-mode` are nucleotide-only.

## Execution and storage

`seqgrep` chooses the backend automatically:

```text
exact nucleotide, jobs=1     native string search
exact amino acid, jobs=1     native string search
IUPAC nucleotide, jobs=1     encoded compatibility scanner
any mode, jobs>1             packed shared-memory scanner
```

Packed target representations:

```text
canonical nucleotide          2 bits/base
exact mixed nucleotide        5 bits/symbol
query-mode mixed target       2-bit bases + validity/gap bitmaps
both-mode mixed target        5-bit IUPAC masks
protein                       5 bits/residue
```

Queries remain unpacked because they are usually short and faster to access as
one code or mask per byte.

## Features

- exact nucleotide search by default
- configurable IUPAC ambiguity: `none`, `query`, or `both`
- exact amino-acid search
- reverse-complement nucleotide search
- circular search, including patterns longer than the record
- FASTA and multiline FASTQ input
- plain text and gzip input
- overlapping matches
- packed shared-memory multiprocessing
- compact worker result transfer for repetitive queries
- TSV output with 1-based inclusive coordinates

## Usage

```text
seqgrep [-h]
        [-t {nucleotide,amino-acid}]
        [--ambig-mode {none,query,both}]
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

Each match is a TSV row:

```text
record  strand  start  end  matched  circular
```

Coordinates are 1-based and inclusive. For a circular match crossing the
record boundary, `end` wraps and `circular` is `true`.

## Human chromosome smoke test

Extract a canonical pattern from a chromosome instead of sampling an assembly
`N` block:

```bash
pattern=$(
  gzip -cd chr21.fa.gz |
  awk '!/^>/{printf "%s", $0}' |
  tr '[:lower:]' '[:upper:]' |
  grep -oE '[ACGT]{32}' |
  head -n 1
)
```

Compare serial and multiprocessing output:

```bash
seqgrep "$pattern" chr21.fa.gz --jobs 1 > serial.tsv
seqgrep "$pattern" chr21.fa.gz --jobs 4 --chunk-size 5000000 > parallel.tsv
diff -u serial.tsv parallel.tsv
```

Test query-side ambiguity while ignoring target `N` regions:

```bash
ambig="${pattern:0:8}NNNN${pattern:12}"
seqgrep "$ambig" chr21.fa.gz \
  --ambig-mode query \
  --jobs 4 \
  --chunk-size 5000000
```

Test symmetric ambiguity explicitly:

```bash
seqgrep AAAAA chr21.fa.gz --ambig-mode both
```

Be aware that this may report many possible matches inside unknown target
regions.

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
exact.py      native exact serial search
window.py     serial IUPAC compatibility search
chunked.py    packed shared-memory multiprocessing search
planner.py    backend selection
```
