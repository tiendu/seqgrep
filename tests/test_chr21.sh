#!/usr/bin/env bash
set -euo pipefail

# Real-data integration test for seqgrep using GRCh38 chromosome 21.
#
# Override the executable when testing a local installation:
#
#   SEQGREP_BIN=target/release/seqgrep tests/test_chr21.sh
#
# The chromosome is cached under tests/data/ and is not removed after the test.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
FASTA="${DATA_DIR}/chr21.fa.gz"

NCBI_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000021.9&rettype=fasta&retmode=text"
SEQGREP_BIN="${SEQGREP_BIN:-seqgrep}"

SEARCH_START="${SEARCH_START:-15000001}"
PATTERN_LENGTH="${PATTERN_LENGTH:-32}"
SCAN_WINDOW="${SCAN_WINDOW:-100000}"
JOBS="${JOBS:-4}"
CHUNK_SIZE="${CHUNK_SIZE:-5000000}"

require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        printf '[x] Required command not found: %s\n' "$1" >&2
        exit 1
    fi
}

for command in curl gzip awk diff grep mktemp mv rm; do
    require_command "$command"
done

if [[ ! -x "$SEQGREP_BIN" ]] && ! command -v "$SEQGREP_BIN" >/dev/null 2>&1; then
    printf '[x] seqgrep executable not found: %s\n' "$SEQGREP_BIN" >&2
    printf '    Set SEQGREP_BIN to an executable path, for example:\n' >&2
    printf '    SEQGREP_BIN=target/release/seqgrep %s\n' "$0" >&2
    exit 1
fi

mkdir -p "$DATA_DIR"

download_chr21() {
    local plain_partial="${FASTA}.txt.partial"
    local gzip_partial="${FASTA}.partial"

    printf '[*] Downloading GRCh38 chromosome 21 from NCBI...\n'
    rm -f "$plain_partial" "$gzip_partial"

    if ! curl \
        --fail \
        --location \
        --retry 3 \
        --retry-delay 2 \
        --connect-timeout 15 \
        --output "$plain_partial" \
        "$NCBI_URL"
    then
        printf '[x] Failed to download chromosome 21 from NCBI.\n' >&2
        rm -f "$plain_partial" "$gzip_partial"
        return 1
    fi

    if ! awk '
        NR == 1 {
            if ($0 !~ /^>/) {
                exit 1
            }

            print ">chr21"
            next
        }

        {
            print
        }
    ' "$plain_partial" | gzip -c >"$gzip_partial"
    then
        printf '[x] NCBI response was not valid FASTA or compression failed.\n' >&2
        rm -f "$plain_partial" "$gzip_partial"
        return 1
    fi

    rm -f "$plain_partial"

    if ! gzip -t "$gzip_partial"; then
        printf '[x] Downloaded chromosome archive failed gzip validation.\n' >&2
        rm -f "$gzip_partial"
        return 1
    fi

    mv "$gzip_partial" "$FASTA"
    printf '[v] Downloaded and validated: %s\n' "$FASTA"
}

if [[ ! -f "$FASTA" ]]; then
    download_chr21
elif ! gzip -t "$FASTA" 2>/dev/null; then
    printf '[!] Cached FASTA is invalid; downloading it again.\n'
    rm -f "$FASTA"
    download_chr21
else
    printf '[v] Using cached chromosome: %s\n' "$FASTA"
fi

WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/seqgrep-chr21.XXXXXX")"
trap 'rm -rf "$WORK_DIR"' EXIT

find_canonical_pattern() {
    gzip -cd "$FASTA" |
        awk \
            -v search_start="$SEARCH_START" \
            -v pattern_length="$PATTERN_LENGTH" \
            -v scan_window="$SCAN_WINDOW" '
                /^>/ {
                    next
                }

                found {
                    next
                }

                {
                    line = toupper($0)
                    line_start = total_length + 1
                    line_end = total_length + length(line)
                    wanted_end = search_start + scan_window + pattern_length - 2

                    if (line_end >= search_start && line_start <= wanted_end) {
                        first = search_start > line_start \
                            ? search_start - line_start + 1 \
                            : 1

                        last = wanted_end < line_end \
                            ? wanted_end - line_start + 1 \
                            : length(line)

                        buffer = buffer substr(line, first, last - first + 1)
                    }

                    total_length = line_end

                    if (total_length >= wanted_end) {
                        for (offset = 1; offset <= length(buffer) - pattern_length + 1; offset++) {
                            candidate = substr(buffer, offset, pattern_length)

                            if (candidate !~ /[^ACGT]/) {
                                print search_start + offset - 1 "\t" candidate
                                found = 1
                                break
                            }
                        }
                    }
                }
            '
}

pattern_record="$(find_canonical_pattern)"

if [[ -z "$pattern_record" ]]; then
    printf \
        '[x] No canonical %s-mer found within %s bases of position %s.\n' \
        "$PATTERN_LENGTH" \
        "$SCAN_WINDOW" \
        "$SEARCH_START" \
        >&2
    exit 1
fi

IFS=$'\t' read -r expected_start pattern <<<"$pattern_record"
expected_end=$((expected_start + PATTERN_LENGTH - 1))

printf '[v] Selected canonical pattern\n'
printf '    chromosome: chr21\n'
printf '    coordinates: %s-%s\n' "$expected_start" "$expected_end"
printf '    pattern: %s\n' "$pattern"

serial_exact="${WORK_DIR}/serial-exact.tsv"
parallel_exact="${WORK_DIR}/parallel-exact.tsv"

printf '[*] Testing exact nucleotide search...\n'

"$SEQGREP_BIN" \
    "$pattern" \
    "$FASTA" \
    --jobs 1 \
    >"$serial_exact"

"$SEQGREP_BIN" \
    "$pattern" \
    "$FASTA" \
    --jobs "$JOBS" \
    --chunk-size "$CHUNK_SIZE" \
    >"$parallel_exact"

diff -u "$serial_exact" "$parallel_exact"

expected_exact="$(
    printf \
        'chr21\t+\t%s\t%s\t%s\tfalse' \
        "$expected_start" \
        "$expected_end" \
        "$pattern"
)"

if ! grep -Fqx "$expected_exact" "$serial_exact"; then
    printf '[x] Exact search did not report the selected chromosome window.\n' >&2
    printf '    Expected: %s\n' "$expected_exact" >&2
    exit 1
fi

printf '[v] Exact serial/parallel outputs match.\n'

ambig="${pattern:0:8}NNNN${pattern:12}"
serial_ambig="${WORK_DIR}/serial-ambig.tsv"
parallel_ambig="${WORK_DIR}/parallel-ambig.tsv"

printf '[*] Testing query-side IUPAC ambiguity...\n'
printf '    pattern: %s\n' "$ambig"

"$SEQGREP_BIN" \
    "$ambig" \
    "$FASTA" \
    --ambig-mode query \
    --jobs 1 \
    >"$serial_ambig"

"$SEQGREP_BIN" \
    "$ambig" \
    "$FASTA" \
    --ambig-mode query \
    --jobs "$JOBS" \
    --chunk-size "$CHUNK_SIZE" \
    >"$parallel_ambig"

diff -u "$serial_ambig" "$parallel_ambig"

if ! grep -Fq "$(
    printf 'chr21\t+\t%s\t%s\t' "$expected_start" "$expected_end"
)" "$serial_ambig"; then
    printf '[x] Ambiguous query did not recover the selected chromosome window.\n' >&2
    exit 1
fi

printf '[v] Query-mode serial/parallel outputs match.\n'

printf '[*] Testing --ambig compatibility alias...\n'

alias_fasta="${WORK_DIR}/alias.fa"
alias_short="${WORK_DIR}/alias-short.tsv"
alias_explicit="${WORK_DIR}/alias-explicit.tsv"

printf '>alias\nAACGTAAA\n' >"$alias_fasta"

"$SEQGREP_BIN" \
    ANGT \
    "$alias_fasta" \
    --ambig \
    >"$alias_short"

"$SEQGREP_BIN" \
    ANGT \
    "$alias_fasta" \
    --ambig-mode query \
    >"$alias_explicit"

diff -u "$alias_short" "$alias_explicit"

expected_alias=$'alias\t+\t2\t5\tACGT\tfalse'

if ! grep -Fqx "$expected_alias" "$alias_short"; then
    printf '[x] --ambig did not produce the expected query-mode match.\n' >&2
    cat "$alias_short" >&2
    exit 1
fi

printf '[v] --ambig remains equivalent to --ambig-mode query.\n'

printf '[*] Testing target ambiguity semantics on a tiny fixture...\n'

tiny_fasta="${WORK_DIR}/ambiguous-target.fa"
printf '>ambiguous\nNNANN\n' >"$tiny_fasta"

query_mode_output="${WORK_DIR}/tiny-query.tsv"
both_mode_output="${WORK_DIR}/tiny-both.tsv"

"$SEQGREP_BIN" \
    AAAAA \
    "$tiny_fasta" \
    --ambig-mode query \
    >"$query_mode_output"

"$SEQGREP_BIN" \
    AAAAA \
    "$tiny_fasta" \
    --ambig-mode both \
    >"$both_mode_output"

if [[ -s "$query_mode_output" ]]; then
    printf '[x] Query mode incorrectly matched ambiguous target bases.\n' >&2
    cat "$query_mode_output" >&2
    exit 1
fi

expected_both=$'ambiguous\t+\t1\t5\tNNANN\tfalse'

if ! grep -Fqx "$expected_both" "$both_mode_output"; then
    printf '[x] Both-side ambiguity did not match AAAAA against NNANN.\n' >&2
    cat "$both_mode_output" >&2
    exit 1
fi

printf '[v] Query mode rejects unknown target bases.\n'
printf '[v] Both mode accepts compatible ambiguous target bases.\n'
printf '[v] Chromosome 21 integration test passed.\n'
