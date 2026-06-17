from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from .models import SequenceType


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="seqgrep",
        description=(
            "Search nucleotide or amino-acid sequences with compact target "
            "encoding. Nucleotide mode is the default."
        ),
    )

    parser.add_argument("pattern", help="Pattern to search, e.g. ATGNNRY or MTEYK")
    parser.add_argument("input", type=Path, help="Input FASTA/FASTQ file, optionally .gz")

    parser.add_argument(
        "-t",
        "--sequence-type",
        choices=tuple(sequence_type.value for sequence_type in SequenceType),
        default=SequenceType.NUCLEOTIDE.value,
        help=(
            "Biological sequence type. Default: nucleotide. Amino-acid mode "
            "uses exact five-bit protein encoding."
        ),
    )
    parser.add_argument(
        "--ambig",
        action="store_true",
        help=("Enable IUPAC ambiguity matching. Valid only with --sequence-type nucleotide."),
    )
    parser.add_argument(
        "--revcomp",
        action="store_true",
        help=(
            "Also search the reverse-complement pattern. Valid only with "
            "--sequence-type nucleotide."
        ),
    )
    parser.add_argument(
        "--circular",
        action="store_true",
        help="Allow matches crossing the sequence boundary.",
    )
    parser.add_argument(
        "--with-header",
        action="store_true",
        help="Print TSV header.",
    )
    parser.add_argument(
        "--format",
        choices=("auto", "fasta", "fastq"),
        default="auto",
        help="Input format. Default: infer from extension.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Number of worker processes for each sequence. Use 1 for serial mode.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Number of candidate start positions per worker job when --jobs > 1.",
    )

    return parser


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.sequence_type == SequenceType.AMINO_ACID.value:
        if args.ambig:
            parser.error("--ambig is only valid for nucleotide sequences")
        if args.revcomp:
            parser.error("--revcomp is only valid for nucleotide sequences")

    if args.jobs < 1:
        parser.error("--jobs must be at least 1")
    if args.chunk_size < 1:
        parser.error("--chunk-size must be at least 1")

    return args
