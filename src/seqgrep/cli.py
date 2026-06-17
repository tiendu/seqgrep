from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .models import SequenceType


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="seqgrep",
        description=(
            "Search nucleotide or amino-acid FASTA/FASTQ records. "
            "The default is exact nucleotide matching."
        ),
    )

    parser.add_argument("pattern", help="Sequence pattern, e.g. ATGNNRY or MTEYK")
    parser.add_argument("input", type=Path, help="FASTA/FASTQ input, optionally gzip")
    parser.add_argument(
        "-t",
        "--sequence-type",
        choices=tuple(item.value for item in SequenceType),
        default=SequenceType.NUCLEOTIDE.value,
        help="Sequence alphabet. Default: nucleotide.",
    )
    parser.add_argument(
        "--ambig",
        action="store_true",
        help="Enable IUPAC nucleotide compatibility matching.",
    )
    parser.add_argument(
        "--revcomp",
        action="store_true",
        help="Also search the reverse complement; nucleotide mode only.",
    )
    parser.add_argument(
        "--circular",
        action="store_true",
        help="Allow matches to cross the sequence boundary.",
    )
    parser.add_argument("--with-header", action="store_true", help="Print a TSV header.")
    parser.add_argument(
        "--format",
        choices=("auto", "fasta", "fastq"),
        default="auto",
        help="Input format. Default: infer from the filename.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=1,
        help="Worker processes per record. Default: 1.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Candidate starts per multiprocessing job.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

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
