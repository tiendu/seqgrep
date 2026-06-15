from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="seqgrep",
        description="Sequence grep with literal search by default and optional IUPAC DNA ambiguity.",
    )

    parser.add_argument("pattern", help="Pattern to search, e.g. ATGNNRY or MTEYK")
    parser.add_argument("input", type=Path, help="Input FASTA/FASTQ file, optionally .gz")

    parser.add_argument(
        "--ambig",
        action="store_true",
        help="Interpret pattern and target using IUPAC nucleotide ambiguity codes.",
    )
    parser.add_argument(
        "--revcomp",
        action="store_true",
        help="Also search reverse-complement pattern. Only valid for nucleotide patterns.",
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
    return parser.parse_args(argv)
