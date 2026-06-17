from __future__ import annotations

import gzip
from collections.abc import Iterable
from pathlib import Path
from typing import Literal, TextIO

from .models import FastaRecord

SequenceFormat = Literal["auto", "fasta", "fastq"]

_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".ffn", ".frn"}
_FASTQ_SUFFIXES = {".fq", ".fastq"}


class FastaFileReader:
    """Read FASTA or FASTQ records from plain text or gzip-compressed files.

    The class keeps the old name for compatibility. Internally, seqgrep only
    needs a record name and sequence, so FASTQ quality strings are validated
    enough to keep parsing sane, then discarded.
    """

    def __init__(self, path: Path | str, fmt: SequenceFormat = "auto") -> None:
        self.path = Path(path)
        self.fmt = fmt

    def read(self) -> Iterable[FastaRecord]:
        fmt = infer_format(self.path, self.fmt)
        with open_text(self.path) as handle:
            if fmt == "fasta":
                yield from parse_fasta(handle)
            elif fmt == "fastq":
                yield from parse_fastq(handle)
            else:
                raise ValueError(f"Unsupported sequence format: {fmt}")


SequenceFileReader = FastaFileReader


def open_text(path: Path | str) -> TextIO:
    path = Path(path)
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def infer_format(path: Path | str, fmt: SequenceFormat = "auto") -> Literal["fasta", "fastq"]:
    if fmt != "auto":
        return fmt

    path = Path(path)
    suffixes = [suffix.lower() for suffix in path.suffixes]

    if suffixes and suffixes[-1] == ".gz":
        suffixes = suffixes[:-1]

    if not suffixes:
        raise ValueError(
            "Could not infer input format from file extension. Use --format fasta or --format fastq."
        )

    suffix = suffixes[-1].lower()
    if suffix in _FASTA_SUFFIXES:
        return "fasta"
    if suffix in _FASTQ_SUFFIXES:
        return "fastq"

    raise ValueError(
        f"Could not infer input format from extension {suffix!r}. "
        "Use --format fasta or --format fastq."
    )


def parse_fasta(handle: TextIO) -> Iterable[FastaRecord]:
    name: str | None = None
    chunks: list[str] = []

    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if name is not None:
                yield FastaRecord(name=name, sequence="".join(chunks))
            name = line[1:].split(maxsplit=1)[0]
            chunks = []
            continue

        if name is None:
            raise ValueError("FASTA sequence appeared before the first header")
        chunks.append(line)

    if name is not None:
        yield FastaRecord(name=name, sequence="".join(chunks))


def parse_fastq(handle: TextIO) -> Iterable[FastaRecord]:
    """Parse FASTQ records, including multiline sequence and quality blocks."""

    line_number = 0

    while True:
        header = _next_nonempty_line(handle)
        if header is None:
            return

        line_number += 1
        if not header.startswith("@"):
            raise ValueError(f"FASTQ record must start with '@' near line {line_number}")

        name = header[1:].split(maxsplit=1)[0]
        seq_lines: list[str] = []

        for raw_line in handle:
            line_number += 1
            line = raw_line.rstrip("\n\r")
            if line.startswith("+"):
                break
            if line:
                seq_lines.append(line.strip())
        else:
            raise ValueError(f"FASTQ record {name!r} is missing '+' quality header")

        sequence = "".join(seq_lines)
        if not sequence:
            raise ValueError(f"FASTQ record {name!r} has an empty sequence")

        quality_length = 0
        while quality_length < len(sequence):
            quality_line = handle.readline()
            if quality_line == "":
                raise ValueError(f"FASTQ record {name!r} ended before quality block was complete")
            line_number += 1
            quality_length += len(quality_line.rstrip("\n\r"))

        if quality_length != len(sequence):
            raise ValueError(
                f"FASTQ record {name!r} has quality length {quality_length}, "
                f"expected {len(sequence)}"
            )

        yield FastaRecord(name=name, sequence=sequence)


def _next_nonempty_line(handle: TextIO) -> str | None:
    for raw_line in handle:
        line = raw_line.rstrip("\n\r")
        if line:
            return line
    return None
