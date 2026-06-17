from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from typing import Protocol


class SequenceType(str, Enum):
    NUCLEOTIDE = "nucleotide"
    AMINO_ACID = "amino-acid"


@dataclass(frozen=True)
class FastaRecord:
    name: str
    sequence: str


@dataclass(frozen=True)
class SearchQuery:
    pattern: str
    revcomp: bool = False
    circular: bool = False
    ambig: bool = False
    sequence_type: SequenceType = SequenceType.NUCLEOTIDE


@dataclass(frozen=True)
class Match:
    record: str
    strand: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive, may wrap when circular=True
    matched: str
    circular: bool = False


class SequenceReader(Protocol):
    def read(self) -> Iterable[FastaRecord]:
        """Yield sequence records."""


class SequenceMatcher(Protocol):
    def search(self, record: FastaRecord, query: SearchQuery) -> Iterable[Match]:
        """Yield matches for a query in one sequence record."""
