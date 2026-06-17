from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from enum import Enum
from typing import Protocol

from .alphabets import normalize_symbols


class SequenceType(str, Enum):
    NUCLEOTIDE = "nucleotide"
    AMINO_ACID = "amino-acid"


class AmbigMode(str, Enum):
    """How IUPAC nucleotide ambiguity is interpreted."""

    NONE = "none"
    QUERY = "query"
    BOTH = "both"


@dataclass(frozen=True, slots=True)
class SequenceRecord:
    name: str
    sequence: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "sequence", normalize_symbols(self.sequence))


@dataclass(frozen=True, slots=True)
class SearchQuery:
    pattern: str
    revcomp: bool = False
    circular: bool = False
    ambig: bool = False
    sequence_type: SequenceType = SequenceType.NUCLEOTIDE
    ambig_mode: AmbigMode = AmbigMode.NONE

    def __post_init__(self) -> None:
        object.__setattr__(self, "pattern", normalize_symbols(self.pattern))

        mode = AmbigMode(self.ambig_mode)
        if self.ambig:
            if mode is AmbigMode.BOTH:
                raise ValueError(
                    "ambig=True is an alias for ambig_mode='query' and cannot "
                    "be combined with ambig_mode='both'"
                )
            mode = AmbigMode.QUERY

        object.__setattr__(self, "ambig_mode", mode)
        # Keep the old boolean attribute meaningful for callers that still
        # inspect it. New code should use ambig_mode when the distinction
        # between query-only and symmetric ambiguity matters.
        object.__setattr__(self, "ambig", mode is not AmbigMode.NONE)


@dataclass(frozen=True, slots=True)
class Match:
    record: str
    strand: str
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive, may wrap when circular=True
    matched: str
    circular: bool = False


class SequenceReader(Protocol):
    def read(self) -> Iterable[SequenceRecord]:
        """Yield sequence records."""


class SequenceMatcher(Protocol):
    def search(self, record: SequenceRecord, query: SearchQuery) -> Iterable[Match]:
        """Yield matches for a query in one sequence record."""
