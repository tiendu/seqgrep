from __future__ import annotations

from dataclasses import dataclass

from .chunked import ChunkedProcessMatcher
from .codecs import IupacNucleotideCodec, NucleotideCodec, ProteinCodec, SequenceCodec
from .models import SearchQuery, SequenceMatcher, SequenceType
from .window import WindowMatcher


@dataclass(frozen=True)
class SearchPlan:
    matcher: SequenceMatcher


class SearchPlanner:
    def plan(self, query: SearchQuery, jobs: int, chunk_size: int) -> SearchPlan:
        codec = self._select_codec(query)

        if jobs > 1:
            return SearchPlan(
                matcher=ChunkedProcessMatcher(
                    codec=codec,
                    workers=jobs,
                    chunk_size=chunk_size,
                )
            )

        return SearchPlan(matcher=WindowMatcher(codec))

    @staticmethod
    def _select_codec(query: SearchQuery) -> SequenceCodec:
        if query.sequence_type == SequenceType.AMINO_ACID:
            if query.ambig:
                raise ValueError("--ambig is only valid for nucleotide sequences")
            if query.revcomp:
                raise ValueError("--revcomp is only valid for nucleotide sequences")
            return ProteinCodec()

        if query.ambig:
            return IupacNucleotideCodec()

        return NucleotideCodec()
