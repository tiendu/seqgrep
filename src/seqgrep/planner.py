from __future__ import annotations

from dataclasses import dataclass

from .chunked import ChunkedProcessMatcher
from .codecs import IupacNucleotideCodec, NucleotideCodec, ProteinCodec
from .exact import ExactMatcher
from .models import SearchQuery, SequenceMatcher, SequenceType
from .window import WindowMatcher


@dataclass(frozen=True, slots=True)
class SearchPlan:
    matcher: SequenceMatcher


class SearchPlanner:
    """Select matching semantics and execution backend."""

    def plan(self, query: SearchQuery, jobs: int, chunk_size: int) -> SearchPlan:
        if jobs < 1:
            raise ValueError("jobs must be at least 1")
        if chunk_size < 1:
            raise ValueError("chunk_size must be at least 1")

        if query.sequence_type is SequenceType.AMINO_ACID:
            if query.ambig:
                raise ValueError("--ambig is only valid for nucleotide sequences")
            if query.revcomp:
                raise ValueError("--revcomp is only valid for nucleotide sequences")

            protein_codec = ProteinCodec()
            if jobs > 1:
                return SearchPlan(ChunkedProcessMatcher(protein_codec, jobs, chunk_size))
            return SearchPlan(ExactMatcher(protein_codec))

        if query.ambig:
            iupac_codec = IupacNucleotideCodec()
            if jobs > 1:
                return SearchPlan(ChunkedProcessMatcher(iupac_codec, jobs, chunk_size))
            return SearchPlan(WindowMatcher(iupac_codec))

        nucleotide_codec = NucleotideCodec()
        if jobs > 1:
            return SearchPlan(ChunkedProcessMatcher(nucleotide_codec, jobs, chunk_size))
        return SearchPlan(ExactMatcher(nucleotide_codec))
