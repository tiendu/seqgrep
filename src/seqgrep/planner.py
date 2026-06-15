from __future__ import annotations

from dataclasses import dataclass

from .chunked import ChunkedProcessMatcher
from .codecs import IupacDnaCodec, LiteralCodec
from .models import SearchQuery, SequenceMatcher
from .window import WindowMatcher


@dataclass(frozen=True)
class SearchPlan:
    matcher: SequenceMatcher


class SearchPlanner:
    def plan(self, query: SearchQuery, jobs: int, chunk_size: int) -> SearchPlan:
        codec = IupacDnaCodec() if query.ambig else LiteralCodec()

        if jobs > 1:
            return SearchPlan(
                matcher=ChunkedProcessMatcher(
                    codec=codec,
                    workers=jobs,
                    chunk_size=chunk_size,
                )
            )

        return SearchPlan(matcher=WindowMatcher(codec))
