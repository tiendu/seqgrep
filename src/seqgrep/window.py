from __future__ import annotations

from collections.abc import Iterable

from .codecs import EncodedTarget, SequenceCodec
from .models import FastaRecord, Match, SearchQuery


class WindowMatcher:
    def __init__(self, codec: SequenceCodec) -> None:
        self.codec = codec

    def search(self, record: FastaRecord, query: SearchQuery) -> Iterable[Match]:
        yield from self._search_one(
            record=record,
            pattern=query.pattern,
            strand="+",
            circular=query.circular,
        )

        if query.revcomp:
            yield from self._search_one(
                record=record,
                pattern=self.codec.reverse_complement(query.pattern),
                strand="-",
                circular=query.circular,
            )

    def _search_one(
        self,
        record: FastaRecord,
        pattern: str,
        strand: str,
        circular: bool,
    ) -> Iterable[Match]:
        seq = self.codec.normalize(record.sequence)
        normalized_pattern = self.codec.normalize(pattern)

        if not normalized_pattern:
            raise ValueError("Pattern must not be empty")

        if not seq:
            return

        query_symbols = self.codec.encode_query(normalized_pattern)
        target_symbols = self.codec.encode_target(seq)

        seq_len = len(target_symbols)
        pattern_len = len(query_symbols)

        if pattern_len > seq_len and not circular:
            return

        max_starts = (
            seq_len
            if circular
            else seq_len - pattern_len + 1
        )

        for zero_start in range(max_starts):
            if not self._window_matches(
                query_symbols=query_symbols,
                target_symbols=target_symbols,
                start=zero_start,
                circular=circular,
            ):
                continue

            zero_end = zero_start + pattern_len - 1
            wraps = circular and zero_end >= seq_len

            yield Match(
                record=record.name,
                strand=strand,
                start=zero_start + 1,
                end=(zero_end % seq_len) + 1 if circular else zero_end + 1,
                matched=self._matched_sequence(
                    seq=seq,
                    start=zero_start,
                    length=pattern_len,
                    circular=circular,
                ),
                circular=wraps,
            )

    def _window_matches(
        self,
        query_symbols: bytes,
        target_symbols: EncodedTarget,
        start: int,
        circular: bool,
    ) -> bool:
        compatible = self.codec.compatible
        symbol_at = target_symbols.symbol_at

        if circular:
            seq_len = len(target_symbols)

            for offset, query_symbol in enumerate(query_symbols):
                target_index = (start + offset) % seq_len

                if not compatible(
                    query_symbol,
                    symbol_at(target_index),
                ):
                    return False

            return True

        for offset, query_symbol in enumerate(query_symbols):
            if not compatible(
                query_symbol,
                symbol_at(start + offset),
            ):
                return False

        return True

    @staticmethod
    def _matched_sequence(
        seq: str,
        start: int,
        length: int,
        circular: bool,
    ) -> str:
        if not circular:
            return seq[start : start + length]

        seq_len = len(seq)

        return "".join(
            seq[(start + offset) % seq_len]
            for offset in range(length)
        )

