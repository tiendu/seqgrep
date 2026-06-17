from __future__ import annotations

from collections.abc import Iterable

from .codecs import EncodedTarget, SequenceCodec
from .models import FastaRecord, Match, SearchQuery


class WindowMatcher:
    """Serial compatibility matcher for non-exact symbol semantics."""

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
        *,
        record: FastaRecord,
        pattern: str,
        strand: str,
        circular: bool,
    ) -> Iterable[Match]:
        sequence = self.codec.normalize(record.sequence)
        normalized_pattern = self.codec.normalize(pattern)

        if not normalized_pattern:
            raise ValueError("Pattern must not be empty")
        if not sequence:
            return

        query_symbols = self.codec.encode_query(normalized_pattern)
        target = self.codec.encode_target(sequence)
        sequence_length = len(target)
        pattern_length = len(query_symbols)

        if pattern_length > sequence_length and not circular:
            return

        total_starts = sequence_length if circular else sequence_length - pattern_length + 1

        for zero_start in range(total_starts):
            if not self._window_matches(
                query_symbols=query_symbols,
                target=target,
                start=zero_start,
                circular=circular,
            ):
                continue

            zero_end = zero_start + pattern_length - 1
            yield Match(
                record=record.name,
                strand=strand,
                start=zero_start + 1,
                end=(zero_end % sequence_length) + 1 if circular else zero_end + 1,
                matched=self._matched_sequence(
                    sequence=sequence,
                    start=zero_start,
                    length=pattern_length,
                    circular=circular,
                ),
                circular=circular and zero_end >= sequence_length,
            )

    def _window_matches(
        self,
        *,
        query_symbols: bytes,
        target: EncodedTarget,
        start: int,
        circular: bool,
    ) -> bool:
        compatible = self.codec.compatible
        symbol_at = target.symbol_at

        if circular:
            sequence_length = len(target)
            for offset, query_symbol in enumerate(query_symbols):
                if not compatible(
                    query_symbol,
                    symbol_at((start + offset) % sequence_length),
                ):
                    return False
            return True

        for offset, query_symbol in enumerate(query_symbols):
            if not compatible(query_symbol, symbol_at(start + offset)):
                return False

        return True

    @staticmethod
    def _matched_sequence(
        *,
        sequence: str,
        start: int,
        length: int,
        circular: bool,
    ) -> str:
        if not circular:
            return sequence[start : start + length]

        sequence_length = len(sequence)
        return "".join(sequence[(start + offset) % sequence_length] for offset in range(length))
