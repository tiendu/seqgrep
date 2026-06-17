from __future__ import annotations

from collections.abc import Iterable, Iterator

from .codecs import ExactSequenceCodec
from .models import FastaRecord, Match, SearchQuery


class ExactMatcher:
    """Fast serial exact matcher backed by Python's native string search."""

    def __init__(self, codec: ExactSequenceCodec) -> None:
        self.codec = codec

    def search(self, record: FastaRecord, query: SearchQuery) -> Iterable[Match]:
        sequence = record.sequence
        pattern = query.pattern

        if not pattern:
            raise ValueError("Pattern must not be empty")
        if not sequence:
            return

        sequence_key = self.codec.comparison_text(sequence)
        pattern_key = self.codec.comparison_text(pattern)

        yield from self._search_one(
            record_name=record.name,
            sequence=sequence,
            sequence_key=sequence_key,
            pattern_key=pattern_key,
            strand="+",
            circular=query.circular,
        )

        if query.revcomp:
            reverse_pattern = self.codec.reverse_complement(pattern)
            yield from self._search_one(
                record_name=record.name,
                sequence=sequence,
                sequence_key=sequence_key,
                pattern_key=self.codec.comparison_text(reverse_pattern),
                strand="-",
                circular=query.circular,
            )

    def _search_one(
        self,
        *,
        record_name: str,
        sequence: str,
        sequence_key: str,
        pattern_key: str,
        strand: str,
        circular: bool,
    ) -> Iterable[Match]:
        sequence_length = len(sequence_key)
        pattern_length = len(pattern_key)

        if pattern_length > sequence_length and not circular:
            return

        for zero_start in self._starts(
            sequence_key=sequence_key,
            pattern_key=pattern_key,
            circular=circular,
        ):
            zero_end = zero_start + pattern_length - 1
            wraps = circular and zero_end >= sequence_length

            yield Match(
                record=record_name,
                strand=strand,
                start=zero_start + 1,
                end=(zero_end % sequence_length) + 1 if circular else zero_end + 1,
                matched=self._matched_sequence(
                    sequence=sequence,
                    start=zero_start,
                    length=pattern_length,
                    circular=circular,
                ),
                circular=wraps,
            )

    @staticmethod
    def _starts(
        *,
        sequence_key: str,
        pattern_key: str,
        circular: bool,
    ) -> Iterator[int]:
        sequence_length = len(sequence_key)

        if circular:
            required_length = sequence_length + len(pattern_key) - 1
            repeat_count = (required_length + sequence_length - 1) // sequence_length
            haystack = (sequence_key * repeat_count)[:required_length]
            maximum_start = sequence_length
        else:
            haystack = sequence_key
            maximum_start = sequence_length - len(pattern_key) + 1

        cursor = 0
        while cursor < maximum_start:
            found = haystack.find(pattern_key, cursor)
            if found < 0 or found >= maximum_start:
                return
            yield found
            cursor = found + 1

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
