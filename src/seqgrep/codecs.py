from __future__ import annotations

from typing import Protocol

from .iupac import encode_iupac, normalize_sequence, reverse_complement


class SequenceCodec(Protocol):
    """Encode symbols and define the matching rule.

    The window scanner does not know whether it is searching DNA, protein,
    literal text, or IUPAC masks. It only asks the codec how to normalize,
    encode, compare, and reverse-complement a pattern.
    """
    def normalize(self, seq: str) -> str:
        ...

    def encode(self, seq: str) -> bytes:
        ...

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        ...

    def reverse_complement(self, pattern: str) -> str:
        ...


class LiteralCodec:
    """Literal character matching.

    This is the default mode. It is safe for proteins because N, R, Y, X, B,
    and Z are treated as normal characters instead of DNA ambiguity symbols.
    """
    def normalize(self, seq: str) -> str:
        return normalize_sequence(seq)

    def encode(self, seq: str) -> bytes:
        normalized = self.normalize(seq)
        try:
            return normalized.encode("ascii")
        except UnicodeEncodeError as exc:
            raise ValueError(
                "Literal search currently supports ASCII sequence symbols only"
            ) from exc

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return query_symbol == target_symbol

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement(pattern)


class IupacDnaCodec:
    """IUPAC nucleotide ambiguity matching using bit masks."""
    def normalize(self, seq: str) -> str:
        return normalize_sequence(seq)

    def encode(self, seq: str) -> bytes:
        return bytes(encode_iupac(seq))

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return (query_symbol & target_symbol) != 0

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement(pattern)
