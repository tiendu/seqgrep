from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Protocol

from .iupac import (
    encode_canonical_2bit,
    encode_iupac_query,
    normalize_sequence,
    reverse_complement,
)


class TargetEncoding(str, Enum):
    """Physical representation used for an encoded target sequence."""

    BYTE = "byte"
    PACKED_CANONICAL_2BIT = "packed-canonical-2bit"


class EncodedTarget(Protocol):
    """Random-access encoded target sequence.

    ``data`` contains the physical bytes suitable for shared memory.
    ``len(target)`` returns the logical number of sequence symbols, which may
    differ from ``len(target.data)`` for packed representations.
    """

    @property
    def data(self) -> bytes:
        ...

    @property
    def encoding(self) -> TargetEncoding:
        ...

    def __len__(self) -> int:
        ...

    def symbol_at(self, index: int) -> int:
        ...


class SequenceCodec(Protocol):
    """Prepare query and target symbols and define their matching rule."""

    def normalize(self, seq: str) -> str:
        ...

    def encode_query(self, seq: str) -> bytes:
        ...

    def encode_target(self, seq: str) -> EncodedTarget:
        ...

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        ...

    def reverse_complement(self, pattern: str) -> str:
        ...


@dataclass(frozen=True, slots=True)
class ByteTarget:
    """One byte per logical target symbol.

    Used for literal ASCII sequences and ambiguous IUPAC targets.
    """

    data: bytes

    @property
    def encoding(self) -> TargetEncoding:
        return TargetEncoding.BYTE

    def __len__(self) -> int:
        return len(self.data)

    def symbol_at(self, index: int) -> int:
        return self.data[index]


@dataclass(frozen=True, slots=True)
class PackedCanonicalMaskTarget:
    """Canonical DNA packed using two bits per base.

    Storage uses A=00, C=01, G=10, and T/U=11. ``symbol_at`` converts the
    two-bit code into the singleton four-bit IUPAC mask expected by
    :class:`IupacDnaCodec`.
    """

    data: bytes
    length: int

    @property
    def encoding(self) -> TargetEncoding:
        return TargetEncoding.PACKED_CANONICAL_2BIT

    def __len__(self) -> int:
        return self.length

    def symbol_at(self, index: int) -> int:
        if index < 0 or index >= self.length:
            raise IndexError(f"Target index out of range: {index}")

        byte_index = index // 4
        shift = 6 - (index % 4) * 2
        code = (self.data[byte_index] >> shift) & 0b11

        return 1 << code


def _encode_ascii(seq: str) -> bytes:
    normalized = normalize_sequence(seq)

    try:
        return normalized.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError(
            "Literal search currently supports ASCII sequence symbols only"
        ) from exc


class LiteralCodec:
    """Exact byte-for-byte matching for DNA, proteins, or literal text."""

    def normalize(self, seq: str) -> str:
        return normalize_sequence(seq)

    def encode_query(self, seq: str) -> bytes:
        return _encode_ascii(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        return ByteTarget(_encode_ascii(seq))

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return query_symbol == target_symbol

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement(pattern)


class IupacDnaCodec:
    """IUPAC matching with adaptive target storage.

    Queries use one four-bit mask per byte. Canonical targets are packed using
    two bits per base; targets containing ambiguity symbols fall back to one
    four-bit mask per byte.
    """

    def normalize(self, seq: str) -> str:
        return normalize_sequence(seq)

    def encode_query(self, seq: str) -> bytes:
        return encode_iupac_query(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        try:
            packed, length = encode_canonical_2bit(seq)
        except ValueError:
            # The canonical encoder also rejects invalid symbols. The fallback
            # validates them against the full IUPAC alphabet and raises a clear
            # error if they are unsupported.
            return ByteTarget(encode_iupac_query(seq))

        return PackedCanonicalMaskTarget(data=packed, length=length)

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return (query_symbol & target_symbol) != 0

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement(pattern)
