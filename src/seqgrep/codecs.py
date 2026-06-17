from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Protocol

from .iupac import (
    encode_canonical_2bit,
    encode_iupac_query,
    encode_nucleotide_query,
    get_2bit_base,
    normalize_sequence,
    reverse_complement,
)
from .protein import (
    encode_protein_5bit,
    encode_protein_query,
    get_protein_5bit,
    normalize_protein,
)


class TargetEncoding(str, Enum):
    """Physical representation used for an encoded target sequence."""

    BYTE = "byte"
    PACKED_NUCLEOTIDE_2BIT = "packed-nucleotide-2bit"
    # Kept under its original name for compatibility. In IUPAC mode the same
    # two-bit storage is exposed as singleton masks rather than raw codes.
    PACKED_CANONICAL_2BIT = "packed-canonical-2bit"
    PACKED_PROTEIN_5BIT = "packed-protein-5bit"


class EncodedTarget(Protocol):
    """Random-access encoded target sequence.

    ``data`` contains physical bytes suitable for shared memory.
    ``len(target)`` returns the logical number of biological symbols.
    """

    @property
    def data(self) -> bytes: ...

    @property
    def encoding(self) -> TargetEncoding: ...

    def __len__(self) -> int: ...

    def symbol_at(self, index: int) -> int: ...


class SequenceCodec(Protocol):
    """Normalize, encode, and compare one biological sequence type."""

    def normalize(self, seq: str) -> str: ...

    def encode_query(self, seq: str) -> bytes: ...

    def encode_target(self, seq: str) -> EncodedTarget: ...

    def compatible(self, query_symbol: int, target_symbol: int) -> bool: ...

    def reverse_complement(self, pattern: str) -> str: ...


def decode_target_symbol(
    storage: bytes | memoryview,
    index: int,
    length: int,
    encoding: TargetEncoding,
) -> int:
    """Decode one logical symbol from bytes or shared-memory storage."""
    if index < 0 or index >= length:
        raise IndexError(f"Target index out of range: {index}")

    if encoding == TargetEncoding.BYTE:
        return storage[index]

    if encoding == TargetEncoding.PACKED_NUCLEOTIDE_2BIT:
        return get_2bit_base(storage, index, length)

    if encoding == TargetEncoding.PACKED_CANONICAL_2BIT:
        return 1 << get_2bit_base(storage, index, length)

    if encoding == TargetEncoding.PACKED_PROTEIN_5BIT:
        return get_protein_5bit(storage, index, length)

    raise ValueError(f"Unsupported target encoding: {encoding}")


@dataclass(frozen=True, slots=True)
class ByteTarget:
    """One byte per logical target symbol."""

    data: bytes

    @property
    def encoding(self) -> TargetEncoding:
        return TargetEncoding.BYTE

    def __len__(self) -> int:
        return len(self.data)

    def symbol_at(self, index: int) -> int:
        return decode_target_symbol(
            self.data,
            index,
            len(self),
            self.encoding,
        )


@dataclass(frozen=True, slots=True)
class PackedCanonicalCodeTarget:
    """Canonical nucleotide target packed at two bits per base.

    ``symbol_at`` returns raw codes 0..3 for exact nucleotide matching.
    """

    data: bytes
    length: int

    @property
    def encoding(self) -> TargetEncoding:
        return TargetEncoding.PACKED_NUCLEOTIDE_2BIT

    def __len__(self) -> int:
        return self.length

    def symbol_at(self, index: int) -> int:
        return decode_target_symbol(
            self.data,
            index,
            self.length,
            self.encoding,
        )


@dataclass(frozen=True, slots=True)
class PackedCanonicalMaskTarget:
    """Canonical nucleotide target packed at two bits per base.

    ``symbol_at`` returns singleton IUPAC masks for ambiguity matching.
    """

    data: bytes
    length: int

    @property
    def encoding(self) -> TargetEncoding:
        return TargetEncoding.PACKED_CANONICAL_2BIT

    def __len__(self) -> int:
        return self.length

    def symbol_at(self, index: int) -> int:
        return decode_target_symbol(
            self.data,
            index,
            self.length,
            self.encoding,
        )


@dataclass(frozen=True, slots=True)
class PackedProteinTarget:
    """Protein target packed at five bits per residue."""

    data: bytes
    length: int

    @property
    def encoding(self) -> TargetEncoding:
        return TargetEncoding.PACKED_PROTEIN_5BIT

    def __len__(self) -> int:
        return self.length

    def symbol_at(self, index: int) -> int:
        return decode_target_symbol(
            self.data,
            index,
            self.length,
            self.encoding,
        )


def _encode_ascii(seq: str) -> bytes:
    normalized = "".join(seq.split()).upper()

    try:
        return normalized.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValueError("Literal search currently supports ASCII sequence symbols only") from exc


class LiteralCodec:
    """Backward-compatible generic ASCII codec.

    The CLI no longer selects this codec automatically. Use nucleotide or
    amino-acid mode for compact biological encodings.
    """

    def normalize(self, seq: str) -> str:
        return "".join(seq.split()).upper()

    def encode_query(self, seq: str) -> bytes:
        return _encode_ascii(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        return ByteTarget(_encode_ascii(seq))

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return query_symbol == target_symbol

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement(pattern)


class NucleotideCodec:
    """Exact nucleotide matching with adaptive two-bit target storage.

    A, C, G, T, and U targets use two bits per base. Targets containing IUPAC
    symbols such as N or R fall back to one exact code per byte. U and T are
    equivalent in nucleotide mode.
    """

    def normalize(self, seq: str) -> str:
        return normalize_sequence(seq)

    def encode_query(self, seq: str) -> bytes:
        return encode_nucleotide_query(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        try:
            packed, length = encode_canonical_2bit(seq)
        except ValueError:
            return ByteTarget(encode_nucleotide_query(seq))

        return PackedCanonicalCodeTarget(data=packed, length=length)

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return query_symbol == target_symbol

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement(pattern)


class IupacNucleotideCodec:
    """IUPAC nucleotide ambiguity matching.

    Canonical targets use two bits per base. Ambiguous targets use one
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
            return ByteTarget(encode_iupac_query(seq))

        return PackedCanonicalMaskTarget(data=packed, length=length)

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return (query_symbol & target_symbol) != 0

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement(pattern)


class ProteinCodec:
    """Exact amino-acid matching with five-bit target storage.

    B, J, Z, X, U, O, ``*``, ``-``, and ``.`` are supported as literal
    symbols. Protein ambiguity matching is intentionally not enabled.
    """

    def normalize(self, seq: str) -> str:
        return normalize_protein(seq)

    def encode_query(self, seq: str) -> bytes:
        return encode_protein_query(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        packed, length = encode_protein_5bit(seq)
        return PackedProteinTarget(data=packed, length=length)

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return query_symbol == target_symbol

    def reverse_complement(self, pattern: str) -> str:
        raise ValueError("Reverse-complement search is not valid for amino acids")


# Backward-compatible name used by earlier versions.
IupacDnaCodec = IupacNucleotideCodec
