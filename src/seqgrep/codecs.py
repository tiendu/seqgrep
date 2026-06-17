from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Protocol

from .alphabets import (
    encode_nucleotide_exact,
    encode_nucleotide_iupac,
    encode_protein_exact,
    normalize_symbols,
    nucleotide_exact_text,
    pack_nucleotide_2bit,
    pack_protein_5bit,
    protein_exact_text,
    reverse_complement_nucleotide,
    unpack_nucleotide_2bit,
    unpack_protein_5bit,
)


class TargetEncoding(str, Enum):
    """Physical representation used for an encoded target sequence."""

    BYTE = "byte"
    NUCLEOTIDE_2BIT_CODE = "nucleotide-2bit-code"
    NUCLEOTIDE_2BIT_MASK = "nucleotide-2bit-mask"
    PROTEIN_5BIT = "protein-5bit"


class EncodedTarget(Protocol):
    """Random-access target with physical bytes and logical symbol length."""

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


class ExactSequenceCodec(SequenceCodec, Protocol):
    """Codec that can provide a text key for exact serial searching."""

    def comparison_text(self, seq: str) -> str: ...


def decode_target_symbol(
    storage: bytes | memoryview,
    index: int,
    length: int,
    encoding: TargetEncoding,
) -> int:
    """Decode one logical target symbol from physical storage."""
    if index < 0 or index >= length:
        raise IndexError(f"Target index out of range: {index}")

    if encoding is TargetEncoding.BYTE:
        return storage[index]

    if encoding is TargetEncoding.NUCLEOTIDE_2BIT_CODE:
        return unpack_nucleotide_2bit(storage, index, length)

    if encoding is TargetEncoding.NUCLEOTIDE_2BIT_MASK:
        return 1 << unpack_nucleotide_2bit(storage, index, length)

    if encoding is TargetEncoding.PROTEIN_5BIT:
        return unpack_protein_5bit(storage, index, length)

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
        return decode_target_symbol(self.data, index, len(self), self.encoding)


@dataclass(frozen=True, slots=True)
class PackedTarget:
    """Packed biological target with an explicit logical length."""

    data: bytes
    length: int
    encoding: TargetEncoding

    def __len__(self) -> int:
        return self.length

    def symbol_at(self, index: int) -> int:
        return decode_target_symbol(self.data, index, self.length, self.encoding)


class NucleotideCodec:
    """Exact nucleotide matching with adaptive two-bit target storage."""

    def normalize(self, seq: str) -> str:
        return normalize_symbols(seq)

    def comparison_text(self, seq: str) -> str:
        return nucleotide_exact_text(seq)

    def encode_query(self, seq: str) -> bytes:
        return encode_nucleotide_exact(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        try:
            packed, length = pack_nucleotide_2bit(seq)
        except ValueError:
            return ByteTarget(encode_nucleotide_exact(seq))

        return PackedTarget(
            data=packed,
            length=length,
            encoding=TargetEncoding.NUCLEOTIDE_2BIT_CODE,
        )

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return query_symbol == target_symbol

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement_nucleotide(pattern)


class IupacNucleotideCodec:
    """IUPAC nucleotide compatibility matching with adaptive storage."""

    def normalize(self, seq: str) -> str:
        return normalize_symbols(seq)

    def encode_query(self, seq: str) -> bytes:
        return encode_nucleotide_iupac(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        try:
            packed, length = pack_nucleotide_2bit(seq)
        except ValueError:
            return ByteTarget(encode_nucleotide_iupac(seq))

        return PackedTarget(
            data=packed,
            length=length,
            encoding=TargetEncoding.NUCLEOTIDE_2BIT_MASK,
        )

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return (query_symbol & target_symbol) != 0

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement_nucleotide(pattern)


class ProteinCodec:
    """Exact amino-acid matching with five-bit target storage."""

    def normalize(self, seq: str) -> str:
        return normalize_symbols(seq)

    def comparison_text(self, seq: str) -> str:
        return protein_exact_text(seq)

    def encode_query(self, seq: str) -> bytes:
        return encode_protein_exact(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        packed, length = pack_protein_5bit(seq)
        return PackedTarget(
            data=packed,
            length=length,
            encoding=TargetEncoding.PROTEIN_5BIT,
        )

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return query_symbol == target_symbol

    def reverse_complement(self, pattern: str) -> str:
        raise ValueError("Reverse-complement search is not valid for amino acids")
