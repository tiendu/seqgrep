from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Protocol

from .alphabets import (
    encode_nucleotide_exact,
    encode_nucleotide_iupac_query,
    encode_protein_exact,
    normalize_symbols,
    nucleotide_exact_text,
    pack_nucleotide_2bit,
    pack_nucleotide_exact_5bit,
    pack_nucleotide_iupac_masks_5bit,
    pack_nucleotide_iupac_target,
    pack_protein_5bit,
    protein_exact_text,
    reverse_complement_nucleotide,
    unpack_nucleotide_2bit,
    unpack_nucleotide_exact_5bit,
    unpack_nucleotide_iupac_mask_5bit,
    unpack_nucleotide_iupac_target,
    unpack_protein_5bit,
)


class TargetEncoding(str, Enum):
    """Physical representation used for an encoded target sequence."""

    BYTE = "byte"
    NUCLEOTIDE_2BIT_CODE = "nucleotide-2bit-code"
    NUCLEOTIDE_EXACT_5BIT = "nucleotide-exact-5bit"
    NUCLEOTIDE_2BIT_MASK = "nucleotide-2bit-mask"
    NUCLEOTIDE_2BIT_VALID_MASK = "nucleotide-2bit-valid-mask"
    NUCLEOTIDE_2BIT_VALID_GAP_MASK = "nucleotide-2bit-valid-gap-mask"
    NUCLEOTIDE_IUPAC_5BIT_MASK = "nucleotide-iupac-5bit-mask"
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

    if encoding is TargetEncoding.NUCLEOTIDE_EXACT_5BIT:
        return unpack_nucleotide_exact_5bit(storage, index, length)

    if encoding is TargetEncoding.NUCLEOTIDE_2BIT_MASK:
        return 1 << unpack_nucleotide_2bit(storage, index, length)

    if encoding is TargetEncoding.NUCLEOTIDE_2BIT_VALID_MASK:
        return unpack_nucleotide_iupac_target(
            storage,
            index,
            length,
            has_gap_bitmap=False,
        )

    if encoding is TargetEncoding.NUCLEOTIDE_2BIT_VALID_GAP_MASK:
        return unpack_nucleotide_iupac_target(
            storage,
            index,
            length,
            has_gap_bitmap=True,
        )

    if encoding is TargetEncoding.NUCLEOTIDE_IUPAC_5BIT_MASK:
        return unpack_nucleotide_iupac_mask_5bit(storage, index, length)

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
    """Exact nucleotide matching with adaptive packed target storage."""

    def normalize(self, seq: str) -> str:
        return normalize_symbols(seq)

    def comparison_text(self, seq: str) -> str:
        return nucleotide_exact_text(seq)

    def encode_query(self, seq: str) -> bytes:
        return encode_nucleotide_exact(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        try:
            packed, length = pack_nucleotide_2bit(seq)
            encoding = TargetEncoding.NUCLEOTIDE_2BIT_CODE
        except ValueError:
            packed, length = pack_nucleotide_exact_5bit(seq)
            encoding = TargetEncoding.NUCLEOTIDE_EXACT_5BIT

        return PackedTarget(
            data=packed,
            length=length,
            encoding=encoding,
        )

    def compatible(self, query_symbol: int, target_symbol: int) -> bool:
        return query_symbol == target_symbol

    def reverse_complement(self, pattern: str) -> str:
        return reverse_complement_nucleotide(pattern)


class IupacNucleotideCodec:
    """IUPAC matching with selectable target ambiguity semantics.

    By default, ambiguity applies only to the query: ambiguous target symbols
    decode to zero and match nothing. With allow_target_ambiguity=True, target
    symbols use their full IUPAC masks and matching becomes symmetric.
    """

    def __init__(self, *, allow_target_ambiguity: bool = False) -> None:
        self.allow_target_ambiguity = allow_target_ambiguity

    def normalize(self, seq: str) -> str:
        return normalize_symbols(seq)

    def encode_query(self, seq: str) -> bytes:
        return encode_nucleotide_iupac_query(seq)

    def encode_target(self, seq: str) -> EncodedTarget:
        try:
            packed, length = pack_nucleotide_2bit(seq)
            encoding = TargetEncoding.NUCLEOTIDE_2BIT_MASK
        except ValueError:
            if self.allow_target_ambiguity:
                packed, length = pack_nucleotide_iupac_masks_5bit(seq)
                encoding = TargetEncoding.NUCLEOTIDE_IUPAC_5BIT_MASK
            else:
                packed, length, has_gaps = pack_nucleotide_iupac_target(seq)
                encoding = (
                    TargetEncoding.NUCLEOTIDE_2BIT_VALID_GAP_MASK
                    if has_gaps
                    else TargetEncoding.NUCLEOTIDE_2BIT_VALID_MASK
                )

        return PackedTarget(
            data=packed,
            length=length,
            encoding=encoding,
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
