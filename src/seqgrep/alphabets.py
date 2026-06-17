from __future__ import annotations

from collections.abc import Mapping

# ---------------------------------------------------------------------------
# Shared normalization
# ---------------------------------------------------------------------------


def normalize_symbols(seq: str) -> str:
    """Remove whitespace and normalize biological symbols to uppercase."""
    return "".join(seq.split()).upper()


# ---------------------------------------------------------------------------
# Nucleotide alphabet
# ---------------------------------------------------------------------------

NUCLEOTIDE_2BIT_CODES: dict[str, int] = {
    "A": 0b00,
    "C": 0b01,
    "G": 0b10,
    "T": 0b11,
    "U": 0b11,
}

# Exact nucleotide mode keeps ambiguity symbols literal. T and U deliberately
# share a code because they are equivalent nucleotide symbols in this tool.
NUCLEOTIDE_EXACT_CODES: dict[str, int] = {
    "A": 0,
    "C": 1,
    "G": 2,
    "T": 3,
    "U": 3,
    "R": 4,
    "Y": 5,
    "S": 6,
    "W": 7,
    "K": 8,
    "M": 9,
    "B": 10,
    "D": 11,
    "H": 12,
    "V": 13,
    "N": 14,
    "-": 15,
    ".": 16,
}

A_MASK = 1 << NUCLEOTIDE_2BIT_CODES["A"]
C_MASK = 1 << NUCLEOTIDE_2BIT_CODES["C"]
G_MASK = 1 << NUCLEOTIDE_2BIT_CODES["G"]
T_MASK = 1 << NUCLEOTIDE_2BIT_CODES["T"]
GAP_MASK = 1 << 4

NUCLEOTIDE_IUPAC_MASKS: dict[str, int] = {
    "A": A_MASK,
    "C": C_MASK,
    "G": G_MASK,
    "T": T_MASK,
    "U": T_MASK,
    "R": A_MASK | G_MASK,
    "Y": C_MASK | T_MASK,
    "S": G_MASK | C_MASK,
    "W": A_MASK | T_MASK,
    "K": G_MASK | T_MASK,
    "M": A_MASK | C_MASK,
    "B": C_MASK | G_MASK | T_MASK,
    "D": A_MASK | G_MASK | T_MASK,
    "H": A_MASK | C_MASK | T_MASK,
    "V": A_MASK | C_MASK | G_MASK,
    "N": A_MASK | C_MASK | G_MASK | T_MASK,
    "-": GAP_MASK,
    ".": GAP_MASK,
}

NUCLEOTIDE_COMPLEMENT: dict[str, str] = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "U": "A",
    "R": "Y",
    "Y": "R",
    "S": "S",
    "W": "W",
    "K": "M",
    "M": "K",
    "B": "V",
    "D": "H",
    "H": "D",
    "V": "B",
    "N": "N",
    "-": "-",
    ".": ".",
}


def _encode_symbols(
    seq: str,
    codes: Mapping[str, int],
    *,
    alphabet_name: str,
) -> bytes:
    normalized = normalize_symbols(seq)
    encoded = bytearray()

    for position, symbol in enumerate(normalized, start=1):
        try:
            encoded.append(codes[symbol])
        except KeyError as exc:
            raise ValueError(
                f"Unsupported {alphabet_name} symbol {symbol!r} at position {position}"
            ) from exc

    return bytes(encoded)


def encode_nucleotide_exact(seq: str) -> bytes:
    """Encode nucleotide symbols for exact matching, one code per byte."""
    return _encode_symbols(
        seq,
        NUCLEOTIDE_EXACT_CODES,
        alphabet_name="nucleotide",
    )


def encode_nucleotide_iupac(seq: str) -> bytes:
    """Encode nucleotide symbols as IUPAC compatibility masks."""
    return _encode_symbols(
        seq,
        NUCLEOTIDE_IUPAC_MASKS,
        alphabet_name="IUPAC nucleotide",
    )


def nucleotide_exact_text(seq: str) -> str:
    """Return validated text used by the exact serial matcher.

    U is converted to T so DNA and RNA spelling compare as equivalent.
    Ambiguity symbols and the two gap characters remain literal.
    """
    normalized = normalize_symbols(seq)
    encode_nucleotide_exact(normalized)  # validation
    return normalized.replace("U", "T")


def pack_nucleotide_2bit(seq: str) -> tuple[bytes, int]:
    """Pack a canonical A/C/G/T/U sequence using two bits per base."""
    normalized = normalize_symbols(seq)
    packed = bytearray((len(normalized) + 3) // 4)

    for index, symbol in enumerate(normalized):
        try:
            code = NUCLEOTIDE_2BIT_CODES[symbol]
        except KeyError as exc:
            raise ValueError(
                "2-bit nucleotide encoding only supports A, C, G, T, or U; "
                f"found {symbol!r} at position {index + 1}"
            ) from exc

        shift = 6 - (index % 4) * 2
        packed[index // 4] |= code << shift

    return bytes(packed), len(normalized)


def unpack_nucleotide_2bit(
    packed: bytes | memoryview,
    index: int,
    length: int,
) -> int:
    """Return one canonical two-bit code from packed storage."""
    if index < 0 or index >= length:
        raise IndexError(f"Nucleotide index out of range: {index}")

    shift = 6 - (index % 4) * 2
    return (packed[index // 4] >> shift) & 0b11


def reverse_complement_nucleotide(seq: str) -> str:
    """Return the reverse complement of an IUPAC nucleotide sequence."""
    normalized = normalize_symbols(seq)

    try:
        return "".join(NUCLEOTIDE_COMPLEMENT[symbol] for symbol in reversed(normalized))
    except KeyError as exc:
        raise ValueError(
            f"Unsupported nucleotide symbol {exc.args[0]!r} in reverse-complement pattern"
        ) from exc


# ---------------------------------------------------------------------------
# Amino-acid alphabet
# ---------------------------------------------------------------------------

# Twenty canonical residues plus common extended sequence symbols. B, J, Z,
# and X remain literal in the current CLI; no protein ambiguity mode is
# enabled. Twenty-nine symbols fit comfortably in five bits.
PROTEIN_SYMBOLS: tuple[str, ...] = tuple("ACDEFGHIKLMNPQRSTVWYBJZXUO*-.")
PROTEIN_5BIT_CODES: dict[str, int] = {symbol: code for code, symbol in enumerate(PROTEIN_SYMBOLS)}

if len(PROTEIN_5BIT_CODES) > 32:
    raise RuntimeError("The protein alphabet does not fit in five bits")


def encode_protein_exact(seq: str) -> bytes:
    """Encode a protein query as one five-bit code per byte."""
    return _encode_symbols(
        seq,
        PROTEIN_5BIT_CODES,
        alphabet_name="amino-acid",
    )


def protein_exact_text(seq: str) -> str:
    """Return validated normalized text for exact serial protein search."""
    normalized = normalize_symbols(seq)
    encode_protein_exact(normalized)  # validation
    return normalized


def pack_protein_5bit(seq: str) -> tuple[bytes, int]:
    """Pack supported amino-acid symbols using five bits per residue."""
    normalized = normalize_symbols(seq)
    packed = bytearray((len(normalized) * 5 + 7) // 8)

    for index, symbol in enumerate(normalized):
        try:
            code = PROTEIN_5BIT_CODES[symbol]
        except KeyError as exc:
            raise ValueError(
                f"Unsupported amino-acid symbol {symbol!r} at position {index + 1}"
            ) from exc

        bit_offset = index * 5
        byte_index = bit_offset // 8
        offset_in_byte = bit_offset % 8
        first_bits = min(5, 8 - offset_in_byte)
        remaining_bits = 5 - first_bits

        if remaining_bits == 0:
            packed[byte_index] |= code << (8 - offset_in_byte - 5)
        else:
            packed[byte_index] |= code >> remaining_bits
            remainder_mask = (1 << remaining_bits) - 1
            packed[byte_index + 1] |= (code & remainder_mask) << (8 - remaining_bits)

    return bytes(packed), len(normalized)


def unpack_protein_5bit(
    packed: bytes | memoryview,
    index: int,
    length: int,
) -> int:
    """Return one five-bit amino-acid code from packed storage."""
    if index < 0 or index >= length:
        raise IndexError(f"Amino-acid index out of range: {index}")

    bit_offset = index * 5
    byte_index = bit_offset // 8
    offset_in_byte = bit_offset % 8

    value = packed[byte_index] << 8
    if byte_index + 1 < len(packed):
        value |= packed[byte_index + 1]

    shift = 16 - offset_in_byte - 5
    return (value >> shift) & 0b11111
