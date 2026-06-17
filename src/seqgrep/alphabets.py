from __future__ import annotations

from collections.abc import Mapping

# ---------------------------------------------------------------------------
# Shared normalization and validation
# ---------------------------------------------------------------------------


def normalize_symbols(seq: str) -> str:
    """Remove whitespace and normalize biological symbols to uppercase."""
    return "".join(seq.split()).upper()


def _validate_symbols(
    normalized: str,
    codes: Mapping[str, int],
    *,
    alphabet_name: str,
) -> None:
    """Validate symbols without allocating a full encoded sequence."""
    invalid = set(normalized).difference(codes)
    if not invalid:
        return

    for position, symbol in enumerate(normalized, start=1):
        if symbol in invalid:
            raise ValueError(
                f"Unsupported {alphabet_name} symbol {symbol!r} at position {position}"
            )

    raise AssertionError("Invalid symbol set was unexpectedly empty")


def _encode_symbols(
    seq: str,
    codes: Mapping[str, int],
    *,
    alphabet_name: str,
) -> bytes:
    normalized = normalize_symbols(seq)
    _validate_symbols(normalized, codes, alphabet_name=alphabet_name)
    return bytes(codes[symbol] for symbol in normalized)


def _pack_5bit_symbols(
    seq: str,
    codes: Mapping[str, int],
    *,
    alphabet_name: str,
) -> tuple[bytes, int]:
    normalized = normalize_symbols(seq)
    _validate_symbols(normalized, codes, alphabet_name=alphabet_name)
    packed = bytearray((len(normalized) * 5 + 7) // 8)

    for index, symbol in enumerate(normalized):
        code = codes[symbol]
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


def _unpack_5bit_symbol(
    packed: bytes | memoryview,
    index: int,
    length: int,
    *,
    alphabet_name: str,
) -> int:
    if index < 0 or index >= length:
        raise IndexError(f"{alphabet_name} index out of range: {index}")

    bit_offset = index * 5
    byte_index = bit_offset // 8
    offset_in_byte = bit_offset % 8

    value = packed[byte_index] << 8
    if byte_index + 1 < len(packed):
        value |= packed[byte_index + 1]

    shift = 16 - offset_in_byte - 5
    return (value >> shift) & 0b11111


def _set_bitmap_bit(bitmap: bytearray, index: int) -> None:
    bitmap[index // 8] |= 1 << (7 - index % 8)


def _bitmap_bit_is_set(
    storage: bytes | memoryview,
    offset: int,
    index: int,
) -> bool:
    return bool(storage[offset + index // 8] & (1 << (7 - index % 8)))


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

# Exact mode keeps ambiguity symbols literal. T and U deliberately share a
# code because this tool treats them as equivalent nucleotide spellings.
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

# Query-side ambiguity expands IUPAC symbols into the canonical bases they may
# represent. This is used by --ambig and --ambig-mode query/both.
NUCLEOTIDE_IUPAC_QUERY_MASKS: dict[str, int] = {
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

# Query-only target semantics are deliberately asymmetric. Canonical bases
# produce singleton masks, while ambiguous bases produce zero and match
# nothing. Symmetric ambiguity uses NUCLEOTIDE_IUPAC_QUERY_MASKS for the target
# as well.
NUCLEOTIDE_IUPAC_TARGET_MASKS: dict[str, int] = {
    "A": A_MASK,
    "C": C_MASK,
    "G": G_MASK,
    "T": T_MASK,
    "U": T_MASK,
    "R": 0,
    "Y": 0,
    "S": 0,
    "W": 0,
    "K": 0,
    "M": 0,
    "B": 0,
    "D": 0,
    "H": 0,
    "V": 0,
    "N": 0,
    "-": GAP_MASK,
    ".": GAP_MASK,
}

# Backward-compatible name. New code should use the explicit query or target
# table instead.
NUCLEOTIDE_IUPAC_MASKS = NUCLEOTIDE_IUPAC_QUERY_MASKS

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


def encode_nucleotide_exact(seq: str) -> bytes:
    """Encode nucleotide symbols for exact matching, one code per byte."""
    return _encode_symbols(
        seq,
        NUCLEOTIDE_EXACT_CODES,
        alphabet_name="nucleotide",
    )


def encode_nucleotide_iupac_query(seq: str) -> bytes:
    """Encode a nucleotide query as IUPAC compatibility masks."""
    return _encode_symbols(
        seq,
        NUCLEOTIDE_IUPAC_QUERY_MASKS,
        alphabet_name="IUPAC nucleotide query",
    )


def encode_nucleotide_iupac_target(seq: str) -> bytes:
    """Encode a target as query-side IUPAC compatibility masks."""
    return _encode_symbols(
        seq,
        NUCLEOTIDE_IUPAC_TARGET_MASKS,
        alphabet_name="IUPAC nucleotide target",
    )


def encode_nucleotide_iupac(seq: str) -> bytes:
    """Backward-compatible alias for query-side IUPAC encoding."""
    return encode_nucleotide_iupac_query(seq)


def nucleotide_exact_text(seq: str) -> str:
    """Return validated text used by the exact serial matcher."""
    normalized = normalize_symbols(seq)
    _validate_symbols(
        normalized,
        NUCLEOTIDE_EXACT_CODES,
        alphabet_name="nucleotide",
    )
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


def pack_nucleotide_exact_5bit(seq: str) -> tuple[bytes, int]:
    """Pack the complete exact nucleotide alphabet using five bits."""
    return _pack_5bit_symbols(
        seq,
        NUCLEOTIDE_EXACT_CODES,
        alphabet_name="nucleotide",
    )


def unpack_nucleotide_exact_5bit(
    packed: bytes | memoryview,
    index: int,
    length: int,
) -> int:
    """Return one exact nucleotide code from five-bit storage."""
    return _unpack_5bit_symbol(
        packed,
        index,
        length,
        alphabet_name="Nucleotide",
    )


def pack_nucleotide_iupac_masks_5bit(seq: str) -> tuple[bytes, int]:
    """Pack full IUPAC compatibility masks using five bits per symbol.

    This representation is used by symmetric ambiguity mode, where ambiguous
    symbols in both the query and target participate in mask intersection.
    Five bits are required because the gap mask uses bit 4.
    """
    return _pack_5bit_symbols(
        seq,
        NUCLEOTIDE_IUPAC_QUERY_MASKS,
        alphabet_name="IUPAC nucleotide",
    )


def unpack_nucleotide_iupac_mask_5bit(
    packed: bytes | memoryview,
    index: int,
    length: int,
) -> int:
    """Return one full IUPAC mask from five-bit storage."""
    return _unpack_5bit_symbol(
        packed,
        index,
        length,
        alphabet_name="IUPAC nucleotide",
    )


def pack_nucleotide_iupac_target(seq: str) -> tuple[bytes, int, bool]:
    """Pack a mixed IUPAC target using bases plus validity bitmaps.

    Layout:

    1. two-bit base codes for every position,
    2. one validity bit per position,
    3. optionally, one gap bit per position.

    Ambiguous target bases have neither bit set and decode to zero. Targets
    without gaps need only three bits per position on average; targets with
    gaps need four.
    """
    normalized = normalize_symbols(seq)
    _validate_symbols(
        normalized,
        NUCLEOTIDE_IUPAC_TARGET_MASKS,
        alphabet_name="IUPAC nucleotide target",
    )

    length = len(normalized)
    codes = bytearray((length + 3) // 4)
    valid = bytearray((length + 7) // 8)
    gaps: bytearray | None = None

    for index, symbol in enumerate(normalized):
        code = NUCLEOTIDE_2BIT_CODES.get(symbol)
        if code is not None:
            shift = 6 - (index % 4) * 2
            codes[index // 4] |= code << shift
            _set_bitmap_bit(valid, index)
            continue

        if symbol in {"-", "."}:
            if gaps is None:
                gaps = bytearray((length + 7) // 8)
            _set_bitmap_bit(gaps, index)

    return bytes(codes + valid + (gaps or b"")), length, gaps is not None


def unpack_nucleotide_iupac_target(
    packed: bytes | memoryview,
    index: int,
    length: int,
    *,
    has_gap_bitmap: bool,
) -> int:
    """Decode one target mask from mixed packed IUPAC storage."""
    if index < 0 or index >= length:
        raise IndexError(f"Nucleotide index out of range: {index}")

    code_size = (length + 3) // 4
    bitmap_size = (length + 7) // 8
    valid_offset = code_size

    if _bitmap_bit_is_set(packed, valid_offset, index):
        return 1 << unpack_nucleotide_2bit(packed, index, length)

    if has_gap_bitmap:
        gap_offset = code_size + bitmap_size
        if _bitmap_bit_is_set(packed, gap_offset, index):
            return GAP_MASK

    return 0


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
# and X remain literal; protein ambiguity expansion is not enabled.
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
    _validate_symbols(
        normalized,
        PROTEIN_5BIT_CODES,
        alphabet_name="amino-acid",
    )
    return normalized


def pack_protein_5bit(seq: str) -> tuple[bytes, int]:
    """Pack supported amino-acid symbols using five bits per residue."""
    return _pack_5bit_symbols(
        seq,
        PROTEIN_5BIT_CODES,
        alphabet_name="amino-acid",
    )


def unpack_protein_5bit(
    packed: bytes | memoryview,
    index: int,
    length: int,
) -> int:
    """Return one five-bit amino-acid code from packed storage."""
    return _unpack_5bit_symbol(
        packed,
        index,
        length,
        alphabet_name="Amino-acid",
    )
