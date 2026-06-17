from __future__ import annotations

PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYBJZXUO*-."
PROTEIN_SYMBOLS: tuple[str, ...] = tuple(PROTEIN_ALPHABET)

PROTEIN_CODES: dict[str, int] = {symbol: code for code, symbol in enumerate(PROTEIN_SYMBOLS)}

if len(PROTEIN_CODES) > 32:
    raise RuntimeError("The protein alphabet does not fit in five bits")


def normalize_protein(seq: str) -> str:
    """Remove whitespace and normalize protein symbols to uppercase."""
    return "".join(seq.split()).upper()


def encode_protein_query(seq: str) -> bytes:
    """Encode a protein query as one five-bit code per byte.

    Query sequences stay unpacked because they are usually short and direct
    byte iteration is faster than repeatedly unpacking both sides.
    """
    normalized = normalize_protein(seq)
    result = bytearray()

    for position, symbol in enumerate(normalized, start=1):
        try:
            result.append(PROTEIN_CODES[symbol])
        except KeyError as exc:
            raise ValueError(
                f"Unsupported amino-acid symbol {symbol!r} at position {position}"
            ) from exc

    return bytes(result)


def encode_protein_5bit(seq: str) -> tuple[bytes, int]:
    """Pack a protein sequence using five bits per residue.

    The supported exact-match alphabet is the 20 canonical amino acids plus
    B, J, Z, X, U, O, stop ``*``, and alignment gaps ``-`` and ``.``.
    Ambiguous-looking protein symbols remain literal symbols.
    """
    normalized = normalize_protein(seq)
    packed = bytearray((len(normalized) * 5 + 7) // 8)

    for index, symbol in enumerate(normalized):
        try:
            code = PROTEIN_CODES[symbol]
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


def get_protein_5bit(packed: bytes | memoryview, index: int, length: int) -> int:
    """Return one five-bit protein code from packed storage."""
    if index < 0 or index >= length:
        raise IndexError(f"Protein index out of range: {index}")

    bit_offset = index * 5
    byte_index = bit_offset // 8
    offset_in_byte = bit_offset % 8

    value = packed[byte_index] << 8
    if byte_index + 1 < len(packed):
        value |= packed[byte_index + 1]

    shift = 16 - offset_in_byte - 5
    return (value >> shift) & 0b11111
