from __future__ import annotations

BASE_CODES: dict[str, int] = {
    "A": 0b00,
    "C": 0b01,
    "G": 0b10,
    "T": 0b11,
    "U": 0b11,
}

A_MASK = 1 << BASE_CODES["A"]
C_MASK = 1 << BASE_CODES["C"]
G_MASK = 1 << BASE_CODES["G"]
T_MASK = 1 << BASE_CODES["T"]

IUPAC_MASKS: dict[str, int] = {
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
}

# Exact nucleotide codes. Canonical symbols retain their two-bit values so a
# canonical target can use the same packed storage. U is treated as T.
NUCLEOTIDE_CODES: dict[str, int] = {
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
}

COMPLEMENT: dict[str, str] = {
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
}


def normalize_sequence(seq: str) -> str:
    """Remove whitespace and normalize nucleotide symbols to uppercase."""
    return "".join(seq.split()).upper()


def encode_nucleotide_query(seq: str) -> bytes:
    """Encode a nucleotide query for exact symbol matching."""
    result = bytearray()

    for position, base in enumerate(normalize_sequence(seq), start=1):
        try:
            result.append(NUCLEOTIDE_CODES[base])
        except KeyError as exc:
            raise ValueError(f"Unsupported nucleotide {base!r} at position {position}") from exc

    return bytes(result)


def encode_iupac_query(seq: str) -> bytes:
    """Encode an IUPAC query as one four-bit compatibility mask per byte."""
    result = bytearray()

    for position, base in enumerate(normalize_sequence(seq), start=1):
        try:
            result.append(IUPAC_MASKS[base])
        except KeyError as exc:
            raise ValueError(
                f"Unsupported IUPAC nucleotide {base!r} at position {position}"
            ) from exc

    return bytes(result)


def encode_iupac(seq: str) -> bytes:
    """Backward-compatible alias for :func:`encode_iupac_query`."""
    return encode_iupac_query(seq)


def encode_canonical_2bit(seq: str) -> tuple[bytes, int]:
    """Pack A, C, G, T, or U using two bits per nucleotide.

    U uses the same code as T. The returned length is the number of logical
    symbols rather than the number of physical bytes.
    """
    normalized = normalize_sequence(seq)
    packed = bytearray((len(normalized) + 3) // 4)

    for index, base in enumerate(normalized):
        try:
            code = BASE_CODES[base]
        except KeyError as exc:
            raise ValueError(
                "2-bit encoding only supports A, C, G, T, or U; "
                f"found {base!r} at position {index + 1}"
            ) from exc

        byte_index = index // 4
        shift = 6 - (index % 4) * 2
        packed[byte_index] |= code << shift

    return bytes(packed), len(normalized)


def get_2bit_base(
    packed: bytes | memoryview,
    index: int,
    length: int,
) -> int:
    """Return one raw two-bit nucleotide code from packed storage."""
    if index < 0 or index >= length:
        raise IndexError(f"Base index out of range: {index}")

    byte_index = index // 4
    shift = 6 - (index % 4) * 2
    return (packed[byte_index] >> shift) & 0b11


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of an IUPAC nucleotide sequence."""
    normalized = normalize_sequence(seq)

    try:
        return "".join(COMPLEMENT[base] for base in reversed(normalized))
    except KeyError as exc:
        raise ValueError(
            f"Unsupported base in reverse-complement pattern: {exc.args[0]!r}"
        ) from exc
