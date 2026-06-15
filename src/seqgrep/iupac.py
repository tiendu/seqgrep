from __future__ import annotations

A = 0b0001
C = 0b0010
G = 0b0100
T = 0b1000

IUPAC_MASKS: dict[str, int] = {
    "A": A,
    "C": C,
    "G": G,
    "T": T,
    "U": T,
    "R": A | G,
    "Y": C | T,
    "S": G | C,
    "W": A | T,
    "K": G | T,
    "M": A | C,
    "B": C | G | T,
    "D": A | G | T,
    "H": A | C | T,
    "V": A | C | G,
    "N": A | C | G | T,
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
    return "".join(seq.split()).upper()


def encode_iupac(seq: str) -> bytes:
    masks = bytearray()
    for pos, base in enumerate(normalize_sequence(seq), start=1):
        try:
            masks.append(IUPAC_MASKS[base])
        except KeyError as exc:
            raise ValueError(f"Unsupported IUPAC nucleotide {base!r} at position {pos}") from exc
    return bytes(masks)


def reverse_complement(seq: str) -> str:
    normalized = normalize_sequence(seq)
    try:
        return "".join(COMPLEMENT[base] for base in reversed(normalized))
    except KeyError as exc:
        raise ValueError(f"Unsupported base in reverse-complement pattern: {exc.args[0]!r}") from exc

