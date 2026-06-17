from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from seqgrep.cli import parse_args
from seqgrep.codecs import (
    ByteTarget,
    IupacDnaCodec,
    IupacNucleotideCodec,
    LiteralCodec,
    NucleotideCodec,
    PackedCanonicalCodeTarget,
    PackedCanonicalMaskTarget,
    PackedProteinTarget,
    ProteinCodec,
    TargetEncoding,
)
from seqgrep.fasta import FastaFileReader
from seqgrep.iupac import (
    encode_canonical_2bit,
    encode_iupac_query,
    encode_nucleotide_query,
    reverse_complement,
)
from seqgrep.models import FastaRecord, Match, SearchQuery, SequenceType
from seqgrep.planner import SearchPlanner
from seqgrep.protein import (
    PROTEIN_CODES,
    PROTEIN_SYMBOLS,
    encode_protein_5bit,
    encode_protein_query,
    get_protein_5bit,
)


def search_with_planner(
    record: FastaRecord,
    query: SearchQuery,
    *,
    jobs: int = 1,
    chunk_size: int = 1_000_000,
) -> list[Match]:
    plan = SearchPlanner().plan(
        query=query,
        jobs=jobs,
        chunk_size=chunk_size,
    )
    return list(plan.matcher.search(record, query))


def hit_values(hits: list[Match]) -> list[tuple[int, int, str]]:
    return [(hit.start, hit.end, hit.matched) for hit in hits]


# ---------------------------------------------------------------------------
# Nucleotide encodings
# ---------------------------------------------------------------------------


def test_exact_nucleotide_query_codes() -> None:
    assert encode_nucleotide_query("ACGTURYN") == bytes([0, 1, 2, 3, 3, 4, 5, 14])


def test_iupac_query_bitmasks() -> None:
    assert encode_iupac_query("ACGT") == bytes([0b0001, 0b0010, 0b0100, 0b1000])
    assert encode_iupac_query("RYN") == bytes([0b0101, 0b1010, 0b1111])


def test_nucleotide_encoders_normalize_case_and_whitespace() -> None:
    assert encode_nucleotide_query("a c\ng\tt") == bytes([0, 1, 2, 3])
    assert encode_iupac_query("a c\ng\tt") == bytes([0b0001, 0b0010, 0b0100, 0b1000])


def test_nucleotide_encoders_reject_invalid_symbol() -> None:
    with pytest.raises(ValueError, match="Unsupported nucleotide 'X' at position 3"):
        encode_nucleotide_query("ATX")

    with pytest.raises(
        ValueError,
        match="Unsupported IUPAC nucleotide 'X' at position 3",
    ):
        encode_iupac_query("ATX")


def test_reverse_complement_ambiguous() -> None:
    assert reverse_complement("ARYN") == "NRYT"


def test_encode_canonical_2bit() -> None:
    packed, length = encode_canonical_2bit("ACGT")

    assert length == 4
    assert packed == bytes([0b00011011])


def test_encode_canonical_2bit_partial_final_byte() -> None:
    packed, length = encode_canonical_2bit("ACGTA")

    assert length == 5
    assert packed == bytes([0b00011011, 0b00000000])


def test_encode_canonical_2bit_treats_u_as_t() -> None:
    assert encode_canonical_2bit("T") == encode_canonical_2bit("U")


def test_encode_canonical_2bit_rejects_ambiguity() -> None:
    with pytest.raises(
        ValueError,
        match="2-bit encoding only supports A, C, G, T, or U",
    ):
        encode_canonical_2bit("ATGN")


# ---------------------------------------------------------------------------
# Protein encodings
# ---------------------------------------------------------------------------


def test_protein_alphabet_fits_in_five_bits() -> None:
    assert len(PROTEIN_CODES) == len(PROTEIN_SYMBOLS)
    assert len(PROTEIN_CODES) <= 32
    assert max(PROTEIN_CODES.values()) <= 0b11111


def test_protein_query_uses_one_code_per_symbol() -> None:
    query = "MTEYKLVVVG"
    encoded = encode_protein_query(query)

    assert len(encoded) == len(query)
    assert encoded == bytes(PROTEIN_CODES[symbol] for symbol in query)


def test_protein_5bit_round_trip_for_full_alphabet() -> None:
    sequence = "".join(PROTEIN_SYMBOLS)
    packed, length = encode_protein_5bit(sequence)

    assert length == len(sequence)
    assert len(packed) == (length * 5 + 7) // 8
    assert [get_protein_5bit(packed, index, length) for index in range(length)] == [
        PROTEIN_CODES[symbol] for symbol in sequence
    ]


def test_protein_5bit_round_trip_across_byte_boundaries() -> None:
    for length in range(1, 65):
        sequence = "".join(PROTEIN_SYMBOLS[index % len(PROTEIN_SYMBOLS)] for index in range(length))
        packed, encoded_length = encode_protein_5bit(sequence)

        assert encoded_length == length
        assert [
            get_protein_5bit(packed, index, encoded_length) for index in range(encoded_length)
        ] == [PROTEIN_CODES[symbol] for symbol in sequence]


def test_protein_encoding_normalizes_case_and_whitespace() -> None:
    assert encode_protein_query("m t\ne\tyk") == encode_protein_query("MTEYK")


def test_protein_encoding_rejects_invalid_symbol() -> None:
    with pytest.raises(
        ValueError,
        match="Unsupported amino-acid symbol '@' at position 3",
    ):
        encode_protein_5bit("MT@K")


# ---------------------------------------------------------------------------
# Codecs and target representations
# ---------------------------------------------------------------------------


def test_legacy_literal_codec_remains_available() -> None:
    codec = LiteralCodec()
    target = codec.encode_target("ARNX")

    assert isinstance(target, ByteTarget)
    assert target.data == b"ARNX"
    assert codec.compatible(ord("A"), ord("A")) is True
    assert codec.compatible(ord("A"), ord("N")) is False


def test_default_nucleotide_codec_packs_canonical_target() -> None:
    codec = NucleotideCodec()
    target = codec.encode_target("ACGTU")

    assert isinstance(target, PackedCanonicalCodeTarget)
    assert target.encoding == TargetEncoding.PACKED_NUCLEOTIDE_2BIT
    assert len(target) == 5
    assert len(target.data) == 2
    assert [target.symbol_at(index) for index in range(len(target))] == [
        0,
        1,
        2,
        3,
        3,
    ]


def test_default_nucleotide_codec_falls_back_for_iupac_target() -> None:
    codec = NucleotideCodec()
    target = codec.encode_target("ANRY")

    assert isinstance(target, ByteTarget)
    assert [target.symbol_at(index) for index in range(len(target))] == [
        0,
        14,
        4,
        5,
    ]


def test_iupac_codec_packs_canonical_target_as_masks() -> None:
    codec = IupacNucleotideCodec()
    target = codec.encode_target("ACGTU")

    assert isinstance(target, PackedCanonicalMaskTarget)
    assert target.encoding == TargetEncoding.PACKED_CANONICAL_2BIT
    assert [target.symbol_at(index) for index in range(len(target))] == [
        0b0001,
        0b0010,
        0b0100,
        0b1000,
        0b1000,
    ]


def test_iupac_codec_falls_back_for_ambiguous_target() -> None:
    codec = IupacNucleotideCodec()
    target = codec.encode_target("ANRY")

    assert isinstance(target, ByteTarget)
    assert [target.symbol_at(index) for index in range(len(target))] == [
        0b0001,
        0b1111,
        0b0101,
        0b1010,
    ]


def test_old_iupac_codec_name_is_preserved() -> None:
    assert IupacDnaCodec is IupacNucleotideCodec


def test_protein_codec_always_uses_packed_5bit_target() -> None:
    codec = ProteinCodec()
    sequence = "MTEYKLVVVGAGGVGKSAL"
    target = codec.encode_target(sequence)

    assert isinstance(target, PackedProteinTarget)
    assert target.encoding == TargetEncoding.PACKED_PROTEIN_5BIT
    assert len(target) == len(sequence)
    assert len(target.data) == (len(sequence) * 5 + 7) // 8
    assert [target.symbol_at(index) for index in range(len(target))] == [
        PROTEIN_CODES[symbol] for symbol in sequence
    ]


def test_packed_targets_reject_padding_indexes() -> None:
    nucleotide_target = NucleotideCodec().encode_target("ACG")
    protein_target = ProteinCodec().encode_target("MTE")

    with pytest.raises(IndexError):
        nucleotide_target.symbol_at(3)

    with pytest.raises(IndexError):
        protein_target.symbol_at(3)


# ---------------------------------------------------------------------------
# Planner and serial behavior
# ---------------------------------------------------------------------------


def test_default_mode_is_exact_nucleotide() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN")

    assert hit_values(search_with_planner(record, query)) == [
        (1, 4, "ATGN"),
    ]


def test_default_nucleotide_mode_treats_u_as_t() -> None:
    record = FastaRecord("rna", "AUGATG")
    query = SearchQuery("ATG")

    assert hit_values(search_with_planner(record, query)) == [
        (1, 3, "AUG"),
        (4, 6, "ATG"),
    ]


def test_iupac_ambig_treats_n_as_any_base() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN", ambig=True)

    assert hit_values(search_with_planner(record, query)) == [
        (1, 4, "ATGN"),
        (5, 8, "ATGA"),
    ]


def test_protein_mode_matches_exact_amino_acids() -> None:
    record = FastaRecord("ras", "XXMTEYKLVVVGAGGVGKSALXX")
    query = SearchQuery(
        "MTEYK",
        sequence_type=SequenceType.AMINO_ACID,
    )

    assert hit_values(search_with_planner(record, query)) == [
        (3, 7, "MTEYK"),
    ]


def test_protein_x_is_literal_not_wildcard() -> None:
    record = FastaRecord("protein", "MAXMAB")
    query = SearchQuery(
        "MAX",
        sequence_type=SequenceType.AMINO_ACID,
    )

    assert hit_values(search_with_planner(record, query)) == [
        (1, 3, "MAX"),
    ]


def test_protein_mode_rejects_ambig() -> None:
    query = SearchQuery(
        "MTEYK",
        ambig=True,
        sequence_type=SequenceType.AMINO_ACID,
    )

    with pytest.raises(
        ValueError,
        match="--ambig is only valid for nucleotide sequences",
    ):
        SearchPlanner().plan(query, jobs=1, chunk_size=100)


def test_protein_mode_rejects_reverse_complement() -> None:
    query = SearchQuery(
        "MTEYK",
        revcomp=True,
        sequence_type=SequenceType.AMINO_ACID,
    )

    with pytest.raises(
        ValueError,
        match="--revcomp is only valid for nucleotide sequences",
    ):
        SearchPlanner().plan(query, jobs=1, chunk_size=100)


def test_nucleotide_reverse_complement_match() -> None:
    record = FastaRecord("seq1", "CCCATGAAGTCCC")
    query = SearchQuery("ACTTCAT", revcomp=True)

    hits = search_with_planner(record, query)

    assert any(
        hit.strand == "-" and hit.start == 4 and hit.end == 10 and hit.matched == "ATGAAGT"
        for hit in hits
    )


def test_circular_nucleotide_match() -> None:
    record = FastaRecord("plasmid", "CCCAAATTT")
    query = SearchQuery("TTTCCC", circular=True)

    hits = search_with_planner(record, query)

    assert len(hits) == 1
    assert hits[0] == Match(
        record="plasmid",
        strand="+",
        start=7,
        end=3,
        matched="TTTCCC",
        circular=True,
    )


def test_circular_protein_match() -> None:
    record = FastaRecord("ring", "CDEMA")
    query = SearchQuery(
        "MAC",
        circular=True,
        sequence_type=SequenceType.AMINO_ACID,
    )

    assert hit_values(search_with_planner(record, query)) == [
        (4, 1, "MAC"),
    ]


# ---------------------------------------------------------------------------
# Chunked matcher parity for every target encoding
# ---------------------------------------------------------------------------


def assert_chunked_matches_serial(
    record: FastaRecord,
    query: SearchQuery,
    *,
    chunk_size: int,
) -> list[Match]:
    serial_hits = search_with_planner(record, query, jobs=1)
    chunked_hits = search_with_planner(
        record,
        query,
        jobs=2,
        chunk_size=chunk_size,
    )

    assert chunked_hits == serial_hits
    return chunked_hits


def test_chunked_exact_nucleotide_packed_target() -> None:
    hits = assert_chunked_matches_serial(
        FastaRecord("seq1", "AAAAATGCAAAAA"),
        SearchQuery("TGCA"),
        chunk_size=6,
    )

    assert hit_values(hits) == [(6, 9, "TGCA")]


def test_chunked_exact_nucleotide_byte_fallback() -> None:
    assert_chunked_matches_serial(
        FastaRecord("seq1", "ATGNATGA"),
        SearchQuery("ATGN"),
        chunk_size=3,
    )


def test_chunked_iupac_packed_target() -> None:
    assert_chunked_matches_serial(
        FastaRecord("seq1", "ATGAATGT"),
        SearchQuery("ATGN", ambig=True),
        chunk_size=3,
    )


def test_chunked_iupac_byte_fallback() -> None:
    assert_chunked_matches_serial(
        FastaRecord("seq1", "ATGNATGA"),
        SearchQuery("ATGN", ambig=True),
        chunk_size=3,
    )


def test_chunked_protein_5bit_target_across_boundary() -> None:
    hits = assert_chunked_matches_serial(
        FastaRecord("protein", "AAAAAMTEYKAAAAA"),
        SearchQuery(
            "MTEYK",
            sequence_type=SequenceType.AMINO_ACID,
        ),
        chunk_size=6,
    )

    assert hit_values(hits) == [(6, 10, "MTEYK")]


def test_chunked_protein_circular_target() -> None:
    assert_chunked_matches_serial(
        FastaRecord("ring", "CDEMA"),
        SearchQuery(
            "MAC",
            circular=True,
            sequence_type=SequenceType.AMINO_ACID,
        ),
        chunk_size=2,
    )


def test_chunked_iupac_revcomp_circular_match() -> None:
    assert_chunked_matches_serial(
        FastaRecord("plasmid", "CCCAAATTT"),
        SearchQuery(
            "GGGAAA",
            revcomp=True,
            circular=True,
            ambig=True,
        ),
        chunk_size=4,
    )


# ---------------------------------------------------------------------------
# CLI validation
# ---------------------------------------------------------------------------


def test_cli_defaults_to_nucleotide(tmp_path: Path) -> None:
    input_path = tmp_path / "seq.fa"
    args = parse_args(["ATG", str(input_path)])

    assert args.sequence_type == SequenceType.NUCLEOTIDE.value
    assert args.ambig is False


def test_cli_accepts_amino_acid_mode(tmp_path: Path) -> None:
    input_path = tmp_path / "protein.fa"
    args = parse_args(["MTEYK", str(input_path), "--sequence-type", "amino-acid"])

    assert args.sequence_type == SequenceType.AMINO_ACID.value


def test_cli_rejects_ambig_for_amino_acid(tmp_path: Path) -> None:
    input_path = tmp_path / "protein.fa"

    with pytest.raises(SystemExit) as exc_info:
        parse_args(
            [
                "MTEYK",
                str(input_path),
                "--sequence-type",
                "amino-acid",
                "--ambig",
            ]
        )

    assert exc_info.value.code == 2


def test_cli_rejects_revcomp_for_amino_acid(tmp_path: Path) -> None:
    input_path = tmp_path / "protein.fa"

    with pytest.raises(SystemExit) as exc_info:
        parse_args(
            [
                "MTEYK",
                str(input_path),
                "--sequence-type",
                "amino-acid",
                "--revcomp",
            ]
        )

    assert exc_info.value.code == 2


# ---------------------------------------------------------------------------
# FASTA / FASTQ input
# ---------------------------------------------------------------------------


def test_read_fasta_plain_and_gzip(tmp_path: Path) -> None:
    fasta = tmp_path / "example.fa"
    fasta.write_text(">seq1\nATG\nCCC\n", encoding="utf-8")

    gz = tmp_path / "example.fa.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(">seq1\nATG\nCCC\n")

    expected = [FastaRecord("seq1", "ATGCCC")]

    assert list(FastaFileReader(fasta).read()) == expected
    assert list(FastaFileReader(gz).read()) == expected


def test_read_fastq_plain_and_gzip(tmp_path: Path) -> None:
    fastq_text = "@read1\nATGC\n+\n!!!!\n"

    fastq = tmp_path / "reads.fastq"
    fastq.write_text(fastq_text, encoding="utf-8")

    gz = tmp_path / "reads.fastq.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(fastq_text)

    expected = [FastaRecord("read1", "ATGC")]

    assert list(FastaFileReader(fastq).read()) == expected
    assert list(FastaFileReader(gz).read()) == expected
