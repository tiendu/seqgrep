from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from seqgrep.codecs import (
    ByteTarget,
    IupacDnaCodec,
    LiteralCodec,
    PackedCanonicalMaskTarget,
)
from seqgrep.fasta import FastaFileReader
from seqgrep.iupac import (
    encode_canonical_2bit,
    encode_iupac_query,
    reverse_complement,
)
from seqgrep.models import FastaRecord, Match, SearchQuery
from seqgrep.planner import SearchPlanner


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
    return [
        (hit.start, hit.end, hit.matched)
        for hit in hits
    ]


# ---------------------------------------------------------------------------
# IUPAC encoding
# ---------------------------------------------------------------------------


def test_iupac_query_bitmasks() -> None:
    assert encode_iupac_query("A") == bytes([0b0001])
    assert encode_iupac_query("C") == bytes([0b0010])
    assert encode_iupac_query("G") == bytes([0b0100])
    assert encode_iupac_query("T") == bytes([0b1000])

    assert encode_iupac_query("R") == bytes([0b0101])
    assert encode_iupac_query("Y") == bytes([0b1010])
    assert encode_iupac_query("N") == bytes([0b1111])


def test_iupac_query_normalizes_case_and_whitespace() -> None:
    assert encode_iupac_query("a c\ng\tt") == bytes(
        [
            0b0001,
            0b0010,
            0b0100,
            0b1000,
        ]
    )


def test_iupac_query_rejects_invalid_symbol() -> None:
    with pytest.raises(
        ValueError,
        match="Unsupported IUPAC nucleotide 'X' at position 3",
    ):
        encode_iupac_query("ATX")


def test_reverse_complement_ambiguous() -> None:
    assert reverse_complement("ARYN") == "NRYT"


# ---------------------------------------------------------------------------
# Canonical two-bit encoding
# ---------------------------------------------------------------------------


def test_encode_canonical_2bit() -> None:
    packed, length = encode_canonical_2bit("ACGT")

    assert length == 4
    assert packed == bytes([0b00011011])


def test_encode_canonical_2bit_partial_final_byte() -> None:
    packed, length = encode_canonical_2bit("ACGTA")

    assert length == 5
    assert packed == bytes(
        [
            0b00011011,
            0b00000000,
        ]
    )


def test_encode_canonical_2bit_treats_u_as_t() -> None:
    packed_t, length_t = encode_canonical_2bit("T")
    packed_u, length_u = encode_canonical_2bit("U")

    assert length_t == length_u == 1
    assert packed_t == packed_u


def test_encode_canonical_2bit_rejects_ambiguity() -> None:
    with pytest.raises(
        ValueError,
        match="2-bit encoding only supports A, C, G, T, or U",
    ):
        encode_canonical_2bit("ATGN")


# ---------------------------------------------------------------------------
# Encoded targets and codecs
# ---------------------------------------------------------------------------


def test_literal_codec_uses_ascii_bytes() -> None:
    codec = LiteralCodec()

    query = codec.encode_query("ARNX")
    target = codec.encode_target("ARNX")

    assert query == b"ARNX"
    assert isinstance(target, ByteTarget)
    assert len(target) == 4

    assert [
        target.symbol_at(index)
        for index in range(len(target))
    ] == list(b"ARNX")


def test_literal_codec_compares_symbols_exactly() -> None:
    codec = LiteralCodec()

    assert codec.compatible(ord("A"), ord("A")) is True
    assert codec.compatible(ord("A"), ord("N")) is False


def test_iupac_codec_packs_canonical_target() -> None:
    codec = IupacDnaCodec()
    target = codec.encode_target("ACGTU")

    assert isinstance(target, PackedCanonicalMaskTarget)
    assert len(target) == 5

    assert [
        target.symbol_at(index)
        for index in range(len(target))
    ] == [
        0b0001,
        0b0010,
        0b0100,
        0b1000,
        0b1000,
    ]


def test_iupac_codec_falls_back_for_ambiguous_target() -> None:
    codec = IupacDnaCodec()
    target = codec.encode_target("ANRY")

    assert isinstance(target, ByteTarget)
    assert len(target) == 4

    assert [
        target.symbol_at(index)
        for index in range(len(target))
    ] == [
        0b0001,
        0b1111,
        0b0101,
        0b1010,
    ]


def test_packed_target_rejects_padding_index() -> None:
    codec = IupacDnaCodec()
    target = codec.encode_target("ACG")

    assert isinstance(target, PackedCanonicalMaskTarget)

    with pytest.raises(IndexError):
        target.symbol_at(3)


def test_iupac_codec_compatibility() -> None:
    codec = IupacDnaCodec()

    r_mask = encode_iupac_query("R")[0]
    y_mask = encode_iupac_query("Y")[0]

    a_mask = encode_iupac_query("A")[0]
    c_mask = encode_iupac_query("C")[0]
    g_mask = encode_iupac_query("G")[0]
    n_mask = encode_iupac_query("N")[0]

    assert codec.compatible(r_mask, a_mask) is True
    assert codec.compatible(r_mask, g_mask) is True
    assert codec.compatible(r_mask, c_mask) is False

    assert codec.compatible(y_mask, c_mask) is True
    assert codec.compatible(r_mask, n_mask) is True


# ---------------------------------------------------------------------------
# Planner and serial matcher behavior
# ---------------------------------------------------------------------------


def test_literal_default_treats_n_as_literal() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN")

    hits = search_with_planner(record, query)

    assert hit_values(hits) == [
        (1, 4, "ATGN"),
    ]


def test_iupac_ambig_treats_n_as_any_base() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN", ambig=True)

    hits = search_with_planner(record, query)

    assert hit_values(hits) == [
        (1, 4, "ATGN"),
        (5, 8, "ATGA"),
    ]


def test_iupac_search_uses_packed_canonical_target() -> None:
    record = FastaRecord("seq1", "ATGAATGT")
    query = SearchQuery("ATGN", ambig=True)

    hits = search_with_planner(record, query)

    assert hit_values(hits) == [
        (1, 4, "ATGA"),
        (5, 8, "ATGT"),
    ]


def test_planner_default_behavior_is_literal() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN")

    hits = search_with_planner(record, query)

    assert hit_values(hits) == [
        (1, 4, "ATGN"),
    ]


def test_planner_ambig_behavior_is_iupac() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN", ambig=True)

    hits = search_with_planner(record, query)

    assert hit_values(hits) == [
        (1, 4, "ATGN"),
        (5, 8, "ATGA"),
    ]


def test_revcomp_match() -> None:
    record = FastaRecord("seq1", "CCCATGAAGTCCC")
    query = SearchQuery("ACTTCAT", revcomp=True)

    hits = search_with_planner(record, query)

    assert any(
        hit.strand == "-"
        and hit.start == 4
        and hit.end == 10
        and hit.matched == "ATGAAGT"
        for hit in hits
    )


def test_circular_match() -> None:
    record = FastaRecord("plasmid", "CCCAAATTT")
    query = SearchQuery("TTTCCC", circular=True)

    hits = search_with_planner(record, query)

    assert len(hits) == 1
    assert hits[0].record == "plasmid"
    assert hits[0].strand == "+"
    assert hits[0].start == 7
    assert hits[0].end == 3
    assert hits[0].matched == "TTTCCC"
    assert hits[0].circular is True


def test_circular_pattern_longer_than_sequence() -> None:
    record = FastaRecord("tiny", "ATG")
    query = SearchQuery("ATGAT", circular=True)

    hits = search_with_planner(record, query)

    assert [
        (
            hit.start,
            hit.end,
            hit.matched,
            hit.circular,
        )
        for hit in hits
    ] == [
        (1, 2, "ATGAT", True),
    ]


def test_revcomp_circular_match() -> None:
    record = FastaRecord("plasmid", "CCCAAATTT")
    query = SearchQuery(
        "GGGAAA",
        revcomp=True,
        circular=True,
    )

    hits = search_with_planner(record, query)

    assert len(hits) == 1
    assert hits[0].strand == "-"
    assert hits[0].start == 7
    assert hits[0].end == 3
    assert hits[0].matched == "TTTCCC"
    assert hits[0].circular is True


# ---------------------------------------------------------------------------
# Chunked matcher parity
# ---------------------------------------------------------------------------


def test_chunked_literal_matches_serial_across_internal_boundary() -> None:
    record = FastaRecord("seq1", "AAAAATGCAAAAA")
    query = SearchQuery("TGCA")

    serial_hits = search_with_planner(
        record,
        query,
        jobs=1,
    )
    chunked_hits = search_with_planner(
        record,
        query,
        jobs=2,
        chunk_size=6,
    )

    assert chunked_hits == serial_hits
    assert hit_values(chunked_hits) == [
        (6, 9, "TGCA"),
    ]


def test_chunked_iupac_ambiguous_target_matches_serial() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN", ambig=True)

    serial_hits = search_with_planner(
        record,
        query,
        jobs=1,
    )
    chunked_hits = search_with_planner(
        record,
        query,
        jobs=2,
        chunk_size=3,
    )

    assert chunked_hits == serial_hits


def test_chunked_iupac_packed_target_matches_serial() -> None:
    record = FastaRecord("seq1", "ATGAATGT")
    query = SearchQuery("ATGN", ambig=True)

    serial_hits = search_with_planner(
        record,
        query,
        jobs=1,
    )
    chunked_hits = search_with_planner(
        record,
        query,
        jobs=2,
        chunk_size=3,
    )

    assert chunked_hits == serial_hits


def test_chunked_revcomp_circular_match() -> None:
    record = FastaRecord("plasmid", "CCCAAATTT")
    query = SearchQuery(
        "GGGAAA",
        revcomp=True,
        circular=True,
    )

    serial_hits = search_with_planner(
        record,
        query,
        jobs=1,
    )
    chunked_hits = search_with_planner(
        record,
        query,
        jobs=2,
        chunk_size=4,
    )

    assert chunked_hits == serial_hits


def test_chunked_iupac_revcomp_circular_match() -> None:
    record = FastaRecord("plasmid", "CCCAAATTT")
    query = SearchQuery(
        "GGGAAA",
        revcomp=True,
        circular=True,
        ambig=True,
    )

    serial_hits = search_with_planner(
        record,
        query,
        jobs=1,
    )
    chunked_hits = search_with_planner(
        record,
        query,
        jobs=2,
        chunk_size=4,
    )

    assert chunked_hits == serial_hits


# ---------------------------------------------------------------------------
# FASTA / FASTQ input
# ---------------------------------------------------------------------------


def test_read_fasta_plain_and_gzip(tmp_path: Path) -> None:
    fasta = tmp_path / "example.fa"
    fasta.write_text(
        ">seq1\nATG\nCCC\n",
        encoding="utf-8",
    )

    gz = tmp_path / "example.fa.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(">seq1\nATG\nCCC\n")

    expected = [
        FastaRecord("seq1", "ATGCCC"),
    ]

    assert list(FastaFileReader(fasta).read()) == expected
    assert list(FastaFileReader(gz).read()) == expected


def test_read_fastq_plain_and_gzip(tmp_path: Path) -> None:
    fastq_text = "@read1\nATGC\n+\n!!!!\n"

    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        fastq_text,
        encoding="utf-8",
    )

    gz = tmp_path / "reads.fastq.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(fastq_text)

    expected = [
        FastaRecord("read1", "ATGC"),
    ]

    assert list(FastaFileReader(fastq).read()) == expected
    assert list(FastaFileReader(gz).read()) == expected

