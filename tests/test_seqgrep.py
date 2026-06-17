from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from seqgrep.alphabets import (
    NUCLEOTIDE_EXACT_CODES,
    NUCLEOTIDE_IUPAC_MASKS,
    NUCLEOTIDE_IUPAC_QUERY_MASKS,
    NUCLEOTIDE_IUPAC_TARGET_MASKS,
    PROTEIN_5BIT_CODES,
    PROTEIN_SYMBOLS,
    encode_nucleotide_exact,
    encode_nucleotide_iupac,
    encode_nucleotide_iupac_query,
    encode_nucleotide_iupac_target,
    encode_protein_exact,
    nucleotide_exact_text,
    pack_nucleotide_2bit,
    pack_nucleotide_iupac_masks_5bit,
    pack_protein_5bit,
    protein_exact_text,
    reverse_complement_nucleotide,
    unpack_nucleotide_2bit,
    unpack_nucleotide_iupac_mask_5bit,
    unpack_protein_5bit,
)
from seqgrep.chunked import ChunkedProcessMatcher
from seqgrep.cli import parse_args
from seqgrep.codecs import (
    IupacNucleotideCodec,
    NucleotideCodec,
    PackedTarget,
    ProteinCodec,
    TargetEncoding,
)
from seqgrep.exact import ExactMatcher
from seqgrep.fastx import FastxReader
from seqgrep.models import AmbigMode, SequenceRecord, Match, SearchQuery, SequenceType
from seqgrep.planner import SearchPlanner
from seqgrep.window import WindowMatcher


def search_with_planner(
    record: SequenceRecord,
    query: SearchQuery,
    *,
    jobs: int = 1,
    chunk_size: int = 1_000_000,
) -> list[Match]:
    plan = SearchPlanner().plan(query, jobs, chunk_size)
    return list(plan.matcher.search(record, query))


def hit_values(hits: list[Match]) -> list[tuple[int, int, str]]:
    return [(hit.start, hit.end, hit.matched) for hit in hits]


# ---------------------------------------------------------------------------
# Unified alphabet module
# ---------------------------------------------------------------------------


def test_exact_nucleotide_codes() -> None:
    expected = bytes(NUCLEOTIDE_EXACT_CODES[symbol] for symbol in "ACGTURYN-.")
    assert encode_nucleotide_exact("ACGTURYN-.") == expected


def test_iupac_nucleotide_query_masks() -> None:
    expected = bytes(NUCLEOTIDE_IUPAC_QUERY_MASKS[symbol] for symbol in "ACGTRYN-.")
    assert encode_nucleotide_iupac_query("ACGTRYN-.") == expected
    assert encode_nucleotide_iupac("ACGTRYN-.") == expected
    assert NUCLEOTIDE_IUPAC_MASKS is NUCLEOTIDE_IUPAC_QUERY_MASKS


def test_iupac_nucleotide_target_masks_are_asymmetric() -> None:
    assert encode_nucleotide_iupac_target("ACGTUNRY-.") == bytes(
        NUCLEOTIDE_IUPAC_TARGET_MASKS[symbol] for symbol in "ACGTUNRY-."
    )
    assert encode_nucleotide_iupac_target("NRY") == bytes([0, 0, 0])


def test_iupac_five_bit_masks_round_trip_for_symmetric_targets() -> None:
    sequence = "ACGTRYSWKMBDHVN-."
    packed, length = pack_nucleotide_iupac_masks_5bit(sequence)

    assert len(packed) == (length * 5 + 7) // 8
    assert [
        unpack_nucleotide_iupac_mask_5bit(packed, index, length) for index in range(length)
    ] == [NUCLEOTIDE_IUPAC_QUERY_MASKS[symbol] for symbol in sequence]


def test_nucleotide_encoders_normalize_case_and_whitespace() -> None:
    assert encode_nucleotide_exact("a c\ng\tt") == bytes([0, 1, 2, 3])
    assert encode_nucleotide_iupac("a c\ng\tt") == bytes([1, 2, 4, 8])


def test_nucleotide_encoders_reject_invalid_symbol() -> None:
    with pytest.raises(ValueError, match="Unsupported nucleotide symbol 'X' at position 3"):
        encode_nucleotide_exact("ATX")

    with pytest.raises(
        ValueError,
        match="Unsupported IUPAC nucleotide query symbol 'X' at position 3",
    ):
        encode_nucleotide_iupac_query("ATX")

    with pytest.raises(
        ValueError,
        match="Unsupported IUPAC nucleotide target symbol 'X' at position 3",
    ):
        encode_nucleotide_iupac_target("ATX")


def test_nucleotide_exact_text_treats_u_as_t() -> None:
    assert nucleotide_exact_text("AUG") == "ATG"


def test_reverse_complement_supports_ambiguity_and_gaps() -> None:
    assert reverse_complement_nucleotide("ARYN-.") == ".-NRYT"


def test_nucleotide_2bit_round_trip() -> None:
    sequence = "ACGTUACGT"
    packed, length = pack_nucleotide_2bit(sequence)

    assert length == len(sequence)
    assert len(packed) == (length + 3) // 4
    assert [unpack_nucleotide_2bit(packed, index, length) for index in range(length)] == [
        NUCLEOTIDE_EXACT_CODES[symbol] for symbol in sequence
    ]


def test_nucleotide_2bit_rejects_noncanonical_symbols() -> None:
    with pytest.raises(
        ValueError,
        match="2-bit nucleotide encoding only supports A, C, G, T, or U",
    ):
        pack_nucleotide_2bit("ATGN")


def test_protein_alphabet_fits_five_bits() -> None:
    assert len(PROTEIN_5BIT_CODES) == len(PROTEIN_SYMBOLS)
    assert len(PROTEIN_5BIT_CODES) <= 32
    assert max(PROTEIN_5BIT_CODES.values()) <= 0b11111


def test_protein_query_codes() -> None:
    sequence = "MTEYKLVVVG"
    assert encode_protein_exact(sequence) == bytes(
        PROTEIN_5BIT_CODES[symbol] for symbol in sequence
    )


def test_protein_5bit_round_trip_full_alphabet() -> None:
    sequence = "".join(PROTEIN_SYMBOLS)
    packed, length = pack_protein_5bit(sequence)

    assert len(packed) == (length * 5 + 7) // 8
    assert [unpack_protein_5bit(packed, index, length) for index in range(length)] == [
        PROTEIN_5BIT_CODES[symbol] for symbol in sequence
    ]


def test_protein_5bit_round_trip_all_boundaries() -> None:
    for length in range(1, 65):
        sequence = "".join(PROTEIN_SYMBOLS[index % len(PROTEIN_SYMBOLS)] for index in range(length))
        packed, encoded_length = pack_protein_5bit(sequence)
        assert [
            unpack_protein_5bit(packed, index, encoded_length) for index in range(encoded_length)
        ] == [PROTEIN_5BIT_CODES[symbol] for symbol in sequence]


def test_protein_validation() -> None:
    assert protein_exact_text("m t\ne\tyk") == "MTEYK"

    with pytest.raises(
        ValueError,
        match="Unsupported amino-acid symbol '@' at position 3",
    ):
        pack_protein_5bit("MT@K")


# ---------------------------------------------------------------------------
# Codecs and storage
# ---------------------------------------------------------------------------


def test_exact_nucleotide_codec_packs_canonical_target() -> None:
    target = NucleotideCodec().encode_target("ACGTU")

    assert isinstance(target, PackedTarget)
    assert target.encoding is TargetEncoding.NUCLEOTIDE_2BIT_CODE
    assert len(target) == 5
    assert len(target.data) == 2
    assert [target.symbol_at(index) for index in range(len(target))] == [0, 1, 2, 3, 3]


def test_exact_nucleotide_codec_uses_five_bits_for_iupac_target() -> None:
    sequence = "ANRY-."
    target = NucleotideCodec().encode_target(sequence)

    assert isinstance(target, PackedTarget)
    assert target.encoding is TargetEncoding.NUCLEOTIDE_EXACT_5BIT
    assert len(target.data) == (len(sequence) * 5 + 7) // 8
    assert [target.symbol_at(index) for index in range(len(target))] == [
        NUCLEOTIDE_EXACT_CODES[symbol] for symbol in sequence
    ]


def test_iupac_codec_packs_canonical_target_as_masks() -> None:
    target = IupacNucleotideCodec().encode_target("ACGTU")

    assert isinstance(target, PackedTarget)
    assert target.encoding is TargetEncoding.NUCLEOTIDE_2BIT_MASK
    assert [target.symbol_at(index) for index in range(len(target))] == [1, 2, 4, 8, 8]


def test_iupac_codec_uses_validity_and_gap_bitmaps_for_mixed_target() -> None:
    sequence = "ANRY-."
    target = IupacNucleotideCodec().encode_target(sequence)

    assert isinstance(target, PackedTarget)
    assert target.encoding is TargetEncoding.NUCLEOTIDE_2BIT_VALID_GAP_MASK
    assert [target.symbol_at(index) for index in range(len(target))] == [
        NUCLEOTIDE_IUPAC_TARGET_MASKS[symbol] for symbol in sequence
    ]


def test_iupac_codec_uses_full_masks_for_symmetric_target() -> None:
    sequence = "ANRY-."
    target = IupacNucleotideCodec(allow_target_ambiguity=True).encode_target(sequence)

    assert isinstance(target, PackedTarget)
    assert target.encoding is TargetEncoding.NUCLEOTIDE_IUPAC_5BIT_MASK
    assert [target.symbol_at(index) for index in range(len(target))] == [
        NUCLEOTIDE_IUPAC_QUERY_MASKS[symbol] for symbol in sequence
    ]


def test_iupac_codec_uses_three_bits_per_base_without_gaps() -> None:
    sequence = "ACGT" + "N" * 12
    target = IupacNucleotideCodec().encode_target(sequence)

    assert target.encoding is TargetEncoding.NUCLEOTIDE_2BIT_VALID_MASK
    assert len(target.data) == (len(sequence) + 3) // 4 + (len(sequence) + 7) // 8


def test_protein_codec_uses_packed_5bit_target() -> None:
    sequence = "MTEYKLVVVGAGGVGKSAL"
    target = ProteinCodec().encode_target(sequence)

    assert isinstance(target, PackedTarget)
    assert target.encoding is TargetEncoding.PROTEIN_5BIT
    assert len(target.data) == (len(sequence) * 5 + 7) // 8
    assert [target.symbol_at(index) for index in range(len(target))] == [
        PROTEIN_5BIT_CODES[symbol] for symbol in sequence
    ]


def test_packed_targets_reject_padding_indexes() -> None:
    nucleotide_target = NucleotideCodec().encode_target("ACG")
    protein_target = ProteinCodec().encode_target("MTE")

    with pytest.raises(IndexError):
        nucleotide_target.symbol_at(3)
    with pytest.raises(IndexError):
        protein_target.symbol_at(3)


# ---------------------------------------------------------------------------
# Planner and exact serial backend
# ---------------------------------------------------------------------------


def test_planner_uses_native_exact_matcher_for_serial_exact_modes() -> None:
    nucleotide = SearchPlanner().plan(SearchQuery("ATG"), jobs=1, chunk_size=100)
    protein = SearchPlanner().plan(
        SearchQuery("MTE", sequence_type=SequenceType.AMINO_ACID),
        jobs=1,
        chunk_size=100,
    )

    assert isinstance(nucleotide.matcher, ExactMatcher)
    assert isinstance(protein.matcher, ExactMatcher)


def test_planner_uses_window_matcher_for_serial_iupac() -> None:
    plan = SearchPlanner().plan(SearchQuery("ATN", ambig=True), jobs=1, chunk_size=100)
    assert isinstance(plan.matcher, WindowMatcher)


def test_search_query_ambig_alias_selects_query_mode() -> None:
    query = SearchQuery("ATN", ambig=True)

    assert query.ambig is True
    assert query.ambig_mode is AmbigMode.QUERY


def test_planner_uses_symmetric_iupac_codec_for_both_mode() -> None:
    query = SearchQuery("AAA", ambig_mode=AmbigMode.BOTH)
    plan = SearchPlanner().plan(query, jobs=1, chunk_size=100)

    assert isinstance(plan.matcher, WindowMatcher)
    assert isinstance(plan.matcher.codec, IupacNucleotideCodec)
    assert plan.matcher.codec.allow_target_ambiguity is True


def test_planner_uses_chunked_matcher_when_jobs_exceed_one() -> None:
    plan = SearchPlanner().plan(SearchQuery("ATG"), jobs=2, chunk_size=100)
    assert isinstance(plan.matcher, ChunkedProcessMatcher)


def test_default_mode_is_exact_nucleotide() -> None:
    hits = search_with_planner(SequenceRecord("seq1", "ATGNATGA"), SearchQuery("ATGN"))
    assert hit_values(hits) == [(1, 4, "ATGN")]


def test_exact_matcher_finds_overlapping_hits() -> None:
    hits = search_with_planner(SequenceRecord("seq1", "AAAA"), SearchQuery("AA"))
    assert hit_values(hits) == [(1, 2, "AA"), (2, 3, "AA"), (3, 4, "AA")]


def test_default_nucleotide_mode_treats_u_as_t() -> None:
    hits = search_with_planner(SequenceRecord("rna", "AUGATG"), SearchQuery("ATG"))
    assert hit_values(hits) == [(1, 3, "AUG"), (4, 6, "ATG")]


def test_exact_nucleotide_gaps_are_literal() -> None:
    record = SequenceRecord("aligned", "A-C.A.C")

    assert hit_values(search_with_planner(record, SearchQuery("A-C"))) == [(1, 3, "A-C")]
    assert hit_values(search_with_planner(record, SearchQuery("A.C"))) == [
        (5, 7, "A.C"),
    ]


def test_iupac_gaps_are_equivalent() -> None:
    record = SequenceRecord("aligned", "A-C A.C")
    query = SearchQuery("A-C", ambig=True)

    assert hit_values(search_with_planner(record, query)) == [
        (1, 3, "A-C"),
        (4, 6, "A.C"),
    ]


def test_iupac_query_n_matches_any_canonical_base() -> None:
    hits = search_with_planner(
        SequenceRecord("seq1", "ACGT"),
        SearchQuery("N", ambig=True),
    )
    assert hit_values(hits) == [
        (1, 1, "A"),
        (2, 2, "C"),
        (3, 3, "G"),
        (4, 4, "T"),
    ]


def test_iupac_target_n_does_not_act_as_wildcard() -> None:
    hits = search_with_planner(
        SequenceRecord("seq1", "ATGNATGA"),
        SearchQuery("ATGN", ambig=True),
    )
    assert hit_values(hits) == [(5, 8, "ATGA")]


def test_iupac_query_does_not_match_unknown_target_region() -> None:
    record = SequenceRecord("assembly_gap", "N" * 256)

    assert search_with_planner(record, SearchQuery("ACGT", ambig=True)) == []
    assert search_with_planner(record, SearchQuery("NNNN", ambig=True)) == []


def test_both_mode_matches_ambiguous_target_symbols() -> None:
    hits = search_with_planner(
        SequenceRecord("uncertain", "NNANN"),
        SearchQuery("AAAAA", ambig_mode=AmbigMode.BOTH),
    )

    assert hit_values(hits) == [(1, 5, "NNANN")]


def test_both_mode_matches_n_query_to_n_target() -> None:
    hits = search_with_planner(
        SequenceRecord("uncertain", "NNNNN"),
        SearchQuery("NNNNN", ambig_mode=AmbigMode.BOTH),
    )

    assert hit_values(hits) == [(1, 5, "NNNNN")]


def test_query_mode_still_rejects_ambiguous_target_symbols() -> None:
    hits = search_with_planner(
        SequenceRecord("uncertain", "NNANN"),
        SearchQuery("AAAAA", ambig_mode=AmbigMode.QUERY),
    )

    assert hits == []


def test_exact_mode_still_matches_target_n_literally() -> None:
    hits = search_with_planner(SequenceRecord("gap", "NNNN"), SearchQuery("NN"))
    assert hit_values(hits) == [(1, 2, "NN"), (2, 3, "NN"), (3, 4, "NN")]


def test_protein_mode_matches_exact_amino_acids() -> None:
    hits = search_with_planner(
        SequenceRecord("ras", "XXMTEYKLVVVGAGGVGKSALXX"),
        SearchQuery("MTEYK", sequence_type=SequenceType.AMINO_ACID),
    )
    assert hit_values(hits) == [(3, 7, "MTEYK")]


def test_protein_extended_symbols_are_literal() -> None:
    hits = search_with_planner(
        SequenceRecord("protein", "MAXMAB"),
        SearchQuery("MAX", sequence_type=SequenceType.AMINO_ACID),
    )
    assert hit_values(hits) == [(1, 3, "MAX")]


def test_protein_mode_rejects_ambig_and_revcomp() -> None:
    with pytest.raises(ValueError, match="--ambig-mode is only valid"):
        SearchPlanner().plan(
            SearchQuery("MTE", ambig=True, sequence_type=SequenceType.AMINO_ACID),
            jobs=1,
            chunk_size=100,
        )

    with pytest.raises(ValueError, match="--revcomp is only valid"):
        SearchPlanner().plan(
            SearchQuery("MTE", revcomp=True, sequence_type=SequenceType.AMINO_ACID),
            jobs=1,
            chunk_size=100,
        )


def test_reverse_complement_match() -> None:
    hits = search_with_planner(
        SequenceRecord("seq1", "CCCATGAAGTCCC"),
        SearchQuery("ACTTCAT", revcomp=True),
    )
    assert any(
        hit.strand == "-" and hit.start == 4 and hit.end == 10 and hit.matched == "ATGAAGT"
        for hit in hits
    )


def test_circular_exact_match() -> None:
    hits = search_with_planner(
        SequenceRecord("plasmid", "CCCAAATTT"),
        SearchQuery("TTTCCC", circular=True),
    )
    assert hits == [
        Match(
            record="plasmid",
            strand="+",
            start=7,
            end=3,
            matched="TTTCCC",
            circular=True,
        )
    ]


def test_circular_pattern_longer_than_sequence() -> None:
    hits = search_with_planner(
        SequenceRecord("tiny", "ATG"),
        SearchQuery("ATGAT", circular=True),
    )
    assert hit_values(hits) == [(1, 2, "ATGAT")]


def test_circular_protein_match() -> None:
    hits = search_with_planner(
        SequenceRecord("ring", "CDEMA"),
        SearchQuery("MAC", circular=True, sequence_type=SequenceType.AMINO_ACID),
    )
    assert hit_values(hits) == [(4, 1, "MAC")]


# ---------------------------------------------------------------------------
# Chunked parity across all storage representations
# ---------------------------------------------------------------------------


def assert_chunked_matches_serial(
    record: SequenceRecord,
    query: SearchQuery,
    *,
    chunk_size: int,
) -> list[Match]:
    serial_hits = search_with_planner(record, query, jobs=1)
    chunked_hits = search_with_planner(record, query, jobs=2, chunk_size=chunk_size)
    assert chunked_hits == serial_hits
    return chunked_hits


def test_chunked_exact_nucleotide_packed_target() -> None:
    hits = assert_chunked_matches_serial(
        SequenceRecord("seq1", "AAAAATGCAAAAA"),
        SearchQuery("TGCA"),
        chunk_size=6,
    )
    assert hit_values(hits) == [(6, 9, "TGCA")]


def test_chunked_exact_nucleotide_five_bit_target() -> None:
    assert_chunked_matches_serial(
        SequenceRecord("seq1", "ATGNATGA"),
        SearchQuery("ATGN"),
        chunk_size=3,
    )


def test_chunked_iupac_packed_target() -> None:
    assert_chunked_matches_serial(
        SequenceRecord("seq1", "ATGAATGT"),
        SearchQuery("ATGN", ambig=True),
        chunk_size=3,
    )


def test_chunked_iupac_validity_and_gap_bitmaps() -> None:
    assert_chunked_matches_serial(
        SequenceRecord("seq1", "A-CATGNA.C"),
        SearchQuery("A.C", ambig=True),
        chunk_size=2,
    )


def test_chunked_iupac_does_not_match_unknown_target_region() -> None:
    hits = assert_chunked_matches_serial(
        SequenceRecord("assembly_gap", "N" * 256),
        SearchQuery("NNNN", ambig=True),
        chunk_size=31,
    )
    assert hits == []


def test_chunked_symmetric_iupac_target_matches_serial() -> None:
    hits = assert_chunked_matches_serial(
        SequenceRecord("uncertain", "NNANN"),
        SearchQuery("AAAAA", ambig_mode=AmbigMode.BOTH),
        chunk_size=2,
    )
    assert hit_values(hits) == [(1, 5, "NNANN")]


def test_chunked_protein_5bit_target() -> None:
    hits = assert_chunked_matches_serial(
        SequenceRecord("protein", "AAAAAMTEYKAAAAA"),
        SearchQuery("MTEYK", sequence_type=SequenceType.AMINO_ACID),
        chunk_size=6,
    )
    assert hit_values(hits) == [(6, 10, "MTEYK")]


def test_chunked_circular_and_reverse_complement() -> None:
    assert_chunked_matches_serial(
        SequenceRecord("plasmid", "CCCAAATTT"),
        SearchQuery("GGGAAA", revcomp=True, circular=True, ambig=True),
        chunk_size=4,
    )


def test_chunked_longer_circular_pattern() -> None:
    assert_chunked_matches_serial(
        SequenceRecord("tiny", "ATG"),
        SearchQuery("ATGAT", circular=True),
        chunk_size=2,
    )


def test_models_normalize_symbols_once_at_the_boundary() -> None:
    assert SequenceRecord("seq", "a c\ng").sequence == "ACG"
    assert SearchQuery("a t\ng").pattern == "ATG"


# ---------------------------------------------------------------------------
# CLI and readers
# ---------------------------------------------------------------------------


def test_cli_defaults_to_nucleotide(tmp_path: Path) -> None:
    args = parse_args(["ATG", str(tmp_path / "seq.fa")])
    assert args.sequence_type == SequenceType.NUCLEOTIDE.value
    assert args.ambig is False
    assert args.ambig_mode == AmbigMode.NONE.value


def test_cli_ambig_alias_selects_query_mode(tmp_path: Path) -> None:
    args = parse_args(["ATN", str(tmp_path / "seq.fa"), "--ambig"])

    assert args.ambig is True
    assert args.ambig_mode == AmbigMode.QUERY.value


def test_cli_accepts_explicit_ambig_modes(tmp_path: Path) -> None:
    path = str(tmp_path / "seq.fa")

    query_args = parse_args(["ATN", path, "--ambig-mode", "query"])
    both_args = parse_args(["AAA", path, "--ambig-mode", "both"])

    assert query_args.ambig_mode == AmbigMode.QUERY.value
    assert both_args.ambig_mode == AmbigMode.BOTH.value


def test_cli_rejects_conflicting_ambig_options(tmp_path: Path) -> None:
    with pytest.raises(SystemExit):
        parse_args(
            [
                "ATN",
                str(tmp_path / "seq.fa"),
                "--ambig",
                "--ambig-mode",
                "both",
            ]
        )


def test_cli_accepts_amino_acid_mode(tmp_path: Path) -> None:
    args = parse_args(["MTEYK", str(tmp_path / "protein.fa"), "--sequence-type", "amino-acid"])
    assert args.sequence_type == SequenceType.AMINO_ACID.value


def test_cli_rejects_invalid_amino_acid_options(tmp_path: Path) -> None:
    path = str(tmp_path / "protein.fa")

    with pytest.raises(SystemExit):
        parse_args(["MTE", path, "-t", "amino-acid", "--ambig"])
    with pytest.raises(SystemExit):
        parse_args(["MTE", path, "-t", "amino-acid", "--ambig-mode", "both"])
    with pytest.raises(SystemExit):
        parse_args(["MTE", path, "-t", "amino-acid", "--revcomp"])


def test_cli_rejects_invalid_worker_values(tmp_path: Path) -> None:
    path = str(tmp_path / "seq.fa")

    with pytest.raises(SystemExit):
        parse_args(["ATG", path, "--jobs", "0"])
    with pytest.raises(SystemExit):
        parse_args(["ATG", path, "--chunk-size", "0"])


def test_read_fasta_plain_and_gzip(tmp_path: Path) -> None:
    fasta = tmp_path / "example.fa"
    fasta.write_text(">seq1 description\nATG\nCCC\n", encoding="utf-8")

    gz = tmp_path / "example.fa.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(">seq1 description\nATG\nCCC\n")

    upper_gz = tmp_path / "example.FA.GZ"
    with gzip.open(upper_gz, "wt", encoding="utf-8") as handle:
        handle.write(">seq1 description\nATG\nCCC\n")

    expected = [SequenceRecord("seq1", "ATGCCC")]
    assert list(FastxReader(fasta).read()) == expected
    assert list(FastxReader(gz).read()) == expected
    assert list(FastxReader(upper_gz).read()) == expected


def test_read_multiline_fastq_plain_and_gzip(tmp_path: Path) -> None:
    fastq_text = "@read1 description\nAT\nGC\n+\n!!\n!!\n"

    fastq = tmp_path / "reads.fastq"
    fastq.write_text(fastq_text, encoding="utf-8")

    gz = tmp_path / "reads.fastq.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(fastq_text)

    expected = [SequenceRecord("read1", "ATGC")]
    assert list(FastxReader(fastq).read()) == expected
    assert list(FastxReader(gz).read()) == expected


def test_fastq_rejects_quality_length_mismatch(tmp_path: Path) -> None:
    fastq = tmp_path / "bad.fastq"
    fastq.write_text("@read1\nATGC\n+\n!!!\n", encoding="utf-8")

    with pytest.raises(ValueError, match="ended before quality block was complete"):
        list(FastxReader(fastq).read())
