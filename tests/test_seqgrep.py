from __future__ import annotations

import gzip
from pathlib import Path

from seqgrep.fasta import FastaFileReader
from seqgrep.iupac import encode_iupac, reverse_complement
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


def test_iupac_bitmasks() -> None:
    assert encode_iupac("R") == bytes([0b0101])
    assert encode_iupac("Y") == bytes([0b1010])
    assert encode_iupac("N") == bytes([0b1111])


def test_reverse_complement_ambiguous() -> None:
    assert reverse_complement("ARYN") == "NRYT"


def test_literal_default_treats_n_as_literal() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN")

    hits = search_with_planner(record, query)

    assert [(hit.start, hit.end, hit.matched) for hit in hits] == [
        (1, 4, "ATGN")
    ]


def test_iupac_ambig_treats_n_as_any_base() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN", ambig=True)

    hits = search_with_planner(record, query)

    assert [(hit.start, hit.end, hit.matched) for hit in hits] == [
        (1, 4, "ATGN"),
        (5, 8, "ATGA"),
    ]


def test_planner_default_behavior_is_literal() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN")

    hits = search_with_planner(record, query)

    assert [(hit.start, hit.end, hit.matched) for hit in hits] == [
        (1, 4, "ATGN")
    ]


def test_planner_ambig_behavior_is_iupac() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN", ambig=True)

    hits = search_with_planner(record, query)

    assert [(hit.start, hit.end, hit.matched) for hit in hits] == [
        (1, 4, "ATGN"),
        (5, 8, "ATGA"),
    ]


def test_revcomp_match() -> None:
    record = FastaRecord("seq1", "CCCATGAAGTCCC")
    query = SearchQuery("ACTTCAT", revcomp=True)

    hits = search_with_planner(record, query)

    assert any(hit.strand == "-" and hit.start == 4 and hit.end == 10 for hit in hits)


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

    assert [(hit.start, hit.end, hit.matched, hit.circular) for hit in hits] == [
        (1, 2, "ATGAT", True)
    ]


def test_revcomp_circular_match() -> None:
    record = FastaRecord("plasmid", "CCCAAATTT")
    query = SearchQuery("GGGAAA", revcomp=True, circular=True)

    hits = search_with_planner(record, query)

    assert len(hits) == 1
    assert hits[0].strand == "-"
    assert hits[0].start == 7
    assert hits[0].end == 3
    assert hits[0].matched == "TTTCCC"
    assert hits[0].circular is True


def test_chunked_literal_matches_serial_across_internal_boundary() -> None:
    record = FastaRecord("seq1", "AAAAATGCAAAAA")
    query = SearchQuery("TGCA")

    serial_hits = search_with_planner(record, query, jobs=1)
    chunked_hits = search_with_planner(record, query, jobs=2, chunk_size=6)

    assert chunked_hits == serial_hits
    assert [(hit.start, hit.end, hit.matched) for hit in chunked_hits] == [
        (6, 9, "TGCA")
    ]


def test_chunked_iupac_matches_serial() -> None:
    record = FastaRecord("seq1", "ATGNATGA")
    query = SearchQuery("ATGN", ambig=True)

    serial_hits = search_with_planner(record, query, jobs=1)
    chunked_hits = search_with_planner(record, query, jobs=2, chunk_size=3)

    assert chunked_hits == serial_hits


def test_chunked_revcomp_circular_match() -> None:
    record = FastaRecord("plasmid", "CCCAAATTT")
    query = SearchQuery("GGGAAA", revcomp=True, circular=True)

    serial_hits = search_with_planner(record, query, jobs=1)
    chunked_hits = search_with_planner(record, query, jobs=2, chunk_size=4)

    assert chunked_hits == serial_hits


def test_chunked_iupac_revcomp_circular_match() -> None:
    record = FastaRecord("plasmid", "CCCAAATTT")
    query = SearchQuery("GGGAAA", revcomp=True, circular=True, ambig=True)

    serial_hits = search_with_planner(record, query, jobs=1)
    chunked_hits = search_with_planner(record, query, jobs=2, chunk_size=4)

    assert chunked_hits == serial_hits


def test_read_fasta_plain_and_gzip(tmp_path: Path) -> None:
    fasta = tmp_path / "example.fa"
    fasta.write_text(">seq1\nATG\nCCC\n", encoding="utf-8")

    gz = tmp_path / "example.fa.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(">seq1\nATG\nCCC\n")

    assert list(FastaFileReader(fasta).read()) == [FastaRecord("seq1", "ATGCCC")]
    assert list(FastaFileReader(gz).read()) == [FastaRecord("seq1", "ATGCCC")]


def test_read_fastq_plain_and_gzip(tmp_path: Path) -> None:
    fastq_text = "@read1\nATGC\n+\n!!!!\n"

    fastq = tmp_path / "reads.fastq"
    fastq.write_text(fastq_text, encoding="utf-8")

    gz = tmp_path / "reads.fastq.gz"
    with gzip.open(gz, "wt", encoding="utf-8") as handle:
        handle.write(fastq_text)

    assert list(FastaFileReader(fastq).read()) == [FastaRecord("read1", "ATGC")]
    assert list(FastaFileReader(gz).read()) == [FastaRecord("read1", "ATGC")]
