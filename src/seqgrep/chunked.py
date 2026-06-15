from __future__ import annotations

from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from multiprocessing import shared_memory

from .codecs import SequenceCodec
from .models import FastaRecord, Match, SearchQuery


@dataclass(frozen=True)
class _ChunkJob:
    shm_name: str
    seq_len: int
    query_symbols: bytes
    start_begin: int
    start_end: int
    strand: str
    circular: bool
    codec: SequenceCodec


@dataclass(frozen=True)
class _RawHit:
    strand: str
    zero_start: int
    zero_end: int
    circular: bool


def _chunk_ranges(total_starts: int, chunk_size: int) -> Iterator[tuple[int, int]]:
    if chunk_size <= 0:
        raise ValueError("chunk_size must be greater than zero")

    for start in range(0, total_starts, chunk_size):
        yield start, min(start + chunk_size, total_starts)


def _window_matches(
    seq_symbols: memoryview,
    seq_len: int,
    query_symbols: bytes,
    zero_start: int,
    circular: bool,
    codec: SequenceCodec,
) -> bool:
    for offset, query_symbol in enumerate(query_symbols):
        target_index = (zero_start + offset) % seq_len if circular else zero_start + offset

        if not codec.compatible(query_symbol, seq_symbols[target_index]):
            return False

    return True


def _search_chunk(job: _ChunkJob) -> list[_RawHit]:
    shm = shared_memory.SharedMemory(name=job.shm_name)
    buf = shm.buf
    assert buf is not None

    seq_symbols = buf.cast("B")

    try:
        hits: list[_RawHit] = []
        query_len = len(job.query_symbols)

        for zero_start in range(job.start_begin, job.start_end):
            if _window_matches(
                seq_symbols=seq_symbols,
                seq_len=job.seq_len,
                query_symbols=job.query_symbols,
                zero_start=zero_start,
                circular=job.circular,
                codec=job.codec,
            ):
                zero_end = zero_start + query_len - 1

                hits.append(
                    _RawHit(
                        strand=job.strand,
                        zero_start=zero_start,
                        zero_end=zero_end,
                        circular=job.circular and zero_end >= job.seq_len,
                    )
                )

        return hits
    finally:
        seq_symbols.release()
        buf.release()
        shm.close()


class ChunkedProcessMatcher:
    """Search one long sequence by splitting candidate start positions.

    Each worker owns a range of possible match starts.

    The worker may read past an internal chunk boundary via shared memory, but
    it only reports matches whose start belongs to its assigned range.
    """

    def __init__(
        self,
        codec: SequenceCodec,
        workers: int | None = None,
        chunk_size: int = 1_000_000,
    ) -> None:
        if workers is not None and workers < 1:
            raise ValueError("workers must be at least 1")
        if chunk_size < 1:
            raise ValueError("chunk_size must be at least 1")

        self.codec = codec
        self.workers = workers
        self.chunk_size = chunk_size

    def search(self, record: FastaRecord, query: SearchQuery) -> Iterable[Match]:
        seq = self.codec.normalize(record.sequence)
        seq_len = len(seq)

        if len(query.pattern) == 0:
            raise ValueError("Pattern must not be empty")
        if seq_len == 0:
            return

        encoded = self.codec.encode(seq)

        shm = shared_memory.SharedMemory(create=True, size=len(encoded))
        buf = shm.buf
        assert buf is not None

        buf[: len(encoded)] = encoded
        buf.release()

        try:
            yield from self._search_shared(
                record_name=record.name,
                seq=seq,
                query=query,
                shm_name=shm.name,
                seq_len=seq_len,
            )
        finally:
            shm.close()
            shm.unlink()

    def _search_shared(
        self,
        record_name: str,
        seq: str,
        query: SearchQuery,
        shm_name: str,
        seq_len: int,
    ) -> Iterable[Match]:
        patterns = [(query.pattern, "+")]

        if query.revcomp:
            patterns.append((self.codec.reverse_complement(query.pattern), "-"))

        with ProcessPoolExecutor(max_workers=self.workers) as pool:
            for pattern, strand in patterns:
                normalized_pattern = self.codec.normalize(pattern)
                query_symbols = self.codec.encode(normalized_pattern)
                query_len = len(query_symbols)

                if query_len == 0:
                    raise ValueError("Pattern must not be empty")
                if query_len > seq_len and not query.circular:
                    continue

                total_starts = seq_len if query.circular else seq_len - query_len + 1
                if total_starts <= 0:
                    continue

                jobs = [
                    _ChunkJob(
                        shm_name=shm_name,
                        seq_len=seq_len,
                        query_symbols=query_symbols,
                        start_begin=start_begin,
                        start_end=start_end,
                        strand=strand,
                        circular=query.circular,
                        codec=self.codec,
                    )
                    for start_begin, start_end in _chunk_ranges(total_starts, self.chunk_size)
                ]

                for raw_hits in pool.map(_search_chunk, jobs):
                    for raw_hit in raw_hits:
                        yield Match(
                            record=record_name,
                            strand=raw_hit.strand,
                            start=raw_hit.zero_start + 1,
                            end=(raw_hit.zero_end % seq_len) + 1
                            if query.circular
                            else raw_hit.zero_end + 1,
                            matched=self._matched_sequence(
                                seq=seq,
                                start=raw_hit.zero_start,
                                length=query_len,
                                circular=query.circular,
                            ),
                            circular=raw_hit.circular,
                        )

    @staticmethod
    def _matched_sequence(seq: str, start: int, length: int, circular: bool) -> str:
        if not circular:
            return seq[start : start + length]

        seq_len = len(seq)
        return "".join(seq[(start + offset) % seq_len] for offset in range(length))
