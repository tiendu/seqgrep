from __future__ import annotations

from array import array
from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from multiprocessing import get_context, shared_memory

from .codecs import SequenceCodec, TargetEncoding, decode_target_symbol
from .models import FastaRecord, Match, SearchQuery


@dataclass(frozen=True, slots=True)
class _ChunkJob:
    shm_name: str
    stored_size: int
    sequence_length: int
    target_encoding: TargetEncoding
    query_symbols: bytes
    start_begin: int
    start_end: int
    circular: bool
    codec: SequenceCodec


def _chunk_ranges(total_starts: int, chunk_size: int) -> Iterator[tuple[int, int]]:
    if chunk_size < 1:
        raise ValueError("chunk_size must be at least 1")

    for start in range(0, total_starts, chunk_size):
        yield start, min(start + chunk_size, total_starts)


def _window_matches(
    *,
    storage: memoryview,
    sequence_length: int,
    target_encoding: TargetEncoding,
    query_symbols: bytes,
    zero_start: int,
    circular: bool,
    codec: SequenceCodec,
) -> bool:
    compatible = codec.compatible

    if circular:
        for offset, query_symbol in enumerate(query_symbols):
            target_symbol = decode_target_symbol(
                storage,
                (zero_start + offset) % sequence_length,
                sequence_length,
                target_encoding,
            )
            if not compatible(query_symbol, target_symbol):
                return False
        return True

    for offset, query_symbol in enumerate(query_symbols):
        target_symbol = decode_target_symbol(
            storage,
            zero_start + offset,
            sequence_length,
            target_encoding,
        )
        if not compatible(query_symbol, target_symbol):
            return False

    return True


def _search_chunk(job: _ChunkJob) -> array[int]:
    shm = shared_memory.SharedMemory(name=job.shm_name)
    buffer = shm.buf
    assert buffer is not None
    storage = buffer[: job.stored_size]

    try:
        hits = array("Q")
        for zero_start in range(job.start_begin, job.start_end):
            if not _window_matches(
                storage=storage,
                sequence_length=job.sequence_length,
                target_encoding=job.target_encoding,
                query_symbols=job.query_symbols,
                zero_start=zero_start,
                circular=job.circular,
                codec=job.codec,
            ):
                continue

            hits.append(zero_start)

        return hits
    finally:
        storage.release()
        buffer.release()
        shm.close()


class ChunkedProcessMatcher:
    """Search one long record using packed shared memory and worker processes."""

    def __init__(
        self,
        codec: SequenceCodec,
        workers: int,
        chunk_size: int = 1_000_000,
    ) -> None:
        if workers < 1:
            raise ValueError("workers must be at least 1")
        if chunk_size < 1:
            raise ValueError("chunk_size must be at least 1")

        self.codec = codec
        self.workers = workers
        self.chunk_size = chunk_size

    def search(self, record: FastaRecord, query: SearchQuery) -> Iterable[Match]:
        sequence = record.sequence
        normalized_pattern = query.pattern

        if not normalized_pattern:
            raise ValueError("Pattern must not be empty")
        if not sequence:
            return

        target = self.codec.encode_target(sequence)
        storage = target.data
        sequence_length = len(target)

        if sequence_length != len(sequence):
            raise ValueError("Encoded target length does not match sequence length")
        if not storage:
            return

        shm = shared_memory.SharedMemory(create=True, size=len(storage))
        buffer = shm.buf
        assert buffer is not None
        buffer[: len(storage)] = storage
        buffer.release()

        try:
            yield from self._search_shared(
                record_name=record.name,
                sequence=sequence,
                query=query,
                normalized_pattern=normalized_pattern,
                shm_name=shm.name,
                stored_size=len(storage),
                sequence_length=sequence_length,
                target_encoding=target.encoding,
            )
        finally:
            shm.close()
            shm.unlink()

    def _search_shared(
        self,
        *,
        record_name: str,
        sequence: str,
        query: SearchQuery,
        normalized_pattern: str,
        shm_name: str,
        stored_size: int,
        sequence_length: int,
        target_encoding: TargetEncoding,
    ) -> Iterable[Match]:
        patterns = [(normalized_pattern, "+")]
        if query.revcomp:
            patterns.append((self.codec.reverse_complement(normalized_pattern), "-"))

        with ProcessPoolExecutor(
            max_workers=self.workers,
            mp_context=get_context("spawn"),
        ) as pool:
            for pattern, strand in patterns:
                query_symbols = self.codec.encode_query(pattern)
                query_length = len(query_symbols)

                if query_length > sequence_length and not query.circular:
                    continue

                total_starts = (
                    sequence_length if query.circular else sequence_length - query_length + 1
                )
                if total_starts <= 0:
                    continue

                jobs = (
                    _ChunkJob(
                        shm_name=shm_name,
                        stored_size=stored_size,
                        sequence_length=sequence_length,
                        target_encoding=target_encoding,
                        query_symbols=query_symbols,
                        start_begin=start_begin,
                        start_end=start_end,
                        circular=query.circular,
                        codec=self.codec,
                    )
                    for start_begin, start_end in _chunk_ranges(
                        total_starts,
                        self.chunk_size,
                    )
                )

                for zero_starts in pool.map(_search_chunk, jobs):
                    for zero_start in zero_starts:
                        zero_end = zero_start + query_length - 1
                        yield Match(
                            record=record_name,
                            strand=strand,
                            start=zero_start + 1,
                            end=(zero_end % sequence_length) + 1
                            if query.circular
                            else zero_end + 1,
                            matched=self._matched_sequence(
                                sequence=sequence,
                                start=zero_start,
                                length=query_length,
                                circular=query.circular,
                            ),
                            circular=query.circular and zero_end >= sequence_length,
                        )

    @staticmethod
    def _matched_sequence(
        *,
        sequence: str,
        start: int,
        length: int,
        circular: bool,
    ) -> str:
        if not circular:
            return sequence[start : start + length]

        sequence_length = len(sequence)
        return "".join(sequence[(start + offset) % sequence_length] for offset in range(length))
