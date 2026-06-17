from __future__ import annotations

from collections.abc import Iterable, Iterator
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from multiprocessing import shared_memory

from .codecs import SequenceCodec, TargetEncoding
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
    strand: str
    circular: bool
    codec: SequenceCodec


@dataclass(frozen=True, slots=True)
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


def _target_symbol(
    storage: memoryview,
    index: int,
    encoding: TargetEncoding,
) -> int:
    """Read one logical target symbol from shared-memory storage."""
    if encoding == TargetEncoding.BYTE:
        return storage[index]

    if encoding == TargetEncoding.PACKED_CANONICAL_2BIT:
        byte_index = index // 4
        shift = 6 - (index % 4) * 2
        code = (storage[byte_index] >> shift) & 0b11

        # IUPAC compatibility expects singleton masks, not raw two-bit codes.
        return 1 << code

    raise ValueError(f"Unsupported target encoding: {encoding}")


def _window_matches(
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
            target_index = (zero_start + offset) % sequence_length
            target_symbol = _target_symbol(
                storage=storage,
                index=target_index,
                encoding=target_encoding,
            )

            if not compatible(query_symbol, target_symbol):
                return False

        return True

    for offset, query_symbol in enumerate(query_symbols):
        target_index = zero_start + offset
        target_symbol = _target_symbol(
            storage=storage,
            index=target_index,
            encoding=target_encoding,
        )

        if not compatible(query_symbol, target_symbol):
            return False

    return True


def _search_chunk(job: _ChunkJob) -> list[_RawHit]:
    shm = shared_memory.SharedMemory(name=job.shm_name)
    buf = shm.buf
    assert buf is not None

    storage = buf[: job.stored_size]

    try:
        hits: list[_RawHit] = []
        query_len = len(job.query_symbols)

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

            zero_end = zero_start + query_len - 1
            hits.append(
                _RawHit(
                    strand=job.strand,
                    zero_start=zero_start,
                    zero_end=zero_end,
                    circular=job.circular and zero_end >= job.sequence_length,
                )
            )

        return hits
    finally:
        storage.release()
        buf.release()
        shm.close()


class ChunkedProcessMatcher:
    """Search one long sequence by splitting candidate start positions.

    Workers share the physical target storage while candidate ranges are based
    on the target's logical sequence length. This distinction is essential for
    packed two-bit DNA, where four logical bases occupy one stored byte.
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
        normalized_pattern = self.codec.normalize(query.pattern)

        if not normalized_pattern:
            raise ValueError("Pattern must not be empty")
        if not seq:
            return

        target = self.codec.encode_target(seq)
        storage = target.data
        sequence_length = len(target)

        if sequence_length != len(seq):
            raise ValueError(
                "Encoded target length does not match normalized sequence length"
            )
        if not storage:
            return

        shm = shared_memory.SharedMemory(create=True, size=len(storage))
        buf = shm.buf
        assert buf is not None

        buf[: len(storage)] = storage
        buf.release()

        try:
            yield from self._search_shared(
                record_name=record.name,
                seq=seq,
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
        record_name: str,
        seq: str,
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

        with ProcessPoolExecutor(max_workers=self.workers) as pool:
            for pattern, strand in patterns:
                query_symbols = self.codec.encode_query(pattern)
                query_len = len(query_symbols)

                if not query_symbols:
                    raise ValueError("Pattern must not be empty")
                if query_len > sequence_length and not query.circular:
                    continue

                total_starts = (
                    sequence_length
                    if query.circular
                    else sequence_length - query_len + 1
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
                        strand=strand,
                        circular=query.circular,
                        codec=self.codec,
                    )
                    for start_begin, start_end in _chunk_ranges(
                        total_starts,
                        self.chunk_size,
                    )
                )

                for raw_hits in pool.map(_search_chunk, jobs):
                    for raw_hit in raw_hits:
                        yield Match(
                            record=record_name,
                            strand=raw_hit.strand,
                            start=raw_hit.zero_start + 1,
                            end=(raw_hit.zero_end % sequence_length) + 1
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
    def _matched_sequence(
        seq: str,
        start: int,
        length: int,
        circular: bool,
    ) -> str:
        if not circular:
            return seq[start : start + length]

        sequence_length = len(seq)
        return "".join(
            seq[(start + offset) % sequence_length]
            for offset in range(length)
        )
