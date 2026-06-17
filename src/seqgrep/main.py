from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence

from .cli import parse_args
from .fasta import FastaFileReader
from .models import Match, SearchQuery, SequenceType
from .planner import SearchPlanner


logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

logger = logging.getLogger(__name__)


def print_match(hit: Match) -> None:
    print(
        f"{hit.record}\t{hit.strand}\t{hit.start}\t{hit.end}\t"
        f"{hit.matched}\t{str(hit.circular).lower()}"
    )


def run(args: argparse.Namespace) -> int:
    query = SearchQuery(
        pattern=args.pattern,
        revcomp=args.revcomp,
        circular=args.circular,
        ambig=args.ambig,
        sequence_type=SequenceType(args.sequence_type),
    )

    plan = SearchPlanner().plan(
        query=query,
        jobs=args.jobs,
        chunk_size=args.chunk_size,
    )

    matcher = plan.matcher
    reader = FastaFileReader(args.input, fmt=args.format)

    if args.with_header:
        print("record\tstrand\tstart\tend\tmatched\tcircular")

    for record in reader.read():
        for hit in matcher.search(record, query):
            print_match(hit)

    return 0


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)

    try:
        return run(args)
    except Exception:
        logger.exception("seqgrep failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
