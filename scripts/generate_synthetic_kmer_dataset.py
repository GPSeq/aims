#!/usr/bin/env python3
"""Generate deterministic FASTA/truth fixtures for AIMS k-mer retrieval."""

from __future__ import annotations

import argparse
import random
from pathlib import Path


DNA = "ACGT"


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(table)[::-1].upper()


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(DNA) for _ in range(length))


def mutate(seq: str, rng: random.Random, every: int) -> str:
    if every <= 0:
        return seq
    out = list(seq)
    for i in range(every - 1, len(out), every):
      choices = [base for base in DNA if base != out[i]]
      out[i] = rng.choice(choices)
    return "".join(out)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate deterministic AIMS k-mer fixtures.")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--documents", type=int, default=64)
    parser.add_argument("--length", type=int, default=1000)
    parser.add_argument("--queries", type=int, default=64)
    parser.add_argument("--query-length", type=int, default=150)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--mutation-period", type=int, default=0,
                        help="Mutate every Nth query base for noisy-query experiments; 0 keeps queries exact.")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    repeat = "ACGT" * 40 + "AAAAACCCCC" * 8
    references: list[tuple[str, str]] = []
    for doc_id in range(args.documents):
        core = random_dna(rng, max(1, args.length - len(repeat)))
        if doc_id % 8 == 1 and references:
            seq = reverse_complement(references[-1][1])
        elif doc_id % 8 == 2:
            seq = repeat + core
        elif doc_id % 8 == 3:
            seq = core[: args.length // 2] + "N" * 24 + core[args.length // 2 :]
        else:
            seq = core[: args.length // 3] + repeat + core[args.length // 3 :]
        references.append((f"doc{doc_id}", seq[: args.length]))

    queries: list[tuple[str, str]] = []
    truth_rows: list[tuple[str, str, str, str]] = []
    for qid in range(args.queries):
        if qid % 10 == 9:
            name = f"q{qid}_absent"
            seq = "N" * args.query_length
            expected = ""
        else:
            doc_id = qid % args.documents
            ref_name, ref_seq = references[doc_id]
            max_start = max(0, len(ref_seq) - args.query_length)
            start = rng.randint(0, max_start) if max_start > 0 else 0
            seq = ref_seq[start : start + args.query_length].replace("N", "A")
            if qid % 6 == 4:
                seq = reverse_complement(seq)
            if args.mutation_period > 0 and qid % 6 == 5:
                seq = mutate(seq, rng, every=args.mutation_period)
            name = f"q{qid}_{ref_name}"
            expected = str(doc_id)
        queries.append((name, seq))
        truth_rows.append((name, "exact_retrieval", expected, "synthetic"))

    write_fasta(args.out_dir / "refs.fa", references)
    write_fasta(args.out_dir / "queries.fa", queries)
    with (args.out_dir / "truth.tsv").open("w", encoding="utf-8") as handle:
        handle.write("query_id\tcomparison_stage\texpected_document_ids\tnote\n")
        for row in truth_rows:
            handle.write("\t".join(row) + "\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
