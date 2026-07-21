# Benchmark Philosophy

AIMS benchmarks must be stage-matched and reproducible. A lookup-only measurement must not be compared as if it were full alignment. A routing-only result must remain labeled as routing-only.

## Fairness Rules

1. Declare `comparison_stage` for every result row.
2. Declare whether early stopping is disabled, conservative, or heuristic.
3. Report bytes read per query and postings decoded per query.
4. Keep truth inputs and dataset manifests versioned with checksums.
5. Separate cold-cache and warm-cache modes.
6. Report command lines, executable versions, thread counts, and host information.

## Required Metrics

1. build wall time and CPU time
2. peak RSS
3. serialized index size
4. p50/p90/p95/p99 latency
5. bytes read per query
6. postings decoded per query
7. seeds generated and queried per query
8. candidate count
9. top-1/top-5/top-10 recall against a declared truth manifest

## Current K-mer Benchmark Driver

`aims_bench` currently supports exact k-mer retrieval benchmarks:

```sh
aims_bench --ref refs.fa --query queries.fa --truth truth.tsv --k 15,19,23,27,31 --topk 10
```

Budgeted benchmark runs should include the exact options in the saved command. Relevant controls are `--max-seeds`, `--max-postings`, `--max-candidates`, and `--hot-seed-threshold`. The output reports skipped hot seeds and budgeted skips separately.

Loading mode is part of the benchmark configuration. Runs using `--mmap` or `--posting-cache-blocks` should not be mixed with copied-load runs unless the command line is reported.

The truth TSV format used by the bundled fixture is:

```text
query_id	comparison_stage	expected_document_ids	note
q0	exact_retrieval	0,2	example
```

Empty `expected_document_ids` means no exact candidate is expected for that query.

The per-query JSONL contract is documented in `benchmarks/schemas/query_result.schema.json`.
