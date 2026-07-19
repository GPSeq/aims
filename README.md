# AIMS (Adaptive Information-per-Memory-access Seed Index)

AIMS is a C++20 research system for exact and adaptive biological sequence retrieval over large redundant DNA collections.

The current implementation is a stage-matched exact retrieval baseline. It builds a multi-k canonical k-mer index, stores coordinate-aware positional postings, orders query seeds by estimated information per byte accessed, accumulates candidate sequences, and emits per-query instrumentation.

## Current Stage

`comparison_stage=exact_retrieval`

This repository does not currently implement full downstream alignment. The active implementation path is k-mer exact retrieval with multi-k planning, compressed postings, mmap loading, and explicit instrumentation.

## Build

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## CLI Examples

```sh
build/aims_build --ref refs.fa --out index.aims --k 15,19,23,27,31
build/aims_query --index index.aims --query queries.fa --topk 10 --emit jsonl
build/aims_validate --index index.aims
build/aims_dump --index index.aims
build/aims_bench --ref refs.fa --query queries.fa --k 15,19,23,27,31 --topk 10 --truth truth.tsv
```

Useful speed/memory controls:

```sh
build/aims_query --index index.aims --query queries.fa --topk 10 --max-seeds 200 --max-postings 50000 --max-candidates 100000 --hot-seed-threshold 10000
```

These options keep the stage labeled as `exact_retrieval`, but they change the retrieval policy. Skipped hot seeds and budgeted skips are reported explicitly in the metrics.

For large indexes, use the mmap reader and optional decoded posting-block cache:

```sh
build/aims_query --index index.aims --query queries.fa --mmap --posting-cache-blocks 4096 --topk 10
build/aims_bench --ref refs.fa --query queries.fa --truth truth.tsv --k 15,19,23 --mmap --posting-cache-blocks 4096
```

`--mmap` keeps compressed posting blocks file-backed instead of copying them into owned vectors at load time. `--posting-cache-blocks` caches recently decoded posting blocks and should be benchmarked because it trades memory for speed.

See [docs/cli_reference.md](/Users/ahmadlutfi/Downloads/aims/docs/cli_reference.md) for all current k-mer options.

## Synthetic Fixture

```sh
build/aims_build --ref tests/data/kmer_refs.fa --out /tmp/aims_kmer_fixture.aims --k 5,7
build/aims_query --index /tmp/aims_kmer_fixture.aims --query tests/data/kmer_queries.fa --topk 3 --emit jsonl
build/aims_bench --ref tests/data/kmer_refs.fa --query tests/data/kmer_queries.fa --truth tests/data/kmer_truth.tsv --k 5,7 --topk 3
```

`aims_query` emits per-query metrics, per-k accounting, and structured `topk_results`. `aims_bench` emits an aggregate `BenchmarkResult` JSON object with top-1/top-5/top-10 recall when a truth TSV is provided.

Generate a larger deterministic stress fixture:

```sh
python3 scripts/generate_synthetic_kmer_dataset.py --out-dir /tmp/aims_synth --documents 256 --length 5000 --queries 512 --query-length 250
build/aims_bench --ref /tmp/aims_synth/refs.fa --query /tmp/aims_synth/queries.fa --truth /tmp/aims_synth/truth.tsv --k 15,19,23 --topk 10
```
