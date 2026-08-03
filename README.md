# AIMS (Adaptive Information-per-Memory-access Seed Index)

[![CI](https://github.com/GPSeq/aims/actions/workflows/cmake-single-platform.yml/badge.svg)](https://github.com/GPSeq/aims/actions/workflows/cmake-single-platform.yml). 
AIMS is a C++20 research system for exact and adaptive biological sequence retrieval over large redundant DNA collections.

The current implementation is a stage-matched exact retrieval baseline. It builds a multi-k canonical k-mer index, stores coordinate-aware positional postings, orders query seeds by estimated information per byte accessed, accumulates candidate sequences, and emits per-query instrumentation.

## Current Stage

`comparison_stage=exact_retrieval`

This repository does not currently implement full downstream alignment. The active implementation path is k-mer exact retrieval with multi-k planning, compressed postings, mmap loading, and explicit instrumentation.

## Build

Debug builds enable sanitizers by default and are useful for development. Use Release for timing or large index builds.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Release build:

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release -DAIMS_ENABLE_ASAN=OFF
cmake --build build-release --parallel
ctest --test-dir build-release --output-on-failure
```

Optional OpenMP build support parallelizes independent k-mer layers during index construction. `--threads N` only uses more than one worker when the binary was configured with `AIMS_ENABLE_OPENMP=ON`.

```sh
cmake -S . -B build-omp -DCMAKE_BUILD_TYPE=Release -DAIMS_ENABLE_ASAN=OFF -DAIMS_ENABLE_OPENMP=ON
cmake --build build-omp --parallel
```

On macOS with Homebrew `libomp`, CMake may need the runtime path explicitly:

```sh
brew install libomp
cmake -S . -B build-omp -DCMAKE_BUILD_TYPE=Release -DAIMS_ENABLE_ASAN=OFF -DAIMS_ENABLE_OPENMP=ON -DOpenMP_ROOT=/opt/homebrew/opt/libomp
```

## CLI Commands

Build a serialized exact k-mer index:

```sh
build-omp/aims_build --ref refs.fa --out index.aims --k 15,19,23,27,31 --threads 5
```

`aims_build` arguments:

| Argument | Required | Meaning |
| --- | --- | --- |
| `--ref refs.fa` | yes | FASTA or FASTQ reference file to index. |
| `--out index.aims` | yes | Output path for the serialized AIMS index. |
| `--k 15,19,23` | yes | Comma-separated k-mer lengths, each in `1..32`; each k becomes one independent layer. |
| `--chunk-size records` | no | Number of input records processed per streaming chunk, default `1024`. |
| `--frequency-thresholds rare,medium,hot` | no | Collection-frequency class thresholds stored in seed metadata; above hot is `very_hot`. |
| `--threads N` | no | OpenMP worker count for independent k layers, default `1`; requires `AIMS_ENABLE_OPENMP=ON` for `N > 1`. |

Query a serialized index:

```sh
build-release/aims_query --index index.aims --query queries.fa --topk 10 --emit jsonl --threads 4
```

`aims_query` arguments:

| Argument | Required | Meaning |
| --- | --- | --- |
| `--index index.aims` | yes | Serialized index produced by `aims_build` or the benchmark path. |
| `--query queries.fa` | yes | FASTA or FASTQ query records. |
| `--topk N` | no | Number of ranked candidates to emit per query, default `10`. |
| `--emit jsonl\|tsv` | no | Output format, default `jsonl`. |
| `--out path` | no | Write query metrics to a file instead of stdout. |
| `--mmap` | no | Memory-map the index so compressed posting blocks stay file-backed. |
| `--posting-cache-blocks N` | no | Cache recently decoded posting blocks; trades memory for repeated-query speed. |
| `--threads N` | no | Number of async query workers, default `1`. |
| `--max-seeds N` | no | Query seed budget; `0` means unlimited. |
| `--max-postings N` | no | Posting decode budget; `0` means unlimited. |
| `--max-candidates N` | no | Candidate accumulator budget; `0` means unlimited. |
| `--hot-seed-threshold N` | no | Skip or doc-only seeds above this collection-frequency threshold. |
| `--hot-seed-class hot\|very-hot` | no | Apply hot-seed policy using build-time frequency classes. |
| `--hot-mode skip\|doc-only` | no | Hot-seed policy, default `skip`. |

Validate and inspect an index:

```sh
build-release/aims_validate --index index.aims --mmap
build-release/aims_dump --index index.aims --mmap
```

`aims_validate` checks the checksum, offsets, k-mer layers, and reference metadata. `aims_dump` prints index metadata. Both accept `--index index.aims` and optional `--mmap`.

Run an end-to-end benchmark that builds an index, serializes it, reloads it, queries it, and emits one benchmark JSON object:

```sh
build-omp/aims_bench --ref refs.fa --query queries.fa --truth truth.tsv --k 15,19,23 --topk 10 --threads 3
```

`aims_bench` arguments:

| Argument | Required | Meaning |
| --- | --- | --- |
| `--ref refs.fa` | yes | FASTA or FASTQ reference file used to build the benchmark index. |
| `--query queries.fa` | yes | FASTA or FASTQ query records. |
| `--k 15,19,23` | yes | Comma-separated k-mer layers to build and query. |
| `--topk N` | no | Number of ranked candidates used for recall metrics, default `10`. |
| `--truth truth.tsv` | no | Truth table used for top-1/top-5/top-10 recall. |
| `--dataset name` | no | Dataset label in the benchmark JSON, default `synthetic_kmer_fixture`. |
| `--mmap` | no | Reload the temporary serialized index through the mmap reader. |
| `--posting-cache-blocks N` | no | Cache decoded posting blocks after loading. |
| `--threads N` | no | OpenMP worker count for the internal index build, default `1`; requires OpenMP for `N > 1`. |
| `--repeats N` | no | Repeat the query workload, default `1`. |
| `--query-metrics-out path` | no | Write per-query metrics JSONL. |
| `--max-seeds N` | no | Query seed budget; `0` means unlimited. |
| `--max-postings N` | no | Posting decode budget; `0` means unlimited. |
| `--max-candidates N` | no | Candidate accumulator budget; `0` means unlimited. |
| `--hot-seed-threshold N` | no | Skip or doc-only seeds above this collection-frequency threshold. |
| `--hot-seed-class hot\|very-hot` | no | Apply hot-seed policy using build-time frequency classes. |
| `--hot-mode skip\|doc-only` | no | Hot-seed policy, default `skip`. |

Budget and hot-seed options keep the stage labeled as `exact_retrieval`, but they change the retrieval policy. Skipped hot seeds and budgeted skips are reported explicitly in the metrics.

See [docs/cli_reference.md](docs/cli_reference.md) for all current k-mer options.

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
