# AIMS CLI Reference

All current commands are k-mer exact-retrieval tools. They do not perform alignment.

## Build

```sh
aims_build --ref refs.fa --out index.aims --k 15,19,23,27,31
aims_build --ref refs.fa --out index.aims --k 15,19,23 --frequency-thresholds 2,16,128
```

Inputs may be FASTA or FASTQ. The index stores source URI, an FNV-1a file checksum, build command, sequence names, sequence lengths, compressed posting blocks, and per-k dictionaries.

`--frequency-thresholds rare,medium,hot` controls build-time collection-frequency classes. Seeds with `cf <= rare` are rare, `cf <= medium` are medium, `cf <= hot` are hot, and larger seeds are very-hot. These classes are stored in the index metadata.

## Query

```sh
aims_query --index index.aims --query queries.fa --topk 10 --emit jsonl --out results.jsonl
```

Important options:

```text
--threads N
--mmap
--posting-cache-blocks N
--max-seeds N
--max-postings N
--max-candidates N
--hot-seed-threshold N
--hot-seed-class hot|very-hot
--hot-mode skip|doc-only
```

`--threads` preserves input-order output. `--mmap` keeps compressed posting blocks file-backed. `--posting-cache-blocks` caches decoded posting blocks. `--hot-seed-threshold` applies a query-time collection-frequency cutoff. `--hot-seed-class` applies the build-time frequency class policy. `--hot-mode doc-only` uses hot seeds once per document/sequence/strand candidate instead of every coordinate.

## Benchmark

```sh
aims_bench --ref refs.fa --query queries.fa --truth truth.tsv --k 15,19,23 --topk 10
```

Useful benchmark options:

```text
--mmap
--posting-cache-blocks N
--repeats N
--query-metrics-out query_metrics.jsonl
--max-seeds N
--max-postings N
--max-candidates N
--hot-seed-threshold N
--hot-seed-class hot|very-hot
--hot-mode skip|doc-only
```

`aims_bench` emits one `BenchmarkResult` JSON object. When truth is provided, it reports top-1/top-5/top-10 recall. Load time, build time, query time, RSS, bytes read per query, and postings decoded per query are reported separately.

## Inspect And Validate

```sh
aims_validate --index index.aims --mmap
aims_dump --index index.aims --mmap
```

Both commands validate the same checksum-protected binary format. `aims_dump` prints source metadata, sequence metadata, and per-k layer summaries.
