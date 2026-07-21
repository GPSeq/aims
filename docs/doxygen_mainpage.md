# AIMS API Documentation

AIMS is a C++20 k-mer exact-retrieval research system for biological sequence collections.

The current public API centers on:

- `aims::build::build_kmer_exact` for multi-k index construction.
- `aims::serialization::write_kmer_exact_index` and `aims::serialization::read_kmer_exact_index_mmap` for stable index persistence.
- `aims::query::search_kmer_exact` for exact candidate retrieval with compressed posting blocks, mmap loading, query budgets, and hot-kmer policies.
- `aims::instrumentation::QueryMetrics` for stage-labeled per-query metrics.

The system is not an aligner. It reports exact seed-supported candidate retrieval metrics.
