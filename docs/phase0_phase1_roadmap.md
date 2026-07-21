# Phase 0 and Phase 1 Roadmap

## Phase 0: Benchmark Harness

1. Define dataset manifest YAML/JSON with reference/query/truth checksums.
2. Define truth manifest schema with coordinate and stage semantics.
3. Implement `aims_bench --config bench.yaml`.
4. Add timing, CPU time, peak RSS, temp disk, and cache mode collection.
5. Add wrapper schemas for minimap2, Winnowmap, LexicMap, Raptor/HIBF-like routers, and SSHash-like exact dictionaries.
6. Require wrapper output to annotate `comparison_stage`.
7. Add deterministic benchmark result schema tests.
8. Add a small public or synthetic CI fixture with one-command execution.

## Phase 1: Fixed-k Exact Baseline

1. Expand checksum-validated binary serialization into the full directory-based mmap layout.
2. Broaden randomized dictionary/posting correctness tests against a naive map.
3. Add exact coordinate golden tests for reverse-complement duplicates.
4. Add edge tests for all-N inputs and highly repetitive references.
5. Formalize the query result JSON schema for top-k candidate output.
6. Add `aims_dump` structured output for dictionary and posting summaries.
7. Expand `aims_validate` checks for offset ranges, checksums, and metadata consistency.
8. Expand benchmark driver modes for build time, lookup-only, and candidate retrieval.
9. Add CI Debug ASan/UBSan and Release builds.
10. Document the exactness boundary: candidate retrieval only, no alignment.
