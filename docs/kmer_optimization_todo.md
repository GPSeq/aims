# AIMS k-mer Optimization TODO

This list tracks optimization work for the k-mer exact-retrieval tool. It intentionally excludes full alignment, learned indexes, syncmers, and strobemers.

## Step 1: Track the k-mer-only optimization roadmap

Status: done in this document.

Goal: keep optimization work scoped to exact k-mer retrieval, compressed posting access, candidate generation, and reproducible metrics.

## Step 2: Build-time frequency thresholds

Status: done.

Goal: make rare, medium, hot, and very-hot seed classes configurable instead of hard-coded. These classes should be stored in seed metadata and used by query policies.

## Step 3: Class-aware hot-seed policy

Status: done.

Goal: allow query runs to skip or doc-deduplicate hot and very-hot seeds based on build-time frequency classes, while reporting the policy in metrics.

## Step 4: Benchmark load-time and RSS metrics

Status: done.

Goal: separate index build time, index load time, query time, peak RSS, serialized size, bytes read per query, and postings decoded per query.

## Step 5: Scientific labeling in docs and CLI help

Status: done.

Goal: make every optimized mode clearly labeled as exact retrieval with an explicit policy. Budgeted or skipped-seed modes must not be described as full unrestricted retrieval.

## Step 6: Tests

Status: done.

Goal: add focused tests for configurable frequency classes, class-aware hot-seed behavior, deterministic metrics, and benchmark schema stability.

## Step 7: Verification

Status: done.

Goal: rebuild, run the complete test suite, and regenerate the algorithm PDF if documentation changed.

## Step 8: Correct query-relative strand accumulation

Status: done.

Goal: retain canonical orientation counts while deduplicating query seeds and combine query and
reference orientation with XOR so same-strand and reverse-complement evidence is not mislabeled or
split by reference canonicalization alone.

## Step 9: Remove avoidable query overhead

Status: done.

Goal: resolve seed IDs once in the planner, reserve planner output, use partial sorting for top-k,
and reuse a fixed number of query workers instead of creating one asynchronous task per query.

## Step 10: Backing-store-safe postings

Status: done.

Goal: make posting-store copies rebind internal references, make mmap and owned compressed blocks
equally serializable, and keep posting budgets deterministic across cache states.

## Step 11: Compact dictionary lookup table

Status: done.

Goal: replace the node-based seed-ID hash map with a contiguous open-addressing table while
retaining deterministic sorted dictionary keys.

## Step 12: Lower-memory index finalization

Status: done.

Goal: release hash-map nodes and raw posting vectors progressively while sorted posting lists are
encoded. This lowers finalization peak memory; disk-backed partitioning remains a later option for
collections that cannot fit their raw postings in memory at all.

## Step 13: Atomic streaming index writer

Status: done.

Goal: write index sections directly to a same-directory temporary file, backpatch the final header,
calculate the checksum without constructing a second full index image in memory, and atomically
rename the completed file.

## Step 14: Bounded query pipeline

Status: done.

Goal: read and execute query records in configurable chunks while preserving deterministic input
order and using a bounded worker set.

## Step 15: Byte-bounded decoded posting cache

Status: done.

Goal: support cache limits in decoded bytes as well as block count and report cache methodology
accurately in benchmarks.

## Step 16: Position-consistent candidate evidence

Status: done.

Goal: retain query positions and group evidence by compatible forward/reverse coordinate diagonals
so repetitive seed Cartesian products do not dominate candidate scores. Keep this as exact
retrieval rather than claiming downstream alignment.

## Step 17: Stress and concurrency verification

Status: done.

Goal: add malformed-input properties, large generated fixtures, threaded deterministic-output
checks, and optional thread-sanitizer coverage where the toolchain supports it.
