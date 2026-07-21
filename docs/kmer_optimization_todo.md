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
