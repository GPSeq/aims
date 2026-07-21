# Initial Exact Fixed-k Baseline Implementation Plan

This baseline is an exact retrieval stage, not a full aligner and not a pseudoaligner. Its purpose is to produce coordinate-aware candidate generation with first-class instrumentation before adding probabilistic routing, compression, and adaptive seed families.

## Scope

1. Parse DNA FASTA references and queries.
2. Split seed extraction at ambiguous bases.
3. Generate canonical k-mers for one configured `k`.
4. Build an exact static dictionary from distinct canonical seed keys.
5. Store plain sorted positional postings for each seed.
6. Query by planning seed order, fetching postings, and accumulating candidate document/sequence/strand scores.
7. Emit per-query metrics including exact bytes read and postings decoded.

## Non-goals

1. No downstream alignment in the core library.
2. No claims about end-to-end sensitivity or superiority.
3. No probabilistic router results mixed into exact retrieval metrics.
4. No compression-dependent conclusions until codecs are implemented and tested.

## Data Flow

`FASTA -> SeedGenerator -> SortedArrayDictionary + PostingStore -> QueryPlanner -> CandidateAccumulator -> metrics`

The first CLI implementation writes a checksum-validated `KmerExactIndex` binary containing seed metadata and plain positional postings for one or more k-mer layers. It is intentionally narrow, but it exercises the same stage boundary expected from later memory-mapped index files.

## Correctness Checks

1. Canonicalization property checks for forward and reverse-complement k-mers.
2. N-handling tests that confirm ambiguous bases split seedable runs.
3. Golden toy FASTA tests for exact seed positions and document frequencies.
4. Randomized lookup tests in Phase 1 to compare dictionary/postings against a naive map.
5. Deterministic metrics schema tests for TSV and JSONL output.

## Instrumentation Required Before Optimization

Every exact retrieval query must report:

1. `comparison_stage=exact_retrieval`
2. `seeds_generated`
3. `seeds_queried`
4. `exact_bytes_read`
5. `postings_decoded`
6. `candidates_created`
7. `candidates_touched`
8. `peak_accumulator_bytes`
9. `early_stop_fired`
10. `elapsed_ns`

Router fields remain present and zero in exact-only runs.
