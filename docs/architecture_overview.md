# AIMS Architecture Overview

AIMS separates seed extraction, routing, exact lookup, posting retrieval, candidate accumulation, and optional downstream alignment handoff.

## Modules

`seeds`
: DNA encoding, reverse-complement handling, and canonical k-mer generation.

`router`
: Probabilistic front-end interfaces. Router outputs are false-positive-prone and must be reported separately from exact lookups.

`index`
: Exact seed dictionaries, seed metadata, frequency classes, and posting stores.

`codecs`
: Compression primitives for compact posting storage. The implemented codecs are monotone delta-varint for integer lists and compressed posting blocks for exact coordinate postings.

`query`
: Planner, global multi-k seed ordering, and candidate accumulator. The accumulator is hash-based and supports candidate budgets.

`serialization`
: Stable binary index header and validation utilities. The complete mmap format will include directories, offsets, checksums, and metadata blocks.

`instrumentation`
: TSV/JSONL metric writers used by CLIs and benchmarks.

## Stage Labels

Benchmark and CLI output must label the comparison stage:

1. `exact_retrieval`
2. `probabilistic_routing`
3. `routing_only`
4. `pseudoalignment`
5. `end_to_end`

Results from different stages are not interchangeable.

## Query Planning And Candidate Orientation

The planner resolves each generated canonical k-mer to a dictionary seed ID once and carries that
ID into execution. Repeated occurrences of the same seed are deduplicated into one posting lookup,
while their information weight and forward/reverse canonical-orientation counts are retained.

A posting stores the reference occurrence's orientation relative to the canonical k-mer. Candidate
orientation is relative to the query and is computed as:

```text
candidate_reverse = query_canonical_reverse XOR reference_canonical_reverse
```

Consequently, an exact same-strand match accumulates under `forward`, and a reverse-complement
match accumulates under `reverse`. Doc-only hot-seed processing applies the same transformation
before deduplicating document/sequence/strand candidates.

## Ownership And Concurrency

Compressed posting blocks may be owned vectors or references into a mapped file. `PostingStore`
exposes the same zero-copy encoded-block view for both forms. Copies rebuild owned block pointers,
retain shared ownership for mapped blocks, and start with an empty decoded cache. Query workers
share an immutable index; only the bounded decoded-block LRU is synchronized.

The query CLI starts at most `--threads` long-lived asynchronous workers and writes completed
results in original input order. Posting budgets use logical posting counts, so cache warmness does
not alter which seeds are admitted.
