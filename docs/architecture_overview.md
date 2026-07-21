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
