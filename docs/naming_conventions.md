# Naming Conventions

Use names that describe the research stage and data contract precisely.

## Preferred Terms

1. `KmerExactIndex` for the current exact retrieval baseline index.
2. `comparison_stage` for benchmark stage labels.
3. `QueryMetrics` for per-query instrumentation.
4. `BenchmarkResult` for aggregated benchmark output.
5. `DatasetManifest` for input dataset declarations.
6. `TruthManifest` for ground-truth declarations.
7. `CandidateAccumulator` for seed-supported candidate scoring.

## Avoid

1. “best” or “better than” in code, docs, or output.
2. “prototype” in benchmark-facing output.
3. ambiguous labels such as “search result” when the stage is exact retrieval, routing-only, or pseudoalignment.
4. mixing heuristic and conservative early stopping under one label.
