# AIMS File Format Sketch

The current CLI uses a checksum-validated binary `KmerExactIndex` for Phase 1 exact retrieval. The target production index format remains binary, memory-mappable, and organized into explicit directories for router, dictionary, posting, and hot-seed payloads.

## Target Binary Layout

```text
[magic][version][endianness][manifest_crc]
[global metadata]
[router directory]
[router payloads]
[dictionary directory]
[dictionary payloads]
[posting directory]
[posting payloads]
[hot-seed tables]
[checksums footer]
```

## Header Fields

1. `magic`
2. `version`
3. `endian_marker`
4. `manifest_crc32`
5. `metadata_offset`
6. `router_directory_offset`
7. `dictionary_directory_offset`
8. `posting_directory_offset`
9. `checksum_footer_offset`

Readers must reject unknown magic values, unsupported versions, inconsistent offsets, invalid enum values, and invalid checksum footers.

## Current KmerExactIndex Payload

The current baseline writer stores:

1. global metadata: document count and number of k-mer layers
2. reference sequence metadata: document ID, sequence ID, name, length
3. per-layer metadata: `k`, document count, distinct seed count
4. dictionary records: seed key, seed family, k, df, cf, estimated posting bytes, frequency class
5. compressed posting blocks: posting count, encoded byte count, delta-varint-like row encoding for document ID, sequence ID, position, strand
6. checksum footer over the header and payload

This format is exact retrieval only. It does not contain router summaries or alignment output.

## Loading Modes

The standard reader copies compressed posting blocks into owned memory. The mmap reader keeps posting blocks as references into a mapped index file and decodes blocks lazily during query execution. Both readers validate the same checksum footer before exposing the index.
