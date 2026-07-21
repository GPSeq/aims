#include "aims/build/index_builder.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>

#include "aims/codecs/posting_block_codec.hpp"

namespace aims::build {
namespace {

// Convert raw collection frequency into the frequency class used by query policies.
FrequencyClass classify_frequency(std::uint64_t cf, const FrequencyThresholds& thresholds) {
  // Seeds occurring at most twice are highly selective and should usually be queried first.
  if (cf <= thresholds.rare_max) {
    return FrequencyClass::Rare;
  }
  // Medium seeds are still useful exact evidence but cost more postings.
  if (cf <= thresholds.medium_max) {
    return FrequencyClass::Medium;
  }
  // Hot seeds may still help recall, but they are expensive enough to make policy visible.
  if (cf <= thresholds.hot_max) {
    return FrequencyClass::Hot;
  }
  // Very hot seeds are usually repetitive sequence evidence and should be budgeted carefully.
  return FrequencyClass::VeryHot;
}

} // namespace

FixedKIndex build_fixed_k_exact(std::span<const io::FastaRecord> records, std::uint16_t k) {
  return build_fixed_k_exact(records, k, KmerBuildOptions{});
}

FixedKIndex build_fixed_k_exact(std::span<const io::FastaRecord> records,
                                std::uint16_t k,
                                const KmerBuildOptions& options) {
  // Group raw positional postings by canonical seed key before building compact arrays.
  std::map<SeedKey, std::vector<index::Posting>> postings_by_key;
  // KmerGenerator owns the DNA canonicalization and ambiguous-base splitting logic.
  seeds::KmerGenerator generator;
  // Phase 1 is k-mer only, canonicalized by default for exact DNA retrieval.
  seeds::SeedGenerationParams params{
      .family = SeedFamily::Kmer,
      .k = k,
      .canonical = true,
  };

  // Visit each reference sequence and emit its seed occurrences.
  for (const auto& record : records) {
    // Generate all valid canonical k-mers for the current reference sequence.
    const auto occurrences = generator.generate(SequenceView{
        .bases = record.sequence,
        .document_id = record.document_id,
        .sequence_id = record.sequence_id,
        .name = record.name,
    }, params);
    // Convert seed occurrences into coordinate-aware posting entries.
    for (const auto& occurrence : occurrences) {
      // The map groups identical canonical seeds into one posting list.
      postings_by_key[occurrence.key].push_back(index::Posting{
          .document_id = occurrence.document_id,
          .sequence_id = occurrence.sequence_id,
          .position = occurrence.position,
          .strand = occurrence.strand,
      });
    }
  }

  // Dictionary keys are stored separately from metadata for fast exact lookup.
  std::vector<SeedKey> keys;
  // Metadata records hold df/cf/cost estimates used by the planner.
  std::vector<index::SeedMetadata> metadata;
  // Plain posting vectors are retained here before the PostingStore compresses them.
  std::vector<std::vector<index::Posting>> postings;
  // The codec is used at build time to estimate compressed bytes per seed.
  codecs::PostingBlockCodec posting_codec;
  // Reserve final arrays because the number of unique seeds is already known.
  keys.reserve(postings_by_key.size());
  metadata.reserve(postings_by_key.size());
  postings.reserve(postings_by_key.size());

  // Convert the ordered map into deterministic dictionary and posting arrays.
  for (auto& [key, posting_list] : postings_by_key) {
    // Sorting is required for deterministic output and delta compression.
    std::sort(posting_list.begin(), posting_list.end());
    // Document frequency is the number of unique documents containing the seed.
    std::set<DocumentId> docs;
    // Collect unique document IDs from the posting list.
    for (const auto& posting : posting_list) {
      docs.insert(posting.document_id);
    }
    // Store this seed key in sorted order.
    keys.push_back(key);
    // Store planner-visible metadata for this seed.
    metadata.push_back(index::SeedMetadata{
        .document_frequency = docs.size(),
        .collection_frequency = posting_list.size(),
        .estimated_posting_bytes = posting_codec.encode(posting_list).size(),
        .frequency_class = classify_frequency(posting_list.size(), options.frequency_thresholds),
    });
    // Move the sorted posting list into the store input.
    postings.push_back(std::move(posting_list));
  }

  // Assemble the fixed-k index layer.
  FixedKIndex index;
  // Record the k value so query metrics and serialization can label this layer.
  index.k = k;
  // In this prototype each input record is treated as one document.
  index.document_count = records.size();
  // Build the exact dictionary with deterministic keys and metadata.
  index.dictionary.build(std::move(keys), std::move(metadata));
  // Build and compress the posting store.
  index.postings.build(std::move(postings));
  // Return a complete exact retrieval layer.
  return index;
}

KmerExactIndex build_kmer_exact(std::span<const io::FastaRecord> records,
                                std::span<const std::uint16_t> k_values) {
  return build_kmer_exact(records, k_values, KmerBuildOptions{});
}

KmerExactIndex build_kmer_exact(std::span<const io::FastaRecord> records,
                                std::span<const std::uint16_t> k_values,
                                const KmerBuildOptions& options) {
  // Multi-k indexes need at least one layer to be meaningful.
  if (k_values.empty()) {
    throw std::invalid_argument("at least one k value is required");
  }
  // The top-level index stores shared sequence metadata plus one layer per k.
  KmerExactIndex index;
  // Record document count once for planner IDF calculations.
  index.document_count = records.size();
  // Sequence metadata is copied from input records for output result names and validation.
  index.sequences.reserve(records.size());
  // Preserve every sequence ID/name/length in deterministic input order.
  for (const auto& record : records) {
    index.sequences.push_back(SequenceMetadata{
        .document_id = record.document_id,
        .sequence_id = record.sequence_id,
        .name = record.name,
        .length = record.sequence.size(),
    });
  }
  // Allocate all requested k layers.
  index.layers.reserve(k_values.size());
  // Build each fixed-k layer independently so ablation and per-k metrics remain clean.
  for (const auto k : k_values) {
    index.layers.push_back(build_fixed_k_exact(records, k, options));
  }
  // Return the complete multi-k exact index.
  return index;
}

} // namespace aims::build
