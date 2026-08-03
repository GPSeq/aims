#include "aims/build/index_builder.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>

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

using PostingMap = std::unordered_map<SeedKey, std::vector<index::Posting>>;

void validate_thread_count(const KmerBuildOptions& options) {
  if (options.thread_count == 0) {
    throw std::invalid_argument("thread_count must be greater than zero");
  }
  if (options.thread_count > static_cast<std::uint32_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("thread_count is too large for OpenMP");
  }
#ifndef AIMS_HAS_OPENMP
  if (options.thread_count > 1) {
    throw std::invalid_argument("--threads requires configuring with AIMS_ENABLE_OPENMP=ON");
  }
#endif
}

std::int64_t checked_parallel_layer_count(std::size_t size) {
  if (size > static_cast<std::size_t>(std::numeric_limits<std::int64_t>::max())) {
    throw std::runtime_error("too many index layers");
  }
  return static_cast<std::int64_t>(size);
}

void collect_fixed_k_postings(PostingMap& postings_by_key,
                              std::span<const io::FastaRecord> records,
                              std::uint16_t k) {
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
}

std::uint64_t document_frequency_for_sorted_postings(
    std::span<const index::Posting> posting_list) {
  std::uint64_t document_frequency = 0;
  DocumentId previous_document = 0;
  bool has_previous = false;
  for (const auto& posting : posting_list) {
    if (!has_previous || posting.document_id != previous_document) {
      ++document_frequency;
      previous_document = posting.document_id;
      has_previous = true;
    }
  }
  return document_frequency;
}

FixedKIndex build_fixed_k_from_postings(PostingMap postings_by_key,
                                        std::uint16_t k,
                                        std::uint64_t document_count,
                                        const KmerBuildOptions& options) {
  // Dictionary keys are stored separately from metadata for fast exact lookup.
  std::vector<SeedKey> keys;
  // Metadata records hold df/cf/cost estimates used by the planner.
  std::vector<index::SeedMetadata> metadata;
  // Posting blocks are encoded directly as each seed list is finalized.
  std::vector<std::vector<std::uint8_t>> encoded_blocks;
  std::vector<std::uint64_t> posting_counts;
  // The codec is used at build time to estimate compressed bytes per seed.
  codecs::PostingBlockCodec posting_codec;
  // Reserve final arrays because the number of unique seeds is already known.
  keys.reserve(postings_by_key.size());
  metadata.reserve(postings_by_key.size());
  encoded_blocks.reserve(postings_by_key.size());
  posting_counts.reserve(postings_by_key.size());

  std::vector<std::pair<SeedKey, std::vector<index::Posting>*>> ordered_postings;
  ordered_postings.reserve(postings_by_key.size());
  for (auto& [key, posting_list] : postings_by_key) {
    ordered_postings.push_back({key, &posting_list});
  }
  std::sort(ordered_postings.begin(), ordered_postings.end(),
            [](const auto& lhs, const auto& rhs) {
              return lhs.first < rhs.first;
            });

  // Convert the hash accumulator into deterministic dictionary and posting arrays.
  for (auto& [key, posting_list_ptr] : ordered_postings) {
    auto& posting_list = *posting_list_ptr;
    // Sorting is required for deterministic output and delta compression.
    if (!std::is_sorted(posting_list.begin(), posting_list.end())) {
      std::sort(posting_list.begin(), posting_list.end());
    }
    const auto document_frequency = document_frequency_for_sorted_postings(posting_list);
    const auto collection_frequency = posting_list.size();
    auto encoded = posting_codec.encode(posting_list);
    const auto encoded_size = encoded.size();
    // Store this seed key in sorted order.
    keys.push_back(key);
    // Store planner-visible metadata for this seed.
    metadata.push_back(index::SeedMetadata{
        .document_frequency = document_frequency,
        .collection_frequency = collection_frequency,
        .estimated_posting_bytes = encoded_size,
        .frequency_class = classify_frequency(collection_frequency, options.frequency_thresholds),
    });
    posting_counts.push_back(collection_frequency);
    encoded_blocks.push_back(std::move(encoded));
  }

  // Assemble the fixed-k index layer.
  FixedKIndex index;
  // Record the k value so query metrics and serialization can label this layer.
  index.k = k;
  // In this prototype each input record is treated as one document.
  index.document_count = document_count;
  // Build the exact dictionary with deterministic keys and metadata.
  index.dictionary.build(std::move(keys), std::move(metadata));
  // Build the posting store from already compressed blocks.
  index.postings.build_encoded(std::move(encoded_blocks), std::move(posting_counts));
  // Return a complete exact retrieval layer.
  return index;
}

} // namespace

FixedKIndex build_fixed_k_exact(std::span<const io::FastaRecord> records, std::uint16_t k) {
  return build_fixed_k_exact(records, k, KmerBuildOptions{});
}

FixedKIndex build_fixed_k_exact(std::span<const io::FastaRecord> records,
                                std::uint16_t k,
                                const KmerBuildOptions& options) {
  validate_thread_count(options);
  // Group raw positional postings by canonical seed key before building compact arrays.
  PostingMap postings_by_key;
  collect_fixed_k_postings(postings_by_key, records, k);
  return build_fixed_k_from_postings(std::move(postings_by_key), k, records.size(), options);
}

KmerExactIndexBuilder::KmerExactIndexBuilder(std::span<const std::uint16_t> k_values,
                                             const KmerBuildOptions& options)
    : options_(options) {
  validate_thread_count(options_);
  // Multi-k indexes need at least one layer to be meaningful.
  if (k_values.empty()) {
    throw std::invalid_argument("at least one k value is required");
  }
  layers_.reserve(k_values.size());
  for (const auto k : k_values) {
    layers_.push_back(LayerAccumulator{.k = k});
  }
}

void KmerExactIndexBuilder::add_records(std::span<const io::FastaRecord> records) {
  if (records.empty()) {
    return;
  }
  sequences_.reserve(sequences_.size() + records.size());
  for (const auto& record : records) {
    sequences_.push_back(SequenceMetadata{
        .document_id = record.document_id,
        .sequence_id = record.sequence_id,
        .name = record.name,
        .length = record.sequence.size(),
    });
  }
  document_count_ += records.size();
  const auto layer_count = checked_parallel_layer_count(layers_.size());
#ifdef AIMS_HAS_OPENMP
#pragma omp parallel for schedule(static) num_threads(static_cast<int>(options_.thread_count)) \
    if (options_.thread_count > 1 && layers_.size() > 1)
#endif
  for (std::int64_t i = 0; i < layer_count; ++i) {
    auto& layer = layers_[static_cast<std::size_t>(i)];
    collect_fixed_k_postings(layer.postings_by_key, records, layer.k);
  }
}

KmerExactIndex KmerExactIndexBuilder::finish() {
  // The top-level index stores shared sequence metadata plus one layer per k.
  KmerExactIndex index;
  // Record document count once for planner IDF calculations.
  index.document_count = document_count_;
  index.sequences = std::move(sequences_);
  // Allocate all requested k layers.
  index.layers.resize(layers_.size());
  // Build each fixed-k layer independently so ablation and per-k metrics remain clean.
  const auto layer_count = checked_parallel_layer_count(layers_.size());
#ifdef AIMS_HAS_OPENMP
#pragma omp parallel for schedule(static) num_threads(static_cast<int>(options_.thread_count)) \
    if (options_.thread_count > 1 && layers_.size() > 1)
#endif
  for (std::int64_t i = 0; i < layer_count; ++i) {
    auto& layer = layers_[static_cast<std::size_t>(i)];
    index.layers[static_cast<std::size_t>(i)] =
        build_fixed_k_from_postings(std::move(layer.postings_by_key), layer.k, document_count_,
                                    options_);
  }
  // Return the complete multi-k exact index.
  return index;
}

KmerExactIndex build_kmer_exact(std::span<const io::FastaRecord> records,
                                std::span<const std::uint16_t> k_values) {
  return build_kmer_exact(records, k_values, KmerBuildOptions{});
}

KmerExactIndex build_kmer_exact(std::span<const io::FastaRecord> records,
                                std::span<const std::uint16_t> k_values,
                                const KmerBuildOptions& options) {
  KmerExactIndexBuilder builder(k_values, options);
  builder.add_records(records);
  return builder.finish();
}

KmerExactIndex build_kmer_exact_from_sequence_file(const std::filesystem::path& path,
                                                  std::span<const std::uint16_t> k_values,
                                                  const KmerBuildOptions& options,
                                                  std::size_t max_records_per_chunk) {
  KmerExactIndexBuilder builder(k_values, options);
  io::read_sequence_chunks(path, max_records_per_chunk,
                           [&](std::span<const io::FastaRecord> records) {
                             builder.add_records(records);
                           });
  return builder.finish();
}

} // namespace aims::build
