#pragma once

#include <cstdint>

#include "aims/build/index_builder.hpp"
#include "aims/instrumentation/metrics.hpp"
#include "aims/io/fasta.hpp"

namespace aims::query {

enum class HotSeedMode {
  Skip,
  DocOnly,
};

struct KmerSearchOptions {
  std::uint32_t topk{10};
  std::uint64_t max_seeds{0};
  std::uint64_t max_postings{0};
  std::uint64_t max_candidates{0};
  std::uint64_t hot_seed_threshold{0};
  bool use_frequency_class_hot_policy{false};
  FrequencyClass hot_seed_min_class{FrequencyClass::VeryHot};
  HotSeedMode hot_seed_mode{HotSeedMode::Skip};
};

/**
 * @fn search_kmer_exact
 * @brief Execute exact k-mer candidate retrieval for one query sequence.
 * @signature instrumentation::QueryMetrics search_kmer_exact(const build::KmerExactIndex& index, const io::FastaRecord& query, KmerSearchOptions options);
 * @param index K-mer exact index, optionally backed by mmap posting blocks.
 * @param query Query record loaded from FASTA or FASTQ input.
 * @param options Top-k, query-budget, candidate-budget, and hot-kmer policy settings.
 * @throws std::runtime_error if a compressed posting block is malformed.
 * @return Per-query metrics containing per-k accounting and structured top-k candidates.
 */
[[nodiscard]] instrumentation::QueryMetrics search_kmer_exact(
    const build::KmerExactIndex& index,
    const io::FastaRecord& query,
    KmerSearchOptions options);

} // namespace aims::query
