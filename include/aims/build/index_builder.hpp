#pragma once

#include <cstdint>
#include <filesystem>
#include <span>
#include <string>
#include <unordered_map>
#include <vector>

#include "aims/index/exact_dictionary.hpp"
#include "aims/index/postings.hpp"
#include "aims/io/fasta.hpp"
#include "aims/seeds/kmer_generator.hpp"

namespace aims::build {

struct FrequencyThresholds {
  std::uint64_t rare_max{2};
  std::uint64_t medium_max{16};
  std::uint64_t hot_max{128};
};

struct KmerBuildOptions {
  FrequencyThresholds frequency_thresholds{};
  std::uint32_t thread_count{1};
};

struct FixedKIndex {
  std::uint16_t k{0};
  std::uint64_t document_count{0};
  index::SortedArrayDictionary dictionary;
  index::PostingStore postings;
};

struct SequenceMetadata {
  DocumentId document_id{0};
  SequenceId sequence_id{0};
  std::string name;
  std::uint64_t length{0};
};

struct KmerExactIndex {
  std::uint64_t document_count{0};
  std::string source_uri;
  std::string source_checksum;
  std::string build_command;
  std::vector<SequenceMetadata> sequences;
  std::vector<FixedKIndex> layers;
};

class KmerExactIndexBuilder {
public:
  KmerExactIndexBuilder(std::span<const std::uint16_t> k_values,
                        const KmerBuildOptions& options = {});

  void add_records(std::span<const io::FastaRecord> records);
  [[nodiscard]] KmerExactIndex finish();

private:
  struct LayerAccumulator {
    std::uint16_t k{0};
    std::unordered_map<SeedKey, std::vector<index::Posting>> postings_by_key;
  };

  KmerBuildOptions options_{};
  std::uint64_t document_count_{0};
  std::vector<SequenceMetadata> sequences_{};
  std::vector<LayerAccumulator> layers_{};
};

/**
 * @fn build_fixed_k_exact
 * @brief Build one exact canonical k-mer index layer for a single k value.
 * @signature FixedKIndex build_fixed_k_exact(std::span<const io::FastaRecord> records, std::uint16_t k, const KmerBuildOptions& options = {});
 * @param records Reference records loaded from FASTA or FASTQ input.
 * @param k K-mer length. The current DNA encoder supports 1 <= k <= 32.
 * @param options Build-time frequency classification thresholds.
 * @throws std::invalid_argument if the generated dictionary/posting inputs are inconsistent.
 * @return A fixed-k exact index layer containing dictionary metadata and compressed posting blocks.
 */
[[nodiscard]] FixedKIndex build_fixed_k_exact(
    std::span<const io::FastaRecord> records,
    std::uint16_t k,
    const KmerBuildOptions& options = {});

/**
 * @fn build_kmer_exact
 * @brief Build a multi-k exact canonical k-mer index.
 * @signature KmerExactIndex build_kmer_exact(std::span<const io::FastaRecord> records, std::span<const std::uint16_t> k_values, const KmerBuildOptions& options = {});
 * @param records Reference records loaded from FASTA or FASTQ input.
 * @param k_values K-mer lengths to index as independent exact retrieval layers.
 * @param options Build-time frequency classification thresholds.
 * @throws std::invalid_argument if k_values is empty or a layer cannot be built.
 * @return A k-mer exact index containing reference metadata and one layer per k value.
 */
[[nodiscard]] KmerExactIndex build_kmer_exact(
    std::span<const io::FastaRecord> records,
    std::span<const std::uint16_t> k_values,
    const KmerBuildOptions& options = {});

[[nodiscard]] KmerExactIndex build_kmer_exact_from_sequence_file(
    const std::filesystem::path& path,
    std::span<const std::uint16_t> k_values,
    const KmerBuildOptions& options = {},
    std::size_t max_records_per_chunk = 1024);

} // namespace aims::build
