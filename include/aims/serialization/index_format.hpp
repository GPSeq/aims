#pragma once

#include <cstdint>
#include <filesystem>
#include <string>

#include "aims/build/index_builder.hpp"

namespace aims::serialization {

constexpr std::uint32_t index_magic = 0x534d4941U; // AIMS, little-endian byte order.
constexpr std::uint16_t index_format_version = 1;
constexpr std::uint64_t checksum_seed = 1469598103934665603ULL;

struct IndexHeader {
  std::uint32_t magic{index_magic};
  std::uint16_t version{index_format_version};
  std::uint16_t endian_marker{0x0102};
  std::uint32_t manifest_crc32{0};
  std::uint64_t metadata_offset{0};
  std::uint64_t router_directory_offset{0};
  std::uint64_t dictionary_directory_offset{0};
  std::uint64_t posting_directory_offset{0};
  std::uint64_t checksum_footer_offset{0};
};

[[nodiscard]] bool is_little_endian() noexcept;
[[nodiscard]] std::string describe_header(const IndexHeader& header);
void validate_header_or_throw(const IndexHeader& header);

void write_fixed_k_exact_index(const std::filesystem::path& path,
                               const build::FixedKIndex& index);
[[nodiscard]] build::FixedKIndex read_fixed_k_exact_index(const std::filesystem::path& path);

/**
 * @fn write_kmer_exact_index
 * @brief Serialize a KmerExactIndex to a checksum-protected binary index file.
 * @signature void write_kmer_exact_index(const std::filesystem::path& path, const build::KmerExactIndex& index);
 * @param path Output index path.
 * @param index K-mer exact index to serialize.
 * @throws std::runtime_error if the index is invalid or the file cannot be written.
 * @return None.
 */
void write_kmer_exact_index(const std::filesystem::path& path,
                            const build::KmerExactIndex& index);

/**
 * @fn read_kmer_exact_index
 * @brief Read a KmerExactIndex by copying compressed posting blocks into owned memory.
 * @signature build::KmerExactIndex read_kmer_exact_index(const std::filesystem::path& path);
 * @param path Input index path.
 * @throws std::runtime_error if the header, offsets, enum values, or checksum are invalid.
 * @return A KmerExactIndex with owned compressed posting blocks.
 */
[[nodiscard]] build::KmerExactIndex read_kmer_exact_index(const std::filesystem::path& path);

/**
 * @fn read_kmer_exact_index_mmap
 * @brief Read a KmerExactIndex with compressed posting blocks backed by an mmap file.
 * @signature build::KmerExactIndex read_kmer_exact_index_mmap(const std::filesystem::path& path, std::uint64_t decoded_block_cache_size);
 * @param path Input index path.
 * @param decoded_block_cache_size Maximum number of decoded posting blocks cached per posting store; zero disables caching.
 * @throws std::runtime_error if mapping, validation, or parsing fails.
 * @return A KmerExactIndex whose posting blocks reference the mapped file.
 */
[[nodiscard]] build::KmerExactIndex read_kmer_exact_index_mmap(
    const std::filesystem::path& path,
    std::uint64_t decoded_block_cache_size = 0);

} // namespace aims::serialization
