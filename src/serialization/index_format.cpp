#include "aims/serialization/index_format.hpp"

#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <span>
#include <type_traits>
#include <vector>

#include "aims/io/mapped_file.hpp"

namespace aims::serialization {
namespace {

constexpr std::uint64_t checksum_prime = 1099511628211ULL;

void checksum_append(std::uint64_t& checksum, std::uint8_t byte) noexcept {
  checksum ^= byte;
  checksum *= checksum_prime;
}

std::uint64_t checksum_bytes(std::span<const std::uint8_t> bytes) noexcept {
  std::uint64_t checksum = checksum_seed;
  for (const auto byte : bytes) {
    checksum_append(checksum, byte);
  }
  return checksum;
}

template <typename UInt>
void append_le(std::vector<std::uint8_t>& bytes, UInt value) {
  static_assert(std::is_unsigned_v<UInt>);
  for (std::size_t i = 0; i < sizeof(UInt); ++i) {
    bytes.push_back(static_cast<std::uint8_t>((value >> (i * 8U)) & 0xffU));
  }
}

template <typename UInt>
UInt read_le(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  static_assert(std::is_unsigned_v<UInt>);
  if (offset + sizeof(UInt) > bytes.size()) {
    throw std::runtime_error("truncated AIMS index");
  }
  UInt value = 0;
  for (std::size_t i = 0; i < sizeof(UInt); ++i) {
    value |= static_cast<UInt>(bytes[offset + i]) << (i * 8U);
  }
  offset += sizeof(UInt);
  return value;
}

void append_header(std::vector<std::uint8_t>& bytes, const IndexHeader& header) {
  append_le<std::uint32_t>(bytes, header.magic);
  append_le<std::uint16_t>(bytes, header.version);
  append_le<std::uint16_t>(bytes, header.endian_marker);
  append_le<std::uint32_t>(bytes, header.manifest_crc32);
  append_le<std::uint64_t>(bytes, header.metadata_offset);
  append_le<std::uint64_t>(bytes, header.router_directory_offset);
  append_le<std::uint64_t>(bytes, header.dictionary_directory_offset);
  append_le<std::uint64_t>(bytes, header.posting_directory_offset);
  append_le<std::uint64_t>(bytes, header.checksum_footer_offset);
}

void append_string(std::vector<std::uint8_t>& bytes, const std::string& value) {
  append_le<std::uint64_t>(bytes, value.size());
  bytes.insert(bytes.end(), value.begin(), value.end());
}

IndexHeader read_header(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  IndexHeader header;
  header.magic = read_le<std::uint32_t>(bytes, offset);
  header.version = read_le<std::uint16_t>(bytes, offset);
  header.endian_marker = read_le<std::uint16_t>(bytes, offset);
  header.manifest_crc32 = read_le<std::uint32_t>(bytes, offset);
  header.metadata_offset = read_le<std::uint64_t>(bytes, offset);
  header.router_directory_offset = read_le<std::uint64_t>(bytes, offset);
  header.dictionary_directory_offset = read_le<std::uint64_t>(bytes, offset);
  header.posting_directory_offset = read_le<std::uint64_t>(bytes, offset);
  header.checksum_footer_offset = read_le<std::uint64_t>(bytes, offset);
  return header;
}

FrequencyClass frequency_class_from_byte(std::uint8_t value) {
  switch (value) {
  case 0:
    return FrequencyClass::Rare;
  case 1:
    return FrequencyClass::Medium;
  case 2:
    return FrequencyClass::Hot;
  case 3:
    return FrequencyClass::VeryHot;
  default:
    throw std::runtime_error("invalid frequency class in AIMS index");
  }
}

SeedFamily seed_family_from_byte(std::uint8_t value) {
  switch (value) {
  case 0:
    return SeedFamily::Kmer;
  case 1:
    return SeedFamily::Syncmer;
  case 2:
    return SeedFamily::Strobemer;
  default:
    throw std::runtime_error("invalid seed family in AIMS index");
  }
}

std::vector<std::uint8_t> read_file(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to open AIMS index: " + path.string());
  }
  return std::vector<std::uint8_t>(std::istreambuf_iterator<char>(in),
                                   std::istreambuf_iterator<char>());
}

std::string read_string(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  const auto size = read_le<std::uint64_t>(bytes, offset);
  if (size > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) ||
      offset + static_cast<std::size_t>(size) > bytes.size()) {
    throw std::runtime_error("truncated string in AIMS index");
  }
  std::string value(reinterpret_cast<const char*>(bytes.data() + offset),
                    static_cast<std::size_t>(size));
  offset += static_cast<std::size_t>(size);
  return value;
}

build::KmerExactIndex parse_kmer_exact_index(std::span<const std::uint8_t> bytes,
                                             std::shared_ptr<const void> block_owner,
                                             bool copy_posting_blocks,
                                             std::uint64_t decoded_block_cache_size) {
  if (bytes.size() < 60) {
    throw std::runtime_error("AIMS index is too small");
  }

  std::size_t offset = 0;
  const auto header = read_header(bytes, offset);
  validate_header_or_throw(header);

  if (header.metadata_offset != offset) {
    throw std::runtime_error("unexpected AIMS metadata offset");
  }
  if (header.checksum_footer_offset + sizeof(std::uint64_t) != bytes.size()) {
    throw std::runtime_error("invalid AIMS checksum footer offset");
  }

  std::size_t footer_offset = static_cast<std::size_t>(header.checksum_footer_offset);
  std::size_t checksum_offset = footer_offset;
  const auto stored_checksum = read_le<std::uint64_t>(bytes, checksum_offset);
  const auto computed_checksum =
      checksum_bytes(std::span<const std::uint8_t>(bytes.data(), footer_offset));
  if (stored_checksum != computed_checksum) {
    throw std::runtime_error("AIMS index checksum validation failed");
  }

  build::KmerExactIndex index;
  index.document_count = read_le<std::uint64_t>(bytes, offset);
  const auto layer_count = read_le<std::uint64_t>(bytes, offset);
  index.source_uri = read_string(bytes, offset);
  index.source_checksum = read_string(bytes, offset);
  index.build_command = read_string(bytes, offset);
  const auto sequence_count = read_le<std::uint64_t>(bytes, offset);
  if (layer_count == 0) {
    throw std::runtime_error("AIMS k-mer index contains no layers");
  }

  if (sequence_count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw std::runtime_error("AIMS sequence metadata is too large for this platform");
  }
  index.sequences.reserve(static_cast<std::size_t>(sequence_count));
  for (std::uint64_t i = 0; i < sequence_count; ++i) {
    build::SequenceMetadata sequence;
    sequence.document_id = read_le<std::uint32_t>(bytes, offset);
    sequence.sequence_id = read_le<std::uint64_t>(bytes, offset);
    sequence.length = read_le<std::uint64_t>(bytes, offset);
    sequence.name = read_string(bytes, offset);
    index.sequences.push_back(std::move(sequence));
  }

  if (header.dictionary_directory_offset != offset) {
    throw std::runtime_error("unexpected AIMS dictionary offset");
  }

  if (layer_count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw std::runtime_error("AIMS layer count is too large for this platform");
  }
  index.layers.reserve(static_cast<std::size_t>(layer_count));

  for (std::uint64_t layer_index = 0; layer_index < layer_count; ++layer_index) {
    build::FixedKIndex layer;
    layer.k = read_le<std::uint16_t>(bytes, offset);
    static_cast<void>(read_le<std::uint16_t>(bytes, offset));
    layer.document_count = read_le<std::uint64_t>(bytes, offset);
    const auto dictionary_count = read_le<std::uint64_t>(bytes, offset);
    if (dictionary_count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
      throw std::runtime_error("AIMS dictionary is too large for this platform");
    }

    std::vector<SeedKey> keys;
    std::vector<index::SeedMetadata> metadata;
    keys.reserve(static_cast<std::size_t>(dictionary_count));
    metadata.reserve(static_cast<std::size_t>(dictionary_count));
    for (std::uint64_t i = 0; i < dictionary_count; ++i) {
      auto key = SeedKey{};
      key.value = read_le<std::uint64_t>(bytes, offset);
      key.family = seed_family_from_byte(read_le<std::uint8_t>(bytes, offset));
      key.k = read_le<std::uint16_t>(bytes, offset);
      keys.push_back(key);

      auto item = index::SeedMetadata{};
      item.document_frequency = read_le<std::uint64_t>(bytes, offset);
      item.collection_frequency = read_le<std::uint64_t>(bytes, offset);
      item.estimated_posting_bytes = read_le<std::uint64_t>(bytes, offset);
      item.frequency_class = frequency_class_from_byte(read_le<std::uint8_t>(bytes, offset));
      metadata.push_back(item);
    }

    const auto posting_list_count = read_le<std::uint64_t>(bytes, offset);
    if (posting_list_count != dictionary_count) {
      throw std::runtime_error("AIMS posting list count mismatch");
    }

    if (copy_posting_blocks) {
      std::vector<std::vector<std::uint8_t>> encoded_blocks;
      std::vector<std::uint64_t> posting_counts;
      encoded_blocks.reserve(static_cast<std::size_t>(posting_list_count));
      posting_counts.reserve(static_cast<std::size_t>(posting_list_count));
      for (std::uint64_t i = 0; i < posting_list_count; ++i) {
        const auto posting_count = read_le<std::uint64_t>(bytes, offset);
        const auto encoded_size = read_le<std::uint64_t>(bytes, offset);
        if (encoded_size > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) ||
            offset + static_cast<std::size_t>(encoded_size) > bytes.size()) {
          throw std::runtime_error("AIMS encoded posting block is too large or truncated");
        }
        std::vector<std::uint8_t> block(
            bytes.begin() + static_cast<std::ptrdiff_t>(offset),
            bytes.begin() + static_cast<std::ptrdiff_t>(offset + encoded_size));
        offset += static_cast<std::size_t>(encoded_size);
        posting_counts.push_back(posting_count);
        encoded_blocks.push_back(std::move(block));
      }
      layer.dictionary.build(std::move(keys), std::move(metadata));
      layer.postings.build_encoded(std::move(encoded_blocks), std::move(posting_counts));
    } else {
      std::vector<std::shared_ptr<const void>> owners;
      std::vector<const std::uint8_t*> block_data;
      std::vector<std::uint64_t> encoded_sizes;
      std::vector<std::uint64_t> posting_counts;
      owners.reserve(static_cast<std::size_t>(posting_list_count));
      block_data.reserve(static_cast<std::size_t>(posting_list_count));
      encoded_sizes.reserve(static_cast<std::size_t>(posting_list_count));
      posting_counts.reserve(static_cast<std::size_t>(posting_list_count));
      for (std::uint64_t i = 0; i < posting_list_count; ++i) {
        const auto posting_count = read_le<std::uint64_t>(bytes, offset);
        const auto encoded_size = read_le<std::uint64_t>(bytes, offset);
        if (encoded_size > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()) ||
            offset + static_cast<std::size_t>(encoded_size) > bytes.size()) {
          throw std::runtime_error("AIMS encoded posting block is too large or truncated");
        }
        owners.push_back(block_owner);
        block_data.push_back(bytes.data() + offset);
        encoded_sizes.push_back(encoded_size);
        posting_counts.push_back(posting_count);
        offset += static_cast<std::size_t>(encoded_size);
      }
      layer.dictionary.build(std::move(keys), std::move(metadata));
      layer.postings.build_encoded_refs(std::move(owners), std::move(block_data),
                                        std::move(encoded_sizes), std::move(posting_counts));
    }
    layer.postings.set_cache_limit(decoded_block_cache_size);
    index.layers.push_back(std::move(layer));
  }
  if (offset != footer_offset) {
    throw std::runtime_error("unexpected trailing payload before AIMS checksum footer");
  }
  return index;
}

} // namespace

bool is_little_endian() noexcept {
  const std::uint16_t value = 0x0001;
  return *reinterpret_cast<const std::uint8_t*>(&value) == 0x01;
}

std::string describe_header(const IndexHeader& header) {
  std::ostringstream out;
  out << "magic=0x" << std::hex << header.magic << std::dec
      << " version=" << header.version
      << " endian_marker=0x" << std::hex << header.endian_marker << std::dec
      << " manifest_crc32=" << header.manifest_crc32;
  return out.str();
}

void validate_header_or_throw(const IndexHeader& header) {
  if (header.magic != index_magic) {
    throw std::runtime_error("invalid AIMS index magic");
  }
  if (header.version != index_format_version) {
    throw std::runtime_error("unsupported AIMS index version");
  }
  if (header.endian_marker != 0x0102) {
    throw std::runtime_error("unexpected AIMS endian marker");
  }
}

void write_fixed_k_exact_index(const std::filesystem::path& path,
                               const build::FixedKIndex& index) {
  build::KmerExactIndex multi_index;
  multi_index.document_count = index.document_count;
  multi_index.layers.push_back(index);
  write_kmer_exact_index(path, multi_index);
}

build::FixedKIndex read_fixed_k_exact_index(const std::filesystem::path& path) {
  auto index = read_kmer_exact_index(path);
  if (index.layers.empty()) {
    throw std::runtime_error("AIMS k-mer index contains no layers");
  }
  return std::move(index.layers.front());
}

void write_kmer_exact_index(const std::filesystem::path& path,
                            const build::KmerExactIndex& index) {
  if (!is_little_endian()) {
    throw std::runtime_error("AIMS index writer currently requires little-endian host order");
  }
  if (index.layers.empty()) {
    throw std::runtime_error("KmerExactIndex requires at least one k-mer layer");
  }
  for (const auto& layer : index.layers) {
    if (layer.dictionary.size() != layer.postings.size()) {
      throw std::runtime_error("dictionary and posting store size mismatch");
    }
  }

  IndexHeader header;
  header.metadata_offset = 52;
  header.router_directory_offset = 0;

  std::vector<std::uint8_t> payload;
  append_le<std::uint64_t>(payload, index.document_count);
  append_le<std::uint64_t>(payload, index.layers.size());
  append_string(payload, index.source_uri);
  append_string(payload, index.source_checksum);
  append_string(payload, index.build_command);
  append_le<std::uint64_t>(payload, index.sequences.size());
  for (const auto& sequence : index.sequences) {
    append_le<std::uint32_t>(payload, sequence.document_id);
    append_le<std::uint64_t>(payload, sequence.sequence_id);
    append_le<std::uint64_t>(payload, sequence.length);
    append_string(payload, sequence.name);
  }

  header.dictionary_directory_offset = header.metadata_offset + payload.size();

  bool posting_offset_set = false;
  for (const auto& layer : index.layers) {
    append_le<std::uint16_t>(payload, layer.k);
    append_le<std::uint16_t>(payload, 0);
    append_le<std::uint64_t>(payload, layer.document_count);
    append_le<std::uint64_t>(payload, layer.dictionary.size());

    const auto keys = layer.dictionary.keys();
    for (std::size_t i = 0; i < keys.size(); ++i) {
      const auto& key = keys[i];
      const auto& metadata = layer.dictionary.metadata(static_cast<SeedId>(i));
      append_le<std::uint64_t>(payload, key.value);
      append_le<std::uint8_t>(payload, static_cast<std::uint8_t>(key.family));
      append_le<std::uint16_t>(payload, key.k);
      append_le<std::uint64_t>(payload, metadata.document_frequency);
      append_le<std::uint64_t>(payload, metadata.collection_frequency);
      append_le<std::uint64_t>(payload, metadata.estimated_posting_bytes);
      append_le<std::uint8_t>(payload, static_cast<std::uint8_t>(metadata.frequency_class));
    }

    if (!posting_offset_set) {
      header.posting_directory_offset = header.metadata_offset + payload.size();
      posting_offset_set = true;
    }
    append_le<std::uint64_t>(payload, layer.postings.size());
    const auto encoded_blocks = layer.postings.encoded_blocks();
    const auto posting_counts = layer.postings.posting_counts();
    for (std::size_t i = 0; i < encoded_blocks.size(); ++i) {
      append_le<std::uint64_t>(payload, posting_counts[i]);
      append_le<std::uint64_t>(payload, encoded_blocks[i].size());
      payload.insert(payload.end(), encoded_blocks[i].begin(), encoded_blocks[i].end());
    }
  }

  header.checksum_footer_offset = header.metadata_offset + payload.size();

  std::vector<std::uint8_t> bytes;
  append_header(bytes, header);
  bytes.insert(bytes.end(), payload.begin(), payload.end());
  append_le<std::uint64_t>(bytes, checksum_bytes(bytes));

  std::ofstream out(path, std::ios::binary);
  if (!out) {
    throw std::runtime_error("failed to write AIMS index: " + path.string());
  }
  out.write(reinterpret_cast<const char*>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
}

build::KmerExactIndex read_kmer_exact_index(const std::filesystem::path& path) {
  const auto bytes = read_file(path);
  return parse_kmer_exact_index(bytes, {}, true, 0);
}

build::KmerExactIndex read_kmer_exact_index_mmap(const std::filesystem::path& path,
                                                 std::uint64_t decoded_block_cache_size) {
  auto mapped = std::make_shared<io::MappedFile>(path);
  std::shared_ptr<const void> owner = mapped;
  return parse_kmer_exact_index(mapped->bytes(), std::move(owner), false,
                                decoded_block_cache_size);
}

} // namespace aims::serialization
