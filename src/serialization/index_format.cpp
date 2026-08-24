#include "aims/serialization/index_format.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
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
void write_le(std::ostream& out, UInt value) {
  static_assert(std::is_unsigned_v<UInt>);
  std::array<std::uint8_t, sizeof(UInt)> bytes{};
  for (std::size_t i = 0; i < sizeof(UInt); ++i) {
    bytes[i] = static_cast<std::uint8_t>((value >> (i * 8U)) & 0xffU);
  }
  out.write(reinterpret_cast<const char*>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
  if (!out) {
    throw std::runtime_error("failed while writing AIMS index");
  }
}

template <typename UInt>
UInt read_le(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  static_assert(std::is_unsigned_v<UInt>);
  if (offset > bytes.size() || sizeof(UInt) > bytes.size() - offset) {
    throw std::runtime_error("truncated AIMS index");
  }
  UInt value = 0;
  for (std::size_t i = 0; i < sizeof(UInt); ++i) {
    value |= static_cast<UInt>(bytes[offset + i]) << (i * 8U);
  }
  offset += sizeof(UInt);
  return value;
}

void write_header(std::ostream& out, const IndexHeader& header) {
  write_le<std::uint32_t>(out, header.magic);
  write_le<std::uint16_t>(out, header.version);
  write_le<std::uint16_t>(out, header.endian_marker);
  write_le<std::uint32_t>(out, header.manifest_crc32);
  write_le<std::uint64_t>(out, header.metadata_offset);
  write_le<std::uint64_t>(out, header.router_directory_offset);
  write_le<std::uint64_t>(out, header.dictionary_directory_offset);
  write_le<std::uint64_t>(out, header.posting_directory_offset);
  write_le<std::uint64_t>(out, header.checksum_footer_offset);
}

void write_bytes(std::ostream& out, std::span<const std::uint8_t> bytes) {
  constexpr std::size_t max_chunk = 1U << 30U;
  while (!bytes.empty()) {
    const auto chunk_size = std::min(bytes.size(), max_chunk);
    out.write(reinterpret_cast<const char*>(bytes.data()),
              static_cast<std::streamsize>(chunk_size));
    if (!out) {
      throw std::runtime_error("failed while writing AIMS index payload");
    }
    bytes = bytes.subspan(chunk_size);
  }
}

void write_string(std::ostream& out, const std::string& value) {
  write_le<std::uint64_t>(out, value.size());
  write_bytes(out, std::span<const std::uint8_t>(
                       reinterpret_cast<const std::uint8_t*>(value.data()), value.size()));
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
      offset > bytes.size() || static_cast<std::size_t>(size) > bytes.size() - offset) {
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
  if (header.checksum_footer_offset >
          static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    throw std::runtime_error("AIMS checksum footer offset is too large for this platform");
  }
  const auto footer_offset = static_cast<std::size_t>(header.checksum_footer_offset);
  if (footer_offset > bytes.size() || bytes.size() - footer_offset != sizeof(std::uint64_t)) {
    throw std::runtime_error("invalid AIMS checksum footer offset");
  }

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
            offset > bytes.size() ||
            static_cast<std::size_t>(encoded_size) > bytes.size() - offset) {
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
            offset > bytes.size() ||
            static_cast<std::size_t>(encoded_size) > bytes.size() - offset) {
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

std::uint64_t output_offset(std::ostream& out) {
  const auto position = out.tellp();
  if (position < 0) {
    throw std::runtime_error("failed to determine AIMS index output offset");
  }
  return static_cast<std::uint64_t>(position);
}

std::filesystem::path temporary_index_path(const std::filesystem::path& destination) {
  static std::atomic<std::uint64_t> counter{0};
  const auto stamp = static_cast<std::uint64_t>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  const auto parent = destination.has_parent_path() ? destination.parent_path()
                                                     : std::filesystem::path{"."};
  const auto base = destination.filename().string();
  for (std::uint64_t attempt = 0; attempt < 100; ++attempt) {
    const auto suffix = stamp + counter.fetch_add(1, std::memory_order_relaxed) + attempt;
    auto candidate = parent / (base + ".tmp." + std::to_string(suffix));
    if (!std::filesystem::exists(candidate)) {
      return candidate;
    }
  }
  throw std::runtime_error("failed to allocate a temporary AIMS index path");
}

class TemporaryIndexGuard {
public:
  explicit TemporaryIndexGuard(std::filesystem::path path) : path_(std::move(path)) {}
  ~TemporaryIndexGuard() {
    if (!committed_) {
      std::error_code ignored;
      std::filesystem::remove(path_, ignored);
    }
  }

  TemporaryIndexGuard(const TemporaryIndexGuard&) = delete;
  TemporaryIndexGuard& operator=(const TemporaryIndexGuard&) = delete;

  [[nodiscard]] const std::filesystem::path& path() const noexcept { return path_; }
  void commit() noexcept { committed_ = true; }

private:
  std::filesystem::path path_;
  bool committed_{false};
};

std::uint64_t checksum_file(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to reopen temporary AIMS index for checksum");
  }
  std::uint64_t checksum = checksum_seed;
  std::array<char, 1U << 20U> buffer{};
  while (in) {
    in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const auto count = in.gcount();
    for (std::streamsize i = 0; i < count; ++i) {
      checksum_append(checksum, static_cast<std::uint8_t>(buffer[static_cast<std::size_t>(i)]));
    }
  }
  if (!in.eof()) {
    throw std::runtime_error("failed while checksumming temporary AIMS index");
  }
  return checksum;
}

void write_index_streaming(const std::filesystem::path& path,
                           std::uint64_t document_count,
                           const std::string& source_uri,
                           const std::string& source_checksum,
                           const std::string& build_command,
                           std::span<const build::SequenceMetadata> sequences,
                           std::span<const build::FixedKIndex> layers) {
  if (!is_little_endian()) {
    throw std::runtime_error("AIMS index writer currently requires little-endian host order");
  }
  if (layers.empty()) {
    throw std::runtime_error("KmerExactIndex requires at least one k-mer layer");
  }
  for (const auto& layer : layers) {
    if (layer.dictionary.size() != layer.postings.size()) {
      throw std::runtime_error("dictionary and posting store size mismatch");
    }
  }

  TemporaryIndexGuard temporary(temporary_index_path(path));
  std::ofstream out(temporary.path(), std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("failed to open temporary AIMS index: " +
                             temporary.path().string());
  }

  constexpr std::uint64_t header_size = 52;
  std::array<std::uint8_t, header_size> placeholder{};
  write_bytes(out, placeholder);

  IndexHeader header;
  header.metadata_offset = header_size;
  header.router_directory_offset = 0;
  write_le<std::uint64_t>(out, document_count);
  write_le<std::uint64_t>(out, layers.size());
  write_string(out, source_uri);
  write_string(out, source_checksum);
  write_string(out, build_command);
  write_le<std::uint64_t>(out, sequences.size());
  for (const auto& sequence : sequences) {
    write_le<std::uint32_t>(out, sequence.document_id);
    write_le<std::uint64_t>(out, sequence.sequence_id);
    write_le<std::uint64_t>(out, sequence.length);
    write_string(out, sequence.name);
  }

  header.dictionary_directory_offset = output_offset(out);
  bool posting_offset_set = false;
  for (const auto& layer : layers) {
    write_le<std::uint16_t>(out, layer.k);
    write_le<std::uint16_t>(out, 0);
    write_le<std::uint64_t>(out, layer.document_count);
    write_le<std::uint64_t>(out, layer.dictionary.size());
    const auto keys = layer.dictionary.keys();
    for (std::size_t i = 0; i < keys.size(); ++i) {
      const auto& key = keys[i];
      const auto& metadata = layer.dictionary.metadata(static_cast<SeedId>(i));
      write_le<std::uint64_t>(out, key.value);
      write_le<std::uint8_t>(out, static_cast<std::uint8_t>(key.family));
      write_le<std::uint16_t>(out, key.k);
      write_le<std::uint64_t>(out, metadata.document_frequency);
      write_le<std::uint64_t>(out, metadata.collection_frequency);
      write_le<std::uint64_t>(out, metadata.estimated_posting_bytes);
      write_le<std::uint8_t>(out, static_cast<std::uint8_t>(metadata.frequency_class));
    }

    if (!posting_offset_set) {
      header.posting_directory_offset = output_offset(out);
      posting_offset_set = true;
    }
    write_le<std::uint64_t>(out, layer.postings.size());
    for (SeedId i = 0; i < layer.postings.size(); ++i) {
      const auto block = layer.postings.encoded_block(i);
      write_le<std::uint64_t>(out, layer.postings.posting_count(i));
      write_le<std::uint64_t>(out, block.size());
      write_bytes(out, block);
    }
  }

  header.checksum_footer_offset = output_offset(out);
  out.seekp(0);
  write_header(out, header);
  out.flush();
  if (!out) {
    throw std::runtime_error("failed to finalize temporary AIMS index");
  }
  out.close();

  const auto checksum = checksum_file(temporary.path());
  std::ofstream footer(temporary.path(), std::ios::binary | std::ios::app);
  if (!footer) {
    throw std::runtime_error("failed to append AIMS index checksum");
  }
  write_le<std::uint64_t>(footer, checksum);
  footer.close();

  std::error_code rename_error;
  std::filesystem::rename(temporary.path(), path, rename_error);
  if (rename_error) {
    throw std::runtime_error("failed to install AIMS index atomically: " +
                             rename_error.message());
  }
  temporary.commit();
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
  write_index_streaming(path, index.document_count, {}, {}, {}, {},
                        std::span<const build::FixedKIndex>(&index, 1));
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
  write_index_streaming(path, index.document_count, index.source_uri, index.source_checksum,
                        index.build_command, index.sequences, index.layers);
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
