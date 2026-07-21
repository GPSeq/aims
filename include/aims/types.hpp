#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <functional>

namespace aims {

using DocumentId = std::uint32_t;
using SequenceId = std::uint64_t;
using Position = std::uint64_t;
using SeedId = std::uint64_t;
using EncodedSeed = std::uint64_t;

enum class Strand : std::uint8_t {
  Forward = 0,
  Reverse = 1,
};

enum class SeedFamily : std::uint8_t {
  Kmer = 0,
  // Reserved for format compatibility; the current implementation builds k-mer indexes only.
  Syncmer = 1,
  // Reserved for format compatibility; the current implementation builds k-mer indexes only.
  Strobemer = 2,
};

enum class FrequencyClass : std::uint8_t {
  Rare = 0,
  Medium = 1,
  Hot = 2,
  VeryHot = 3,
};

struct SequenceView {
  std::string_view bases;
  DocumentId document_id{0};
  SequenceId sequence_id{0};
  std::string_view name;
};

struct SeedKey {
  EncodedSeed value{0};
  SeedFamily family{SeedFamily::Kmer};
  std::uint16_t k{0};

  friend auto operator<=>(const SeedKey&, const SeedKey&) = default;
};

struct SeedOccurrence {
  SeedKey key{};
  DocumentId document_id{0};
  SequenceId sequence_id{0};
  Position position{0};
  Strand strand{Strand::Forward};
};

constexpr EncodedSeed invalid_seed = std::numeric_limits<EncodedSeed>::max();

} // namespace aims

template <>
struct std::hash<aims::SeedKey> {
  [[nodiscard]] std::size_t operator()(const aims::SeedKey& key) const noexcept {
    std::size_t value = std::hash<aims::EncodedSeed>{}(key.value);
    value ^= std::hash<std::uint8_t>{}(static_cast<std::uint8_t>(key.family)) +
             0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
    value ^= std::hash<std::uint16_t>{}(key.k) +
             0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
    return value;
  }
};
