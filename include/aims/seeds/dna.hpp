#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

#include "aims/types.hpp"

namespace aims::seeds {

std::optional<std::uint8_t> encode_base(char base) noexcept;
char decode_base(std::uint8_t bits) noexcept;
std::uint8_t complement_bits(std::uint8_t bits) noexcept;
char complement_base(char base) noexcept;

std::string reverse_complement(std::string_view sequence);

bool can_encode_k(std::uint16_t k) noexcept;
EncodedSeed forward_kmer(std::string_view kmer);
EncodedSeed reverse_complement_kmer(std::string_view kmer);

struct CanonicalSeed {
  EncodedSeed value{invalid_seed};
  Strand strand{Strand::Forward};
};

/**
 * @fn canonical_kmer
 * @brief Compute the canonical two-bit DNA representation for a k-mer.
 * @signature CanonicalSeed canonical_kmer(std::string_view kmer);
 * @param kmer DNA k-mer string containing only A, C, G, and T.
 * @throws None.
 * @return Canonical seed value and strand; invalid_seed is returned when the k-mer cannot be encoded.
 */
CanonicalSeed canonical_kmer(std::string_view kmer);

} // namespace aims::seeds
