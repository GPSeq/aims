#include "aims/seeds/kmer_generator.hpp"

#include <cassert>

#include "aims/seeds/dna.hpp"

namespace aims::seeds {
namespace {

// Emit all valid k-mers from one contiguous unambiguous DNA run.
void flush_run(std::vector<SeedOccurrence>& out,
               std::string_view bases,
               Position run_start,
               SequenceView sequence,
               const SeedGenerationParams& params) {
  // A zero k or a run shorter than k cannot produce seeds.
  if (params.k == 0 || bases.size() < params.k) {
    return;
  }
  // This implementation is intentionally k-mer-only.
  assert(params.family == SeedFamily::Kmer);

  // Slide a window of length k across the valid run.
  for (std::size_t offset = 0; offset + params.k <= bases.size(); ++offset) {
    // The current k-mer view avoids copying sequence characters.
    const auto kmer = bases.substr(offset, params.k);
    // Canonical mode chooses min(forward, reverse-complement) and records strand.
    const auto canonical = params.canonical ? canonical_kmer(kmer)
                                            : CanonicalSeed{forward_kmer(kmer), Strand::Forward};
    // Invalid means an ambiguous base slipped through, so skip defensively.
    if (canonical.value == invalid_seed) {
      continue;
    }
    // Record the seed key, reference IDs, sequence-local coordinate, and strand.
    out.push_back(SeedOccurrence{
        .key = SeedKey{.value = canonical.value, .family = params.family, .k = params.k},
        .document_id = sequence.document_id,
        .sequence_id = sequence.sequence_id,
        .position = run_start + offset,
        .strand = canonical.strand,
    });
  }
}

} // namespace

std::vector<SeedOccurrence>
KmerGenerator::generate(SequenceView sequence, const SeedGenerationParams& params) const {
  // The two-bit representation fits k <= 32 inside a uint64_t.
  assert(params.k <= 32);
  // Store all generated occurrences for this sequence.
  std::vector<SeedOccurrence> out;
  // Reject invalid k values early.
  if (params.k == 0 || !can_encode_k(params.k)) {
    return out;
  }

  // run_begin marks the start of the current A/C/G/T-only segment.
  std::size_t run_begin = 0;
  // in_run tells whether the scan is currently inside a valid DNA segment.
  bool in_run = false;
  // Scan every base so ambiguous characters split seed extraction.
  for (std::size_t i = 0; i < sequence.bases.size(); ++i) {
    // Valid DNA bases extend or start a run.
    if (encode_base(sequence.bases[i]).has_value()) {
      // Starting a new run records its coordinate in the original sequence.
      if (!in_run) {
        run_begin = i;
        in_run = true;
      }
    // An ambiguous base closes the current run.
    } else if (in_run) {
      // Emit k-mers from the run before the ambiguous base.
      flush_run(out, sequence.bases.substr(run_begin, i - run_begin),
                static_cast<Position>(run_begin), sequence, params);
      // The scan is now outside a valid run.
      in_run = false;
    }
  }
  // If the sequence ended inside a valid run, flush the final segment.
  if (in_run) {
    flush_run(out, sequence.bases.substr(run_begin),
              static_cast<Position>(run_begin), sequence, params);
  }
  // Return all exact seed occurrences for this sequence.
  return out;
}

std::optional<std::uint8_t> encode_base(char base) noexcept {
  // Map uppercase and lowercase DNA bases into two-bit values.
  switch (base) {
  case 'A':
  case 'a':
    return 0U;
  case 'C':
  case 'c':
    return 1U;
  case 'G':
  case 'g':
    return 2U;
  case 'T':
  case 't':
    return 3U;
  default:
    // Ambiguous or non-DNA characters are deliberately not encoded.
    return std::nullopt;
  }
}

char decode_base(std::uint8_t bits) noexcept {
  // Lookup table for the inverse two-bit base encoding.
  static constexpr char table[] = {'A', 'C', 'G', 'T'};
  // Masking keeps the function total for any uint8_t input.
  return table[bits & 0x3U];
}

std::uint8_t complement_bits(std::uint8_t bits) noexcept {
  // With A=0,C=1,G=2,T=3, bitwise inversion over two bits gives the complement.
  return static_cast<std::uint8_t>((~static_cast<unsigned>(bits)) & 0x3U);
}

char complement_base(char base) noexcept {
  // Encode first so lowercase and uppercase share one complement implementation.
  const auto bits = encode_base(base);
  // Unknown bases are represented as N in reverse-complement strings.
  return bits.has_value() ? decode_base(complement_bits(*bits)) : 'N';
}

std::string reverse_complement(std::string_view sequence) {
  // Allocate the final reverse-complement string once.
  std::string out;
  // Preserve length exactly.
  out.reserve(sequence.size());
  // Read bases in reverse order and complement each one.
  for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
    out.push_back(complement_base(*it));
  }
  // Return the reverse-complement sequence.
  return out;
}

bool can_encode_k(std::uint16_t k) noexcept {
  // A uint64_t stores at most 32 two-bit bases.
  return k > 0 && k <= 32;
}

EncodedSeed forward_kmer(std::string_view kmer) {
  // Reject k values that cannot fit the selected integer representation.
  if (!can_encode_k(static_cast<std::uint16_t>(kmer.size()))) {
    return invalid_seed;
  }
  // Build the two-bit integer from left to right.
  EncodedSeed value = 0;
  // Encode every base.
  for (const char base : kmer) {
    // Convert A/C/G/T into two bits.
    const auto bits = encode_base(base);
    // Any ambiguous base invalidates the whole seed.
    if (!bits.has_value()) {
      return invalid_seed;
    }
    // Shift the previous bases and append the current base bits.
    value = (value << 2U) | static_cast<EncodedSeed>(*bits);
  }
  // Return the forward-strand encoding.
  return value;
}

EncodedSeed reverse_complement_kmer(std::string_view kmer) {
  // Reject k values that cannot fit the selected integer representation.
  if (!can_encode_k(static_cast<std::uint16_t>(kmer.size()))) {
    return invalid_seed;
  }
  // Build the reverse-complement two-bit integer from right to left.
  EncodedSeed value = 0;
  // Visit the input in reverse order.
  for (auto it = kmer.rbegin(); it != kmer.rend(); ++it) {
    // Convert the original base into two bits.
    const auto bits = encode_base(*it);
    // Any ambiguous base invalidates the whole seed.
    if (!bits.has_value()) {
      return invalid_seed;
    }
    // Append the complemented base bits.
    value = (value << 2U) | static_cast<EncodedSeed>(complement_bits(*bits));
  }
  // Return the reverse-complement encoding.
  return value;
}

CanonicalSeed canonical_kmer(std::string_view kmer) {
  // Encode the k-mer as written.
  const EncodedSeed fwd = forward_kmer(kmer);
  // Encode the reverse-complement k-mer.
  const EncodedSeed rev = reverse_complement_kmer(kmer);
  // If either representation failed, the canonical seed is invalid.
  if (fwd == invalid_seed || rev == invalid_seed) {
    return {};
  }
  // Choose the smaller integer so both DNA strands map to the same seed key.
  if (rev < fwd) {
    return CanonicalSeed{.value = rev, .strand = Strand::Reverse};
  }
  // Forward is used on ties, including palindromic k-mers.
  return CanonicalSeed{.value = fwd, .strand = Strand::Forward};
}

} // namespace aims::seeds
