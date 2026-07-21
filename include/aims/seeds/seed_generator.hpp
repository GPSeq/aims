#pragma once

#include <vector>

#include "aims/types.hpp"

namespace aims::seeds {

struct SeedGenerationParams {
  SeedFamily family{SeedFamily::Kmer};
  std::uint16_t k{0};
  bool canonical{true};
};

class SeedGenerator {
public:
  virtual ~SeedGenerator() = default;
  [[nodiscard]] virtual std::vector<SeedOccurrence>
  generate(SequenceView sequence, const SeedGenerationParams& params) const = 0;
};

} // namespace aims::seeds
