#pragma once

#include "aims/seeds/seed_generator.hpp"

namespace aims::seeds {

class KmerGenerator final : public SeedGenerator {
public:
  [[nodiscard]] std::vector<SeedOccurrence>
  generate(SequenceView sequence, const SeedGenerationParams& params) const override;
};

} // namespace aims::seeds
