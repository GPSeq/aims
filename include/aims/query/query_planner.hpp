#pragma once

#include <string>
#include <vector>

#include "aims/index/exact_dictionary.hpp"
#include "aims/seeds/seed_generator.hpp"

namespace aims::query {

struct PlannedSeed {
  SeedOccurrence occurrence{};
  double information{0.0};
  double estimated_cost{1.0};
  double priority{0.0};
  std::string skip_reason;
};

struct QueryPlan {
  SeedFamily family{SeedFamily::Kmer};
  std::uint16_t k{0};
  std::vector<PlannedSeed> seeds;
  std::uint64_t seeds_generated{0};
  std::string label;
};

class QueryPlanner {
public:
  [[nodiscard]] QueryPlan plan(
      std::vector<SeedOccurrence> generated,
      const index::ExactDictionary& dictionary,
      std::uint64_t document_count) const;
};

} // namespace aims::query
