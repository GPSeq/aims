#pragma once

#include <cstdint>
#include <span>
#include <vector>

#include "aims/instrumentation/metrics.hpp"
#include "aims/types.hpp"

namespace aims::router {

struct CandidateBin {
  std::uint32_t bin_id{0};
  double approximate_score{0.0};
};

struct RouteBudget {
  std::uint64_t max_seeds{0};
  std::uint64_t max_bytes{0};
};

class RouterBackend {
public:
  virtual ~RouterBackend() = default;
  [[nodiscard]] virtual std::vector<CandidateBin>
  route(std::span<const SeedKey> query_seed_set,
        RouteBudget budget,
        instrumentation::QueryMetrics& metrics) const = 0;
};

} // namespace aims::router
