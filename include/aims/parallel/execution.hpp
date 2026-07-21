#pragma once

#include <cstdint>

namespace aims::parallel {

struct ExecutionPolicy {
  std::uint32_t threads{1};
  bool numa_aware{false};
  bool deterministic_order{true};
};

} // namespace aims::parallel
