#pragma once

#include <cstdint>
#include <string>

namespace aims::scoring {

enum class EarlyStopMode : std::uint8_t {
  Off = 0,
  ConservativeBound = 1,
  Heuristic = 2,
};

struct EarlyStopDecision {
  bool stop{false};
  EarlyStopMode mode{EarlyStopMode::Off};
  std::string rule{"none"};
  double best_score{0.0};
  double runner_up_score{0.0};
  double remaining_upper_bound{0.0};
};

[[nodiscard]] inline EarlyStopDecision conservative_bound_decision(
    double best_score,
    double runner_up_score,
    double remaining_upper_bound) {
  return EarlyStopDecision{
      .stop = best_score > runner_up_score + remaining_upper_bound,
      .mode = EarlyStopMode::ConservativeBound,
      .rule = "best_gt_runner_up_plus_remaining_upper_bound",
      .best_score = best_score,
      .runner_up_score = runner_up_score,
      .remaining_upper_bound = remaining_upper_bound,
  };
}

} // namespace aims::scoring
