#include "aims/query/query_planner.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace aims::query {

QueryPlan QueryPlanner::plan(std::vector<SeedOccurrence> generated,
                             const index::ExactDictionary& dictionary,
                             std::uint64_t document_count) const {
  QueryPlan plan;
  plan.seeds_generated = generated.size();
  if (!generated.empty()) {
    plan.family = generated.front().key.family;
    plan.k = generated.front().key.k;
  }

  for (auto& occurrence : generated) {
    PlannedSeed planned;
    planned.occurrence = occurrence;
    const auto seed_id = dictionary.id(occurrence.key);
    if (!seed_id.has_value()) {
      planned.skip_reason = "absent_from_exact_dictionary";
      plan.seeds.push_back(std::move(planned));
      continue;
    }

    const auto& meta = dictionary.metadata(*seed_id);
    const double idf = std::log((static_cast<double>(document_count) + 1.0) /
                                (static_cast<double>(meta.document_frequency) + 1.0));
    const double survival_estimate = 1.0;
    const double length_bonus = std::max(1.0, static_cast<double>(occurrence.key.k) / 15.0);
    const double family_reliability = 1.0;
    planned.information = idf * survival_estimate * length_bonus * family_reliability;
    planned.estimated_cost =
        static_cast<double>(std::max<std::uint64_t>(1, meta.estimated_posting_bytes));
    planned.priority = planned.information / planned.estimated_cost;
    plan.seeds.push_back(std::move(planned));
  }

  std::stable_sort(plan.seeds.begin(), plan.seeds.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs.skip_reason.empty() != rhs.skip_reason.empty()) {
      return lhs.skip_reason.empty();
    }
    return lhs.priority > rhs.priority;
  });

  std::ostringstream label;
  label << "family=kmer"
        << ";k=" << plan.k
        << ";ordering=information_per_estimated_byte";
  plan.label = label.str();
  return plan;
}

} // namespace aims::query
