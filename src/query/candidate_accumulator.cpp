#include "aims/query/candidate_accumulator.hpp"

#include <algorithm>
#include <functional>

namespace aims::query {

CandidateAccumulator::CandidateAccumulator(CandidateAccumulatorOptions options)
    : options_(options) {}

std::size_t CandidateAccumulator::CandidateKeyHash::operator()(const CandidateKey& key) const noexcept {
  std::size_t value = std::hash<std::uint32_t>{}(key.document_id);
  value ^= std::hash<std::uint64_t>{}(key.sequence_id) + 0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
  value ^= std::hash<std::uint8_t>{}(static_cast<std::uint8_t>(key.strand)) +
           0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
  return value;
}

void CandidateAccumulator::add_postings(index::PostingView postings,
                                        double seed_weight,
                                        instrumentation::QueryMetrics& metrics) {
  metrics.postings_decoded += postings.size();
  metrics.exact_bytes_read += postings.bytes();

  for (const auto& posting : postings) {
    add_posting(posting, seed_weight, metrics);
  }
  metrics.peak_accumulator_bytes =
      std::max(metrics.peak_accumulator_bytes, estimated_bytes());
}

void CandidateAccumulator::add_posting(const index::Posting& posting,
                                       double seed_weight,
                                       instrumentation::QueryMetrics& metrics,
                                       std::uint32_t support_increment) {
  const auto key = CandidateKey{
      .document_id = posting.document_id,
      .sequence_id = posting.sequence_id,
      .strand = posting.strand,
  };
  auto it = candidates_.find(key);
  if (it == candidates_.end()) {
    if (options_.max_candidates != 0 && candidates_.size() >= options_.max_candidates) {
      ++metrics.candidates_touched;
      return;
    }
    candidates_.emplace(key, CandidateScore{
        .document_id = posting.document_id,
        .sequence_id = posting.sequence_id,
        .strand = posting.strand,
        .supporting_seeds = support_increment,
        .score = seed_weight,
    });
    ++metrics.candidates_created;
  } else {
    it->second.supporting_seeds += support_increment;
    it->second.score += seed_weight;
  }
  ++metrics.candidates_touched;
  metrics.peak_accumulator_bytes =
      std::max(metrics.peak_accumulator_bytes, estimated_bytes());
}

std::vector<CandidateScore> CandidateAccumulator::top_k(std::uint32_t k) const {
  std::vector<CandidateScore> out;
  out.reserve(candidates_.size());
  for (const auto& [_, candidate] : candidates_) {
    out.push_back(candidate);
  }
  std::sort(out.begin(), out.end(), [](const CandidateScore& lhs, const CandidateScore& rhs) {
    if (lhs.score != rhs.score) {
      return lhs.score > rhs.score;
    }
    if (lhs.supporting_seeds != rhs.supporting_seeds) {
      return lhs.supporting_seeds > rhs.supporting_seeds;
    }
    return lhs < rhs;
  });
  if (out.size() > k) {
    out.resize(k);
  }
  return out;
}

std::uint64_t CandidateAccumulator::candidate_count() const noexcept {
  return candidates_.size();
}

std::uint64_t CandidateAccumulator::estimated_bytes() const noexcept {
  return candidates_.bucket_count() * sizeof(void*) +
         candidates_.size() * (sizeof(CandidateKey) + sizeof(CandidateScore));
}

} // namespace aims::query
