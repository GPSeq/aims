#include "aims/query/candidate_accumulator.hpp"

#include <algorithm>
#include <functional>

namespace aims::query {

CandidateAccumulator::CandidateAccumulator(CandidateAccumulatorOptions options)
    : options_(options) {}

namespace {

bool candidate_ranks_before(const CandidateScore& lhs, const CandidateScore& rhs) {
  if (lhs.score != rhs.score) {
    return lhs.score > rhs.score;
  }
  if (lhs.supporting_seeds != rhs.supporting_seeds) {
    return lhs.supporting_seeds > rhs.supporting_seeds;
  }
  return lhs < rhs;
}

} // namespace

std::size_t CandidateAccumulator::CandidateKeyHash::operator()(const CandidateKey& key) const noexcept {
  std::size_t value = std::hash<std::uint32_t>{}(key.document_id);
  value ^= std::hash<std::uint64_t>{}(key.sequence_id) + 0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
  value ^= std::hash<std::uint8_t>{}(static_cast<std::uint8_t>(key.strand)) +
           0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
  return value;
}

std::size_t CandidateAccumulator::DiagonalKeyHash::operator()(const DiagonalKey& key) const noexcept {
  std::size_t value = CandidateKeyHash{}(key.candidate);
  value ^= std::hash<std::uint64_t>{}(key.diagonal) + 0x9e3779b97f4a7c15ULL + (value << 6U) +
           (value >> 2U);
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
    candidates_.emplace(key, CandidateEntry{
        .score = CandidateScore{
            .document_id = posting.document_id,
            .sequence_id = posting.sequence_id,
            .strand = posting.strand,
            .supporting_seeds = support_increment,
            .score = seed_weight,
        },
        .base_support = support_increment,
        .base_score = seed_weight,
    });
    ++metrics.candidates_created;
  } else {
    it->second.base_support += support_increment;
    it->second.base_score += seed_weight;
    it->second.score.supporting_seeds =
        it->second.base_support + it->second.best_diagonal_support;
    it->second.score.score = it->second.base_score + it->second.best_diagonal_score;
  }
  ++metrics.candidates_touched;
  metrics.peak_accumulator_bytes =
      std::max(metrics.peak_accumulator_bytes, estimated_bytes());
}

void CandidateAccumulator::add_positional_posting(const index::Posting& posting,
                                                  Position query_position,
                                                  std::uint16_t k,
                                                  double seed_weight,
                                                  instrumentation::QueryMetrics& metrics) {
  const auto candidate_key = CandidateKey{
      .document_id = posting.document_id,
      .sequence_id = posting.sequence_id,
      .strand = posting.strand,
  };
  auto candidate = candidates_.find(candidate_key);
  if (candidate == candidates_.end()) {
    if (options_.max_candidates != 0 && candidates_.size() >= options_.max_candidates) {
      ++metrics.candidates_touched;
      return;
    }
    candidate = candidates_.emplace(candidate_key, CandidateEntry{
                                                       .score = CandidateScore{
                                                           .document_id = posting.document_id,
                                                           .sequence_id = posting.sequence_id,
                                                           .strand = posting.strand,
                                                       },
                                                   }).first;
    ++metrics.candidates_created;
  }

  // Forward matches preserve coordinate direction. Reverse matches use ref + query + k so the
  // same reverse-complement placement has one k-independent diagonal across multiple layers.
  const auto diagonal = posting.strand == Strand::Forward
                            ? posting.position - query_position
                            : posting.position + query_position + static_cast<Position>(k);
  auto& bin = diagonals_[DiagonalKey{.candidate = candidate_key, .diagonal = diagonal}];
  ++bin.supporting_seeds;
  bin.score += seed_weight;
  if (bin.score > candidate->second.best_diagonal_score ||
      (bin.score == candidate->second.best_diagonal_score &&
       bin.supporting_seeds > candidate->second.best_diagonal_support)) {
    candidate->second.best_diagonal_score = bin.score;
    candidate->second.best_diagonal_support = bin.supporting_seeds;
    candidate->second.score.score = candidate->second.base_score + bin.score;
    candidate->second.score.supporting_seeds =
        candidate->second.base_support + bin.supporting_seeds;
  }
  ++metrics.candidates_touched;
  metrics.peak_accumulator_bytes =
      std::max(metrics.peak_accumulator_bytes, estimated_bytes());
}

std::vector<CandidateScore> CandidateAccumulator::top_k(std::uint32_t k) const {
  if (k == 0 || candidates_.empty()) {
    return {};
  }
  std::vector<CandidateScore> out;
  out.reserve(candidates_.size());
  for (const auto& [_, candidate] : candidates_) {
    out.push_back(candidate.score);
  }
  const auto result_size = std::min<std::size_t>(out.size(), k);
  std::partial_sort(out.begin(), out.begin() + static_cast<std::ptrdiff_t>(result_size), out.end(),
                    candidate_ranks_before);
  out.resize(result_size);
  return out;
}

std::uint64_t CandidateAccumulator::candidate_count() const noexcept {
  return candidates_.size();
}

std::uint64_t CandidateAccumulator::estimated_bytes() const noexcept {
  return candidates_.bucket_count() * sizeof(void*) +
         candidates_.size() * (sizeof(CandidateKey) + sizeof(CandidateEntry)) +
         diagonals_.bucket_count() * sizeof(void*) +
         diagonals_.size() * (sizeof(DiagonalKey) + sizeof(DiagonalScore));
}

} // namespace aims::query
