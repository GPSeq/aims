#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "aims/index/postings.hpp"
#include "aims/instrumentation/metrics.hpp"

namespace aims::query {

struct CandidateScore {
  DocumentId document_id{0};
  SequenceId sequence_id{0};
  Strand strand{Strand::Forward};
  std::uint32_t supporting_seeds{0};
  double score{0.0};

  friend auto operator<=>(const CandidateScore&, const CandidateScore&) = default;
};

struct CandidateAccumulatorOptions {
  std::uint64_t max_candidates{0};
};

class CandidateAccumulator {
public:
  explicit CandidateAccumulator(CandidateAccumulatorOptions options = {});

  /**
   * @fn add_posting
   * @brief Add one decoded posting to the candidate accumulator.
   * @signature void add_posting(const index::Posting& posting, double seed_weight, instrumentation::QueryMetrics& metrics, std::uint32_t support_increment);
   * @param posting Decoded exact seed occurrence.
   * @param seed_weight Information-derived weight for the queried seed.
   * @param metrics Query metrics updated with candidate counts and memory estimates.
   * @param support_increment Number of deduplicated query occurrences represented by this posting.
   * @throws None.
   * @return None.
   */
  void add_posting(const index::Posting& posting,
                   double seed_weight,
                   instrumentation::QueryMetrics& metrics,
                   std::uint32_t support_increment = 1);

  void add_postings(index::PostingView postings,
                   double seed_weight,
                   instrumentation::QueryMetrics& metrics);

  /**
   * @fn top_k
   * @brief Return the highest-scoring candidate records.
   * @signature std::vector<CandidateScore> top_k(std::uint32_t k) const;
   * @param k Maximum number of candidates to return.
   * @throws None.
   * @return Candidate records sorted by score, supporting seeds, and stable identifiers.
   */
  [[nodiscard]] std::vector<CandidateScore> top_k(std::uint32_t k) const;
  [[nodiscard]] std::uint64_t candidate_count() const noexcept;
  [[nodiscard]] std::uint64_t estimated_bytes() const noexcept;

private:
  struct CandidateKey {
    DocumentId document_id{0};
    SequenceId sequence_id{0};
    Strand strand{Strand::Forward};

    friend bool operator==(const CandidateKey&, const CandidateKey&) = default;
  };

  struct CandidateKeyHash {
    [[nodiscard]] std::size_t operator()(const CandidateKey& key) const noexcept;
  };

  CandidateAccumulatorOptions options_{};
  std::unordered_map<CandidateKey, CandidateScore, CandidateKeyHash> candidates_{};
};

} // namespace aims::query
