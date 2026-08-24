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

  /**
   * Add coordinate-aware evidence. Only the strongest exact coordinate diagonal contributes to
   * each candidate, preventing repetitive Cartesian seed matches from inflating its rank.
   */
  void add_positional_posting(const index::Posting& posting,
                              Position query_position,
                              std::uint16_t k,
                              double seed_weight,
                              instrumentation::QueryMetrics& metrics);

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

  struct CandidateEntry {
    CandidateScore score{};
    std::uint32_t base_support{0};
    double base_score{0.0};
    std::uint32_t best_diagonal_support{0};
    double best_diagonal_score{0.0};
  };

  struct DiagonalKey {
    CandidateKey candidate{};
    std::uint64_t diagonal{0};

    friend bool operator==(const DiagonalKey&, const DiagonalKey&) = default;
  };

  struct DiagonalKeyHash {
    [[nodiscard]] std::size_t operator()(const DiagonalKey& key) const noexcept;
  };

  struct DiagonalScore {
    std::uint32_t supporting_seeds{0};
    double score{0.0};
  };

  CandidateAccumulatorOptions options_{};
  std::unordered_map<CandidateKey, CandidateEntry, CandidateKeyHash> candidates_{};
  std::unordered_map<DiagonalKey, DiagonalScore, DiagonalKeyHash> diagonals_{};
};

} // namespace aims::query
