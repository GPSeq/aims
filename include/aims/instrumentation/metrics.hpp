#pragma once

#include <chrono>
#include <cstdint>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "aims/types.hpp"

namespace aims::instrumentation {

enum class ComparisonStage {
  ExactRetrieval,
  Pseudoalignment,
  ProbabilisticRouting,
  RoutingOnly,
  EndToEnd,
};

struct LayerMetrics {
  std::uint16_t k{0};
  std::uint64_t seeds_generated{0};
  std::uint64_t seeds_queried{0};
  std::uint64_t seeds_skipped_hot{0};
  std::uint64_t seeds_skipped_frequency_class{0};
  std::uint64_t seeds_doc_only_hot{0};
  std::uint64_t seeds_skipped_budget{0};
  std::uint64_t exact_bytes_read{0};
  std::uint64_t postings_decoded{0};
};

struct TopKResult {
  std::uint32_t rank{0};
  DocumentId document_id{0};
  SequenceId sequence_id{0};
  std::string sequence_name;
  Strand strand{Strand::Forward};
  std::uint32_t supporting_seeds{0};
  double score{0.0};
};

struct QueryMetrics {
  std::string query_id;
  std::string plan_label;
  std::uint64_t seeds_generated{0};
  std::uint64_t seeds_queried{0};
  std::uint64_t seeds_skipped_hot{0};
  std::uint64_t seeds_skipped_frequency_class{0};
  std::uint64_t seeds_doc_only_hot{0};
  std::uint64_t seeds_skipped_budget{0};
  std::uint64_t postings_skipped_budget{0};
  std::uint64_t router_seeds_used{0};
  std::uint64_t bins_returned{0};
  std::uint64_t router_bytes_read{0};
  std::uint64_t exact_bytes_read{0};
  std::uint64_t postings_decoded{0};
  std::uint64_t candidates_created{0};
  std::uint64_t candidates_touched{0};
  std::uint64_t peak_accumulator_bytes{0};
  bool early_stop_fired{false};
  std::string early_stop_rule{"none"};
  std::chrono::nanoseconds elapsed{0};
  ComparisonStage stage{ComparisonStage::ExactRetrieval};
  std::vector<LayerMetrics> per_k;
  std::vector<TopKResult> topk_results;
};

std::string stage_name(ComparisonStage stage);
void write_query_metrics_tsv_header(std::ostream& out);

/**
 * @fn write_query_metrics_tsv
 * @brief Write one QueryMetrics record as a TSV row.
 * @signature void write_query_metrics_tsv(std::ostream& out, const QueryMetrics& metrics);
 * @param out Output stream receiving the TSV row.
 * @param metrics Query metrics to serialize.
 * @throws std::ios_base::failure if the stream is configured to throw and writing fails.
 * @return None.
 */
void write_query_metrics_tsv(std::ostream& out, const QueryMetrics& metrics);

/**
 * @fn write_query_metrics_jsonl
 * @brief Write one QueryMetrics record as a JSONL object.
 * @signature void write_query_metrics_jsonl(std::ostream& out, const QueryMetrics& metrics);
 * @param out Output stream receiving the JSON object and trailing newline.
 * @param metrics Query metrics to serialize.
 * @throws std::ios_base::failure if the stream is configured to throw and writing fails.
 * @return None.
 */
void write_query_metrics_jsonl(std::ostream& out, const QueryMetrics& metrics);

} // namespace aims::instrumentation
