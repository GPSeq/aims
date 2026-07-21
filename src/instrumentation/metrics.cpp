#include "aims/instrumentation/metrics.hpp"

#include <iomanip>

namespace aims::instrumentation {
namespace {

std::string escape_json(const std::string& input) {
  std::string out;
  out.reserve(input.size());
  for (const char c : input) {
    switch (c) {
    case '\\':
      out += "\\\\";
      break;
    case '"':
      out += "\\\"";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\t':
      out += "\\t";
      break;
    default:
      out.push_back(c);
      break;
    }
  }
  return out;
}

std::string strand_name(Strand strand) {
  return strand == Strand::Forward ? "forward" : "reverse";
}

void write_per_k_json(std::ostream& out, const std::vector<LayerMetrics>& layers) {
  out << "[";
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const auto& layer = layers[i];
    out << (i == 0 ? "" : ",")
        << "{"
        << "\"k\":" << layer.k << ","
        << "\"seeds_generated\":" << layer.seeds_generated << ","
        << "\"seeds_queried\":" << layer.seeds_queried << ","
        << "\"seeds_skipped_hot\":" << layer.seeds_skipped_hot << ","
        << "\"seeds_skipped_frequency_class\":" << layer.seeds_skipped_frequency_class << ","
        << "\"seeds_doc_only_hot\":" << layer.seeds_doc_only_hot << ","
        << "\"seeds_skipped_budget\":" << layer.seeds_skipped_budget << ","
        << "\"exact_bytes_read\":" << layer.exact_bytes_read << ","
        << "\"postings_decoded\":" << layer.postings_decoded
        << "}";
  }
  out << "]";
}

void write_topk_json(std::ostream& out, const std::vector<TopKResult>& results) {
  out << "[";
  for (std::size_t i = 0; i < results.size(); ++i) {
    const auto& result = results[i];
    out << (i == 0 ? "" : ",")
        << "{"
        << "\"rank\":" << result.rank << ","
        << "\"document_id\":" << result.document_id << ","
        << "\"sequence_id\":" << result.sequence_id << ","
        << "\"sequence_name\":\"" << escape_json(result.sequence_name) << "\","
        << "\"strand\":\"" << strand_name(result.strand) << "\","
        << "\"supporting_seeds\":" << result.supporting_seeds << ","
        << "\"score\":" << result.score
        << "}";
  }
  out << "]";
}

std::string compact_per_k(const std::vector<LayerMetrics>& layers) {
  std::ostringstream out;
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const auto& layer = layers[i];
    out << (i == 0 ? "" : ";")
        << "k=" << layer.k
        << ",seeds_generated=" << layer.seeds_generated
        << ",seeds_queried=" << layer.seeds_queried
        << ",seeds_skipped_hot=" << layer.seeds_skipped_hot
        << ",seeds_skipped_frequency_class=" << layer.seeds_skipped_frequency_class
        << ",seeds_doc_only_hot=" << layer.seeds_doc_only_hot
        << ",seeds_skipped_budget=" << layer.seeds_skipped_budget
        << ",exact_bytes_read=" << layer.exact_bytes_read
        << ",postings_decoded=" << layer.postings_decoded;
  }
  return out.str();
}

std::string compact_topk(const std::vector<TopKResult>& results) {
  std::ostringstream out;
  for (std::size_t i = 0; i < results.size(); ++i) {
    const auto& result = results[i];
    out << (i == 0 ? "" : ";")
        << "rank=" << result.rank
        << ",document_id=" << result.document_id
        << ",sequence_id=" << result.sequence_id
        << ",sequence_name=" << result.sequence_name
        << ",strand=" << strand_name(result.strand)
        << ",supporting_seeds=" << result.supporting_seeds
        << ",score=" << result.score;
  }
  return out.str();
}

} // namespace

std::string stage_name(ComparisonStage stage) {
  switch (stage) {
  case ComparisonStage::ExactRetrieval:
    return "exact_retrieval";
  case ComparisonStage::Pseudoalignment:
    return "pseudoalignment";
  case ComparisonStage::ProbabilisticRouting:
    return "probabilistic_routing";
  case ComparisonStage::RoutingOnly:
    return "routing_only";
  case ComparisonStage::EndToEnd:
    return "end_to_end";
  }
  return "unknown";
}

void write_query_metrics_tsv_header(std::ostream& out) {
  out << "query_id\tcomparison_stage\tplan_label\tseeds_generated\tseeds_queried"
      << "\tseeds_skipped_hot\tseeds_skipped_frequency_class\tseeds_doc_only_hot"
      << "\tseeds_skipped_budget\tpostings_skipped_budget"
      << "\trouter_seeds_used\tbins_returned\trouter_bytes_read\texact_bytes_read"
      << "\tpostings_decoded\tcandidates_created\tcandidates_touched"
      << "\tpeak_accumulator_bytes\tearly_stop_fired\tearly_stop_rule"
      << "\tper_k_metrics\ttopk_results\telapsed_ns\n";
}

void write_query_metrics_tsv(std::ostream& out, const QueryMetrics& metrics) {
  out << metrics.query_id << '\t'
      << stage_name(metrics.stage) << '\t'
      << metrics.plan_label << '\t'
      << metrics.seeds_generated << '\t'
      << metrics.seeds_queried << '\t'
      << metrics.seeds_skipped_hot << '\t'
      << metrics.seeds_skipped_frequency_class << '\t'
      << metrics.seeds_doc_only_hot << '\t'
      << metrics.seeds_skipped_budget << '\t'
      << metrics.postings_skipped_budget << '\t'
      << metrics.router_seeds_used << '\t'
      << metrics.bins_returned << '\t'
      << metrics.router_bytes_read << '\t'
      << metrics.exact_bytes_read << '\t'
      << metrics.postings_decoded << '\t'
      << metrics.candidates_created << '\t'
      << metrics.candidates_touched << '\t'
      << metrics.peak_accumulator_bytes << '\t'
      << (metrics.early_stop_fired ? "true" : "false") << '\t'
      << metrics.early_stop_rule << '\t'
      << compact_per_k(metrics.per_k) << '\t'
      << compact_topk(metrics.topk_results) << '\t'
      << metrics.elapsed.count() << '\n';
}

void write_query_metrics_jsonl(std::ostream& out, const QueryMetrics& metrics) {
  out << "{"
      << "\"query_id\":\"" << escape_json(metrics.query_id) << "\","
      << "\"comparison_stage\":\"" << stage_name(metrics.stage) << "\","
      << "\"plan_label\":\"" << escape_json(metrics.plan_label) << "\","
      << "\"seeds_generated\":" << metrics.seeds_generated << ","
      << "\"seeds_queried\":" << metrics.seeds_queried << ","
      << "\"seeds_skipped_hot\":" << metrics.seeds_skipped_hot << ","
      << "\"seeds_skipped_frequency_class\":" << metrics.seeds_skipped_frequency_class << ","
      << "\"seeds_doc_only_hot\":" << metrics.seeds_doc_only_hot << ","
      << "\"seeds_skipped_budget\":" << metrics.seeds_skipped_budget << ","
      << "\"postings_skipped_budget\":" << metrics.postings_skipped_budget << ","
      << "\"router_seeds_used\":" << metrics.router_seeds_used << ","
      << "\"bins_returned\":" << metrics.bins_returned << ","
      << "\"router_bytes_read\":" << metrics.router_bytes_read << ","
      << "\"exact_bytes_read\":" << metrics.exact_bytes_read << ","
      << "\"postings_decoded\":" << metrics.postings_decoded << ","
      << "\"candidates_created\":" << metrics.candidates_created << ","
      << "\"candidates_touched\":" << metrics.candidates_touched << ","
      << "\"peak_accumulator_bytes\":" << metrics.peak_accumulator_bytes << ","
      << "\"early_stop_fired\":" << (metrics.early_stop_fired ? "true" : "false") << ","
      << "\"early_stop_rule\":\"" << escape_json(metrics.early_stop_rule) << "\","
      << "\"per_k_metrics\":";
  write_per_k_json(out, metrics.per_k);
  out << ","
      << "\"topk_results\":";
  write_topk_json(out, metrics.topk_results);
  out << ","
      << "\"elapsed_ns\":" << metrics.elapsed.count()
      << "}\n";
}

} // namespace aims::instrumentation
