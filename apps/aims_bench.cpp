#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/resource.h>

#include "aims/build/index_builder.hpp"
#include "aims/instrumentation/metrics.hpp"
#include "aims/io/fasta.hpp"
#include "aims/query/kmer_search.hpp"
#include "aims/serialization/index_format.hpp"

namespace {

std::string arg_or(int argc, char** argv, const std::string& name, const std::string& fallback) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == name) {
      return argv[i + 1];
    }
  }
  return fallback;
}

std::string require_arg(int argc, char** argv, const std::string& name) {
  const auto value = arg_or(argc, argv, name, {});
  if (value.empty()) {
    throw std::runtime_error("missing required argument: " + name);
  }
  return value;
}

aims::query::HotSeedMode parse_hot_seed_mode(const std::string& value) {
  if (value == "skip") {
    return aims::query::HotSeedMode::Skip;
  }
  if (value == "doc-only") {
    return aims::query::HotSeedMode::DocOnly;
  }
  throw std::runtime_error("invalid --hot-mode, expected skip or doc-only");
}

aims::FrequencyClass parse_frequency_class(const std::string& value) {
  if (value == "rare") {
    return aims::FrequencyClass::Rare;
  }
  if (value == "medium") {
    return aims::FrequencyClass::Medium;
  }
  if (value == "hot") {
    return aims::FrequencyClass::Hot;
  }
  if (value == "very-hot" || value == "very_hot") {
    return aims::FrequencyClass::VeryHot;
  }
  throw std::runtime_error("invalid --hot-seed-class, expected rare, medium, hot, or very-hot");
}

bool has_flag(int argc, char** argv, const std::string& name) {
  for (int i = 1; i < argc; ++i) {
    if (argv[i] == name) {
      return true;
    }
  }
  return false;
}

void print_help(std::ostream& out) {
  out << "Usage: aims_bench --ref refs.fa --query queries.fa --k 15,19 [--topk 10] [--truth truth.tsv]\n"
      << "       [--mmap] [--posting-cache-blocks N]\n"
      << "       [--repeats N] [--query-metrics-out metrics.jsonl]\n"
      << "       [--max-seeds N] [--max-postings N] [--max-candidates N]\n"
      << "       [--hot-seed-threshold N] [--hot-seed-class hot|very-hot] [--hot-mode skip|doc-only]\n"
      << "Stage: exact_retrieval. Emits one BenchmarkResult JSON object.\n";
}

std::vector<std::uint16_t> parse_k_values(const std::string& value) {
  std::set<std::uint16_t> unique;
  std::stringstream input(value);
  std::string token;
  while (std::getline(input, token, ',')) {
    if (token.empty()) {
      continue;
    }
    const auto parsed = std::stoul(token);
    if (parsed == 0 || parsed > 32) {
      throw std::runtime_error("KmerExactIndex supports 1 <= k <= 32");
    }
    unique.insert(static_cast<std::uint16_t>(parsed));
  }
  if (unique.empty()) {
    throw std::runtime_error("at least one k value is required");
  }
  return {unique.begin(), unique.end()};
}

std::string escape_json(const std::string& value) {
  std::string out;
  out.reserve(value.size());
  for (const char c : value) {
    switch (c) {
    case '\\':
      out += "\\\\";
      break;
    case '"':
      out += "\\\"";
      break;
    default:
      out.push_back(c);
      break;
    }
  }
  return out;
}

std::string command_line(int argc, char** argv) {
  std::ostringstream out;
  for (int i = 0; i < argc; ++i) {
    out << (i == 0 ? "" : " ") << argv[i];
  }
  return out.str();
}

std::string fnv1a_file_checksum(const std::filesystem::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("failed to checksum file: " + path.string());
  }
  std::uint64_t hash = 1469598103934665603ULL;
  char c = 0;
  while (in.get(c)) {
    hash ^= static_cast<unsigned char>(c);
    hash *= 1099511628211ULL;
  }
  std::ostringstream out;
  out << std::hex << std::setw(16) << std::setfill('0') << hash;
  return out.str();
}

std::set<aims::DocumentId> parse_document_ids(const std::string& value) {
  std::set<aims::DocumentId> ids;
  std::stringstream input(value);
  std::string token;
  while (std::getline(input, token, ',')) {
    if (!token.empty()) {
      ids.insert(static_cast<aims::DocumentId>(std::stoul(token)));
    }
  }
  return ids;
}

std::map<std::string, std::set<aims::DocumentId>>
read_truth_tsv(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open truth TSV: " + path.string());
  }
  std::map<std::string, std::set<aims::DocumentId>> truth;
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (first) {
      first = false;
      if (line.rfind("query_id\t", 0) == 0) {
        continue;
      }
    }
    std::stringstream row(line);
    std::string query_id;
    std::string stage;
    std::string expected_ids;
    std::getline(row, query_id, '\t');
    std::getline(row, stage, '\t');
    std::getline(row, expected_ids, '\t');
    if (stage != "exact_retrieval") {
      throw std::runtime_error("truth TSV contains non exact_retrieval stage");
    }
    truth[query_id] = parse_document_ids(expected_ids);
  }
  return truth;
}

bool has_truth_hit(const aims::instrumentation::QueryMetrics& metrics,
                   const std::set<aims::DocumentId>& expected,
                   std::uint32_t topn) {
  if (expected.empty()) {
    return metrics.topk_results.empty();
  }
  for (const auto& result : metrics.topk_results) {
    if (result.rank > topn) {
      continue;
    }
    if (expected.contains(result.document_id)) {
      return true;
    }
  }
  return false;
}

double mean(const std::vector<double>& values) {
  if (values.empty()) {
    return 0.0;
  }
  return std::accumulate(values.begin(), values.end(), 0.0) /
         static_cast<double>(values.size());
}

double percentile(std::vector<double> values, double p) {
  if (values.empty()) {
    return 0.0;
  }
  std::sort(values.begin(), values.end());
  const auto rank = static_cast<std::size_t>(
      std::clamp(p, 0.0, 1.0) * static_cast<double>(values.size() - 1));
  return values[rank];
}

double process_cpu_ms() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) {
    return 0.0;
  }
  const auto user_ms = static_cast<double>(usage.ru_utime.tv_sec) * 1000.0 +
                       static_cast<double>(usage.ru_utime.tv_usec) / 1000.0;
  const auto system_ms = static_cast<double>(usage.ru_stime.tv_sec) * 1000.0 +
                         static_cast<double>(usage.ru_stime.tv_usec) / 1000.0;
  return user_ms + system_ms;
}

std::uint64_t peak_rss_bytes() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0) {
    return 0;
  }
#if defined(__APPLE__)
  return static_cast<std::uint64_t>(usage.ru_maxrss);
#else
  return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024ULL;
#endif
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (has_flag(argc, argv, "--help") || has_flag(argc, argv, "-h")) {
      print_help(std::cout);
      return 0;
    }
    const auto ref_path = std::filesystem::path(require_arg(argc, argv, "--ref"));
    const auto query_path = std::filesystem::path(require_arg(argc, argv, "--query"));
    const auto k_values = parse_k_values(require_arg(argc, argv, "--k"));
    const auto topk = static_cast<std::uint32_t>(std::stoul(arg_or(argc, argv, "--topk", "10")));
    const auto dataset_name = arg_or(argc, argv, "--dataset", "synthetic_kmer_fixture");
    const auto truth_path_arg = arg_or(argc, argv, "--truth", "");
    const auto max_seeds = std::stoull(arg_or(argc, argv, "--max-seeds", "0"));
    const auto max_postings = std::stoull(arg_or(argc, argv, "--max-postings", "0"));
    const auto max_candidates = std::stoull(arg_or(argc, argv, "--max-candidates", "0"));
    const auto hot_seed_threshold = std::stoull(arg_or(argc, argv, "--hot-seed-threshold", "0"));
    const auto hot_seed_class_arg = arg_or(argc, argv, "--hot-seed-class", "");
    const auto hot_seed_mode = parse_hot_seed_mode(arg_or(argc, argv, "--hot-mode", "skip"));
    const auto use_mmap = has_flag(argc, argv, "--mmap");
    const auto posting_cache_blocks =
        std::stoull(arg_or(argc, argv, "--posting-cache-blocks", "0"));
    const auto repeats = std::max<std::uint64_t>(1, std::stoull(arg_or(argc, argv, "--repeats", "1")));
    const auto query_metrics_out = arg_or(argc, argv, "--query-metrics-out", "");
    const bool budgeted_run = max_seeds != 0 || max_postings != 0 || max_candidates != 0 ||
                              hot_seed_threshold != 0 || !hot_seed_class_arg.empty();

    const auto records = aims::io::read_sequences(ref_path);
    const auto queries = aims::io::read_sequences(query_path);
    const auto reference_checksum = "fnv1a64:" + fnv1a_file_checksum(ref_path);

    const auto build_cpu_start = process_cpu_ms();
    const auto build_start = std::chrono::steady_clock::now();
    auto index = aims::build::build_kmer_exact(records, k_values);
    index.source_uri = ref_path.string();
    index.source_checksum = reference_checksum;
    index.build_command = command_line(argc, argv);
    const auto build_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - build_start);
    const auto build_cpu_elapsed = process_cpu_ms() - build_cpu_start;

    const auto temp_index_path = std::filesystem::temp_directory_path() /
                                 ("aims_bench_" + dataset_name + ".aims");
    aims::serialization::write_kmer_exact_index(temp_index_path, index);
    const auto serialized_size = std::filesystem::file_size(temp_index_path);
    const auto load_cpu_start = process_cpu_ms();
    const auto load_start = std::chrono::steady_clock::now();
    const auto loaded_index =
        use_mmap ? aims::serialization::read_kmer_exact_index_mmap(temp_index_path,
                                                                   posting_cache_blocks)
                 : aims::serialization::read_kmer_exact_index(temp_index_path);
    const auto load_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - load_start);
    const auto load_cpu_elapsed = process_cpu_ms() - load_cpu_start;
    if (!use_mmap && posting_cache_blocks != 0) {
      for (const auto& layer : loaded_index.layers) {
        layer.postings.set_cache_limit(posting_cache_blocks);
      }
    }
    std::filesystem::remove(temp_index_path);

    std::vector<double> latencies;
    std::vector<double> bytes_read;
    std::vector<double> postings_decoded;
    std::vector<double> seeds_generated;
    std::vector<double> seeds_used;
    std::vector<double> seeds_skipped_hot;
    std::vector<double> seeds_skipped_frequency_class;
    std::vector<double> seeds_doc_only_hot;
    std::vector<double> seeds_skipped_budget;
    std::vector<double> postings_skipped_budget;
    std::vector<double> candidate_count;
    std::ofstream query_metrics_file;
    if (!query_metrics_out.empty()) {
      query_metrics_file.open(query_metrics_out);
      if (!query_metrics_file) {
        throw std::runtime_error("failed to open query metrics output: " + query_metrics_out);
      }
    }
    const auto truth = truth_path_arg.empty()
                           ? std::map<std::string, std::set<aims::DocumentId>>{}
                           : read_truth_tsv(std::filesystem::path(truth_path_arg));
    std::uint64_t truth_queries = 0;
    std::uint64_t top1_hits = 0;
    std::uint64_t top5_hits = 0;
    std::uint64_t top10_hits = 0;
    latencies.reserve(queries.size());

    const auto query_cpu_start = process_cpu_ms();
    const auto query_start = std::chrono::steady_clock::now();
    for (std::uint64_t repeat = 0; repeat < repeats; ++repeat) {
      for (const auto& query : queries) {
        const auto metrics = aims::query::search_kmer_exact(
            loaded_index, query, aims::query::KmerSearchOptions{
                .topk = topk,
                .max_seeds = max_seeds,
                .max_postings = max_postings,
                .max_candidates = max_candidates,
                .hot_seed_threshold = hot_seed_threshold,
                .use_frequency_class_hot_policy = !hot_seed_class_arg.empty(),
                .hot_seed_min_class = hot_seed_class_arg.empty()
                                          ? aims::FrequencyClass::VeryHot
                                          : parse_frequency_class(hot_seed_class_arg),
                .hot_seed_mode = hot_seed_mode,
            });
        if (query_metrics_file) {
          aims::instrumentation::write_query_metrics_jsonl(query_metrics_file, metrics);
        }
        latencies.push_back(static_cast<double>(metrics.elapsed.count()));
        bytes_read.push_back(static_cast<double>(metrics.router_bytes_read + metrics.exact_bytes_read));
        postings_decoded.push_back(static_cast<double>(metrics.postings_decoded));
        seeds_generated.push_back(static_cast<double>(metrics.seeds_generated));
        seeds_used.push_back(static_cast<double>(metrics.seeds_queried));
        seeds_skipped_hot.push_back(static_cast<double>(metrics.seeds_skipped_hot));
        seeds_skipped_frequency_class.push_back(
            static_cast<double>(metrics.seeds_skipped_frequency_class));
        seeds_doc_only_hot.push_back(static_cast<double>(metrics.seeds_doc_only_hot));
        seeds_skipped_budget.push_back(static_cast<double>(metrics.seeds_skipped_budget));
        postings_skipped_budget.push_back(static_cast<double>(metrics.postings_skipped_budget));
        candidate_count.push_back(static_cast<double>(metrics.candidates_created));
        if (repeat == 0) {
          const auto truth_it = truth.find(query.name);
          if (truth_it != truth.end()) {
            ++truth_queries;
            top1_hits += has_truth_hit(metrics, truth_it->second, 1) ? 1U : 0U;
            top5_hits += has_truth_hit(metrics, truth_it->second, 5) ? 1U : 0U;
            top10_hits += has_truth_hit(metrics, truth_it->second, 10) ? 1U : 0U;
          }
        }
      }
    }
    const auto query_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - query_start);
    const auto query_cpu_elapsed = process_cpu_ms() - query_cpu_start;

    const auto recall = [truth_queries](std::uint64_t hits) {
      return truth_queries == 0 ? 0.0 : static_cast<double>(hits) / static_cast<double>(truth_queries);
    };

    std::cout << "{"
              << "\"schema_version\":\"0.1.0\","
              << "\"run_id\":\"" << dataset_name << "_kmer_exact\","
              << "\"tool\":\"AIMS\","
              << "\"tool_version\":\"0.1.0\","
              << "\"comparison_stage\":\"exact_retrieval\","
              << "\"dataset\":{"
              << "\"name\":\"" << escape_json(dataset_name) << "\","
              << "\"reference_uri\":\"" << escape_json(ref_path.string()) << "\","
              << "\"reference_checksum\":\"" << reference_checksum << "\","
              << "\"query_uri\":\"" << escape_json(query_path.string()) << "\"";
    if (!truth_path_arg.empty()) {
      std::cout << ",\"truth_uri\":\"" << escape_json(truth_path_arg) << "\"";
    }
    std::cout
              << "},"
              << "\"cache_mode\":\"warm\","
              << "\"early_stop_mode\":\"" << (budgeted_run ? "heuristic" : "off") << "\","
              << "\"loading_mode\":\"" << (use_mmap ? "mmap" : "copy") << "\","
              << "\"repeats\":" << repeats << ","
              << "\"index\":{"
              << "\"seed_families\":[\"kmer\"],"
              << "\"k_values\":[";
    for (std::size_t i = 0; i < k_values.size(); ++i) {
      std::cout << (i == 0 ? "" : ",") << k_values[i];
    }
    std::cout << "],"
              << "\"serialized_size_bytes\":" << serialized_size << ","
              << "\"build_wall_ms\":" << build_elapsed.count() << ","
              << "\"build_cpu_ms\":" << build_cpu_elapsed << ","
              << "\"load_wall_ms\":" << load_elapsed.count() << ","
              << "\"load_cpu_ms\":" << load_cpu_elapsed << ","
              << "\"peak_rss_bytes\":" << peak_rss_bytes()
              << "},"
              << "\"metrics\":{"
              << "\"queries\":" << queries.size() << ","
              << "\"query_runs\":" << (queries.size() * repeats) << ","
              << "\"query_wall_ms_total\":" << query_elapsed.count() << ","
              << "\"query_cpu_ms_total\":" << query_cpu_elapsed << ","
              << "\"latency_ns_p50\":" << percentile(latencies, 0.50) << ","
              << "\"latency_ns_p90\":" << percentile(latencies, 0.90) << ","
              << "\"latency_ns_p95\":" << percentile(latencies, 0.95) << ","
              << "\"latency_ns_p99\":" << percentile(latencies, 0.99) << ","
              << "\"bytes_read_per_query\":" << mean(bytes_read) << ","
              << "\"postings_decoded_per_query\":" << mean(postings_decoded) << ","
              << "\"seeds_generated_per_query\":" << mean(seeds_generated) << ","
              << "\"seeds_used_per_query\":" << mean(seeds_used) << ","
              << "\"seeds_skipped_hot_per_query\":" << mean(seeds_skipped_hot) << ","
              << "\"seeds_skipped_frequency_class_per_query\":"
              << mean(seeds_skipped_frequency_class) << ","
              << "\"seeds_doc_only_hot_per_query\":" << mean(seeds_doc_only_hot) << ","
              << "\"seeds_skipped_budget_per_query\":" << mean(seeds_skipped_budget) << ","
              << "\"postings_skipped_budget_per_query\":" << mean(postings_skipped_budget) << ","
              << "\"candidate_count_per_query\":" << mean(candidate_count) << ","
              << "\"truth_queries\":" << truth_queries << ","
              << "\"top1_recall\":" << recall(top1_hits) << ","
              << "\"top5_recall\":" << recall(top5_hits) << ","
              << "\"top10_recall\":" << recall(top10_hits)
              << "}"
              << "}\n";

    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "aims_bench: " << ex.what() << '\n';
    return 1;
  }
}
