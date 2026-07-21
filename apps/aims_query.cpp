#include <chrono>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <future>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

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
  out << "Usage: aims_query --index index.aims --query queries.fa [--topk 10] [--emit jsonl|tsv] [--out results.jsonl]\n"
      << "       [--mmap] [--posting-cache-blocks N]\n"
      << "       [--threads N]\n"
      << "       [--max-seeds N] [--max-postings N] [--max-candidates N]\n"
      << "       [--hot-seed-threshold N] [--hot-seed-class hot|very-hot] [--hot-mode skip|doc-only]\n"
      << "Stage: exact_retrieval. Output includes query metrics, per-k metrics, and structured top-k candidates.\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (has_flag(argc, argv, "--help") || has_flag(argc, argv, "-h")) {
      print_help(std::cout);
      return 0;
    }
    const auto index_path = std::filesystem::path(require_arg(argc, argv, "--index"));
    const auto query_path = std::filesystem::path(require_arg(argc, argv, "--query"));
    const auto topk = static_cast<std::uint32_t>(std::stoul(arg_or(argc, argv, "--topk", "10")));
    const auto emit = arg_or(argc, argv, "--emit", "jsonl");
    const auto max_seeds = std::stoull(arg_or(argc, argv, "--max-seeds", "0"));
    const auto max_postings = std::stoull(arg_or(argc, argv, "--max-postings", "0"));
    const auto max_candidates = std::stoull(arg_or(argc, argv, "--max-candidates", "0"));
    const auto hot_seed_threshold = std::stoull(arg_or(argc, argv, "--hot-seed-threshold", "0"));
    const auto hot_seed_class_arg = arg_or(argc, argv, "--hot-seed-class", "");
    const auto hot_seed_mode = parse_hot_seed_mode(arg_or(argc, argv, "--hot-mode", "skip"));
    const auto out_path = arg_or(argc, argv, "--out", "");
    const auto use_mmap = has_flag(argc, argv, "--mmap");
    const auto posting_cache_blocks =
        std::stoull(arg_or(argc, argv, "--posting-cache-blocks", "0"));
    const auto threads = std::max<std::uint64_t>(1, std::stoull(arg_or(argc, argv, "--threads", "1")));

    const auto index = use_mmap
                           ? aims::serialization::read_kmer_exact_index_mmap(index_path,
                                                                             posting_cache_blocks)
                           : aims::serialization::read_kmer_exact_index(index_path);
    if (!use_mmap && posting_cache_blocks != 0) {
      for (const auto& layer : index.layers) {
        layer.postings.set_cache_limit(posting_cache_blocks);
      }
    }
    const auto queries = aims::io::read_sequences(query_path);

    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!out_path.empty()) {
      out_file.open(out_path);
      if (!out_file) {
        throw std::runtime_error("failed to open output file: " + out_path);
      }
      out = &out_file;
    }

    if (emit == "tsv") {
      aims::instrumentation::write_query_metrics_tsv_header(*out);
    }

    const auto options = aims::query::KmerSearchOptions{
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
    };

    std::vector<aims::instrumentation::QueryMetrics> results(queries.size());
    if (threads == 1 || queries.size() <= 1) {
      for (std::size_t i = 0; i < queries.size(); ++i) {
        results[i] = aims::query::search_kmer_exact(index, queries[i], options);
      }
    } else {
      std::vector<std::future<void>> futures;
      std::size_t next = 0;
      while (next < queries.size()) {
        while (futures.size() < threads && next < queries.size()) {
          const auto query_index = next++;
          futures.push_back(std::async(std::launch::async, [&, query_index] {
            results[query_index] = aims::query::search_kmer_exact(index, queries[query_index], options);
          }));
        }
        futures.front().get();
        futures.erase(futures.begin());
      }
      for (auto& future : futures) {
        future.get();
      }
    }

    for (const auto& metrics : results) {

      if (emit == "tsv") {
        aims::instrumentation::write_query_metrics_tsv(*out, metrics);
      } else {
        aims::instrumentation::write_query_metrics_jsonl(*out, metrics);
      }

      std::cerr << "query=" << metrics.query_id
                << " comparison_stage=exact_retrieval"
                << " index=KmerExactIndex"
                << " topk";
      for (const auto& result : metrics.topk_results) {
        std::cerr << " doc=" << result.document_id
                  << ":seq=" << result.sequence_id
                  << ":name=" << result.sequence_name
                  << ":strand=" << (result.strand == aims::Strand::Forward ? "forward" : "reverse")
                  << ":support=" << result.supporting_seeds
                  << ":score=" << result.score;
      }
      std::cerr << '\n';
    }
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "aims_query: " << ex.what() << '\n';
    return 1;
  }
}
