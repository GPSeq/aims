#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "aims/build/index_builder.hpp"
#include "aims/io/fasta.hpp"
#include "aims/serialization/index_format.hpp"

namespace {

std::string require_arg(int argc, char** argv, const std::string& name) {
  for (int i = 1; i + 1 < argc; ++i) {
    if (argv[i] == name) {
      return argv[i + 1];
    }
  }
  throw std::runtime_error("missing required argument: " + name);
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
  out << "Usage: aims_build --ref refs.fa --out index.aims --k 15,19,23,27,31\n"
      << "Builds a checksum-validated KmerExactIndex for exact_retrieval benchmarks.\n"
      << "Options:\n"
      << "  --frequency-thresholds rare,medium,hot\n"
      << "      Build-time cf thresholds for rare/medium/hot classes; above hot is very_hot.\n";
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

aims::build::FrequencyThresholds parse_frequency_thresholds(const std::string& value) {
  std::vector<std::uint64_t> parsed_values;
  std::stringstream input(value);
  std::string token;
  while (std::getline(input, token, ',')) {
    if (token.empty()) {
      continue;
    }
    parsed_values.push_back(std::stoull(token));
  }
  if (parsed_values.size() != 3) {
    throw std::runtime_error("--frequency-thresholds expects exactly rare,medium,hot");
  }
  if (parsed_values[0] > parsed_values[1] || parsed_values[1] > parsed_values[2]) {
    throw std::runtime_error("--frequency-thresholds must be nondecreasing");
  }
  return aims::build::FrequencyThresholds{
      .rare_max = parsed_values[0],
      .medium_max = parsed_values[1],
      .hot_max = parsed_values[2],
  };
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

} // namespace

int main(int argc, char** argv) {
  try {
    if (has_flag(argc, argv, "--help") || has_flag(argc, argv, "-h")) {
      print_help(std::cout);
      return 0;
    }
    const auto ref_path = std::filesystem::path(require_arg(argc, argv, "--ref"));
    const auto out_path = std::filesystem::path(require_arg(argc, argv, "--out"));
    const auto k_values = parse_k_values(require_arg(argc, argv, "--k"));
    aims::build::KmerBuildOptions build_options;
    if (has_flag(argc, argv, "--frequency-thresholds")) {
      build_options.frequency_thresholds =
          parse_frequency_thresholds(require_arg(argc, argv, "--frequency-thresholds"));
    }

    const auto records = aims::io::read_sequences(ref_path);
    auto index = aims::build::build_kmer_exact(records, k_values, build_options);
    index.source_uri = ref_path.string();
    index.source_checksum = "fnv1a64:" + fnv1a_file_checksum(ref_path);
    index.build_command = command_line(argc, argv);
    aims::serialization::write_kmer_exact_index(out_path, index);

    std::uint64_t distinct_seeds = 0;
    for (const auto& layer : index.layers) {
      distinct_seeds += layer.dictionary.size();
    }
    std::cerr << "built KmerExactIndex comparison_stage=exact_retrieval records=" << records.size()
              << " layers=" << index.layers.size()
              << " distinct_seeds_total=" << distinct_seeds
              << " frequency_thresholds="
              << build_options.frequency_thresholds.rare_max << ","
              << build_options.frequency_thresholds.medium_max << ","
              << build_options.frequency_thresholds.hot_max
              << " k_values=";
    for (std::size_t i = 0; i < k_values.size(); ++i) {
      std::cerr << (i == 0 ? "" : ",") << k_values[i];
    }
    std::cerr << '\n';
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "aims_build: " << ex.what() << '\n';
    return 1;
  }
}
