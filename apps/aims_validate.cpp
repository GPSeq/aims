#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

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
  out << "Usage: aims_validate --index index.aims [--mmap]\n"
      << "Validates checksum, offsets, k-mer layers, and reference metadata.\n";
}

} // namespace

int main(int argc, char** argv) {
  try {
    if (has_flag(argc, argv, "--help") || has_flag(argc, argv, "-h")) {
      print_help(std::cout);
      return 0;
    }
    const auto path = std::filesystem::path(require_arg(argc, argv, "--index"));
    const auto index = has_flag(argc, argv, "--mmap")
                           ? aims::serialization::read_kmer_exact_index_mmap(path)
                           : aims::serialization::read_kmer_exact_index(path);
    if (index.sequences.empty()) {
      throw std::runtime_error("index has no reference sequence metadata");
    }
    for (const auto& layer : index.layers) {
      if (layer.document_count != index.document_count) {
        throw std::runtime_error("layer document count does not match index document count");
      }
      if (layer.dictionary.size() != layer.postings.size()) {
        throw std::runtime_error("layer dictionary and posting list counts differ");
      }
    }
    std::uint64_t distinct_seeds = 0;
    for (const auto& layer : index.layers) {
      distinct_seeds += layer.dictionary.size();
    }
    std::cout << "valid\tcomparison_stage=exact_retrieval\tformat=KmerExactIndex"
              << "\tlayers=" << index.layers.size()
              << "\tsequences=" << index.sequences.size()
              << "\tdistinct_seeds_total=" << distinct_seeds << '\n';
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "aims_validate: " << ex.what() << '\n';
    return 1;
  }
}
