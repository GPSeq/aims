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
  out << "Usage: aims_dump --index index.aims [--mmap]\n"
      << "Prints KmerExactIndex metadata for exact_retrieval inspection.\n";
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
    std::cout << "format\tKmerExactIndex\n";
    std::cout << "comparison_stage\texact_retrieval\n";
    std::cout << "source_uri\t" << index.source_uri << '\n';
    std::cout << "source_checksum\t" << index.source_checksum << '\n';
    std::cout << "build_command\t" << index.build_command << '\n';
    std::cout << "documents\t" << index.document_count << '\n';
    std::cout << "sequences\t" << index.sequences.size() << '\n';
    for (const auto& sequence : index.sequences) {
      std::cout << "sequence\tdocument_id=" << sequence.document_id
                << "\tsequence_id=" << sequence.sequence_id
                << "\tname=" << sequence.name
                << "\tlength=" << sequence.length << '\n';
    }
    std::cout << "layers\t" << index.layers.size() << '\n';
    for (const auto& layer : index.layers) {
      std::cout << "layer\tk=" << layer.k
                << "\tdistinct_seeds=" << layer.dictionary.size()
                << "\tposting_lists=" << layer.postings.size() << '\n';
    }
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "aims_dump: " << ex.what() << '\n';
    return 1;
  }
}
