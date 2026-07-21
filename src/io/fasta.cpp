#include "aims/io/fasta.hpp"

#include <fstream>
#include <stdexcept>

namespace aims::io {

std::vector<FastaRecord> read_fasta(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open FASTA: " + path.string());
  }

  std::vector<FastaRecord> records;
  FastaRecord current;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (line.front() == '>') {
      if (!current.name.empty()) {
        current.sequence_id = records.size();
        current.document_id = static_cast<DocumentId>(records.size());
        records.push_back(std::move(current));
        current = FastaRecord{};
      }
      current.name = line.substr(1);
    } else {
      current.sequence += line;
    }
  }
  if (!current.name.empty()) {
    current.sequence_id = records.size();
    current.document_id = static_cast<DocumentId>(records.size());
    records.push_back(std::move(current));
  }
  return records;
}

std::vector<FastaRecord> read_fastq(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open FASTQ: " + path.string());
  }

  std::vector<FastaRecord> records;
  std::string name;
  std::string sequence;
  std::string plus;
  std::string quality;
  while (std::getline(in, name)) {
    if (name.empty()) {
      continue;
    }
    if (name.front() != '@') {
      throw std::runtime_error("invalid FASTQ record header in: " + path.string());
    }
    if (!std::getline(in, sequence) || !std::getline(in, plus) || !std::getline(in, quality)) {
      throw std::runtime_error("truncated FASTQ record in: " + path.string());
    }
    if (plus.empty() || plus.front() != '+') {
      throw std::runtime_error("invalid FASTQ plus line in: " + path.string());
    }
    if (quality.size() != sequence.size()) {
      throw std::runtime_error("FASTQ quality length mismatch in: " + path.string());
    }
    const auto id = records.size();
    records.push_back(FastaRecord{
        .name = name.substr(1),
        .sequence = sequence,
        .document_id = static_cast<DocumentId>(id),
        .sequence_id = id,
    });
  }
  return records;
}

std::vector<FastaRecord> read_sequences(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open sequence file: " + path.string());
  }
  char first = '\0';
  in >> first;
  if (first == '>') {
    return read_fasta(path);
  }
  if (first == '@') {
    return read_fastq(path);
  }
  throw std::runtime_error("unknown sequence file format: " + path.string());
}

void write_fasta(const std::filesystem::path& path, std::span<const FastaRecord> records) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("failed to write FASTA: " + path.string());
  }
  for (const auto& record : records) {
    out << '>' << record.name << '\n';
    out << record.sequence << '\n';
  }
}

} // namespace aims::io
