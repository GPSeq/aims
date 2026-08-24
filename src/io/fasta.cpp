#include "aims/io/fasta.hpp"

#include <fstream>
#include <limits>
#include <stdexcept>

namespace aims::io {
namespace {

void strip_carriage_return(std::string& line) {
  if (!line.empty() && line.back() == '\r') {
    line.pop_back();
  }
}

FastaRecord make_record(std::string name, std::string sequence, std::uint64_t id) {
  if (id > std::numeric_limits<DocumentId>::max()) {
    throw std::runtime_error("too many sequence records for 32-bit document identifiers");
  }
  return FastaRecord{
      .name = std::move(name),
      .sequence = std::move(sequence),
      .document_id = static_cast<DocumentId>(id),
      .sequence_id = id,
  };
}

void flush_chunk(std::vector<FastaRecord>& chunk, const FastaRecordChunkVisitor& visitor) {
  if (chunk.empty()) {
    return;
  }
  visitor(chunk);
  chunk.clear();
}

void reserve_chunk_capacity(std::vector<FastaRecord>& chunk, std::size_t max_records_per_chunk) {
  constexpr std::size_t max_eager_reserve = 65536;
  if (max_records_per_chunk <= max_eager_reserve) {
    chunk.reserve(max_records_per_chunk);
  }
}

void append_record(std::vector<FastaRecord>& chunk,
                   const FastaRecordChunkVisitor& visitor,
                   std::size_t max_records_per_chunk,
                   std::string name,
                   std::string sequence,
                   std::uint64_t& next_id) {
  chunk.push_back(make_record(std::move(name), std::move(sequence), next_id));
  ++next_id;
  if (chunk.size() >= max_records_per_chunk) {
    flush_chunk(chunk, visitor);
  }
}

void read_fasta_chunks(const std::filesystem::path& path,
                       std::size_t max_records_per_chunk,
                       const FastaRecordChunkVisitor& visitor) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open FASTA: " + path.string());
  }

  std::vector<FastaRecord> chunk;
  reserve_chunk_capacity(chunk, max_records_per_chunk);
  std::uint64_t next_id = 0;
  std::string name;
  std::string sequence;
  std::string line;
  bool has_record = false;
  while (std::getline(in, line)) {
    strip_carriage_return(line);
    if (line.empty()) {
      continue;
    }
    if (line.front() == '>') {
      if (has_record) {
        append_record(chunk, visitor, max_records_per_chunk, std::move(name), std::move(sequence),
                      next_id);
        sequence.clear();
      }
      name = line.substr(1);
      has_record = true;
    } else {
      if (!has_record) {
        throw std::runtime_error("FASTA sequence data before header in: " + path.string());
      }
      sequence += line;
    }
  }
  if (has_record) {
    append_record(chunk, visitor, max_records_per_chunk, std::move(name), std::move(sequence),
                  next_id);
  }
  flush_chunk(chunk, visitor);
}

void read_fastq_chunks(const std::filesystem::path& path,
                       std::size_t max_records_per_chunk,
                       const FastaRecordChunkVisitor& visitor) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open FASTQ: " + path.string());
  }

  std::vector<FastaRecord> chunk;
  reserve_chunk_capacity(chunk, max_records_per_chunk);
  std::uint64_t next_id = 0;
  std::string name;
  std::string sequence;
  std::string plus;
  std::string quality;
  while (std::getline(in, name)) {
    strip_carriage_return(name);
    if (name.empty()) {
      continue;
    }
    if (name.front() != '@') {
      throw std::runtime_error("invalid FASTQ record header in: " + path.string());
    }
    if (!std::getline(in, sequence) || !std::getline(in, plus) || !std::getline(in, quality)) {
      throw std::runtime_error("truncated FASTQ record in: " + path.string());
    }
    strip_carriage_return(sequence);
    strip_carriage_return(plus);
    strip_carriage_return(quality);
    if (plus.empty() || plus.front() != '+') {
      throw std::runtime_error("invalid FASTQ plus line in: " + path.string());
    }
    if (quality.size() != sequence.size()) {
      throw std::runtime_error("FASTQ quality length mismatch in: " + path.string());
    }
    append_record(chunk, visitor, max_records_per_chunk, name.substr(1), std::move(sequence),
                  next_id);
  }
  flush_chunk(chunk, visitor);
}

char detect_sequence_format(const std::filesystem::path& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("failed to open sequence file: " + path.string());
  }
  char first = '\0';
  in >> first;
  return first;
}

} // namespace

std::vector<FastaRecord> read_fasta(const std::filesystem::path& path) {
  std::vector<FastaRecord> records;
  read_fasta_chunks(path, std::numeric_limits<std::size_t>::max(),
                    [&](std::span<const FastaRecord> chunk) {
                      records.insert(records.end(), chunk.begin(), chunk.end());
                    });
  return records;
}

std::vector<FastaRecord> read_fastq(const std::filesystem::path& path) {
  std::vector<FastaRecord> records;
  read_fastq_chunks(path, std::numeric_limits<std::size_t>::max(),
                    [&](std::span<const FastaRecord> chunk) {
                      records.insert(records.end(), chunk.begin(), chunk.end());
                    });
  return records;
}

std::vector<FastaRecord> read_sequences(const std::filesystem::path& path) {
  const char first = detect_sequence_format(path);
  if (first == '>') {
    return read_fasta(path);
  }
  if (first == '@') {
    return read_fastq(path);
  }
  throw std::runtime_error("unknown sequence file format: " + path.string());
}

void read_sequence_chunks(const std::filesystem::path& path,
                          std::size_t max_records_per_chunk,
                          const FastaRecordChunkVisitor& visitor) {
  if (max_records_per_chunk == 0) {
    throw std::invalid_argument("max_records_per_chunk must be greater than zero");
  }
  const char first = detect_sequence_format(path);
  if (first == '>') {
    read_fasta_chunks(path, max_records_per_chunk, visitor);
    return;
  }
  if (first == '@') {
    read_fastq_chunks(path, max_records_per_chunk, visitor);
    return;
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
