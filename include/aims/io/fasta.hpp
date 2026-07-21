#pragma once

#include <filesystem>
#include <span>
#include <string>
#include <vector>

#include "aims/types.hpp"

namespace aims::io {

struct FastaRecord {
  std::string name;
  std::string sequence;
  DocumentId document_id{0};
  SequenceId sequence_id{0};
};

/**
 * @fn read_fasta
 * @brief Read FASTA records into memory.
 * @signature std::vector<FastaRecord> read_fasta(const std::filesystem::path& path);
 * @param path FASTA file path.
 * @throws std::runtime_error if the file cannot be opened.
 * @return Records with assigned document and sequence identifiers.
 */
[[nodiscard]] std::vector<FastaRecord> read_fasta(const std::filesystem::path& path);

/**
 * @fn read_fastq
 * @brief Read FASTQ records into memory while validating record structure.
 * @signature std::vector<FastaRecord> read_fastq(const std::filesystem::path& path);
 * @param path FASTQ file path.
 * @throws std::runtime_error if the file cannot be opened or a record is malformed.
 * @return Records with assigned document and sequence identifiers; quality strings are not retained.
 */
[[nodiscard]] std::vector<FastaRecord> read_fastq(const std::filesystem::path& path);

/**
 * @fn read_sequences
 * @brief Auto-detect FASTA or FASTQ input and read sequence records.
 * @signature std::vector<FastaRecord> read_sequences(const std::filesystem::path& path);
 * @param path FASTA or FASTQ file path.
 * @throws std::runtime_error if the format is unknown or parsing fails.
 * @return Sequence records suitable for indexing or querying.
 */
[[nodiscard]] std::vector<FastaRecord> read_sequences(const std::filesystem::path& path);
void write_fasta(const std::filesystem::path& path, std::span<const FastaRecord> records);

} // namespace aims::io
