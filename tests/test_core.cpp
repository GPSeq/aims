#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>
#include <vector>

#include "aims/build/index_builder.hpp"
#include "aims/codecs/delta_varint.hpp"
#include "aims/codecs/posting_block_codec.hpp"
#include "aims/instrumentation/metrics.hpp"
#include "aims/io/fasta.hpp"
#include "aims/query/candidate_accumulator.hpp"
#include "aims/query/kmer_search.hpp"
#include "aims/query/query_planner.hpp"
#include "aims/seeds/dna.hpp"
#include "aims/seeds/kmer_generator.hpp"
#include "aims/serialization/index_format.hpp"

namespace {

[[noreturn]] void check_failed(const char* expression, const char* file, int line) {
  std::cerr << file << ':' << line << ": check failed: " << expression << '\n';
  std::abort();
}

#define CHECK(condition) \
  ((condition) ? static_cast<void>(0) : check_failed(#condition, __FILE__, __LINE__))

void test_dna_encoding_and_canonicalization() {
  using namespace aims::seeds;
  CHECK(encode_base('A') == 0);
  CHECK(encode_base('c') == 1);
  CHECK(!encode_base('N').has_value());
  CHECK(reverse_complement("ACGTN") == "NACGT");

  const auto fwd = forward_kmer("ACG");
  const auto rc = reverse_complement_kmer("CGT");
  CHECK(fwd == rc);
  CHECK(canonical_kmer("CGT").value == fwd);
}

void test_kmer_generation_splits_ambiguous_runs() {
  aims::seeds::KmerGenerator generator;
  const auto seeds = generator.generate(aims::SequenceView{
      .bases = "AACNAAA",
      .document_id = 7,
      .sequence_id = 9,
      .name = "q",
  }, aims::seeds::SeedGenerationParams{
      .family = aims::SeedFamily::Kmer,
      .k = 3,
      .canonical = true,
  });
  CHECK(seeds.size() == 2);
  CHECK(seeds[0].position == 0);
  CHECK(seeds[1].position == 4);
}

void test_fastq_reader() {
  const auto path = std::filesystem::temp_directory_path() / "aims_fastq_reader.fastq";
  {
    std::ofstream out(path);
    out << "@q0\nAACCGGTT\n+\nFFFFFFFF\n@q1\nTTTT\n+\nIIII\n";
  }
  const auto records = aims::io::read_sequences(path);
  std::filesystem::remove(path);
  CHECK(records.size() == 2);
  CHECK(records[0].name == "q0");
  CHECK(records[0].sequence == "AACCGGTT");
  CHECK(records[1].sequence_id == 1);
}

void test_crlf_sequence_reader() {
  const auto path = std::filesystem::temp_directory_path() / "aims_crlf_reader.fa";
  {
    std::ofstream out(path, std::ios::binary);
    out << ">q0\r\nAACCGG\r\nTT\r\n";
  }
  const auto records = aims::io::read_sequences(path);
  std::filesystem::remove(path);
  CHECK(records.size() == 1);
  CHECK(records[0].name == "q0");
  CHECK(records[0].sequence == "AACCGGTT");
}

void test_chunked_sequence_reader_and_streaming_index_build() {
  const auto path = std::filesystem::temp_directory_path() / "aims_chunked_reader.fa";
  {
    std::ofstream out(path);
    out << ">r0\nAACCGGTT\n>r1\nTTTTAACC\n>r2\nCCCCGGGG\n";
  }

  std::vector<std::size_t> chunk_sizes;
  std::vector<std::string> names;
  aims::io::read_sequence_chunks(path, 2, [&](std::span<const aims::io::FastaRecord> chunk) {
    chunk_sizes.push_back(chunk.size());
    for (const auto& record : chunk) {
      names.push_back(record.name);
      CHECK(record.sequence_id == names.size() - 1);
      CHECK(record.document_id == names.size() - 1);
    }
  });
  CHECK((chunk_sizes == std::vector<std::size_t>{2, 1}));
  CHECK((names == std::vector<std::string>{"r0", "r1", "r2"}));

  const std::vector<std::uint16_t> k_values = {4, 5};
  const auto streamed = aims::build::build_kmer_exact_from_sequence_file(path, k_values, {}, 2);
  const auto records = aims::io::read_sequences(path);
  const auto in_memory = aims::build::build_kmer_exact(records, k_values);
  std::filesystem::remove(path);

  CHECK(streamed.document_count == in_memory.document_count);
  CHECK(streamed.sequences.size() == in_memory.sequences.size());
  CHECK(streamed.layers.size() == in_memory.layers.size());
  for (std::size_t layer_index = 0; layer_index < streamed.layers.size(); ++layer_index) {
    const auto& streamed_layer = streamed.layers[layer_index];
    const auto& in_memory_layer = in_memory.layers[layer_index];
    CHECK(streamed_layer.k == in_memory_layer.k);
    CHECK(streamed_layer.dictionary.size() == in_memory_layer.dictionary.size());
    CHECK(streamed_layer.postings.size() == in_memory_layer.postings.size());
    const auto streamed_keys = streamed_layer.dictionary.keys();
    const auto in_memory_keys = in_memory_layer.dictionary.keys();
    for (std::size_t seed_index = 0; seed_index < streamed_keys.size(); ++seed_index) {
      CHECK(streamed_keys[seed_index] == in_memory_keys[seed_index]);
      CHECK(streamed_layer.dictionary.metadata(seed_index).document_frequency ==
             in_memory_layer.dictionary.metadata(seed_index).document_frequency);
      CHECK(streamed_layer.dictionary.metadata(seed_index).collection_frequency ==
             in_memory_layer.dictionary.metadata(seed_index).collection_frequency);
      CHECK(streamed_layer.postings.fetch(seed_index).size() ==
             in_memory_layer.postings.fetch(seed_index).size());
    }
  }
}

void test_fixed_k_index_and_lookup() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AACCGGTT", .document_id = 0, .sequence_id = 0},
      {.name = "r1", .sequence = "AACCAA", .document_id = 1, .sequence_id = 1},
  };
  const auto index = aims::build::build_fixed_k_exact(records, 4);
  const auto key = aims::SeedKey{
      .value = aims::seeds::canonical_kmer("AACC").value,
      .family = aims::SeedFamily::Kmer,
      .k = 4,
  };
  const auto id = index.dictionary.id(key);
  CHECK(id.has_value());
  CHECK(index.dictionary.metadata(*id).document_frequency == 2);
  CHECK(index.postings.fetch(*id).size() == 3);
}

void test_configurable_frequency_classes() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AAAAACCCCC", .document_id = 0, .sequence_id = 0},
      {.name = "r1", .sequence = "AAAAAGGGGG", .document_id = 1, .sequence_id = 1},
  };
  const auto index = aims::build::build_fixed_k_exact(
      records, 5, aims::build::KmerBuildOptions{
                      .frequency_thresholds = aims::build::FrequencyThresholds{
                          .rare_max = 0,
                          .medium_max = 1,
                          .hot_max = 2,
                      },
                  });
  const auto key = aims::SeedKey{
      .value = aims::seeds::canonical_kmer("AAAAA").value,
      .family = aims::SeedFamily::Kmer,
      .k = 5,
  };
  const auto id = index.dictionary.id(key);
  CHECK(id.has_value());
  CHECK(index.dictionary.metadata(*id).collection_frequency == 2);
  CHECK(index.dictionary.metadata(*id).frequency_class == aims::FrequencyClass::Hot);
}

void test_planner_and_accumulator_metrics() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AACCGGTT", .document_id = 0, .sequence_id = 0},
      {.name = "r1", .sequence = "TTTTGGTT", .document_id = 1, .sequence_id = 1},
  };
  const auto index = aims::build::build_fixed_k_exact(records, 4);
  aims::seeds::KmerGenerator generator;
  auto generated = generator.generate(aims::SequenceView{
      .bases = "AACCGGTT",
      .document_id = 0,
      .sequence_id = 0,
      .name = "q",
  }, aims::seeds::SeedGenerationParams{
      .family = aims::SeedFamily::Kmer,
      .k = 4,
      .canonical = true,
  });

  aims::query::QueryPlanner planner;
  const auto plan = planner.plan(std::move(generated), index.dictionary, index.document_count);
  CHECK(plan.seeds_generated == 5);
  CHECK(!plan.seeds.empty());
  CHECK(plan.seeds.front().skip_reason.empty());

  aims::instrumentation::QueryMetrics metrics;
  aims::query::CandidateAccumulator accumulator;
  for (const auto& planned : plan.seeds) {
    if (!planned.skip_reason.empty()) {
      continue;
    }
    const auto id = index.dictionary.id(planned.occurrence.key);
    CHECK(id.has_value());
    accumulator.add_postings(index.postings.fetch(*id), planned.information, metrics);
  }
  CHECK(metrics.postings_decoded > 0);
  CHECK(metrics.exact_bytes_read > 0);
  CHECK(accumulator.candidate_count() >= 2);
  const auto top = accumulator.top_k(1);
  CHECK(top.size() == 1);
  CHECK(top[0].supporting_seeds >= 1);
}

void test_metrics_serialization_schema_fields() {
  aims::instrumentation::QueryMetrics metrics;
  metrics.query_id = "q0";
  metrics.plan_label = "family=kmer;k=4";
  metrics.seeds_generated = 3;
  metrics.elapsed = std::chrono::nanoseconds{42};

  std::ostringstream header;
  aims::instrumentation::write_query_metrics_tsv_header(header);
  CHECK(header.str().find("postings_decoded") != std::string::npos);
  CHECK(header.str().find("exact_bytes_read") != std::string::npos);
  CHECK(header.str().find("seeds_skipped_frequency_class") != std::string::npos);
  CHECK(header.str().find("seeds_doc_only_hot") != std::string::npos);
  CHECK(header.str().find("topk_results") != std::string::npos);

  std::ostringstream json;
  aims::instrumentation::write_query_metrics_jsonl(json, metrics);
  CHECK(json.str().find("\"comparison_stage\":\"exact_retrieval\"") != std::string::npos);
  CHECK(json.str().find("\"seeds_skipped_frequency_class\":") != std::string::npos);
  CHECK(json.str().find("\"seeds_doc_only_hot\":") != std::string::npos);
  CHECK(json.str().find("\"per_k_metrics\":") != std::string::npos);
  CHECK(json.str().find("\"topk_results\":") != std::string::npos);
  CHECK(json.str().find("\"elapsed_ns\":42") != std::string::npos);
}

void test_fixed_k_index_binary_round_trip() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AACCGGTT", .document_id = 0, .sequence_id = 0},
      {.name = "r1", .sequence = "TTTTAACC", .document_id = 1, .sequence_id = 1},
  };
  const auto index = aims::build::build_fixed_k_exact(records, 4);
  const auto path = std::filesystem::temp_directory_path() / "aims_fixed_k_exact_roundtrip.aims";
  aims::serialization::write_fixed_k_exact_index(path, index);
  const auto loaded = aims::serialization::read_fixed_k_exact_index(path);
  std::filesystem::remove(path);

  CHECK(loaded.k == index.k);
  CHECK(loaded.document_count == index.document_count);
  CHECK(loaded.dictionary.size() == index.dictionary.size());
  CHECK(loaded.postings.size() == index.postings.size());

  const auto key = aims::SeedKey{
      .value = aims::seeds::canonical_kmer("AACC").value,
      .family = aims::SeedFamily::Kmer,
      .k = 4,
  };
  const auto original_id = index.dictionary.id(key);
  const auto loaded_id = loaded.dictionary.id(key);
  CHECK(original_id.has_value());
  CHECK(loaded_id.has_value());
  CHECK(index.postings.fetch(*original_id).size() == loaded.postings.fetch(*loaded_id).size());
}

void test_multi_k_index_binary_round_trip() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "chrA", .sequence = "AACCGGTTAACCGGTT", .document_id = 0, .sequence_id = 0},
      {.name = "chrB", .sequence = "GGGGAAAACCCCTTTT", .document_id = 1, .sequence_id = 1},
  };
  const std::vector<std::uint16_t> k_values = {4, 6};
  const auto index = aims::build::build_kmer_exact(records, k_values);
  auto metadata_index = index;
  metadata_index.source_uri = "refs.fa";
  metadata_index.source_checksum = "fnv1a64:test";
  metadata_index.build_command = "aims_build --ref refs.fa";
  const auto path = std::filesystem::temp_directory_path() / "aims_kmer_exact_roundtrip.aims";
  aims::serialization::write_kmer_exact_index(path, metadata_index);
  const auto loaded = aims::serialization::read_kmer_exact_index(path);
  std::filesystem::remove(path);

  CHECK(loaded.document_count == 2);
  CHECK(loaded.source_uri == "refs.fa");
  CHECK(loaded.source_checksum == "fnv1a64:test");
  CHECK(loaded.sequences.size() == 2);
  CHECK(loaded.sequences[0].name == "chrA");
  CHECK(loaded.layers.size() == 2);
  CHECK(loaded.layers[0].k == 4);
  CHECK(loaded.layers[1].k == 6);
  CHECK(loaded.layers[0].dictionary.size() > 0);
  CHECK(loaded.layers[1].dictionary.size() > 0);
}

void test_mmap_index_reader_with_cache() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "chrA", .sequence = "AACCGGTTAACCGGTT", .document_id = 0, .sequence_id = 0},
      {.name = "chrB", .sequence = "GGGGAAAACCCCTTTT", .document_id = 1, .sequence_id = 1},
  };
  const std::vector<std::uint16_t> k_values = {4, 6};
  const auto index = aims::build::build_kmer_exact(records, k_values);
  const auto path = std::filesystem::temp_directory_path() / "aims_kmer_exact_mmap.aims";
  aims::serialization::write_kmer_exact_index(path, index);
  const auto loaded = aims::serialization::read_kmer_exact_index_mmap(path, 2);
  std::filesystem::remove(path);

  const auto query = aims::io::FastaRecord{
      .name = "q",
      .sequence = "AACCGGTT",
      .document_id = 0,
      .sequence_id = 0,
  };
  const auto metrics = aims::query::search_kmer_exact(
      loaded, query, aims::query::KmerSearchOptions{.topk = 2});
  CHECK(metrics.seeds_queried > 0);
  CHECK(!metrics.topk_results.empty());
  CHECK(metrics.exact_bytes_read > 0);
}

void test_posting_store_copy_rebinds_owned_blocks() {
  aims::index::PostingStore copied;
  aims::SeedId copied_seed_id = 0;
  std::uint64_t expected_count = 0;
  {
    const std::vector<aims::io::FastaRecord> records = {
        {.name = "r0", .sequence = "AACCGGTT", .document_id = 0, .sequence_id = 0},
    };
    const auto index = aims::build::build_fixed_k_exact(records, 4);
    const auto key = aims::SeedKey{
        .value = aims::seeds::canonical_kmer("AACC").value,
        .family = aims::SeedFamily::Kmer,
        .k = 4,
    };
    const auto seed_id = index.dictionary.id(key);
    CHECK(seed_id.has_value());
    copied_seed_id = *seed_id;
    expected_count = index.postings.fetch(*seed_id).size();
    index.postings.set_cache_limit(1);
    static_cast<void>(index.postings.fetch(*seed_id));
    copied = index.postings;
  }
  CHECK(copied.fetch(copied_seed_id).size() == expected_count);
}

void test_mmap_index_can_be_reserialized() {
  const std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AACCGGTTAACC", .document_id = 0, .sequence_id = 0},
      {.name = "r1", .sequence = "TTTTAACCGGTT", .document_id = 1, .sequence_id = 1},
  };
  const std::vector<std::uint16_t> k_values = {4, 6};
  const auto original = aims::build::build_kmer_exact(records, k_values);
  const auto source_path =
      std::filesystem::temp_directory_path() / "aims_mmap_reserialize_source.aims";
  const auto copy_path =
      std::filesystem::temp_directory_path() / "aims_mmap_reserialize_copy.aims";
  aims::serialization::write_kmer_exact_index(source_path, original);
  aims::build::KmerExactIndex mapped_copy;
  {
    const auto mapped = aims::serialization::read_kmer_exact_index_mmap(source_path, 1);
    mapped_copy = mapped;
  }
  std::filesystem::remove(source_path);
  aims::serialization::write_kmer_exact_index(copy_path, mapped_copy);
  const auto copied = aims::serialization::read_kmer_exact_index(copy_path);
  std::filesystem::remove(copy_path);

  CHECK(copied.layers.size() == original.layers.size());
  for (std::size_t layer = 0; layer < original.layers.size(); ++layer) {
    CHECK(copied.layers[layer].postings.size() == original.layers[layer].postings.size());
    for (aims::SeedId seed = 0; seed < original.layers[layer].postings.size(); ++seed) {
      const auto copied_view = copied.layers[layer].postings.fetch(seed);
      const auto original_view = original.layers[layer].postings.fetch(seed);
      const auto copied_postings = copied_view.span();
      const auto original_postings = original_view.span();
      CHECK(copied_postings.size() == original_postings.size());
      CHECK(std::equal(copied_postings.begin(), copied_postings.end(),
                       original_postings.begin()));
    }
  }
}

void test_search_reports_query_relative_strand() {
  const std::vector<aims::io::FastaRecord> records = {
      {.name = "forward", .sequence = "AACCGTAA", .document_id = 0, .sequence_id = 0},
      {.name = "reverse", .sequence = "TTACGGTT", .document_id = 1, .sequence_id = 1},
  };
  const std::vector<std::uint16_t> k_values = {5};
  const auto index = aims::build::build_kmer_exact(records, k_values);
  const auto query = aims::io::FastaRecord{
      .name = "q", .sequence = "CCGTA", .document_id = 0, .sequence_id = 0};
  const auto metrics = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{.topk = 2});

  CHECK(metrics.topk_results.size() == 2);
  CHECK(metrics.topk_results[0].document_id == 0);
  CHECK(metrics.topk_results[0].strand == aims::Strand::Forward);
  CHECK(metrics.topk_results[1].document_id == 1);
  CHECK(metrics.topk_results[1].strand == aims::Strand::Reverse);
}

void test_posting_budget_is_cache_state_independent() {
  const std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AACCGGTTAACCGGTT", .document_id = 0, .sequence_id = 0},
      {.name = "r1", .sequence = "TTTTAACCGGTTAAAA", .document_id = 1, .sequence_id = 1},
  };
  const std::vector<std::uint16_t> k_values = {4};
  const auto index = aims::build::build_kmer_exact(records, k_values);
  index.layers.front().postings.set_cache_limit(64);
  const auto query = aims::io::FastaRecord{
      .name = "q", .sequence = "AACCGGTT", .document_id = 0, .sequence_id = 0};
  const auto options = aims::query::KmerSearchOptions{.topk = 3, .max_postings = 10};
  const auto cold = aims::query::search_kmer_exact(index, query, options);
  const auto warm = aims::query::search_kmer_exact(index, query, options);

  CHECK(cold.seeds_queried == warm.seeds_queried);
  CHECK(cold.seeds_skipped_budget == warm.seeds_skipped_budget);
  CHECK(cold.topk_results.size() == warm.topk_results.size());
  for (std::size_t i = 0; i < cold.topk_results.size(); ++i) {
    CHECK(cold.topk_results[i].document_id == warm.topk_results[i].document_id);
    CHECK(cold.topk_results[i].sequence_id == warm.topk_results[i].sequence_id);
    CHECK(cold.topk_results[i].strand == warm.topk_results[i].strand);
    CHECK(cold.topk_results[i].supporting_seeds == warm.topk_results[i].supporting_seeds);
    CHECK(cold.topk_results[i].score == warm.topk_results[i].score);
  }
}

void test_randomized_lookup_against_naive_map() {
  std::mt19937 rng(7);
  constexpr char bases[] = {'A', 'C', 'G', 'T'};
  std::vector<aims::io::FastaRecord> records;
  for (std::uint32_t i = 0; i < 8; ++i) {
    std::string sequence;
    for (std::uint32_t j = 0; j < 40; ++j) {
      sequence.push_back(bases[rng() % 4]);
    }
    records.push_back(aims::io::FastaRecord{
        .name = "random_" + std::to_string(i),
        .sequence = sequence,
        .document_id = i,
        .sequence_id = i,
    });
  }

  constexpr std::uint16_t k = 5;
  const auto index = aims::build::build_fixed_k_exact(records, k);
  std::map<aims::SeedKey, std::uint64_t> naive_counts;
  aims::seeds::KmerGenerator generator;
  for (const auto& record : records) {
    const auto seeds = generator.generate(aims::SequenceView{
        .bases = record.sequence,
        .document_id = record.document_id,
        .sequence_id = record.sequence_id,
        .name = record.name,
    }, aims::seeds::SeedGenerationParams{
        .family = aims::SeedFamily::Kmer,
        .k = k,
        .canonical = true,
    });
    for (const auto& seed : seeds) {
      ++naive_counts[seed.key];
    }
  }

  CHECK(index.dictionary.size() == naive_counts.size());
  for (const auto& [key, count] : naive_counts) {
    const auto id = index.dictionary.id(key);
    CHECK(id.has_value());
    CHECK(index.postings.fetch(*id).size() == count);
    CHECK(index.dictionary.metadata(*id).collection_frequency == count);
  }
}

void test_corrupt_index_checksum_is_rejected() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AACCGGTT", .document_id = 0, .sequence_id = 0},
  };
  const std::vector<std::uint16_t> k_values = {4};
  const auto index = aims::build::build_kmer_exact(records, k_values);
  const auto path = std::filesystem::temp_directory_path() / "aims_corrupt_checksum.aims";
  aims::serialization::write_kmer_exact_index(path, index);

  {
    std::fstream file(path, std::ios::binary | std::ios::in | std::ios::out);
    CHECK(file);
    file.seekg(-1, std::ios::end);
    char byte = 0;
    file.read(&byte, 1);
    file.seekp(-1, std::ios::end);
    byte = static_cast<char>(byte ^ 0x7f);
    file.write(&byte, 1);
  }

  bool rejected = false;
  try {
    static_cast<void>(aims::serialization::read_kmer_exact_index(path));
  } catch (const std::exception&) {
    rejected = true;
  }
  std::filesystem::remove(path);
  CHECK(rejected);
}

void test_delta_varint_codec_round_trip() {
  const std::vector<std::uint64_t> values = {3, 5, 5, 20, 1024, 65536};
  aims::codecs::DeltaVarintCodec codec;
  const auto encoded = codec.encode(values);
  aims::codecs::DecodeMetrics metrics;
  const auto decoded = codec.decode(encoded, metrics);
  CHECK(decoded == values);
  CHECK(metrics.bytes_read == encoded.size());
  CHECK(metrics.values_decoded == values.size());
  CHECK(encoded.size() < values.size() * sizeof(std::uint64_t));
}

void test_posting_block_codec_round_trip() {
  const std::vector<aims::index::Posting> postings = {
      {.document_id = 0, .sequence_id = 0, .position = 3, .strand = aims::Strand::Forward},
      {.document_id = 0, .sequence_id = 0, .position = 8, .strand = aims::Strand::Reverse},
      {.document_id = 0, .sequence_id = 1, .position = 2, .strand = aims::Strand::Forward},
      {.document_id = 2, .sequence_id = 3, .position = 1, .strand = aims::Strand::Reverse},
  };
  aims::codecs::PostingBlockCodec codec;
  const auto encoded = codec.encode(postings);
  const auto decoded = codec.decode(encoded, postings.size());
  CHECK(decoded == postings);
  CHECK(encoded.size() < postings.size() * sizeof(aims::index::Posting));
}

void test_codecs_reject_out_of_range_varints() {
  std::vector<std::uint8_t> out_of_range(9, 0x80U);
  out_of_range.push_back(0x02U);

  bool delta_rejected = false;
  try {
    aims::codecs::DecodeMetrics metrics;
    aims::codecs::DeltaVarintCodec codec;
    static_cast<void>(codec.decode(out_of_range, metrics));
  } catch (const std::exception&) {
    delta_rejected = true;
  }
  CHECK(delta_rejected);

  bool posting_rejected = false;
  try {
    aims::codecs::PostingBlockCodec codec;
    static_cast<void>(codec.decode(out_of_range, 1));
  } catch (const std::exception&) {
    posting_rejected = true;
  }
  CHECK(posting_rejected);
}

void test_kmer_search_budgets_and_hot_seed_skip() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AAAAACCCCCAAAAACCCCC", .document_id = 0, .sequence_id = 0},
      {.name = "r1", .sequence = "AAAAACCCCC", .document_id = 1, .sequence_id = 1},
  };
  const std::vector<std::uint16_t> k_values = {5};
  const auto index = aims::build::build_kmer_exact(records, k_values);
  const auto query = aims::io::FastaRecord{
      .name = "q",
      .sequence = "AAAAACCCCC",
      .document_id = 0,
      .sequence_id = 0,
  };

  const auto hot_metrics = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{.topk = 3, .hot_seed_threshold = 1});
  CHECK(hot_metrics.seeds_skipped_hot > 0);
  CHECK(hot_metrics.per_k.front().seeds_skipped_hot > 0);

  const auto doc_only_metrics = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{
          .topk = 3,
          .hot_seed_threshold = 1,
          .hot_seed_mode = aims::query::HotSeedMode::DocOnly,
      });
  CHECK(doc_only_metrics.seeds_skipped_hot == 0);
  CHECK(doc_only_metrics.seeds_queried > 0);
  CHECK(!doc_only_metrics.topk_results.empty());

  const auto budget_metrics = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{.topk = 3, .max_seeds = 1});
  CHECK(budget_metrics.seeds_queried == 1);
  CHECK(budget_metrics.seeds_skipped_budget > 0);
  CHECK(budget_metrics.early_stop_fired);
}

void test_kmer_search_frequency_class_hot_policy() {
  std::vector<aims::io::FastaRecord> records = {
      {.name = "r0", .sequence = "AAAAACCCCCAAAAA", .document_id = 0, .sequence_id = 0},
      {.name = "r1", .sequence = "AAAAAGGGGGAAAAA", .document_id = 1, .sequence_id = 1},
  };
  const std::vector<std::uint16_t> k_values = {5};
  const auto index = aims::build::build_kmer_exact(
      records, k_values, aims::build::KmerBuildOptions{
                             .frequency_thresholds = aims::build::FrequencyThresholds{
                                 .rare_max = 0,
                                 .medium_max = 1,
                                 .hot_max = 2,
                             },
                         });
  const auto query = aims::io::FastaRecord{
      .name = "q",
      .sequence = "AAAAACCCCC",
      .document_id = 0,
      .sequence_id = 0,
  };

  const auto skipped = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{
                        .topk = 3,
                        .use_frequency_class_hot_policy = true,
                        .hot_seed_min_class = aims::FrequencyClass::Hot,
                    });
  CHECK(skipped.seeds_skipped_frequency_class > 0);
  CHECK(skipped.seeds_skipped_hot >= skipped.seeds_skipped_frequency_class);

  const auto doc_only = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{
                        .topk = 3,
                        .use_frequency_class_hot_policy = true,
                        .hot_seed_min_class = aims::FrequencyClass::Hot,
                        .hot_seed_mode = aims::query::HotSeedMode::DocOnly,
                    });
  CHECK(doc_only.seeds_doc_only_hot > 0);
  CHECK(doc_only.seeds_skipped_frequency_class == 0);
  CHECK(doc_only.seeds_queried > 0);
}

} // namespace

int main() {
  test_dna_encoding_and_canonicalization();
  test_fastq_reader();
  test_crlf_sequence_reader();
  test_chunked_sequence_reader_and_streaming_index_build();
  test_kmer_generation_splits_ambiguous_runs();
  test_fixed_k_index_and_lookup();
  test_configurable_frequency_classes();
  test_planner_and_accumulator_metrics();
  test_metrics_serialization_schema_fields();
  test_fixed_k_index_binary_round_trip();
  test_multi_k_index_binary_round_trip();
  test_mmap_index_reader_with_cache();
  test_posting_store_copy_rebinds_owned_blocks();
  test_mmap_index_can_be_reserialized();
  test_search_reports_query_relative_strand();
  test_posting_budget_is_cache_state_independent();
  test_randomized_lookup_against_naive_map();
  test_corrupt_index_checksum_is_rejected();
  test_delta_varint_codec_round_trip();
  test_posting_block_codec_round_trip();
  test_codecs_reject_out_of_range_varints();
  test_kmer_search_budgets_and_hot_seed_skip();
  test_kmer_search_frequency_class_hot_policy();
  return 0;
}
