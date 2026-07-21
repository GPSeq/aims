#include <cassert>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <sstream>
#include <vector>

#include "aims/build/index_builder.hpp"
#include "aims/codecs/delta_varint.hpp"
#include "aims/codecs/posting_block_codec.hpp"
#include "aims/instrumentation/metrics.hpp"
#include "aims/query/candidate_accumulator.hpp"
#include "aims/query/kmer_search.hpp"
#include "aims/query/query_planner.hpp"
#include "aims/seeds/dna.hpp"
#include "aims/seeds/kmer_generator.hpp"
#include "aims/serialization/index_format.hpp"

namespace {

void test_dna_encoding_and_canonicalization() {
  using namespace aims::seeds;
  assert(encode_base('A') == 0);
  assert(encode_base('c') == 1);
  assert(!encode_base('N').has_value());
  assert(reverse_complement("ACGTN") == "NACGT");

  const auto fwd = forward_kmer("ACG");
  const auto rc = reverse_complement_kmer("CGT");
  assert(fwd == rc);
  assert(canonical_kmer("CGT").value == fwd);
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
  assert(seeds.size() == 2);
  assert(seeds[0].position == 0);
  assert(seeds[1].position == 4);
}

void test_fastq_reader() {
  const auto path = std::filesystem::temp_directory_path() / "aims_fastq_reader.fastq";
  {
    std::ofstream out(path);
    out << "@q0\nAACCGGTT\n+\nFFFFFFFF\n@q1\nTTTT\n+\nIIII\n";
  }
  const auto records = aims::io::read_sequences(path);
  std::filesystem::remove(path);
  assert(records.size() == 2);
  assert(records[0].name == "q0");
  assert(records[0].sequence == "AACCGGTT");
  assert(records[1].sequence_id == 1);
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
  assert(id.has_value());
  assert(index.dictionary.metadata(*id).document_frequency == 2);
  assert(index.postings.fetch(*id).size() == 3);
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
  assert(id.has_value());
  assert(index.dictionary.metadata(*id).collection_frequency == 2);
  assert(index.dictionary.metadata(*id).frequency_class == aims::FrequencyClass::Hot);
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
  assert(plan.seeds_generated == 5);
  assert(!plan.seeds.empty());
  assert(plan.seeds.front().skip_reason.empty());

  aims::instrumentation::QueryMetrics metrics;
  aims::query::CandidateAccumulator accumulator;
  for (const auto& planned : plan.seeds) {
    if (!planned.skip_reason.empty()) {
      continue;
    }
    const auto id = index.dictionary.id(planned.occurrence.key);
    assert(id.has_value());
    accumulator.add_postings(index.postings.fetch(*id), planned.information, metrics);
  }
  assert(metrics.postings_decoded > 0);
  assert(metrics.exact_bytes_read > 0);
  assert(accumulator.candidate_count() >= 2);
  const auto top = accumulator.top_k(1);
  assert(top.size() == 1);
  assert(top[0].supporting_seeds >= 1);
}

void test_metrics_serialization_schema_fields() {
  aims::instrumentation::QueryMetrics metrics;
  metrics.query_id = "q0";
  metrics.plan_label = "family=kmer;k=4";
  metrics.seeds_generated = 3;
  metrics.elapsed = std::chrono::nanoseconds{42};

  std::ostringstream header;
  aims::instrumentation::write_query_metrics_tsv_header(header);
  assert(header.str().find("postings_decoded") != std::string::npos);
  assert(header.str().find("exact_bytes_read") != std::string::npos);
  assert(header.str().find("seeds_skipped_frequency_class") != std::string::npos);
  assert(header.str().find("seeds_doc_only_hot") != std::string::npos);
  assert(header.str().find("topk_results") != std::string::npos);

  std::ostringstream json;
  aims::instrumentation::write_query_metrics_jsonl(json, metrics);
  assert(json.str().find("\"comparison_stage\":\"exact_retrieval\"") != std::string::npos);
  assert(json.str().find("\"seeds_skipped_frequency_class\":") != std::string::npos);
  assert(json.str().find("\"seeds_doc_only_hot\":") != std::string::npos);
  assert(json.str().find("\"per_k_metrics\":") != std::string::npos);
  assert(json.str().find("\"topk_results\":") != std::string::npos);
  assert(json.str().find("\"elapsed_ns\":42") != std::string::npos);
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

  assert(loaded.k == index.k);
  assert(loaded.document_count == index.document_count);
  assert(loaded.dictionary.size() == index.dictionary.size());
  assert(loaded.postings.size() == index.postings.size());

  const auto key = aims::SeedKey{
      .value = aims::seeds::canonical_kmer("AACC").value,
      .family = aims::SeedFamily::Kmer,
      .k = 4,
  };
  const auto original_id = index.dictionary.id(key);
  const auto loaded_id = loaded.dictionary.id(key);
  assert(original_id.has_value());
  assert(loaded_id.has_value());
  assert(index.postings.fetch(*original_id).size() == loaded.postings.fetch(*loaded_id).size());
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

  assert(loaded.document_count == 2);
  assert(loaded.source_uri == "refs.fa");
  assert(loaded.source_checksum == "fnv1a64:test");
  assert(loaded.sequences.size() == 2);
  assert(loaded.sequences[0].name == "chrA");
  assert(loaded.layers.size() == 2);
  assert(loaded.layers[0].k == 4);
  assert(loaded.layers[1].k == 6);
  assert(loaded.layers[0].dictionary.size() > 0);
  assert(loaded.layers[1].dictionary.size() > 0);
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
  assert(metrics.seeds_queried > 0);
  assert(!metrics.topk_results.empty());
  assert(metrics.exact_bytes_read > 0);
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

  assert(index.dictionary.size() == naive_counts.size());
  for (const auto& [key, count] : naive_counts) {
    const auto id = index.dictionary.id(key);
    assert(id.has_value());
    assert(index.postings.fetch(*id).size() == count);
    assert(index.dictionary.metadata(*id).collection_frequency == count);
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
    assert(file);
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
  assert(rejected);
}

void test_delta_varint_codec_round_trip() {
  const std::vector<std::uint64_t> values = {3, 5, 5, 20, 1024, 65536};
  aims::codecs::DeltaVarintCodec codec;
  const auto encoded = codec.encode(values);
  aims::codecs::DecodeMetrics metrics;
  const auto decoded = codec.decode(encoded, metrics);
  assert(decoded == values);
  assert(metrics.bytes_read == encoded.size());
  assert(metrics.values_decoded == values.size());
  assert(encoded.size() < values.size() * sizeof(std::uint64_t));
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
  assert(decoded == postings);
  assert(encoded.size() < postings.size() * sizeof(aims::index::Posting));
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
  assert(hot_metrics.seeds_skipped_hot > 0);
  assert(hot_metrics.per_k.front().seeds_skipped_hot > 0);

  const auto doc_only_metrics = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{
          .topk = 3,
          .hot_seed_threshold = 1,
          .hot_seed_mode = aims::query::HotSeedMode::DocOnly,
      });
  assert(doc_only_metrics.seeds_skipped_hot == 0);
  assert(doc_only_metrics.seeds_queried > 0);
  assert(!doc_only_metrics.topk_results.empty());

  const auto budget_metrics = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{.topk = 3, .max_seeds = 1});
  assert(budget_metrics.seeds_queried == 1);
  assert(budget_metrics.seeds_skipped_budget > 0);
  assert(budget_metrics.early_stop_fired);
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
  assert(skipped.seeds_skipped_frequency_class > 0);
  assert(skipped.seeds_skipped_hot >= skipped.seeds_skipped_frequency_class);

  const auto doc_only = aims::query::search_kmer_exact(
      index, query, aims::query::KmerSearchOptions{
                        .topk = 3,
                        .use_frequency_class_hot_policy = true,
                        .hot_seed_min_class = aims::FrequencyClass::Hot,
                        .hot_seed_mode = aims::query::HotSeedMode::DocOnly,
                    });
  assert(doc_only.seeds_doc_only_hot > 0);
  assert(doc_only.seeds_skipped_frequency_class == 0);
  assert(doc_only.seeds_queried > 0);
}

} // namespace

int main() {
  test_dna_encoding_and_canonicalization();
  test_fastq_reader();
  test_kmer_generation_splits_ambiguous_runs();
  test_fixed_k_index_and_lookup();
  test_configurable_frequency_classes();
  test_planner_and_accumulator_metrics();
  test_metrics_serialization_schema_fields();
  test_fixed_k_index_binary_round_trip();
  test_multi_k_index_binary_round_trip();
  test_mmap_index_reader_with_cache();
  test_randomized_lookup_against_naive_map();
  test_corrupt_index_checksum_is_rejected();
  test_delta_varint_codec_round_trip();
  test_posting_block_codec_round_trip();
  test_kmer_search_budgets_and_hot_seed_skip();
  test_kmer_search_frequency_class_hot_policy();
  return 0;
}
