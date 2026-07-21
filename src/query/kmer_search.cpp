#include "aims/query/kmer_search.hpp"

#include <algorithm>
#include <chrono>
#include <functional>
#include <sstream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "aims/query/candidate_accumulator.hpp"
#include "aims/query/query_planner.hpp"
#include "aims/seeds/kmer_generator.hpp"

namespace aims::query {
namespace {

// A SeedTask is the unit of work after query seeds have been planned and deduplicated.
struct SeedTask {
  // The k-mer layer that owns the seed dictionary and posting store for this seed.
  std::size_t layer_index{0};
  // The exact dictionary identifier for the seed within that k layer.
  SeedId seed_id{0};
  // The accumulated information value for this seed, including repeated query occurrences.
  double information{0.0};
  // The estimated compressed byte cost used to rank information per memory access.
  double estimated_cost{1.0};
  // The final priority used for global cross-k seed ordering.
  double priority{0.0};
  // The number of repeated query occurrences represented by this one posting lookup.
  std::uint32_t occurrence_count{0};
};

// Deduplication key: the same seed ID in different k layers is not the same lookup.
struct SeedTaskKey {
  // Identifies which k layer owns the seed ID.
  std::size_t layer_index{0};
  // Identifies the seed inside the selected layer.
  SeedId seed_id{0};

  // Default equality is enough because both fields define identity.
  friend bool operator==(const SeedTaskKey&, const SeedTaskKey&) = default;
};

// Hash function for SeedTaskKey so repeated query seeds can be merged in an unordered_map.
struct SeedTaskKeyHash {
  [[nodiscard]] std::size_t operator()(const SeedTaskKey& key) const noexcept {
    // Start with the layer hash so identical seed IDs in different layers remain separate.
    std::size_t value = std::hash<std::size_t>{}(key.layer_index);
    // Mix the seed ID into the hash using a standard 64-bit combine constant.
    value ^= std::hash<SeedId>{}(key.seed_id) + 0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
    // Return the combined hash for unordered_map lookup.
    return value;
  }
};

// Used only by doc-only hot-seed mode to ensure one contribution per candidate.
struct CandidateSeenKey {
  // Document identifier for the candidate.
  DocumentId document_id{0};
  // Sequence identifier for the candidate.
  SequenceId sequence_id{0};
  // Strand is part of the candidate key because AIMS preserves strand evidence.
  Strand strand{Strand::Forward};

  // Default equality is correct for this plain identity tuple.
  friend bool operator==(const CandidateSeenKey&, const CandidateSeenKey&) = default;
};

// Hash function for the doc-only hot-seed candidate identity.
struct CandidateSeenKeyHash {
  [[nodiscard]] std::size_t operator()(const CandidateSeenKey& key) const noexcept {
    // Start with the document identifier.
    std::size_t value = std::hash<DocumentId>{}(key.document_id);
    // Mix in the sequence identifier.
    value ^= std::hash<SequenceId>{}(key.sequence_id) + 0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
    // Mix in the strand bit.
    value ^= std::hash<std::uint8_t>{}(static_cast<std::uint8_t>(key.strand)) + 0x9e3779b97f4a7c15ULL + (value << 6U) + (value >> 2U);
    // Return the combined candidate hash.
    return value;
  }
};

// Build a compact comma-separated label such as "15,19,23" for metrics output.
std::string k_values_label(const build::KmerExactIndex& index) {
  // Use a stream to avoid manual string edge cases.
  std::ostringstream out;
  // Visit every indexed k layer in stored order.
  for (std::size_t i = 0; i < index.layers.size(); ++i) {
    // Add commas only between values.
    out << (i == 0 ? "" : ",") << index.layers[i].k;
  }
  // Return the metric-ready label.
  return out.str();
}

// Convert frequency class to a stable label for plan and benchmark output.
std::string_view frequency_class_name(FrequencyClass frequency_class) {
  // Keep labels lowercase and underscore-separated for JSON/TSV friendliness.
  switch (frequency_class) {
  case FrequencyClass::Rare:
    return "rare";
  case FrequencyClass::Medium:
    return "medium";
  case FrequencyClass::Hot:
    return "hot";
  case FrequencyClass::VeryHot:
    return "very_hot";
  }
  // Defensive fallback for future enum values.
  return "unknown";
}

// Frequency classes are ordered from most selective to least selective.
bool frequency_class_at_least(FrequencyClass observed, FrequencyClass minimum) {
  // Numeric enum order is part of the public type definition.
  return static_cast<std::uint8_t>(observed) >= static_cast<std::uint8_t>(minimum);
}

// Resolve a sequence ID back to its stored name for structured top-k output.
std::string_view sequence_name_for(const build::KmerExactIndex& index, SequenceId id) {
  // Sequence metadata is small compared with postings, so a linear scan is acceptable here.
  for (const auto& sequence : index.sequences) {
    // Return the matching name once the sequence ID is found.
    if (sequence.sequence_id == id) {
      return sequence.name;
    }
  }
  // Missing metadata is not fatal for retrieval; return an empty display name.
  return {};
}

} // namespace

instrumentation::QueryMetrics search_kmer_exact(const build::KmerExactIndex& index,
                                                const io::FastaRecord& query,
                                                KmerSearchOptions options) {
  // The generator turns query DNA into canonical k-mer occurrences.
  seeds::KmerGenerator generator;
  // The planner estimates information per compressed byte for each seed.
  QueryPlanner planner;
  // The accumulator stores candidate scores keyed by document, sequence, and strand.
  CandidateAccumulator accumulator(CandidateAccumulatorOptions{
      // A zero value means no candidate cap; otherwise new candidates stop after this limit.
      .max_candidates = options.max_candidates,
  });

  // Metrics are the public contract for every query.
  auto metrics = instrumentation::QueryMetrics{};
  // Preserve the input query name in all output records.
  metrics.query_id = query.name;
  // This tool stage is exact seed retrieval, not alignment.
  metrics.stage = instrumentation::ComparisonStage::ExactRetrieval;
  // Describe the query plan in a machine-visible but compact string.
  metrics.plan_label = "family=kmer;k=" + k_values_label(index) +
                       ";ordering=global_information_per_estimated_byte";
  // Record configured seed budget in the plan label.
  if (options.max_seeds != 0) {
    metrics.plan_label += ";max_seeds=" + std::to_string(options.max_seeds);
  }
  // Record configured posting budget in the plan label.
  if (options.max_postings != 0) {
    metrics.plan_label += ";max_postings=" + std::to_string(options.max_postings);
  }
  // Record configured candidate budget in the plan label.
  if (options.max_candidates != 0) {
    metrics.plan_label += ";max_candidates=" + std::to_string(options.max_candidates);
  }
  // Record hot-seed threshold and policy when enabled.
  if (options.hot_seed_threshold != 0) {
    metrics.plan_label += ";hot_seed_threshold=" + std::to_string(options.hot_seed_threshold);
    metrics.plan_label += options.hot_seed_mode == HotSeedMode::DocOnly ? ";hot_seed_mode=doc_only"
                                                                        : ";hot_seed_mode=skip";
  }
  // Record class-based hot policy when enabled.
  if (options.use_frequency_class_hot_policy) {
    metrics.plan_label += ";hot_seed_min_class=" +
                          std::string(frequency_class_name(options.hot_seed_min_class));
    metrics.plan_label += options.hot_seed_mode == HotSeedMode::DocOnly ? ";hot_seed_mode=doc_only"
                                                                        : ";hot_seed_mode=skip";
  }
  // Preallocate per-k metric rows so later layer_index lookup is stable.
  metrics.per_k.reserve(index.layers.size());

  // Start the wall-clock timer after setup and before query work.
  const auto start = std::chrono::steady_clock::now();
  // seed_tasks will contain one globally ordered lookup per unique query seed.
  std::vector<SeedTask> seed_tasks;
  // seed_task_index maps (layer, seed_id) to the existing task for deduplication.
  std::unordered_map<SeedTaskKey, std::size_t, SeedTaskKeyHash> seed_task_index;

  // Generate and plan seeds independently for every indexed k layer.
  for (std::size_t layer_index = 0; layer_index < index.layers.size(); ++layer_index) {
    // Fetch the layer so dictionary and postings stay paired.
    const auto& layer = index.layers[layer_index];
    // Generate canonical k-mer occurrences for this query and this k value.
    auto occurrences = generator.generate(SequenceView{
        .bases = query.sequence,
        .document_id = 0,
        .sequence_id = 0,
        .name = query.name,
    }, seeds::SeedGenerationParams{
        .family = SeedFamily::Kmer,
        .k = layer.k,
        .canonical = true,
    });

    // Score seeds within this layer using dictionary metadata and estimated posting costs.
    auto plan = planner.plan(std::move(occurrences), layer.dictionary, index.document_count);
    // Initialize the metrics row for this k layer.
    instrumentation::LayerMetrics layer_metrics;
    // Store the k value so output can be interpreted without external context.
    layer_metrics.k = layer.k;
    // Count all generated seeds, including absent seeds that will not be looked up.
    layer_metrics.seeds_generated = plan.seeds_generated;
    // Add this layer's generated seed count to the query total.
    metrics.seeds_generated += plan.seeds_generated;

    // Convert planned seeds into deduplicated global lookup tasks.
    for (auto& seed : plan.seeds) {
      // Seeds with a skip reason are absent from the exact dictionary or otherwise unusable.
      if (seed.skip_reason.empty()) {
        // Look up the exact seed ID for this layer.
        const auto seed_id = layer.dictionary.id(seed.occurrence.key);
        // A missing ID should be rare after planning, but skip defensively.
        if (!seed_id.has_value()) {
          continue;
        }
        // Build the deduplication key.
        const auto key = SeedTaskKey{.layer_index = layer_index, .seed_id = *seed_id};
        // Check whether this query has already requested the same seed in this layer.
        const auto found = seed_task_index.find(key);
        // First occurrence creates a new posting lookup task.
        if (found == seed_task_index.end()) {
          // Remember where the task is stored so repeated query seeds can merge into it.
          seed_task_index.emplace(key, seed_tasks.size());
          // Store the seed information and cost from the planner.
          seed_tasks.push_back(SeedTask{
              .layer_index = layer_index,
              .seed_id = *seed_id,
              .information = seed.information,
              .estimated_cost = seed.estimated_cost,
              .priority = seed.priority,
              .occurrence_count = 1,
          });
        } else {
          // Repeated query seeds reuse one posting lookup but add evidence weight.
          auto& task = seed_tasks[found->second];
          // Accumulate information so repeated seeds contribute to scoring.
          task.information += seed.information;
          // Count repeated occurrences so candidate support reflects multiplicity.
          ++task.occurrence_count;
          // Recompute priority after changing the information value.
          task.priority = task.information / task.estimated_cost;
        }
      }
    }
    // Store this k layer's metrics row after planning.
    metrics.per_k.push_back(layer_metrics);
  }

  // Globally order seed lookups across all k values by information per estimated byte.
  std::stable_sort(seed_tasks.begin(), seed_tasks.end(), [](const auto& lhs, const auto& rhs) {
    return lhs.priority > rhs.priority;
  });

  // Execute posting lookups in global priority order.
  for (const auto& task : seed_tasks) {
    // Resolve the layer for this task.
    const auto& layer = index.layers[task.layer_index];
    // Fetch the mutable metrics row for this layer.
    auto& layer_metrics = metrics.per_k[task.layer_index];
    // Fetch metadata for budget checks and hot-seed policy decisions.
    const auto& metadata = layer.dictionary.metadata(task.seed_id);

    // Enforce a maximum number of unique seed lookups when configured.
    if (options.max_seeds != 0 && metrics.seeds_queried >= options.max_seeds) {
      // Count skipped query occurrences represented by this deduplicated task.
      metrics.seeds_skipped_budget += task.occurrence_count;
      // Count the same skip in the relevant k layer.
      layer_metrics.seeds_skipped_budget += task.occurrence_count;
      // Mark that a budget rule affected the query result.
      metrics.early_stop_fired = true;
      // Record the specific rule for benchmark transparency.
      metrics.early_stop_rule = "max_seeds_budget";
      // Do not fetch this posting block.
      continue;
    }

    // Enforce a maximum decoded-posting budget when configured.
    if (options.max_postings != 0 &&
        metrics.postings_decoded + metadata.collection_frequency > options.max_postings) {
      // Count the skipped task at query level.
      metrics.postings_skipped_budget += task.occurrence_count;
      // Count the skipped task at layer level.
      layer_metrics.seeds_skipped_budget += task.occurrence_count;
      // Mark budget-driven early stop.
      metrics.early_stop_fired = true;
      // Record which budget stopped this task.
      metrics.early_stop_rule = "max_postings_budget";
      // Do not fetch this posting block.
      continue;
    }

    // Numeric threshold policy is active only when a nonzero threshold is configured.
    const bool hot_by_threshold =
        options.hot_seed_threshold != 0 &&
        metadata.collection_frequency > options.hot_seed_threshold;
    // Class policy is active only when explicitly enabled by the caller.
    const bool hot_by_frequency_class =
        options.use_frequency_class_hot_policy &&
        frequency_class_at_least(metadata.frequency_class, options.hot_seed_min_class);
    // Apply hot-kmer policy when either visible hot condition is true.
    if (hot_by_threshold || hot_by_frequency_class) {
      // Skip mode avoids reading this high-frequency posting block.
      if (options.hot_seed_mode == HotSeedMode::Skip) {
        // Count skipped repeated query occurrences globally.
        metrics.seeds_skipped_hot += task.occurrence_count;
        // Count skipped repeated query occurrences for this k layer.
        layer_metrics.seeds_skipped_hot += task.occurrence_count;
        // Count skips specifically caused by build-time frequency classes.
        if (hot_by_frequency_class) {
          metrics.seeds_skipped_frequency_class += task.occurrence_count;
          layer_metrics.seeds_skipped_frequency_class += task.occurrence_count;
        }
        // Continue to the next seed task without touching postings.
        continue;
      }
      // Count seeds that are handled with doc-only evidence instead of positional contribution.
      metrics.seeds_doc_only_hot += task.occurrence_count;
      layer_metrics.seeds_doc_only_hot += task.occurrence_count;
      // Save counters so layer deltas can be computed after streaming decode.
      const auto before_bytes = metrics.exact_bytes_read;
      const auto before_postings = metrics.postings_decoded;
      // Doc-only mode counts each document/sequence/strand once for this hot seed.
      std::unordered_set<CandidateSeenKey, CandidateSeenKeyHash> seen;
      // Stream the compressed posting block without materializing it as a vector.
      const auto read_stats = layer.postings.visit_postings(
          task.seed_id,
          [&](const index::Posting& posting) {
            // Candidate identity for doc-only hot-seed contribution.
            const auto key = CandidateSeenKey{
                .document_id = posting.document_id,
                .sequence_id = posting.sequence_id,
                .strand = posting.strand,
            };
            // Only the first coordinate per candidate receives this hot seed's evidence.
            if (seen.insert(key).second) {
              accumulator.add_posting(posting, task.information, metrics, task.occurrence_count);
            }
          });
      // Charge compressed bytes read for this posting block.
      metrics.exact_bytes_read += read_stats.encoded_bytes_read;
      // Charge logical decoded postings for this posting block.
      metrics.postings_decoded += read_stats.postings_decoded;
      // Count one unique seed lookup.
      ++metrics.seeds_queried;
      // Count one unique seed lookup in this k layer.
      ++layer_metrics.seeds_queried;
      // Add layer-local compressed byte delta.
      layer_metrics.exact_bytes_read += metrics.exact_bytes_read - before_bytes;
      // Add layer-local decoded posting delta.
      layer_metrics.postings_decoded += metrics.postings_decoded - before_postings;
      // Hot doc-only mode has fully handled this task.
      continue;
    }

    // Save counters so layer deltas can be computed after normal streaming decode.
    const auto before_bytes = metrics.exact_bytes_read;
    const auto before_postings = metrics.postings_decoded;
    // Stream this seed's compressed posting block directly into the accumulator.
    const auto read_stats = layer.postings.visit_postings(
        task.seed_id,
        [&](const index::Posting& posting) {
          // Add every exact coordinate occurrence for non-hot or unrestricted seeds.
          accumulator.add_posting(posting, task.information, metrics, task.occurrence_count);
        });
    // Charge compressed bytes read for this posting block.
    metrics.exact_bytes_read += read_stats.encoded_bytes_read;
    // Charge logical decoded postings for this posting block.
    metrics.postings_decoded += read_stats.postings_decoded;
    // Count one unique seed lookup.
    ++metrics.seeds_queried;

    // Count the lookup in this k layer.
    ++layer_metrics.seeds_queried;
    // Add layer-local compressed byte delta.
    layer_metrics.exact_bytes_read += metrics.exact_bytes_read - before_bytes;
    // Add layer-local decoded posting delta.
    layer_metrics.postings_decoded += metrics.postings_decoded - before_postings;
  }

  // Extract the final ranked candidate set.
  const auto top = accumulator.top_k(options.topk);
  // Store the total number of candidates created during accumulation.
  metrics.candidates_created = accumulator.candidate_count();
  // Reserve output rows for structured JSON/TSV top-k results.
  metrics.topk_results.reserve(top.size());
  // Convert internal candidate scores into public metric result records.
  for (std::size_t i = 0; i < top.size(); ++i) {
    // Read the candidate chosen for this rank.
    const auto& candidate = top[i];
    // Push a stable, self-describing top-k result row.
    metrics.topk_results.push_back(instrumentation::TopKResult{
        .rank = static_cast<std::uint32_t>(i + 1),
        .document_id = candidate.document_id,
        .sequence_id = candidate.sequence_id,
        .sequence_name = std::string(sequence_name_for(index, candidate.sequence_id)),
        .strand = candidate.strand,
        .supporting_seeds = candidate.supporting_seeds,
        .score = candidate.score,
    });
  }

  // Stop timing after all candidate reporting data has been prepared.
  metrics.elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::steady_clock::now() - start);
  // Return the complete per-query metrics object to CLI and benchmark callers.
  return metrics;
}

} // namespace aims::query
