#pragma once

#include <memory>
#include <span>
#include <functional>
#include <list>
#include <mutex>
#include <unordered_map>
#include <vector>

#include "aims/types.hpp"

namespace aims::index {

struct Posting {
  DocumentId document_id{0};
  SequenceId sequence_id{0};
  Position position{0};
  Strand strand{Strand::Forward};

  friend auto operator<=>(const Posting&, const Posting&) = default;
};

class PostingView {
public:
  explicit PostingView(std::span<const Posting> postings) noexcept;
  PostingView(std::shared_ptr<const std::vector<Posting>> owned_postings,
              std::uint64_t encoded_bytes) noexcept;
  [[nodiscard]] auto begin() const noexcept { return postings_.begin(); }
  [[nodiscard]] auto end() const noexcept { return postings_.end(); }
  [[nodiscard]] std::uint64_t size() const noexcept;
  [[nodiscard]] std::uint64_t bytes() const noexcept;
  [[nodiscard]] std::span<const Posting> span() const noexcept;

private:
  std::shared_ptr<const std::vector<Posting>> owned_postings_{};
  std::span<const Posting> postings_;
  std::uint64_t encoded_bytes_{0};
};

class PostingStore {
public:
  struct ReadStats {
    std::uint64_t postings_decoded{0};
    std::uint64_t encoded_bytes_read{0};
  };

  void build(std::vector<std::vector<Posting>> postings_by_seed);
  void build_encoded(std::vector<std::vector<std::uint8_t>> encoded_blocks,
                     std::vector<std::uint64_t> posting_counts);

  /**
   * @fn build_encoded_refs
   * @brief Initialize the posting store from externally owned compressed block references.
   * @signature void build_encoded_refs(std::vector<std::shared_ptr<const void>> owners, std::vector<const std::uint8_t*> block_data, std::vector<std::uint64_t> encoded_sizes, std::vector<std::uint64_t> posting_counts);
   * @param owners Shared lifetime owners for mapped or externally stored block memory.
   * @param block_data Pointers to compressed posting block bytes.
   * @param encoded_sizes Encoded byte length for each posting block.
   * @param posting_counts Logical posting count for each block.
   * @throws std::invalid_argument if vector sizes differ.
   * @return None.
   */
  void build_encoded_refs(std::vector<std::shared_ptr<const void>> owners,
                          std::vector<const std::uint8_t*> block_data,
                          std::vector<std::uint64_t> encoded_sizes,
                          std::vector<std::uint64_t> posting_counts);
  void set_cache_limit(std::uint64_t max_decoded_blocks) const;
  [[nodiscard]] PostingView fetch(SeedId id) const;
  [[nodiscard]] PostingView fetch_doc_only(SeedId id) const;

  /**
   * @fn visit_postings
   * @brief Stream postings from one compressed block into a visitor callback.
   * @signature ReadStats visit_postings(SeedId id, const std::function<void(const Posting&)>& visitor) const;
   * @param id Seed identifier whose posting block should be decoded.
   * @param visitor Callback invoked once per decoded posting.
   * @throws std::runtime_error if the compressed posting block is malformed.
   * @return ReadStats with logical postings decoded and compressed bytes read.
   */
  [[nodiscard]] ReadStats visit_postings(
      SeedId id,
      const std::function<void(const Posting&)>& visitor) const;
  [[nodiscard]] std::uint64_t size() const noexcept;
  [[nodiscard]] std::span<const std::vector<std::uint8_t>> encoded_blocks() const noexcept;
  [[nodiscard]] std::span<const std::uint64_t> posting_counts() const noexcept;

private:
  struct EncodedBlockRef {
    std::shared_ptr<const void> owner{};
    const std::uint8_t* data{nullptr};
    std::uint64_t encoded_size{0};
    std::uint64_t posting_count{0};
  };

  struct CacheEntry {
    SeedId id{0};
    std::shared_ptr<const std::vector<Posting>> postings;
  };

  [[nodiscard]] std::span<const std::uint8_t> block_bytes(std::size_t index) const noexcept;
  [[nodiscard]] std::shared_ptr<const std::vector<Posting>> decode_block(std::size_t index) const;
  void cache_put(SeedId id, std::shared_ptr<const std::vector<Posting>> postings) const;
  [[nodiscard]] std::shared_ptr<const std::vector<Posting>> cache_get(SeedId id) const;

  std::vector<std::vector<std::uint8_t>> encoded_blocks_{};
  std::vector<std::uint64_t> posting_counts_{};
  std::vector<EncodedBlockRef> block_refs_{};
  mutable std::uint64_t max_cached_blocks_{0};
  mutable std::shared_ptr<std::mutex> cache_mutex_{std::make_shared<std::mutex>()};
  mutable std::list<CacheEntry> cache_lru_{};
  mutable std::unordered_map<SeedId, std::list<CacheEntry>::iterator> cache_index_{};
};

} // namespace aims::index
