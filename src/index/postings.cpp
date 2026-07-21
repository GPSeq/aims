#include "aims/index/postings.hpp"

#include <aims/codecs/posting_block_codec.hpp>

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace aims::index {

PostingView::PostingView(std::span<const Posting> postings) noexcept : postings_(postings) {}

PostingView::PostingView(std::shared_ptr<const std::vector<Posting>> owned_postings,
                         std::uint64_t encoded_bytes) noexcept
    : owned_postings_(std::move(owned_postings)),
      postings_(*owned_postings_),
      encoded_bytes_(encoded_bytes) {}

std::uint64_t PostingView::size() const noexcept {
  return postings_.size();
}

std::uint64_t PostingView::bytes() const noexcept {
  return encoded_bytes_ == 0 ? postings_.size_bytes() : encoded_bytes_;
}

std::span<const Posting> PostingView::span() const noexcept {
  return postings_;
}

void PostingStore::build(std::vector<std::vector<Posting>> postings_by_seed) {
  codecs::PostingBlockCodec codec;
  encoded_blocks_.clear();
  posting_counts_.clear();
  block_refs_.clear();
  cache_lru_.clear();
  cache_index_.clear();
  encoded_blocks_.reserve(postings_by_seed.size());
  posting_counts_.reserve(postings_by_seed.size());
  for (auto& postings : postings_by_seed) {
    std::sort(postings.begin(), postings.end());
    posting_counts_.push_back(postings.size());
    encoded_blocks_.push_back(codec.encode(postings));
  }
  block_refs_.reserve(encoded_blocks_.size());
  for (std::size_t i = 0; i < encoded_blocks_.size(); ++i) {
    block_refs_.push_back(EncodedBlockRef{
        .data = encoded_blocks_[i].data(),
        .encoded_size = encoded_blocks_[i].size(),
        .posting_count = posting_counts_[i],
    });
  }
}

void PostingStore::build_encoded(std::vector<std::vector<std::uint8_t>> encoded_blocks,
                                 std::vector<std::uint64_t> posting_counts) {
  if (encoded_blocks.size() != posting_counts.size()) {
    throw std::invalid_argument("encoded posting blocks and counts must have equal length");
  }
  encoded_blocks_ = std::move(encoded_blocks);
  posting_counts_ = std::move(posting_counts);
  block_refs_.clear();
  cache_lru_.clear();
  cache_index_.clear();
  block_refs_.reserve(encoded_blocks_.size());
  for (std::size_t i = 0; i < encoded_blocks_.size(); ++i) {
    block_refs_.push_back(EncodedBlockRef{
        .data = encoded_blocks_[i].data(),
        .encoded_size = encoded_blocks_[i].size(),
        .posting_count = posting_counts_[i],
    });
  }
}

void PostingStore::build_encoded_refs(std::vector<std::shared_ptr<const void>> owners,
                                      std::vector<const std::uint8_t*> block_data,
                                      std::vector<std::uint64_t> encoded_sizes,
                                      std::vector<std::uint64_t> posting_counts) {
  if (owners.size() != block_data.size() || block_data.size() != encoded_sizes.size() ||
      encoded_sizes.size() != posting_counts.size()) {
    throw std::invalid_argument("encoded posting references must have equal length");
  }
  encoded_blocks_.clear();
  posting_counts_ = std::move(posting_counts);
  block_refs_.clear();
  cache_lru_.clear();
  cache_index_.clear();
  block_refs_.reserve(block_data.size());
  for (std::size_t i = 0; i < block_data.size(); ++i) {
    block_refs_.push_back(EncodedBlockRef{
        .owner = std::move(owners[i]),
        .data = block_data[i],
        .encoded_size = encoded_sizes[i],
        .posting_count = posting_counts_[i],
    });
  }
}

void PostingStore::set_cache_limit(std::uint64_t max_decoded_blocks) const {
  std::lock_guard<std::mutex> lock(*cache_mutex_);
  max_cached_blocks_ = max_decoded_blocks;
  while (cache_lru_.size() > max_cached_blocks_) {
    cache_index_.erase(cache_lru_.back().id);
    cache_lru_.pop_back();
  }
}

PostingView PostingStore::fetch(SeedId id) const {
  const auto index = static_cast<std::size_t>(id);
  assert(index < block_refs_.size());
  auto decoded = cache_get(id);
  if (decoded == nullptr) {
    decoded = decode_block(index);
    cache_put(id, decoded);
  }
  return PostingView{std::move(decoded), block_refs_.at(index).encoded_size};
}

PostingView PostingStore::fetch_doc_only(SeedId id) const {
  return fetch(id);
}

PostingStore::ReadStats PostingStore::visit_postings(
    SeedId id,
    const std::function<void(const Posting&)>& visitor) const {
  const auto index = static_cast<std::size_t>(id);
  assert(index < block_refs_.size());
  if (auto cached = cache_get(id); cached != nullptr) {
    for (const auto& posting : *cached) {
      visitor(posting);
    }
    return ReadStats{
        .postings_decoded = 0,
        .encoded_bytes_read = 0,
    };
  } else if (max_cached_blocks_ != 0) {
    auto decoded = decode_block(index);
    for (const auto& posting : *decoded) {
      visitor(posting);
    }
    cache_put(id, std::move(decoded));
  } else {
    codecs::PostingBlockCodec codec;
    codec.visit(block_bytes(index), block_refs_.at(index).posting_count, visitor);
  }
  return ReadStats{
      .postings_decoded = block_refs_.at(index).posting_count,
      .encoded_bytes_read = block_refs_.at(index).encoded_size,
  };
}

std::uint64_t PostingStore::size() const noexcept {
  return block_refs_.size();
}

std::span<const std::vector<std::uint8_t>> PostingStore::encoded_blocks() const noexcept {
  return encoded_blocks_;
}

std::span<const std::uint64_t> PostingStore::posting_counts() const noexcept {
  return posting_counts_;
}

std::span<const std::uint8_t> PostingStore::block_bytes(std::size_t index) const noexcept {
  const auto& block = block_refs_[index];
  return std::span<const std::uint8_t>(block.data, static_cast<std::size_t>(block.encoded_size));
}

std::shared_ptr<const std::vector<Posting>> PostingStore::decode_block(std::size_t index) const {
  codecs::PostingBlockCodec codec;
  return std::make_shared<const std::vector<Posting>>(
      codec.decode(block_bytes(index), block_refs_.at(index).posting_count));
}

void PostingStore::cache_put(SeedId id, std::shared_ptr<const std::vector<Posting>> postings) const {
  std::lock_guard<std::mutex> lock(*cache_mutex_);
  if (max_cached_blocks_ == 0) {
    return;
  }
  if (auto found = cache_index_.find(id); found != cache_index_.end()) {
    cache_lru_.erase(found->second);
    cache_index_.erase(found);
  }
  cache_lru_.push_front(CacheEntry{.id = id, .postings = std::move(postings)});
  cache_index_[id] = cache_lru_.begin();
  while (cache_lru_.size() > max_cached_blocks_) {
    cache_index_.erase(cache_lru_.back().id);
    cache_lru_.pop_back();
  }
}

std::shared_ptr<const std::vector<Posting>> PostingStore::cache_get(SeedId id) const {
  std::lock_guard<std::mutex> lock(*cache_mutex_);
  const auto found = cache_index_.find(id);
  if (found == cache_index_.end()) {
    return nullptr;
  }
  cache_lru_.splice(cache_lru_.begin(), cache_lru_, found->second);
  return found->second->postings;
}

} // namespace aims::index
