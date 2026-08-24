#include "aims/index/postings.hpp"

#include <aims/codecs/posting_block_codec.hpp>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <stdexcept>

namespace aims::index {

PostingStore::PostingStore(const PostingStore& other) {
  std::lock_guard<std::mutex> lock(*other.cache_mutex_);
  encoded_blocks_ = other.encoded_blocks_;
  posting_counts_ = other.posting_counts_;
  max_cached_blocks_ = other.max_cached_blocks_;
  max_cached_bytes_ = other.max_cached_bytes_;

  block_refs_.reserve(other.block_refs_.size());
  const bool owns_all_blocks = other.encoded_blocks_.size() == other.block_refs_.size();
  for (std::size_t i = 0; i < other.block_refs_.size(); ++i) {
    const auto& source = other.block_refs_[i];
    block_refs_.push_back(EncodedBlockRef{
        .owner = owns_all_blocks ? std::shared_ptr<const void>{} : source.owner,
        .data = owns_all_blocks ? encoded_blocks_[i].data() : source.data,
        .encoded_size = source.encoded_size,
        .posting_count = source.posting_count,
    });
  }
}

PostingStore& PostingStore::operator=(const PostingStore& other) {
  if (this != &other) {
    PostingStore copy(other);
    *this = std::move(copy);
  }
  return *this;
}

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
  decoding_ids_.clear();
  active_cache_readers_.clear();
  cached_bytes_ = 0;
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
  decoding_ids_.clear();
  active_cache_readers_.clear();
  cached_bytes_ = 0;
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
  decoding_ids_.clear();
  active_cache_readers_.clear();
  cached_bytes_ = 0;
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
  trim_cache_locked();
}

void PostingStore::set_cache_byte_limit(std::uint64_t max_decoded_bytes) const {
  std::lock_guard<std::mutex> lock(*cache_mutex_);
  max_cached_bytes_ = max_decoded_bytes;
  trim_cache_locked();
}

void PostingStore::clear_cache() const {
  std::lock_guard<std::mutex> lock(*cache_mutex_);
  cache_lru_.clear();
  cache_index_.clear();
  cached_bytes_ = 0;
}

PostingView PostingStore::fetch(SeedId id) const {
  const auto index = static_cast<std::size_t>(id);
  assert(index < block_refs_.size());
  bool decoded_here = false;
  auto decoded = cache_get_or_decode(id, decoded_here);
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
  if (cache_enabled()) {
    bool decoded_here = false;
    auto decoded = cache_get_or_decode(id, decoded_here, true);
    try {
      for (const auto& posting : *decoded) {
        visitor(posting);
      }
    } catch (...) {
      release_cache_reader(id);
      throw;
    }
    release_cache_reader(id);
    return ReadStats{
        .postings_decoded = decoded_here ? block_refs_.at(index).posting_count : 0,
        .encoded_bytes_read = decoded_here ? block_refs_.at(index).encoded_size : 0,
    };
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

std::span<const std::uint8_t> PostingStore::encoded_block(SeedId id) const {
  const auto index = static_cast<std::size_t>(id);
  const auto& block = block_refs_.at(index);
  return std::span<const std::uint8_t>(block.data, static_cast<std::size_t>(block.encoded_size));
}

std::uint64_t PostingStore::posting_count(SeedId id) const {
  return block_refs_.at(static_cast<std::size_t>(id)).posting_count;
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

void PostingStore::cache_put_locked(
    SeedId id,
    std::shared_ptr<const std::vector<Posting>> postings) const {
  if (max_cached_blocks_ == 0 && max_cached_bytes_ == 0) {
    return;
  }
  if (auto found = cache_index_.find(id); found != cache_index_.end()) {
    cached_bytes_ -= found->second->decoded_bytes;
    cache_lru_.erase(found->second);
    cache_index_.erase(found);
  }
  const auto decoded_bytes =
      static_cast<std::uint64_t>(postings->size() * sizeof(Posting));
  // An entry that cannot fit on its own would make a byte-bounded cache
  // permanently exceed its contract. Serve it to the current caller without
  // retaining it instead.
  if (max_cached_bytes_ != 0 && decoded_bytes > max_cached_bytes_) {
    return;
  }
  cache_lru_.push_front(CacheEntry{
      .id = id,
      .postings = std::move(postings),
      .decoded_bytes = decoded_bytes,
  });
  cache_index_[id] = cache_lru_.begin();
  cached_bytes_ += decoded_bytes;
  trim_cache_locked();
}

bool PostingStore::cache_enabled() const {
  std::lock_guard<std::mutex> lock(*cache_mutex_);
  return max_cached_blocks_ != 0 || max_cached_bytes_ != 0;
}

std::shared_ptr<const std::vector<Posting>> PostingStore::cache_get_or_decode(
    SeedId id,
    bool& decoded_here,
    bool pin_reader) const {
  const auto index = static_cast<std::size_t>(id);
  std::unique_lock<std::mutex> lock(*cache_mutex_);
  while (true) {
    if (const auto found = cache_index_.find(id); found != cache_index_.end()) {
      cache_lru_.splice(cache_lru_.begin(), cache_lru_, found->second);
      auto postings = found->second->postings;
      if (pin_reader) {
        ++active_cache_readers_[id];
      }
      decoded_here = false;
      return postings;
    }
    if (max_cached_blocks_ == 0 && max_cached_bytes_ == 0) {
      lock.unlock();
      decoded_here = true;
      return decode_block(index);
    }
    if (decoding_ids_.insert(id).second) {
      break;
    }
    cache_cv_->wait(lock);
  }

  lock.unlock();
  std::shared_ptr<const std::vector<Posting>> decoded;
  try {
    decoded = decode_block(index);
  } catch (...) {
    lock.lock();
    decoding_ids_.erase(id);
    lock.unlock();
    cache_cv_->notify_all();
    throw;
  }

  lock.lock();
  cache_put_locked(id, decoded);
  if (pin_reader && cache_index_.contains(id)) {
    ++active_cache_readers_[id];
  }
  decoding_ids_.erase(id);
  lock.unlock();
  cache_cv_->notify_all();
  decoded_here = true;
  return decoded;
}

void PostingStore::release_cache_reader(SeedId id) const {
  std::lock_guard<std::mutex> lock(*cache_mutex_);
  const auto found = active_cache_readers_.find(id);
  if (found == active_cache_readers_.end()) {
    return;
  }
  if (--found->second == 0) {
    active_cache_readers_.erase(found);
  }
  trim_cache_locked();
}

void PostingStore::trim_cache_locked() const {
  const auto exceeds_limit = [&] {
    const bool disabled = max_cached_blocks_ == 0 && max_cached_bytes_ == 0;
    return disabled ||
           (max_cached_blocks_ != 0 && cache_lru_.size() > max_cached_blocks_) ||
           (max_cached_bytes_ != 0 && cached_bytes_ > max_cached_bytes_);
  };
  while (!cache_lru_.empty() && exceeds_limit()) {
    auto evict = std::prev(cache_lru_.end());
    while (active_cache_readers_.contains(evict->id)) {
      if (evict == cache_lru_.begin()) {
        return;
      }
      --evict;
    }
    cached_bytes_ -= evict->decoded_bytes;
    cache_index_.erase(evict->id);
    cache_lru_.erase(evict);
  }
}

} // namespace aims::index
