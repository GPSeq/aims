#include "aims/index/exact_dictionary.hpp"

#include <algorithm>
#include <cassert>
#include <stdexcept>

namespace aims::index {

void SortedArrayDictionary::build(std::vector<SeedKey> keys, std::vector<SeedMetadata> metadata) {
  if (keys.size() != metadata.size()) {
    throw std::invalid_argument("dictionary keys and metadata must have equal length");
  }

  std::vector<std::size_t> order(keys.size());
  for (std::size_t i = 0; i < order.size(); ++i) {
    order[i] = i;
  }
  std::sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return keys[lhs] < keys[rhs];
  });

  std::vector<SeedKey> sorted_keys;
  std::vector<SeedMetadata> sorted_metadata;
  sorted_keys.reserve(keys.size());
  sorted_metadata.reserve(metadata.size());
  for (const std::size_t idx : order) {
    sorted_keys.push_back(keys[idx]);
    sorted_metadata.push_back(metadata[idx]);
  }
  keys_ = std::move(sorted_keys);
  metadata_ = std::move(sorted_metadata);
  id_index_.clear();
  id_index_.reserve(keys_.size());
  for (std::size_t i = 0; i < keys_.size(); ++i) {
    id_index_.emplace(keys_[i], static_cast<SeedId>(i));
  }
}

bool SortedArrayDictionary::contains(const SeedKey& key) const {
  return id(key).has_value();
}

std::optional<SeedId> SortedArrayDictionary::id(const SeedKey& key) const {
  const auto found = id_index_.find(key);
  if (found != id_index_.end()) {
    return found->second;
  }
  const auto it = std::lower_bound(keys_.begin(), keys_.end(), key);
  if (it == keys_.end() || *it != key) {
    return std::nullopt;
  }
  return static_cast<SeedId>(std::distance(keys_.begin(), it));
}

const SeedMetadata& SortedArrayDictionary::metadata(SeedId id) const {
  assert(id < metadata_.size());
  return metadata_.at(static_cast<std::size_t>(id));
}

std::uint64_t SortedArrayDictionary::size() const noexcept {
  return keys_.size();
}

std::span<const SeedKey> SortedArrayDictionary::keys() const noexcept {
  return keys_;
}

} // namespace aims::index
