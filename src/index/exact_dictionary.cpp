#include "aims/index/exact_dictionary.hpp"

#include <algorithm>
#include <cassert>
#include <limits>
#include <stdexcept>

namespace aims::index {

void SortedArrayDictionary::build(std::vector<SeedKey> keys, std::vector<SeedMetadata> metadata) {
  if (keys.size() != metadata.size()) {
    throw std::invalid_argument("dictionary keys and metadata must have equal length");
  }

  if (std::is_sorted(keys.begin(), keys.end())) {
    keys_ = std::move(keys);
    metadata_ = std::move(metadata);
  } else {
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
  }
  if (std::adjacent_find(keys_.begin(), keys_.end()) != keys_.end()) {
    throw std::invalid_argument("dictionary keys must be unique");
  }
  id_slots_.clear();
  id_slot_mask_ = 0;
  if (keys_.empty()) {
    return;
  }
  if (keys_.size() > std::numeric_limits<std::size_t>::max() / 2U) {
    throw std::length_error("dictionary is too large for the lookup table");
  }
  const auto required_slots = keys_.size() * 2U;
  std::size_t slot_count = 1;
  while (slot_count < required_slots) {
    if (slot_count > std::numeric_limits<std::size_t>::max() / 2U) {
      throw std::length_error("dictionary lookup table size overflow");
    }
    slot_count *= 2U;
  }
  constexpr SeedId empty_slot = std::numeric_limits<SeedId>::max();
  id_slots_.assign(slot_count, empty_slot);
  id_slot_mask_ = slot_count - 1U;
  for (std::size_t i = 0; i < keys_.size(); ++i) {
    auto slot = std::hash<SeedKey>{}(keys_[i]) & id_slot_mask_;
    while (id_slots_[slot] != empty_slot) {
      slot = (slot + 1U) & id_slot_mask_;
    }
    id_slots_[slot] = static_cast<SeedId>(i);
  }
}

bool SortedArrayDictionary::contains(const SeedKey& key) const {
  return id(key).has_value();
}

std::optional<SeedId> SortedArrayDictionary::id(const SeedKey& key) const {
  if (id_slots_.empty()) {
    return std::nullopt;
  }
  constexpr SeedId empty_slot = std::numeric_limits<SeedId>::max();
  auto slot = std::hash<SeedKey>{}(key) & id_slot_mask_;
  while (id_slots_[slot] != empty_slot) {
    const auto id = id_slots_[slot];
    if (keys_[static_cast<std::size_t>(id)] == key) {
      return id;
    }
    slot = (slot + 1U) & id_slot_mask_;
  }
  return std::nullopt;
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
