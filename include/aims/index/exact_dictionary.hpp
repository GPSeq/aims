#pragma once

#include <optional>
#include <span>
#include <unordered_map>
#include <vector>

#include "aims/types.hpp"

namespace aims::index {

struct SeedMetadata {
  std::uint64_t document_frequency{0};
  std::uint64_t collection_frequency{0};
  std::uint64_t estimated_posting_bytes{0};
  FrequencyClass frequency_class{FrequencyClass::Rare};
};

class ExactDictionary {
public:
  virtual ~ExactDictionary() = default;
  [[nodiscard]] virtual bool contains(const SeedKey& key) const = 0;
  [[nodiscard]] virtual std::optional<SeedId> id(const SeedKey& key) const = 0;
  [[nodiscard]] virtual const SeedMetadata& metadata(SeedId id) const = 0;
  [[nodiscard]] virtual std::uint64_t size() const noexcept = 0;
};

class SortedArrayDictionary final : public ExactDictionary {
public:
  void build(std::vector<SeedKey> keys, std::vector<SeedMetadata> metadata);

  [[nodiscard]] bool contains(const SeedKey& key) const override;
  [[nodiscard]] std::optional<SeedId> id(const SeedKey& key) const override;
  [[nodiscard]] const SeedMetadata& metadata(SeedId id) const override;
  [[nodiscard]] std::uint64_t size() const noexcept override;
  [[nodiscard]] std::span<const SeedKey> keys() const noexcept;

private:
  std::vector<SeedKey> keys_{};
  std::vector<SeedMetadata> metadata_{};
  std::unordered_map<SeedKey, SeedId> id_index_{};
};

} // namespace aims::index
