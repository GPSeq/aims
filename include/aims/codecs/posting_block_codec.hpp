#pragma once

#include <cstdint>
#include <functional>
#include <span>
#include <vector>

#include "aims/index/postings.hpp"

namespace aims::codecs {

class PostingBlockCodec {
public:
  /**
   * @fn encode
   * @brief Compress a sorted coordinate-aware posting block.
   * @signature std::vector<std::uint8_t> encode(std::span<const index::Posting> postings) const;
   * @param postings Sorted postings for one seed identifier.
   * @throws std::invalid_argument if postings are not sorted by document, sequence, and position.
   * @return Encoded posting block bytes.
   */
  [[nodiscard]] std::vector<std::uint8_t>
  encode(std::span<const index::Posting> postings) const;

  /**
   * @fn decode
   * @brief Decode a compressed posting block into an owned vector.
   * @signature std::vector<index::Posting> decode(std::span<const std::uint8_t> bytes, std::uint64_t expected_count) const;
   * @param bytes Encoded posting block bytes.
   * @param expected_count Number of postings expected in the block.
   * @throws std::runtime_error if the block is truncated or malformed.
   * @return Decoded postings.
   */
  [[nodiscard]] std::vector<index::Posting>
  decode(std::span<const std::uint8_t> bytes, std::uint64_t expected_count) const;

  /**
   * @fn visit
   * @brief Decode a compressed posting block and stream postings into a visitor callback.
   * @signature void visit(std::span<const std::uint8_t> bytes, std::uint64_t expected_count, const std::function<void(const index::Posting&)>& visitor) const;
   * @param bytes Encoded posting block bytes.
   * @param expected_count Number of postings expected in the block.
   * @param visitor Callback invoked once per decoded posting.
   * @throws std::runtime_error if the block is truncated or malformed.
   * @return None.
   */
  void visit(std::span<const std::uint8_t> bytes,
             std::uint64_t expected_count,
             const std::function<void(const index::Posting&)>& visitor) const;
};

} // namespace aims::codecs
