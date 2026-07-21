#include "aims/codecs/posting_block_codec.hpp"

#include <stdexcept>

namespace aims::codecs {
namespace {

// Write an unsigned integer using little-endian base-128 variable-length encoding.
void append_varint(std::vector<std::uint8_t>& out, std::uint64_t value) {
  // Continue while more than seven payload bits remain.
  while (value >= 0x80U) {
    // Store the low seven bits and set the continuation bit.
    out.push_back(static_cast<std::uint8_t>((value & 0x7fU) | 0x80U));
    // Drop the seven bits that were just written.
    value >>= 7U;
  }
  // Store the final byte without the continuation bit.
  out.push_back(static_cast<std::uint8_t>(value));
}

// Read one unsigned varint and advance offset to the next encoded field.
std::uint64_t read_varint(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  // Accumulate decoded bits here.
  std::uint64_t value = 0;
  // Track how far the next seven payload bits should be shifted.
  std::uint32_t shift = 0;
  // Read until a byte without the continuation bit is found.
  while (offset < bytes.size()) {
    // Consume the next byte from the compressed block.
    const auto byte = bytes[offset++];
    // Merge the seven payload bits into the decoded integer.
    value |= static_cast<std::uint64_t>(byte & 0x7fU) << shift;
    // A clear continuation bit means this integer is complete.
    if ((byte & 0x80U) == 0) {
      return value;
    }
    // The next byte contributes seven higher-order bits.
    shift += 7U;
    // More than 64 shifted bits means the block is malformed.
    if (shift >= 64) {
      throw std::runtime_error("invalid posting block varint");
    }
  }
  // Reaching the end before a terminating byte means the block is truncated.
  throw std::runtime_error("truncated posting block");
}

} // namespace

std::vector<std::uint8_t>
PostingBlockCodec::encode(std::span<const index::Posting> postings) const {
  // Encoded bytes are appended to this buffer.
  std::vector<std::uint8_t> out;
  // Reserve a practical estimate to avoid repeated reallocations for typical small deltas.
  out.reserve(postings.size() * 8U);

  // Previous document ID is needed for delta encoding.
  std::uint32_t previous_doc = 0;
  // Previous sequence ID is needed for sequence-local delta encoding.
  std::uint64_t previous_seq = 0;
  // Previous position is needed for coordinate delta encoding.
  std::uint64_t previous_pos = 0;
  // The first posting is stored as absolute coordinates.
  bool first = true;

  // Encode every sorted posting in order.
  for (const auto& posting : postings) {
    // Delta encoding assumes nondecreasing document identifiers.
    if (!first && posting.document_id < previous_doc) {
      throw std::invalid_argument("posting block must be sorted by document");
    }

    // A same-document posting can delta encode the sequence ID.
    const bool same_doc = !first && posting.document_id == previous_doc;
    // Sequence IDs must be sorted within a document.
    if (same_doc && posting.sequence_id < previous_seq) {
      throw std::invalid_argument("posting block must be sorted by sequence within document");
    }
    // A same-sequence posting can delta encode the position.
    const bool same_seq = same_doc && posting.sequence_id == previous_seq;
    // Positions must be sorted within a sequence.
    if (same_seq && posting.position < previous_pos) {
      throw std::invalid_argument("posting block must be sorted by position within sequence");
    }

    // Store document as an absolute value for the first posting, otherwise as a delta.
    append_varint(out, first ? posting.document_id : posting.document_id - previous_doc);
    // Store sequence as a delta only when the document did not change.
    append_varint(out, same_doc ? posting.sequence_id - previous_seq : posting.sequence_id);
    // Store position as a delta only when the sequence did not change.
    append_varint(out, same_seq ? posting.position - previous_pos : posting.position);
    // Store strand as a compact numeric field.
    append_varint(out, static_cast<std::uint8_t>(posting.strand));

    // Remember this posting so the next one can be encoded as deltas.
    previous_doc = posting.document_id;
    previous_seq = posting.sequence_id;
    previous_pos = posting.position;
    // All following records are no longer the first posting.
    first = false;
  }
  // Return the compressed posting block.
  return out;
}

std::vector<index::Posting>
PostingBlockCodec::decode(std::span<const std::uint8_t> bytes,
                          std::uint64_t expected_count) const {
  // Materialized decoded postings are stored here.
  std::vector<index::Posting> postings;
  // The expected count is known from dictionary metadata.
  postings.reserve(static_cast<std::size_t>(expected_count));
  // Reuse the streaming decoder so materialized and streaming behavior stay identical.
  visit(bytes, expected_count, [&](const index::Posting& posting) {
    // Copy each decoded posting into the output vector.
    postings.push_back(posting);
  });
  // Return all decoded postings.
  return postings;
}

void PostingBlockCodec::visit(std::span<const std::uint8_t> bytes,
                              std::uint64_t expected_count,
                              const std::function<void(const index::Posting&)>& visitor) const {

  // Offset tracks the next unread byte inside the compressed block.
  std::size_t offset = 0;
  // Previous document ID is needed to reconstruct deltas.
  std::uint32_t previous_doc = 0;
  // Previous sequence ID is needed to reconstruct sequence deltas.
  std::uint64_t previous_seq = 0;
  // Previous position is needed to reconstruct position deltas.
  std::uint64_t previous_pos = 0;

  // Decode exactly the posting count declared by metadata.
  for (std::uint64_t i = 0; i < expected_count; ++i) {
    // Decode the document absolute value or delta.
    const auto doc_delta = read_varint(bytes, offset);
    // Decode the sequence absolute value or delta.
    const auto seq_delta = read_varint(bytes, offset);
    // Decode the position absolute value or delta.
    const auto pos_delta = read_varint(bytes, offset);
    // Decode the strand value.
    const auto strand = read_varint(bytes, offset);

    // Reconstruct the document ID from the first absolute value or later delta.
    const auto document_id = i == 0 ? static_cast<DocumentId>(doc_delta)
                                    : static_cast<DocumentId>(previous_doc + doc_delta);
    // Reconstruct sequence ID as a delta only if the document stayed the same.
    const auto sequence_id =
        i != 0 && document_id == previous_doc ? previous_seq + seq_delta : seq_delta;
    // Reconstruct position as a delta only if the sequence stayed the same.
    const auto position =
        i != 0 && document_id == previous_doc && sequence_id == previous_seq
            ? previous_pos + pos_delta
            : pos_delta;
    // Strand has only two legal values in this DNA index.
    if (strand > 1) {
      throw std::runtime_error("invalid strand in posting block");
    }

    // Build a public posting object for the caller.
    const auto posting = index::Posting{
        .document_id = document_id,
        .sequence_id = sequence_id,
        .position = position,
        .strand = static_cast<Strand>(strand),
    };
    // Stream the posting to the caller without requiring full materialization.
    visitor(posting);
    // Save reconstructed values for the next delta.
    previous_doc = document_id;
    previous_seq = sequence_id;
    previous_pos = position;
  }

  // All bytes must be consumed; otherwise the block does not match its metadata.
  if (offset != bytes.size()) {
    throw std::runtime_error("trailing bytes in posting block");
  }
}

} // namespace aims::codecs
