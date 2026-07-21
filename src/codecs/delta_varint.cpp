#include "aims/codecs/delta_varint.hpp"

#include <stdexcept>

namespace aims::codecs {
namespace {

void append_varint(std::vector<std::uint8_t>& out, std::uint64_t value) {
  while (value >= 0x80U) {
    out.push_back(static_cast<std::uint8_t>((value & 0x7fU) | 0x80U));
    value >>= 7U;
  }
  out.push_back(static_cast<std::uint8_t>(value));
}

std::uint64_t read_varint(std::span<const std::uint8_t> bytes, std::size_t& offset) {
  std::uint64_t value = 0;
  std::uint32_t shift = 0;
  while (offset < bytes.size()) {
    const auto byte = bytes[offset++];
    value |= static_cast<std::uint64_t>(byte & 0x7fU) << shift;
    if ((byte & 0x80U) == 0) {
      return value;
    }
    shift += 7U;
    if (shift >= 64) {
      throw std::runtime_error("invalid varint in delta-varint stream");
    }
  }
  throw std::runtime_error("truncated delta-varint stream");
}

} // namespace

CodecKind DeltaVarintCodec::kind() const noexcept {
  return CodecKind::Delta;
}

std::vector<std::uint8_t> DeltaVarintCodec::encode(std::span<const std::uint64_t> values) const {
  std::vector<std::uint8_t> out;
  append_varint(out, values.size());
  std::uint64_t previous = 0;
  bool first = true;
  for (const auto value : values) {
    if (!first && value < previous) {
      throw std::invalid_argument("DeltaVarintCodec requires monotone nondecreasing values");
    }
    append_varint(out, first ? value : value - previous);
    previous = value;
    first = false;
  }
  return out;
}

std::vector<std::uint64_t> DeltaVarintCodec::decode(std::span<const std::uint8_t> bytes,
                                                    DecodeMetrics& metrics) const {
  std::size_t offset = 0;
  const auto count = read_varint(bytes, offset);
  std::vector<std::uint64_t> values;
  values.reserve(static_cast<std::size_t>(count));
  std::uint64_t previous = 0;
  for (std::uint64_t i = 0; i < count; ++i) {
    const auto delta = read_varint(bytes, offset);
    const auto value = i == 0 ? delta : previous + delta;
    values.push_back(value);
    previous = value;
  }
  if (offset != bytes.size()) {
    throw std::runtime_error("trailing bytes in delta-varint stream");
  }
  metrics.bytes_read += bytes.size();
  metrics.values_decoded += values.size();
  return values;
}

} // namespace aims::codecs
