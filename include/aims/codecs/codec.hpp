#pragma once

#include <cstdint>
#include <span>
#include <vector>

namespace aims::codecs {

enum class CodecKind : std::uint8_t {
  Plain = 0,
  Delta = 1,
  EliasFano = 2,
  BlockPacked = 3,
  RoaringDocSet = 4,
};

struct DecodeMetrics {
  std::uint64_t bytes_read{0};
  std::uint64_t values_decoded{0};
};

class IntegerListCodec {
public:
  virtual ~IntegerListCodec() = default;
  [[nodiscard]] virtual CodecKind kind() const noexcept = 0;
  [[nodiscard]] virtual std::vector<std::uint8_t>
  encode(std::span<const std::uint64_t> values) const = 0;
  [[nodiscard]] virtual std::vector<std::uint64_t>
  decode(std::span<const std::uint8_t> bytes, DecodeMetrics& metrics) const = 0;
};

} // namespace aims::codecs
