#pragma once

#include "aims/codecs/codec.hpp"

namespace aims::codecs {

class DeltaVarintCodec final : public IntegerListCodec {
public:
  [[nodiscard]] CodecKind kind() const noexcept override;
  [[nodiscard]] std::vector<std::uint8_t>
  encode(std::span<const std::uint64_t> values) const override;
  [[nodiscard]] std::vector<std::uint64_t>
  decode(std::span<const std::uint8_t> bytes, DecodeMetrics& metrics) const override;
};

} // namespace aims::codecs
