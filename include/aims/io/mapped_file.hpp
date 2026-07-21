#pragma once

#include <cstdint>
#include <filesystem>
#include <span>

namespace aims::io {

class MappedFile {
public:
  explicit MappedFile(const std::filesystem::path& path);
  ~MappedFile();

  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;
  MappedFile(MappedFile&& other) noexcept;
  MappedFile& operator=(MappedFile&& other) noexcept;

  [[nodiscard]] std::span<const std::uint8_t> bytes() const noexcept;
  [[nodiscard]] std::uint64_t size() const noexcept;

private:
  void close() noexcept;

  int fd_{-1};
  const std::uint8_t* data_{nullptr};
  std::uint64_t size_{0};
};

} // namespace aims::io
