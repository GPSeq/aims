#include "aims/io/mapped_file.hpp"

#include <cerrno>
#include <cstring>
#include <stdexcept>
#include <string>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace aims::io {

MappedFile::MappedFile(const std::filesystem::path& path) {
  fd_ = ::open(path.c_str(), O_RDONLY);
  if (fd_ < 0) {
    throw std::runtime_error("failed to open mapped file: " + path.string() +
                             ": " + std::strerror(errno));
  }

  struct stat st {};
  if (::fstat(fd_, &st) != 0) {
    close();
    throw std::runtime_error("failed to stat mapped file: " + path.string() +
                             ": " + std::strerror(errno));
  }
  if (st.st_size <= 0) {
    close();
    throw std::runtime_error("cannot mmap empty file: " + path.string());
  }
  size_ = static_cast<std::uint64_t>(st.st_size);

  void* mapped = ::mmap(nullptr, static_cast<std::size_t>(size_), PROT_READ,
                        MAP_PRIVATE, fd_, 0);
  if (mapped == MAP_FAILED) {
    close();
    throw std::runtime_error("failed to mmap file: " + path.string() +
                             ": " + std::strerror(errno));
  }
  data_ = static_cast<const std::uint8_t*>(mapped);
}

MappedFile::~MappedFile() {
  close();
}

MappedFile::MappedFile(MappedFile&& other) noexcept
    : fd_(other.fd_), data_(other.data_), size_(other.size_) {
  other.fd_ = -1;
  other.data_ = nullptr;
  other.size_ = 0;
}

MappedFile& MappedFile::operator=(MappedFile&& other) noexcept {
  if (this != &other) {
    close();
    fd_ = other.fd_;
    data_ = other.data_;
    size_ = other.size_;
    other.fd_ = -1;
    other.data_ = nullptr;
    other.size_ = 0;
  }
  return *this;
}

std::span<const std::uint8_t> MappedFile::bytes() const noexcept {
  return std::span<const std::uint8_t>(data_, static_cast<std::size_t>(size_));
}

std::uint64_t MappedFile::size() const noexcept {
  return size_;
}

void MappedFile::close() noexcept {
  if (data_ != nullptr) {
    ::munmap(const_cast<std::uint8_t*>(data_), static_cast<std::size_t>(size_));
    data_ = nullptr;
  }
  if (fd_ >= 0) {
    ::close(fd_);
    fd_ = -1;
  }
  size_ = 0;
}

} // namespace aims::io
