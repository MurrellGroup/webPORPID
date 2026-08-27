#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <span>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

namespace webporpid::binary {

class Reader {
 public:
  explicit Reader(std::span<const std::uint8_t> source) : source_(source) {}

  template <class T> bool number(T& value) {
    static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>);
    if (offset_ + sizeof(T) > source_.size()) return false;
    std::array<std::uint8_t, sizeof(T)> raw{};
    std::memcpy(raw.data(), source_.data() + offset_, sizeof(T));
    if constexpr (std::endian::native == std::endian::big) std::reverse(raw.begin(), raw.end());
    std::memcpy(&value, raw.data(), sizeof(T));
    offset_ += sizeof(T);
    return true;
  }

  bool string(std::string& value) {
    std::uint32_t length = 0;
    if (!number(length) || offset_ + length > source_.size()) return false;
    value.assign(reinterpret_cast<const char*>(source_.data() + offset_), length);
    offset_ += length;
    return true;
  }

  bool magic(std::string_view expected) {
    if (offset_ + expected.size() > source_.size()) return false;
    const auto actual = std::string_view(reinterpret_cast<const char*>(source_.data() + offset_), expected.size());
    if (actual != expected) return false;
    offset_ += expected.size();
    return true;
  }

  std::span<const std::uint8_t> bytes(std::size_t length) {
    if (offset_ + length > source_.size()) return {};
    auto result = source_.subspan(offset_, length);
    offset_ += length;
    return result;
  }

  [[nodiscard]] bool done() const { return offset_ == source_.size(); }
  [[nodiscard]] std::size_t remaining() const { return source_.size() - offset_; }

 private:
  std::span<const std::uint8_t> source_;
  std::size_t offset_ = 0;
};

template <class T> void number(std::vector<std::uint8_t>& output, T value) {
  static_assert(std::is_integral_v<T> || std::is_floating_point_v<T>);
  std::array<std::uint8_t, sizeof(T)> raw{};
  std::memcpy(raw.data(), &value, sizeof(T));
  if constexpr (std::endian::native == std::endian::big) std::reverse(raw.begin(), raw.end());
  output.insert(output.end(), raw.begin(), raw.end());
}

inline void string(std::vector<std::uint8_t>& output, std::string_view value) {
  number(output, static_cast<std::uint32_t>(value.size()));
  output.insert(output.end(), value.begin(), value.end());
}

inline void magic(std::vector<std::uint8_t>& output, std::string_view value) {
  output.insert(output.end(), value.begin(), value.end());
}

} // namespace webporpid::binary
