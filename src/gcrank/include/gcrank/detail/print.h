#pragma once
//
// Minimal std::print replacement built on std::format + std::fwrite.
//
// std::print is C++23 and needs a very recent standard library (libstdc++ 14+
// or recent MSVC); std::format (C++20) plus fwrite keeps the toolchain
// requirement lower without changing any output.
//
#include <cstdio>
#include <format>
#include <string>
#include <utility>

namespace gcrank::detail {

// Writes the formatted text to `out` (no trailing newline, no flush).
template <typename... Args>
void print(std::FILE* out, std::format_string<Args...> fmt, Args&&... args) {
  const std::string text = std::format(fmt, std::forward<Args>(args)...);
  std::fwrite(text.data(), 1, text.size(), out);
}

// Same, to stdout.
template <typename... Args>
void print(std::format_string<Args...> fmt, Args&&... args) {
  print(stdout, fmt, std::forward<Args>(args)...);
}

}  // namespace gcrank::detail
