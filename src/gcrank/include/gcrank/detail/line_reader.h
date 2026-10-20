#pragma once
//
// Buffered line reader over an InputStream.
//
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "gcrank/detail/input_stream.h"

namespace gcrank {

// Splits an input stream into lines without ever copying a line.
//
// A line ends at '\n'; a single '\r' immediately before it is stripped, so
// both LF and CRLF files are read correctly. A line that spans the internal
// buffer is returned as consecutive *pieces* instead of being materialized:
// every piece but the last reports `Status::Partial`, the last `Status::Line`
// (which may be empty for an empty line). The reported pointer stays valid
// until the next call to next().
class LineReader {
public:
  enum class Status { Partial, Line, End };

  explicit LineReader(InputStream& in, std::size_t buffer_bytes = 1u << 16);

  // Advances to the next piece of input. `*data` may be nullptr when a piece
  // is empty (an empty line, or the remainder of a line that ended exactly at
  // a buffer boundary).
  Status next(const char** data, std::size_t* len);

  // 1-based number of the line the last returned piece belongs to (the line
  // currently being read, for a Partial piece).
  uint64_t line_number() const { return line_number_; }

private:
  bool fill();

  InputStream& in_;
  std::vector<char> buf_;
  std::size_t begin_ = 0;
  std::size_t end_ = 0;
  bool eof_ = false;
  bool in_line_ = false;  // a Partial piece of the current line was returned
  uint64_t line_number_ = 0;
};

}  // namespace gcrank
