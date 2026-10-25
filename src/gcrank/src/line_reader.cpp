#include "gcrank/detail/line_reader.h"

#include <cstring>
#include <string_view>

namespace gcrank {

namespace {
// Length of `text` without a trailing '\r', so CRLF input reads like LF.
constexpr std::size_t without_cr(std::string_view text) {
  return !text.empty() && text.back() == '\r' ? text.size() - 1 : text.size();
}
}  // namespace

LineReader::LineReader(InputStream& in, std::size_t buffer_bytes) : in_(in) {
  if (buffer_bytes < 64) buffer_bytes = 64;
  buf_.resize(buffer_bytes);
}

bool LineReader::fill() {
  if (eof_) return false;
  begin_ = 0;
  end_ = in_.read(buf_.data(), buf_.size());
  if (end_ == 0) {
    eof_ = true;
    in_line_ = false;
    return false;
  }
  return true;
}

LineReader::Status LineReader::next(const char** data, std::size_t* len) {
  for (;;) {
    if (begin_ < end_) {
      const char* base = buf_.data();
      const char* nl =
          static_cast<const char*>(std::memchr(base + begin_, '\n', end_ - begin_));
      if (nl != nullptr) {
        const char* line = base + begin_;
        std::size_t line_len = static_cast<std::size_t>(nl - line);
        begin_ = static_cast<std::size_t>(nl - base) + 1;
        if (!in_line_) ++line_number_;
        in_line_ = false;
        line_len = without_cr(std::string_view(line, line_len));
        *data = line_len > 0 ? line : nullptr;
        *len = line_len;
        return Status::Line;
      }
      // No terminator in the rest of the buffer: hand out what is there as a
      // piece of the current line, with no copy. The pointer stays valid until
      // the next call; the caller must treat the line as unterminated.
      if (!in_line_) {
        ++line_number_;
        in_line_ = true;
      }
      *data = base + begin_;
      *len = end_ - begin_;
      begin_ = end_;
      return Status::Partial;
    }
    if (!fill()) return Status::End;
  }
}

}  // namespace gcrank
