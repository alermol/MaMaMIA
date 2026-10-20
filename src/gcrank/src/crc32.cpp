#include "gcrank/detail/crc32.h"

#include <algorithm>
#include <format>
#include <fstream>
#include <vector>

#include <zlib.h>

#include "gcrank/error.h"

namespace gcrank {

namespace {

constexpr std::size_t kBufferBytes = 1u << 18;

}  // namespace

void CRC32::update(const void* data, std::size_t n) {
  // zlib's crc32() takes and returns the *final* (complemented) CRC, while
  // value_ keeps the inverted running state: convert on the way in and out.
#if ZLIB_VERNUM >= 0x1290
  value_ = static_cast<uint32_t>(::crc32_z(static_cast<uLong>(value_ ^ 0xFFFFFFFFu),
                                           static_cast<const Bytef*>(data),
                                           static_cast<z_size_t>(n))) ^
           0xFFFFFFFFu;
#else
  // crc32_z() needs zlib 1.2.9 or newer; otherwise use the classic crc32(),
  // whose length argument is 32 bits, in chunks.
  uLong crc = static_cast<uLong>(value_ ^ 0xFFFFFFFFu);
  const uint8_t* p = static_cast<const uint8_t*>(data);
  while (n > 0) {
    const uInt chunk = static_cast<uInt>(std::min<std::size_t>(n, 0xFFFFFFFFu));
    crc = ::crc32(crc, p, chunk);
    p += chunk;
    n -= chunk;
  }
  value_ = static_cast<uint32_t>(crc) ^ 0xFFFFFFFFu;
#endif
}

uint32_t crc32_file(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path + "' for reading");
  CRC32 crc;
  // Heap allocated: a 256 KiB array on the stack can overflow the small stacks
  // some platforms give secondary threads.
  std::vector<char> buf(kBufferBytes);
  while (in) {
    in.read(buf.data(), static_cast<std::streamsize>(buf.size()));
    const std::streamsize got = in.gcount();
    if (got > 0) crc.update(buf.data(), static_cast<std::size_t>(got));
  }
  if (in.bad()) throw IoError("read error on '" + path + "'");
  return crc.value();
}

std::string crc32_hex(uint32_t crc) {
  return std::format("{:08x}", crc);
}

}  // namespace gcrank
