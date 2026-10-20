#pragma once
//
// CRC-32 (IEEE 802.3, reflected polynomial 0xEDB88320).
//
// A thin wrapper over zlib's crc32(), which is a required dependency of the
// library anyway. The result is bit-identical to the table-driven version this
// replaced (and to common `cksum`-style tools); zlib's hand-tuned variant is
// roughly 5x faster.
//
#include <cstddef>
#include <cstdint>
#include <string>

namespace gcrank {

// Streaming CRC-32. The default-constructed object is ready for updates.
class CRC32 {
public:
  CRC32() = default;

  // Mixes `n` bytes into the running checksum.
  void update(const void* data, std::size_t n);

  // Final CRC-32 of everything fed so far (after the final XOR).
  uint32_t value() const { return value_ ^ 0xFFFFFFFFu; }

  void reset() { value_ = 0xFFFFFFFFu; }

  // One-shot convenience.
  static uint32_t compute(const void* data, std::size_t n) {
    CRC32 c;
    c.update(data, n);
    return c.value();
  }

private:
  uint32_t value_ = 0xFFFFFFFFu;
};

// CRC-32 of the exact bytes of a file. Throws IoError on failure to read.
uint32_t crc32_file(const std::string& path);

// Formats a checksum as 8 lowercase hex digits.
std::string crc32_hex(uint32_t crc);

}  // namespace gcrank
