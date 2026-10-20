#pragma once
//
// SHA-256 (FIPS 180-4), implemented from scratch: the library depends only on
// the C++23 standard library and zlib, and zlib has no hash function.
//
#include <cstddef>
#include <cstdint>
#include <string>

namespace gcrank {

class SHA256 {
public:
  SHA256() { reset(); }

  void reset();
  void update(const void* data, std::size_t n);
  // Writes the 32-byte digest of all data fed so far.
  void finish(uint8_t out[32]);

private:
  void transform(const uint8_t block[64]);

  uint32_t h_[8];
  uint8_t buf_[64];
  std::size_t buf_len_ = 0;
  uint64_t total_bits_ = 0;
};

// Lowercase hex digest of a memory buffer.
std::string sha256_hex(const void* data, std::size_t n);

// Lowercase hex encoding of an already-computed 32-byte digest.
std::string sha256_hex_from_digest(const uint8_t digest[32]);

// Lowercase hex digest of the exact bytes of a file. Throws IoError.
std::string sha256_file_hex(const std::string& path);

}  // namespace gcrank
