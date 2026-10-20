#pragma once
//
// Byte-oriented input streams (plain and gzip) and a read-only memory mapping.
//
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>

#include "gcrank/detail/sha256.h"

namespace gcrank {

// Abstract byte source.
class InputStream {
public:
  virtual ~InputStream();

  // Reads up to `n` bytes into `dst`; returns the number of bytes read.
  // Returns 0 only at end of input. Throws InputError on a read failure.
  virtual std::size_t read(void* dst, std::size_t n) = 0;

  // Source path, for diagnostics.
  virtual const std::string& name() const = 0;
};

// Opens `path` and picks the decompressor by sniffing the gzip magic bytes
// (0x1f 0x8b), never by looking at the file extension. The path "-" reads
// standard input (plain or gzip; zlib detects which, and a dup'd descriptor is
// used so the stream can be closed safely).
// Throws InputError when the file cannot be opened.
std::unique_ptr<InputStream> open_input(const std::string& path);

// Uncompressed file reader.
class PlainInputStream final : public InputStream {
public:
  explicit PlainInputStream(std::string path);
  // Reads from an already-open stream; does not take ownership of `file`.
  PlainInputStream(std::FILE* file, std::string name);
  ~PlainInputStream() override;

  std::size_t read(void* dst, std::size_t n) override;
  const std::string& name() const override { return name_; }

private:
  std::FILE* file_ = nullptr;
  bool owned_ = false;
  std::string name_;
};

// gzip reader (handles concatenated gzip members; transparently passes through
// data that is not gzip, which is what zlib's gzread does).
class GzipInputStream final : public InputStream {
public:
  explicit GzipInputStream(std::string path);
  // Adopts an open gzFile (closes it on destruction).
  GzipInputStream(void* gz, std::string name);
  ~GzipInputStream() override;

  std::size_t read(void* dst, std::size_t n) override;
  const std::string& name() const override { return name_; }

private:
  void* gz_ = nullptr;  // gzFile, kept opaque to avoid a zlib.h include here
  std::string name_;
};

// Decorator that hashes everything read through the wrapped stream and counts
// the bytes, so a build can hash and parse the input in a single pass (no
// second read, no time-of-check/time-of-use window between the two).
class HashingInputStream final : public InputStream {
public:
  explicit HashingInputStream(std::unique_ptr<InputStream> inner);
  ~HashingInputStream() override;

  std::size_t read(void* dst, std::size_t n) override;
  const std::string& name() const override { return inner_->name(); }

  // Bytes read so far.
  uint64_t bytes() const { return bytes_; }
  // Lowercase hex SHA-256 of everything read so far; safe to call once the
  // stream is exhausted (the digest is of the bytes actually consumed).
  std::string sha256_hex() const;

private:
  std::unique_ptr<InputStream> inner_;
  SHA256 hash_;
  uint64_t bytes_ = 0;
};

// RAII read-only memory mapping of a whole file.
//
// A default-constructed / reset MappedFile is invalid; `map()` throws IoError
// when the mapping cannot be established, so callers can fall back to reading.
// Mapping an empty file succeeds and yields data() == nullptr, size() == 0.
//
// The mapping is advisory (MADV_RANDOM): the OS reads in the pages the caller
// touches. Truncating or rewriting the underlying file while it is mapped
// raises SIGBUS on the next access to a removed page, so a mapped index must
// not be modified in place (gcrank's own save path builds a new directory and
// renames it into place, which is safe; `cp` over a mapped index is not).
class MappedFile {
public:
  MappedFile() = default;
  ~MappedFile();

  MappedFile(MappedFile&& other) noexcept;
  MappedFile& operator=(MappedFile&& other) noexcept;
  MappedFile(const MappedFile&) = delete;
  MappedFile& operator=(const MappedFile&) = delete;

  // Maps the whole file read-only. Throws IoError on failure.
  static MappedFile map(const std::string& path);

  const uint8_t* data() const { return data_; }
  std::size_t size() const { return size_; }
  // True when the object holds a usable mapping (including an empty file).
  bool valid() const { return mapped_; }
  const std::string& name() const { return name_; }

  void reset();

private:
  const uint8_t* data_ = nullptr;
  std::size_t size_ = 0;
  std::string name_;
  bool mapped_ = false;
#ifdef _WIN32
  void* file_handle_ = nullptr;  // HANDLE
  void* mapping_handle_ = nullptr;
#endif
};

}  // namespace gcrank
