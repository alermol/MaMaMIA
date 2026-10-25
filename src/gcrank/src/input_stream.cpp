#include "gcrank/detail/input_stream.h"

#include <zlib.h>

#include <array>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <utility>

#include "gcrank/error.h"

#ifndef _WIN32
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#else
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace gcrank {

namespace {
#if defined(_WIN32)
// UTF-8 to UTF-16 for the wide-character file APIs, so non-ASCII paths work
// on Windows (fopen/gzopen take the system ANSI codepage instead).
std::wstring utf8_to_wide(const std::string& path) {
  const int wide_len = ::MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
  if (wide_len <= 0) throw IoError("cannot convert path '" + path + "' to UTF-16");
  std::wstring wide(static_cast<std::size_t>(wide_len), L'\0');
  ::MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, wide.data(), wide_len);
  return wide;
}
#endif
}  // namespace

InputStream::~InputStream() = default;

// ---- PlainInputStream ----------------------------------------------------

PlainInputStream::PlainInputStream(std::string path) : name_(std::move(path)) {
#if defined(_WIN32)
  file_ = ::_wfopen(utf8_to_wide(name_).c_str(), L"rb");
#else
  file_ = std::fopen(name_.c_str(), "rb");
#endif
  if (file_ == nullptr)
    throw InputError("cannot open '" + name_ + "' for reading: " + std::strerror(errno));
}

PlainInputStream::~PlainInputStream() {
  if (file_ != nullptr) std::fclose(file_);
}

std::size_t PlainInputStream::read(void* dst, std::size_t n) {
  if (n == 0) return 0;
  const std::size_t got = std::fread(dst, 1, n, file_);
  if (got < n && std::ferror(file_))
    throw InputError("read error in '" + name_ + "'");
  return got;
}

// ---- GzipInputStream -----------------------------------------------------

GzipInputStream::GzipInputStream(std::string path) : name_(std::move(path)) {
#if defined(_WIN32)
  gzFile z = gzopen_w(utf8_to_wide(name_).c_str(), L"rb");
#else
  gzFile z = gzopen(name_.c_str(), "rb");
#endif
  if (z == nullptr)
    throw InputError("cannot open gzip file '" + name_ + "'");
  gz_ = z;
  // The default input buffer is 8 KiB; a megabyte keeps the inflate loop fed.
  gzbuffer(z, 1u << 20);
}

GzipInputStream::GzipInputStream(void* gz, std::string name)
    : gz_(gz), name_(std::move(name)) {
  if (gz_ == nullptr) throw UsageError("GzipInputStream: null gzFile");
  gzbuffer(static_cast<gzFile>(gz_), 1u << 20);
}

GzipInputStream::~GzipInputStream() {
  if (gz_ != nullptr) gzclose(static_cast<gzFile>(gz_));
}

std::size_t GzipInputStream::read(void* dst, std::size_t n) {
  if (n == 0) return 0;
  const std::size_t kMaxChunk = 1u << 30;  // gzread takes an unsigned length
  const unsigned want = static_cast<unsigned>(n < kMaxChunk ? n : kMaxChunk);
  gzFile gz = static_cast<gzFile>(gz_);
  const int got = gzread(gz, dst, want);
  if (got < 0) {
    int errnum = 0;
    const char* msg = gzerror(gz, &errnum);
    if (errnum == Z_ERRNO) msg = std::strerror(errno);
    throw InputError("gzip read error in '" + name_ + "': " + (msg != nullptr ? msg : "unknown"));
  }
  if (got == 0) {
    // A zero return means "end of stream", not necessarily a clean one:
    // zlib reports a truncated member through gzerror() (Z_BUF_ERROR,
    // "unexpected end of file") instead of failing the gzread() itself, so a
    // partial download would otherwise look like a shorter genome.
    int errnum = 0;
    const char* msg = gzerror(gz, &errnum);
    if (errnum != Z_OK && errnum != Z_STREAM_END) {
      if (errnum == Z_ERRNO) msg = std::strerror(errno);
      throw InputError("gzip read error in '" + name_ + "': " +
                       (msg != nullptr ? msg : "unexpected end of file"));
    }
    return 0;
  }
  return static_cast<std::size_t>(got);
}

// ---- HashingInputStream --------------------------------------------------

HashingInputStream::HashingInputStream(std::unique_ptr<InputStream> inner)
    : inner_(std::move(inner)) {
  if (inner_ == nullptr) throw UsageError("HashingInputStream: null stream");
}

HashingInputStream::~HashingInputStream() = default;

std::size_t HashingInputStream::read(void* dst, std::size_t n) {
  const std::size_t got = inner_->read(dst, n);
  if (got > 0) {
    hash_.update(dst, got);
    bytes_ += got;
  }
  return got;
}

std::string HashingInputStream::sha256_hex() const {
  // SHA256::finish() mutates the state (it pads and finalizes), so finish a
  // copy; this keeps the accessor const and leaves the running hash usable.
  SHA256 copy = hash_;
  uint8_t digest[32];
  copy.finish(digest);
  return sha256_hex_from_digest(digest);
}

// ---- open_input ----------------------------------------------------------

std::unique_ptr<InputStream> open_input(const std::string& path) {
  if (path == "-") {
    // Standard input: hand a duplicate of the descriptor to zlib, which
    // auto-detects gzip and otherwise passes the bytes through unchanged.
    const int duplicated =
#if defined(_WIN32)
        ::_dup(::_fileno(stdin));
#else
        ::dup(::fileno(stdin));
#endif
    if (duplicated < 0)
      throw InputError("cannot duplicate stdin: " + std::string(std::strerror(errno)));
    gzFile gz = gzdopen(duplicated, "rb");
    if (gz == nullptr) {
#if defined(_WIN32)
      ::_close(duplicated);
#else
      ::close(duplicated);
#endif
      throw InputError("cannot open stdin for reading");
    }
    return std::make_unique<GzipInputStream>(static_cast<void*>(gz), "<stdin>");
  }
  std::array<unsigned char, 2> magic{};
  std::size_t got = 0;
#if defined(_WIN32)
  std::FILE* sniff = ::_wfopen(utf8_to_wide(path).c_str(), L"rb");
#else
  std::FILE* sniff = std::fopen(path.c_str(), "rb");
#endif
  if (sniff != nullptr) {
    got = std::fread(magic.data(), 1, magic.size(), sniff);
    std::fclose(sniff);
  } else {
    throw InputError("cannot open '" + path + "' for reading: " + std::strerror(errno));
  }
  if (got == magic.size() && magic[0] == 0x1f && magic[1] == 0x8b)
    return std::make_unique<GzipInputStream>(path);
  return std::make_unique<PlainInputStream>(path);
}

// ---- MappedFile ----------------------------------------------------------

MappedFile::MappedFile(MappedFile&& other) noexcept { *this = std::move(other); }

MappedFile& MappedFile::operator=(MappedFile&& other) noexcept {
  if (this != &other) {
    reset();
    data_ = std::exchange(other.data_, nullptr);
    size_ = std::exchange(other.size_, 0);
    name_ = std::move(other.name_);
    mapped_ = std::exchange(other.mapped_, false);
#ifdef _WIN32
    file_handle_ = std::exchange(other.file_handle_, nullptr);
    mapping_handle_ = std::exchange(other.mapping_handle_, nullptr);
#endif
  }
  return *this;
}

void MappedFile::reset() {
  if (!mapped_) return;
#ifdef _WIN32
  if (data_ != nullptr) ::UnmapViewOfFile(data_);
  if (mapping_handle_ != nullptr) ::CloseHandle(static_cast<HANDLE>(mapping_handle_));
  if (file_handle_ != nullptr) ::CloseHandle(static_cast<HANDLE>(file_handle_));
  file_handle_ = nullptr;
  mapping_handle_ = nullptr;
#else
  if (data_ != nullptr) ::munmap(const_cast<uint8_t*>(data_), size_);
#endif
  data_ = nullptr;
  size_ = 0;
  mapped_ = false;
}

MappedFile MappedFile::map(const std::string& path) {
  MappedFile mf;
  mf.name_ = path;
#ifndef _WIN32
  const int fd = ::open(path.c_str(), O_RDONLY);
  if (fd < 0)
    throw IoError("cannot open '" + path + "' for mmap: " + std::strerror(errno));
  struct stat st;
  std::memset(&st, 0, sizeof(st));
  if (::fstat(fd, &st) != 0) {
    const int err = errno;
    ::close(fd);
    throw IoError("cannot stat '" + path + "': " + std::strerror(err));
  }
  if (st.st_size <= 0) {
    ::close(fd);
    mf.mapped_ = true;  // empty file: valid, no pages to map
    return mf;
  }
  const std::size_t length = static_cast<std::size_t>(st.st_size);
  void* p = ::mmap(nullptr, length, PROT_READ, MAP_PRIVATE, fd, 0);
  const int err = errno;
  ::close(fd);  // the mapping stays valid after the descriptor is closed
  if (p == MAP_FAILED)
    throw IoError("mmap failed for '" + path + "': " + std::strerror(err));
  mf.data_ = static_cast<const uint8_t*>(p);
  mf.size_ = length;
  mf.mapped_ = true;
  // Queries are random-access and the OS should not read ahead; ignoring a
  // failure is fine (some filesystems or kernels may not support it).
  ::madvise(p, length, MADV_RANDOM);
  return mf;
#else
  // Convert the UTF-8 path to UTF-16 and map read-only.
  const std::wstring wide = utf8_to_wide(path);
  HANDLE fh = ::CreateFileW(wide.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr,
                            OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
  if (fh == INVALID_HANDLE_VALUE) throw IoError("cannot open '" + path + "' for mmap");
  LARGE_INTEGER size;
  if (!::GetFileSizeEx(fh, &size)) {
    ::CloseHandle(fh);
    throw IoError("cannot stat '" + path + "'");
  }
  if (size.QuadPart == 0) {
    ::CloseHandle(fh);
    mf.mapped_ = true;
    return mf;
  }
  HANDLE mh = ::CreateFileMappingW(fh, nullptr, PAGE_READONLY, 0, 0, nullptr);
  if (mh == nullptr) {
    ::CloseHandle(fh);
    throw IoError("CreateFileMapping failed for '" + path + "'");
  }
  void* p = ::MapViewOfFile(mh, FILE_MAP_READ, 0, 0, 0);
  if (p == nullptr) {
    ::CloseHandle(mh);
    ::CloseHandle(fh);
    throw IoError("MapViewOfFile failed for '" + path + "'");
  }
  mf.data_ = static_cast<const uint8_t*>(p);
  mf.size_ = static_cast<std::size_t>(size.QuadPart);
  mf.file_handle_ = fh;
  mf.mapping_handle_ = mh;
  mf.mapped_ = true;
  return mf;
#endif
}

MappedFile::~MappedFile() { reset(); }

}  // namespace gcrank
