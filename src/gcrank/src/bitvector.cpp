#include "gcrank/detail/bitvector.h"

#include <algorithm>
#include <bit>
#include <cstring>
#include <fstream>

#include "gcrank/detail/crc32.h"

namespace gcrank {

namespace {
// The format stores 64-bit words little-endian byte first, independent of the
// host byte order.
inline uint64_t load_le64(const uint8_t* p) {
  uint64_t v = 0;
  std::memcpy(&v, p, sizeof(v));
  if constexpr (std::endian::native == std::endian::big) v = std::byteswap(v);
  return v;
}
inline void store_le64(uint8_t* p, uint64_t v) {
  if constexpr (std::endian::native == std::endian::big) v = std::byteswap(v);
  std::memcpy(p, &v, sizeof(v));
}
// Rank scalars are always written little-endian (the manifest records the
// byte order, so a reader on a big-endian host swaps).
inline void store_le16(uint8_t* p, uint16_t v) {
  p[0] = static_cast<uint8_t>(v & 0xFFu);
  p[1] = static_cast<uint8_t>(v >> 8);
}
}  // namespace

// ---- build -----------------------------------------------------------------

void BitVector::clear() {
  words_.clear();
  words_view_ = {};
  num_words_ = 0;
  size_ = 0;
  clear_rank();
  maps_.clear();
  mapped_ = false;
}

void BitVector::clear_rank() {
  sb_owned_.clear();
  bb_owned_.clear();
  sb_view_ = {};
  bb_view_ = {};
  num_sb_ = 0;
  num_bb_ = 0;
}

void BitVector::bind_rank_views() {
  sb_view_ = sb_owned_;
  bb_view_ = bb_owned_;
}

void BitVector::push_bit(bool value) {
  if (mapped_) throw Error("BitVector::push_bit on a mapped vector");
  if ((size_ & 63u) == 0) {
    words_.push_back(0);
    num_words_ = words_.size();
  }
  if (value) words_.back() |= (uint64_t(1) << (63 - (size_ & 63u)));
  ++size_;
}

void BitVector::append_bits(std::span<const uint64_t> src, uint64_t nbits) {
  if (nbits == 0) return;
  if (mapped_) throw Error("BitVector::append_bits on a mapped vector");
  const uint64_t nsrc = (nbits + 63) / 64;
  if (src.size() < nsrc)
    throw UsageError("BitVector::append_bits: source holds fewer than nbits bits");
  const uint64_t rem = size_ & 63u;
  if (rem == 0) {
    // Whole words: copy the block directly instead of inserting word by word.
    const std::size_t old_words = words_.size();
    words_.resize(old_words + static_cast<std::size_t>(nsrc));
    if (nsrc > 0)
      std::memcpy(words_.data() + old_words, src.data(),
                  static_cast<std::size_t>(nsrc) * sizeof(uint64_t));
    num_words_ = words_.size();
    size_ += nbits;
    return;
  }
  const unsigned shift = static_cast<unsigned>(64 - rem);
  const std::size_t old_words = words_.size();
  words_.back() |= src[0] >> rem;
  const std::size_t appended = static_cast<std::size_t>(nsrc - 1);
  words_.resize(old_words + appended);
  uint64_t carry = src[0] << shift;
  for (std::size_t i = 1; i < static_cast<std::size_t>(nsrc); ++i) {
    words_[old_words + i - 1] = carry | (src[i] >> rem);
    carry = src[i] << shift;
  }
  size_ += nbits;
  const uint64_t new_words = (size_ + 63) / 64;
  if (new_words > old_words + appended) words_.push_back(carry);
  num_words_ = words_.size();
}

void BitVector::finalize() {
  num_sb_ = rank_superblock_entries(size_);
  num_bb_ = rank_basic_block_entries(size_);
  sb_owned_.assign(static_cast<std::size_t>(num_sb_), 0);
  bb_owned_.assign(static_cast<std::size_t>(num_bb_), 0);
  const std::span<const uint64_t> words = words_data();
  uint64_t total = 0;
  uint64_t next_sb = 0;
  uint64_t next_bb = 0;
  for (uint64_t i = 0; i < num_words_; ++i) {
    const uint64_t pos = i << 6;
    if ((pos & (kSuperblockBits - 1)) == 0) {
      const uint64_t b = pos >> kSuperblockShift;
      while (next_sb <= b) sb_owned_[static_cast<std::size_t>(next_sb++)] = total;
    }
    if ((pos & (kBasicBlockBits - 1)) == 0) {
      const uint64_t b = pos >> kBasicBlockShift;
      const uint64_t base = sb_owned_[pos >> kSuperblockShift];
      while (next_bb <= b)
        bb_owned_[static_cast<std::size_t>(next_bb++)] = static_cast<uint16_t>(total - base);
    }
    total += static_cast<uint64_t>(std::popcount(words[static_cast<std::size_t>(i)]));
  }
  while (next_sb <= num_sb_ - 1) sb_owned_[static_cast<std::size_t>(next_sb++)] = total;
  while (next_bb <= num_bb_ - 1) {
    const uint64_t base = sb_owned_[(next_bb * kBasicBlockBits) >> kSuperblockShift];
    bb_owned_[static_cast<std::size_t>(next_bb++)] = static_cast<uint16_t>(total - base);
  }
  bind_rank_views();
}

bool BitVector::bit(uint64_t i) const {
  if (i >= size_) throw QueryError("BitVector::bit: index out of range");
  return ((words_data()[static_cast<std::size_t>(i >> 6)] >> (63 - (i & 63u))) & 1u) != 0;
}

uint64_t BitVector::rank1(uint64_t pos) const {
  if (sb_view_.empty())
    throw Error("BitVector::rank1: rank index not built (call finalize() or load_rank())");
  if (pos > size_) throw QueryError("BitVector::rank1: position out of range");
  const uint64_t sb = pos >> kSuperblockShift;
  const uint64_t bb = pos >> kBasicBlockShift;
  uint64_t r = sb_view_[static_cast<std::size_t>(sb)] + bb_view_[static_cast<std::size_t>(bb)];
  const unsigned bit_in_bb = static_cast<unsigned>(pos & (kBasicBlockBits - 1));
  const uint64_t base_word = bb * (kBasicBlockBits / 64);
  const uint64_t nwords = bit_in_bb >> 6;
  const std::span<const uint64_t> words = words_data();
  for (uint64_t i = 0; i < nwords; ++i)
    r += static_cast<uint64_t>(std::popcount(words[static_cast<std::size_t>(base_word + i)]));
  const unsigned tail = bit_in_bb & 63u;
  if (tail != 0)
    r += static_cast<uint64_t>(
        std::popcount(words[static_cast<std::size_t>(base_word + nwords)] & mask_high_bits(tail)));
  return r;
}

// ---- persistence -----------------------------------------------------------

BitVector& BitVector::operator=(const BitVector& other) {
  if (this == &other) return *this;
  clear();
  size_ = other.size_;
  num_words_ = other.num_words_;
  if (num_words_ > 0) {
    const std::span<const uint64_t> src = other.words_data();
    words_.assign(src.begin(), src.end());
  }
  if (!other.sb_view_.empty()) {
    sb_owned_.assign(other.sb_view_.begin(), other.sb_view_.end());
    bb_owned_.assign(other.bb_view_.begin(), other.bb_view_.end());
    num_sb_ = other.num_sb_;
    num_bb_ = other.num_bb_;
    bind_rank_views();
  }
  return *this;
}

void BitVector::read_file_exact(const std::string& path, std::vector<char>& out,
                                uint64_t expect) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path + "' for reading");
  in.seekg(0, std::ios::end);
  const std::streamoff fsize = in.tellg();
  if (fsize < 0 || static_cast<uint64_t>(fsize) != expect)
    throw FormatError("'" + path + "': file size " + std::to_string(fsize) +
                      " does not match the expected " + std::to_string(expect) +
                      " bytes implied by the manifest");
  in.seekg(0, std::ios::beg);
  out.resize(static_cast<std::size_t>(expect));
  if (expect > 0) {
    in.read(out.data(), static_cast<std::streamsize>(expect));
    if (static_cast<uint64_t>(in.gcount()) != expect || in.bad())
      throw IoError("short read on '" + path + "'");
  }
}

FileInfo BitVector::save_bv(const std::string& path) const {
  const uint64_t total_words = bv_bytes(size_) / 8;
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path + "' for writing");
  CRC32 crc;
  const bool le = host_is_little_endian();
  const std::span<const uint64_t> words = words_data();
  constexpr uint64_t kChunkWords = 8192;  // 64 KiB buffer
  std::vector<uint8_t> buf(static_cast<std::size_t>(kChunkWords) * 8);
  for (uint64_t i = 0; i < total_words;) {
    const uint64_t chunk = std::min<uint64_t>(total_words - i, kChunkWords);
    const std::size_t bytes = static_cast<std::size_t>(chunk) * 8;
    if (le && i + chunk <= num_words_) {
      std::memcpy(buf.data(), words.data() + static_cast<std::size_t>(i), bytes);
    } else {
      uint8_t* p = buf.data();
      for (uint64_t j = 0; j < chunk; ++j, p += 8)
        store_le64(p, (i + j) < num_words_ ? words[static_cast<std::size_t>(i + j)] : 0);
    }
    out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(bytes));
    if (!out) throw OutputError("write error on '" + path + "'");
    crc.update(buf.data(), bytes);
    i += chunk;
  }
  out.flush();
  if (!out) throw OutputError("write error on '" + path + "'");
  return FileInfo{total_words * 8, crc.value()};
}

FileInfo BitVector::save_rank(const std::string& path) const {
  if (sb_view_.empty())
    throw Error("BitVector::save_rank: no rank index (call finalize())");
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path + "' for writing");
  CRC32 crc;
  // Scalars are always stored little-endian; the manifest's `endianness`
  // records that, so readers on big-endian hosts swap.
  const std::size_t sb_bytes = sb_view_.size() * sizeof(uint64_t);
  const std::size_t bb_bytes = bb_view_.size() * sizeof(uint16_t);
  std::vector<uint8_t> buf(sb_bytes + bb_bytes);
  uint8_t* p = buf.data();
  for (const uint64_t rank : sb_view_) {
    store_le64(p, rank);
    p += sizeof(uint64_t);
  }
  for (const uint16_t rank : bb_view_) {
    store_le16(p, rank);
    p += sizeof(uint16_t);
  }
  out.write(reinterpret_cast<const char*>(buf.data()),
            static_cast<std::streamsize>(buf.size()));
  out.flush();
  if (!out) throw OutputError("write error on '" + path + "'");
  crc.update(buf.data(), buf.size());
  return FileInfo{sb_bytes + bb_bytes, crc.value()};
}

void BitVector::load_bv(const std::string& path, uint64_t nbits) {
  const uint64_t expect = bv_bytes(nbits);
  std::vector<char> raw;
  read_file_exact(path, raw, expect);
  words_view_ = {};
  maps_.clear();
  mapped_ = false;
  size_ = nbits;
  num_words_ = expect / 8;
  words_.assign(static_cast<std::size_t>(num_words_), 0);
  if (host_is_little_endian()) {
    if (expect > 0) std::memcpy(words_.data(), raw.data(), static_cast<std::size_t>(expect));
  } else {
    const uint8_t* p = reinterpret_cast<const uint8_t*>(raw.data());
    for (uint64_t i = 0; i < num_words_; ++i)
      words_[static_cast<std::size_t>(i)] = load_le64(p + 8 * i);
  }
}

void BitVector::load_rank(const std::string& path, uint64_t nbits, bool stored_little_endian) {
  const uint64_t nsb = rank_superblock_entries(nbits);
  const uint64_t nbb = rank_basic_block_entries(nbits);
  std::vector<char> raw;
  read_file_exact(path, raw, rank_bytes(nbits));
  const uint8_t* base = reinterpret_cast<const uint8_t*>(raw.data());
  const bool swap = (stored_little_endian != host_is_little_endian());
  clear_rank();
  num_sb_ = nsb;
  num_bb_ = nbb;
  sb_owned_.resize(static_cast<std::size_t>(nsb));
  bb_owned_.resize(static_cast<std::size_t>(nbb));
  for (uint64_t i = 0; i < nsb; ++i) {
    uint64_t v = 0;
    std::memcpy(&v, base + 8 * i, sizeof(v));
    sb_owned_[static_cast<std::size_t>(i)] = swap ? std::byteswap(v) : v;
  }
  for (uint64_t j = 0; j < nbb; ++j) {
    uint16_t v = 0;
    std::memcpy(&v, base + 8 * nsb + 2 * j, sizeof(v));
    bb_owned_[static_cast<std::size_t>(j)] = swap ? std::byteswap(v) : v;
  }
  bind_rank_views();
  validate_rank_index(path);
}

void BitVector::map_bv(MappedFile file, uint64_t nbits) {
  const uint64_t expect = bv_bytes(nbits);
  if (file.size() != expect)
    throw FormatError("'" + file.name() + "': mapped size " + std::to_string(file.size()) +
                      " does not match the expected " + std::to_string(expect) + " bytes");
  words_view_ = {};
  words_.clear();
  size_ = nbits;
  num_words_ = expect / 8;
  if (host_is_little_endian() && file.data() != nullptr) {
    words_view_ = std::span<const uint64_t>(
        reinterpret_cast<const uint64_t*>(file.data()), static_cast<std::size_t>(num_words_));
    maps_.push_back(std::move(file));  // keeps the view alive
    mapped_ = true;
  } else {
    // A big-endian host (or an empty vector): materialize into heap storage.
    // Nothing stays mapped, so is_mapped() must say so (the file is released).
    if (num_words_ > 0) {
      words_.assign(static_cast<std::size_t>(num_words_), 0);
      const uint8_t* p = file.data();
      for (uint64_t i = 0; i < num_words_; ++i)
        words_[static_cast<std::size_t>(i)] = load_le64(p + 8 * i);
    }
    mapped_ = false;
  }
}

void BitVector::map_rank(MappedFile file, uint64_t nbits, bool stored_little_endian) {
  const uint64_t nsb = rank_superblock_entries(nbits);
  const uint64_t nbb = rank_basic_block_entries(nbits);
  const uint64_t expect = rank_bytes(nbits);
  if (file.size() != expect)
    throw FormatError("'" + file.name() + "': mapped size " + std::to_string(file.size()) +
                      " does not match the expected " + std::to_string(expect) + " bytes");
  const std::string name = file.name();
  clear_rank();
  num_sb_ = nsb;
  num_bb_ = nbb;
  if (stored_little_endian == host_is_little_endian()) {
    sb_view_ = std::span<const uint64_t>(reinterpret_cast<const uint64_t*>(file.data()),
                                         static_cast<std::size_t>(nsb));
    bb_view_ = std::span<const uint16_t>(
        reinterpret_cast<const uint16_t*>(file.data() + nsb * 8), static_cast<std::size_t>(nbb));
    maps_.push_back(std::move(file));  // keeps the views alive
    mapped_ = true;
  } else {
    // Byte order differs: materialize swapped copies on the heap instead of
    // swapping in every rank1() call.
    const uint8_t* base = file.data();
    sb_owned_.resize(static_cast<std::size_t>(nsb));
    bb_owned_.resize(static_cast<std::size_t>(nbb));
    for (uint64_t i = 0; i < nsb; ++i) {
      uint64_t v = 0;
      std::memcpy(&v, base + 8 * i, sizeof(v));
      sb_owned_[static_cast<std::size_t>(i)] = std::byteswap(v);
    }
    for (uint64_t j = 0; j < nbb; ++j) {
      uint16_t v = 0;
      std::memcpy(&v, base + 8 * nsb + 2 * j, sizeof(v));
      bb_owned_[static_cast<std::size_t>(j)] = std::byteswap(v);
    }
    bind_rank_views();
  }
  validate_rank_index(name);
}

void BitVector::validate_rank_index(const std::string& path) const {
  uint64_t prev = 0;
  for (const uint64_t rank : sb_view_) {
    if (rank < prev)
      throw FormatError("'" + path + "': rank index is not monotonic");
    prev = rank;
  }
  // Within a superblock, basic-block ranks are relative to the superblock start
  // and must therefore restart at zero and never decrease. Padding entries past
  // the end of the vector repeat the end-of-vector rank, which also satisfies
  // this.
  constexpr std::size_t kBasicBlocksPerSuperblock = kSuperblockBits / kBasicBlockBits;
  for (std::size_t j = 0; j < bb_view_.size(); ++j) {
    const uint16_t rank = bb_view_[j];
    if (rank > kSuperblockBits)
      throw FormatError("'" + path + "': basic-block rank entry out of range");
    if (j % kBasicBlocksPerSuperblock == 0) {
      if (rank != 0)
        throw FormatError("'" + path +
                          "': basic-block rank does not restart at a superblock boundary");
    } else if (rank < bb_view_[j - 1]) {
      throw FormatError("'" + path +
                        "': basic-block ranks are not monotonic within a superblock");
    }
  }
}

}  // namespace gcrank
