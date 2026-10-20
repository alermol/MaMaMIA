#pragma once
//
// Bitvector with a two-level rank index.
//
// Bit i lives in word i/64 at bit position (63 - i%64): the file is a stream
// of 64-bit little-endian words, so the last byte of each 8-byte group holds
// the most significant bits of its word.
//
// Rank index (per bitvector):
//   superblock = 4096 bits, sb_rank[i] = absolute rank at i*4096   (uint64)
//   basic block = 512 bits, bb_rank[j] = rank from the superblock
//                                        start at j*512            (uint16)
// Both arrays carry one extra trailing entry (the rank at the end of the
// vector), so rank1(pos) is defined for every pos in [0, size]. Rank support
// costs ~4.7% of the bitvector's size.
//
#include <bit>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

#include "gcrank/error.h"
#include "gcrank/format.h"
#include "gcrank/detail/input_stream.h"

namespace gcrank {

// Returns a mask with the top `n` bits set (n <= 64).
inline constexpr uint64_t mask_high_bits(unsigned n) {
  if (n == 0) return 0;
  if (n >= 64) return ~uint64_t(0);
  return ~uint64_t(0) << (64 - n);
}

class BitVector {
public:
  BitVector() = default;
  ~BitVector() = default;

  // Copying materializes mapped data; moving keeps views valid.
  BitVector(const BitVector& other) { *this = other; }
  BitVector& operator=(const BitVector& other);
  BitVector(BitVector&&) noexcept = default;
  BitVector& operator=(BitVector&&) noexcept = default;

  // ---- build ---------------------------------------------------------------

  void clear();

  // Appends one bit.
  void push_bit(bool value);

  // Appends `nbits` bits from `src`, a packed block in the same MSB-first
  // order produced by the FASTA packer; padding bits after `nbits` in the last
  // word must be zero. `src` must hold at least ceil(nbits / 64) words.
  void append_bits(std::span<const uint64_t> src, uint64_t nbits);

  // Builds the rank index over the current contents. Must be called after all
  // appends and before any rank1() query.
  void finalize();
  bool has_rank() const { return !sb_view_.empty(); }

  // ---- access --------------------------------------------------------------

  uint64_t size() const { return size_; }  // logical bits
  uint64_t num_words() const { return num_words_; }
  bool bit(uint64_t i) const;

  // Number of set bits in [0, pos). Requires finalize()/load_rank() and
  // pos <= size(); throws QueryError otherwise.
  uint64_t rank1(uint64_t pos) const;

  // Raw word access (little-endian grouping); safe for i < num_words().
  uint64_t word(uint64_t i) const { return words_data()[static_cast<std::size_t>(i)]; }

  // ---- persistence ---------------------------------------------------------

  // Writes the .bv file (zero padded to a whole number of superblocks) and
  // returns its size and the CRC-32 of the bytes as written.
  FileInfo save_bv(const std::string& path) const;

  // Writes the .rb file, always little-endian (the manifest records the byte
  // order, so readers on big-endian hosts swap).
  FileInfo save_rank(const std::string& path) const;

  // Loads a .bv / .rb file into heap storage. `nbits` is the logical size.
  void load_bv(const std::string& path, uint64_t nbits);
  void load_rank(const std::string& path, uint64_t nbits, bool stored_little_endian);

  // Same, but using a read-only mmap view instead of copying. When the stored
  // byte order differs from the host's, the rank arrays are copied to be
  // swapped (the bitvector bytes never need swapping).
  void map_bv(MappedFile file, uint64_t nbits);
  void map_rank(MappedFile file, uint64_t nbits, bool stored_little_endian);

  // True when at least part of the data is a view into a mapped file rather
  // than heap storage (a big-endian host materializes even mapped input).
  bool is_mapped() const { return mapped_; }

private:
  std::span<const uint64_t> words_data() const {
    return words_view_.empty() ? std::span<const uint64_t>(words_) : words_view_;
  }
  void clear_rank();
  void bind_rank_views();
  void validate_rank_index(const std::string& path) const;
  static void read_file_exact(const std::string& path, std::vector<char>& out, uint64_t expect);

  std::vector<uint64_t> words_;             // owned words (build / heap load)
  std::span<const uint64_t> words_view_{};  // mapped words
  uint64_t num_words_ = 0;
  uint64_t size_ = 0;                       // logical bits, excluding padding

  std::vector<uint64_t> sb_owned_;
  std::vector<uint16_t> bb_owned_;
  std::span<const uint64_t> sb_view_{};
  std::span<const uint16_t> bb_view_{};
  uint64_t num_sb_ = 0;
  uint64_t num_bb_ = 0;

  std::vector<MappedFile> maps_;  // keeps mapped files alive
  bool mapped_ = false;
};

}  // namespace gcrank
