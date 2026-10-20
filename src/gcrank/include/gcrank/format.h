#pragma once
//
// On-disk format contract for `.gcrank` indexes, format_version "0.1.0".
//
// Layout of `<NAME>.gcrank/`:
//
//   manifest.json     canonical JSON (sorted keys, 2-space indent, LF)
//   contigs.tsv       <name>\t<length>\t<offset>  (LF, no header, 0-based)
//   genome.bv         GC bitvector (1 = unambiguous G/C)
//   genome.rb         rank index for genome.bv
//   valid.bv          valid-base bitvector (1 = unambiguous A/C/G/T)
//   valid.rb          rank index for valid.bv
//   <file>.crc32      sidecar checksums, "<8 lowercase hex digits>\n"
//
// Bit order (bit_order = "msb0"): the file is a sequence of 64-bit words
// stored little-endian (the first byte of each 8-byte group is the least
// significant byte of the word). Bit `i` of a bitvector is bit `63 - (i % 64)`
// of word `i / 64`, so in terms of file bytes it is bit `(7 - i % 8)` of byte
// `8 * (i / 64) + 7 - (i % 64) / 8`. This grouping is fixed by the format and
// is host-independent.
//
// genome.bv / valid.bv: `bv_bytes(nbits)` bytes, zero padded to a whole number
// of superblocks (4096 bits). The logical size (total_bases) lives in the
// manifest. Ambiguous input (N, R, Y, S, W, K, M, B, D, H, V) sets neither
// bit: those positions are neither GC nor valid. U is not ambiguous: it is the
// RNA spelling of T, so it counts as valid and never as GC.
//
// genome.rb / valid.rb: multi-byte scalars stored in the byte order recorded
// by the manifest's `endianness` field. Writers store little-endian; readers
// byte-swap when the recorded order differs from the host's, so an index loads
// on any host and an index written in the opposite order is read correctly.
// The layout is
//   (num_superblocks(nbits) + 1) x uint64   absolute rank at i*4096
//   (num_basic_blocks(nbits) + 1) x uint16  rank within the superblock at j*512
// The extra trailing entry of each array is the rank at the end of the vector,
// so rank1(pos) is defined for every pos in [0, nbits].
//
// The manifest's `index` object carries the population counts of the two
// bitvectors (`total_gc`, `total_valid`); loading cross-checks them against the
// end-of-vector ranks, so a corrupt rank index is rejected even when the full
// data-file checksums are skipped.
//
#include <array>
#include <bit>
#include <charconv>
#include <cstdint>
#include <string>
#include <string_view>
#include <system_error>

#include "gcrank/error.h"

namespace gcrank {

// ---- versions ------------------------------------------------------------

// Version of the on-disk format (the public contract).
inline constexpr const char* kFormatVersion = "0.1.0";
// Version of the writer (library/builder).
inline constexpr const char* kBuilderVersion = "0.1.0";

// ---- index layout --------------------------------------------------------

inline constexpr uint32_t kSuperblockBits = 4096;
inline constexpr uint32_t kBasicBlockBits = 512;
// Bit-position shifts; guarded so the two stay in step.
inline constexpr uint32_t kSuperblockShift = 12;  // 1 << 12 == 4096
inline constexpr uint32_t kBasicBlockShift = 9;   // 1 << 9 == 512
static_assert((uint64_t{1} << kSuperblockShift) == kSuperblockBits,
              "kSuperblockShift must match kSuperblockBits");
static_assert((uint64_t{1} << kBasicBlockShift) == kBasicBlockBits,
              "kBasicBlockShift must match kBasicBlockBits");
inline constexpr const char* kBitOrderMsb0 = "msb0";

inline constexpr const char* kManifestFileName = "manifest.json";
inline constexpr const char* kContigsFileName = "contigs.tsv";
inline constexpr const char* kGenomeBvFileName = "genome.bv";
inline constexpr const char* kGenomeRbFileName = "genome.rb";
inline constexpr const char* kValidBvFileName = "valid.bv";
inline constexpr const char* kValidRbFileName = "valid.rb";
inline constexpr const char* kCrcSuffix = ".crc32";

inline constexpr const char* kEndiannessLittle = "little";
inline constexpr const char* kEndiannessBig = "big";

// Data files covered by checksums (manifest `files` and .crc32 sidecars).
inline constexpr const char* kDataFiles[] = {
    kContigsFileName, kGenomeBvFileName, kGenomeRbFileName,
    kValidBvFileName, kValidRbFileName,
};

// ---- sizes ---------------------------------------------------------------

inline constexpr uint64_t num_superblocks(uint64_t nbits) {
  return (nbits + kSuperblockBits - 1) / kSuperblockBits;
}
inline constexpr uint64_t num_basic_blocks(uint64_t nbits) {
  return (nbits + kBasicBlockBits - 1) / kBasicBlockBits;
}
// Physical size of a .bv file for `nbits` logical bits (zero padded).
inline constexpr uint64_t bv_bytes(uint64_t nbits) {
  return num_superblocks(nbits) * (kSuperblockBits / 8);
}
inline constexpr uint64_t rank_superblock_entries(uint64_t nbits) {
  return num_superblocks(nbits) + 1;
}
inline constexpr uint64_t rank_basic_block_entries(uint64_t nbits) {
  return num_basic_blocks(nbits) + 1;
}
// Physical size of a .rb file for `nbits` logical bits.
inline constexpr uint64_t rank_bytes(uint64_t nbits) {
  return rank_superblock_entries(nbits) * sizeof(uint64_t) +
         rank_basic_block_entries(nbits) * sizeof(uint16_t);
}

// Size and checksum of one data file (manifest `files` entry).
struct FileInfo {
  uint64_t size = 0;
  uint32_t crc32 = 0;
};

// ---- host properties -----------------------------------------------------

// True when the host stores multi-byte integers least significant byte first.
inline constexpr bool host_is_little_endian() {
  return std::endian::native == std::endian::little;
}

inline const char* host_endianness_name() {
  return host_is_little_endian() ? kEndiannessLittle : kEndiannessBig;
}

// ---- base semantics ------------------------------------------------------
//
// Semantics of `genome.bv` and `valid.bv`, and of the FASTA parser:
//   genome.bv bit = 1 for an unambiguous G/C
//   valid.bv  bit = 1 for an unambiguous A/C/G/T (U counts as T)
//
// Accepted sequence bytes (upper and lower case are equivalent; the parser
// normalizes to upper case):
//
//   A C G T   unambiguous DNA. A/C/G/T sets the valid bit; G/C also sets the
//             GC bit.
//   U         uracil, the RNA spelling of T. It is normalized to T and is
//             therefore valid but never GC (this is what `bedtools nuc` does
//             as well).
//   N         unknown base.
//   R Y S W K M B D H V
//             IUPAC ambiguity codes (R = A/G, Y = C/T, S = G/C, W = A/T,
//             K = G/T, M = A/C, B = C/G/T, D = A/G/T, H = A/C/T, V = A/C/G).
//
// Ambiguous input sets neither bit: it is neither GC nor valid, so it never
// contributes to gc_count and gc_fraction(remove_ambiguous = true) removes it
// from the denominator. Note that a code which merely *may* be G/C (S, B, V)
// is not counted as GC, so gc_count is a lower bound for windows that contain
// ambiguity codes.
//
// Whitespace inside sequence lines is ignored; any other byte (including gaps
// such as '-', '.' and '*') is an InputError.

enum class BaseClass : uint8_t {
  Invalid = 0,
  A = 1,
  C = 2,
  G = 3,
  T = 4,
  Ambiguous = 5,  // N and the IUPAC ambiguity codes
  Whitespace = 6,
};

// Whitespace bytes skipped inside sequence lines and header lines.
inline constexpr std::string_view kWhitespaceBytes = " \t\v\f\r";

// The IUPAC ambiguity codes and N; upper case (the lower case letters are
// mapped to the same class).
inline constexpr const char* kAmbiguousBasesUpper = "NRYSWKMBDHV";
inline constexpr const char* kAmbiguousBasesLower = "nryswkmbdhv";

struct BaseTable {
  std::array<BaseClass, 256> v{};
  constexpr BaseTable() {
    v[static_cast<unsigned>('A')] = BaseClass::A;
    v[static_cast<unsigned>('C')] = BaseClass::C;
    v[static_cast<unsigned>('G')] = BaseClass::G;
    v[static_cast<unsigned>('T')] = BaseClass::T;
    v[static_cast<unsigned>('a')] = BaseClass::A;
    v[static_cast<unsigned>('c')] = BaseClass::C;
    v[static_cast<unsigned>('g')] = BaseClass::G;
    v[static_cast<unsigned>('t')] = BaseClass::T;
    // Uracil is the RNA spelling of thymine.
    v[static_cast<unsigned>('U')] = BaseClass::T;
    v[static_cast<unsigned>('u')] = BaseClass::T;
    for (int i = 0; kAmbiguousBasesUpper[i] != '\0'; ++i)
      v[static_cast<unsigned char>(kAmbiguousBasesUpper[i])] = BaseClass::Ambiguous;
    for (int i = 0; kAmbiguousBasesLower[i] != '\0'; ++i)
      v[static_cast<unsigned char>(kAmbiguousBasesLower[i])] = BaseClass::Ambiguous;
    v[static_cast<unsigned>(' ')] = BaseClass::Whitespace;
    v[static_cast<unsigned>('\t')] = BaseClass::Whitespace;
    v[static_cast<unsigned>('\v')] = BaseClass::Whitespace;
    v[static_cast<unsigned>('\f')] = BaseClass::Whitespace;
    v[static_cast<unsigned>('\r')] = BaseClass::Whitespace;
  }
};

inline constexpr BaseTable kBaseTable{};

// Upper-case normalization of accepted sequence bytes, with U (uracil) folded
// onto T. Applied to every base the parser hands out.
struct UpperTable {
  std::array<char, 256> v{};
  constexpr UpperTable() {
    for (int c = 0; c < 256; ++c) v[static_cast<std::size_t>(c)] = static_cast<char>(c);
    for (int c = 'a'; c <= 'z'; ++c)
      v[static_cast<std::size_t>(c)] = static_cast<char>(c - 32);
    v[static_cast<unsigned>('U')] = 'T';
    v[static_cast<unsigned>('u')] = 'T';
  }
};

inline constexpr UpperTable kUpperTable{};

// Classification lookup for one byte of sequence data.
inline constexpr BaseClass base_class(char c) {
  return kBaseTable.v[static_cast<unsigned char>(c)];
}
// Upper-case (and U -> T) normalization for the parsed sequence.
inline constexpr char normalized_base(char c) {
  return kUpperTable.v[static_cast<unsigned char>(c)];
}
// True for an unambiguous G/C (the genome.bv bit).
inline constexpr bool base_is_gc(BaseClass cls) {
  return cls == BaseClass::G || cls == BaseClass::C;
}
// True for an unambiguous A/C/G/T/U (the valid.bv bit).
inline constexpr bool base_is_valid(BaseClass cls) {
  return cls >= BaseClass::A && cls <= BaseClass::T;
}
// True for N and the IUPAC ambiguity codes.
inline constexpr bool base_is_ambiguous(BaseClass cls) { return cls == BaseClass::Ambiguous; }

// ---- version parsing -----------------------------------------------------

// Returns the major version of a "MAJOR.MINOR.PATCH" string.
// Throws FormatError when the string is not a version number.
inline int parse_major_version(std::string_view version) {
  int major = 0;
  const char* first = version.data();
  const char* last = first + version.size();
  const std::from_chars_result parsed = std::from_chars(first, last, major);
  if (parsed.ec != std::errc{} || parsed.ptr == first || parsed.ptr == last ||
      *parsed.ptr != '.' || major < 0 || major > 1000000) {
    throw FormatError("malformed version string '" + std::string(version) + "'");
  }
  return major;
}

}  // namespace gcrank
