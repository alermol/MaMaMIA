#pragma once
//
// The index manifest (manifest.json) and a minimal JSON codec for exactly
// this schema.
//
// The writer is canonical: object keys are sorted, indentation is two spaces,
// lines end with LF, and the document ends with a single newline. The reader
// is strict about the required keys and the major format version, and ignores
// unknown keys so that minor revisions stay loadable.
//
#include <cstdint>
#include <map>
#include <string>
#include <string_view>

#include "gcrank/format.h"

namespace gcrank {

// Where the index came from.
struct SourceInfo {
  std::string path;         // path as passed to the builder ("-" for stdin)
  // Hex SHA-256 of the parsed byte stream: the decompressed bytes for gzip
  // input, so a file and its gzip form hash identically. The source is hashed
  // while it is parsed, so this always describes exactly the indexed bytes.
  std::string sha256;
  uint64_t file_bytes = 0;  // bytes parsed (decompressed size for gzip input)
  uint64_t total_bases = 0; // bases over all contigs (including N)
  uint64_t num_contigs = 0;
};

struct Manifest {
  // Identity
  std::string format_version = kFormatVersion;
  std::string builder_version = kBuilderVersion;
  std::string created;  // ISO-8601 UTC, honours SOURCE_DATE_EPOCH

  // Source FASTA
  SourceInfo source;

  // Index description
  std::string bit_order = kBitOrderMsb0;
  // Byte order of the .rb scalars. Writers store little-endian; readers swap
  // when this says "big", so an index in the opposite byte order is read
  // correctly.
  std::string endianness = kEndiannessLittle;
  uint32_t superblock_bits = kSuperblockBits;
  uint32_t basic_block_bits = kBasicBlockBits;
  bool remove_ambiguous_supported = true;

  // Population counts of the two bitvectors (`total_gc` and `total_valid` in
  // the JSON `index` object). Optional on read; on load they are cross-checked
  // against the end-of-vector ranks, which catches a corrupt rank index even
  // when the full data-file CRCs are not validated.
  bool has_totals = false;
  uint64_t total_gc = 0;     // set bits of genome.bv
  uint64_t total_valid = 0;  // set bits of valid.bv

  // Data file sizes and checksums
  std::map<std::string, FileInfo> files;

  // Informative description of what the bits mean
  std::map<std::string, std::string> semantics;

  // Canonical JSON text (sorted keys, 2-space indent, LF, trailing newline).
  std::string to_json() const;

  // Parses and validates. `origin` only appears in error messages.
  // Throws FormatError on malformed JSON, a missing required key, an
  // unsupported bit order/endianness, or a different major format version.
  static Manifest from_json(std::string_view text, const std::string& origin = "<memory>");

  void save(const std::string& path) const;
  static Manifest load(const std::string& path);
};

// "YYYY-MM-DDTHH:MM:SSZ" from SOURCE_DATE_EPOCH, or the current UTC time when
// the variable is unset. This is the only wall-clock read in the build path.
std::string created_from_env();

// Informative `semantics` block written by the builder.
std::map<std::string, std::string> default_semantics();

}  // namespace gcrank
