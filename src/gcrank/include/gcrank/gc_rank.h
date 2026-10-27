#pragma once
//
// GcRank: build, store and query a GC-content index over a FASTA file.
//
// Queries are O(1) in the window length: they are answered from rank indices
// over two bitvectors (GC and valid bases), so a window can cover a whole
// chromosome without scanning it. Every query exists in two flavours: global
// coordinates over the concatenated contigs, and contig-relative coordinates
// (`gc_count("chr1", start, end)`); the latter reject windows that leave the
// named contig.
//
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "gcrank/manifest.h"

namespace gcrank {

// One contig of the indexed genome.
struct ContigInfo {
  std::string name;
  uint64_t length = 0;  // bases, including N
  uint64_t offset = 0;  // 0-based offset of the first base
};

// Counts and GC fraction for one window, computed with one rank lookup pair
// per bitvector. `gc_fraction` is NaN when the fraction is undefined: an empty
// window, or remove_ambiguous with no unambiguous A/C/G/T base.
struct WindowStats {
  uint64_t length = 0;
  uint64_t gc = 0;
  uint64_t valid = 0;
  double gc_fraction = 0.0;
};

// Severity of a message handed to a user-installed log sink.
enum class LogLevel { Info, Warning };

// Optional log callback. The library never writes to stderr by itself: install
// a sink to receive warnings (always) and progress messages (with
// BuildOptions::verbose). It may be called from any thread, so a sink that
// writes to a shared stream must synchronize.
using LogSink = std::function<void(LogLevel, std::string_view)>;

struct BuildOptions {
  // Bases per parallel work unit. Blocks never span contigs and the number of
  // blocks only affects scheduling, never the bytes written. Each in-flight
  // block costs roughly 1.25x chunk_bits bytes (raw bases plus the two packed
  // bitvector slices), across threads + 2 slots, so lower this for large
  // thread counts on memory-constrained machines.
  uint64_t chunk_bits = 16000000ull;
  unsigned threads = 0;  // worker threads; 0 means auto = default_threads(),
                         // the producer, the hash helper and the assembling
                         // caller are extra
  bool verbose = false;  // progress messages through the log sink
  LogSink log;           // optional; silent when empty
};

struct OpenOptions {
  // Re-read every data file and verify it against the manifest checksum. This
  // costs one full pass over the index (~12 s for a 4 GB index), which defeats
  // the point of mmap for large indexes; pass false to rely on the cheap
  // structural checks (file sizes, sidecar checksums, rank totals) instead.
  bool validate_crc = true;
  bool mmap = false;  // map data files instead of reading them
  LogSink log;        // receives warnings such as the mmap fallback
};

// Default worker count: min(8, hardware_concurrency()), never below 1.
[[nodiscard]] unsigned default_threads();

class GcRank {
public:
  GcRank();
  ~GcRank();
  GcRank(GcRank&&) noexcept;
  GcRank& operator=(GcRank&&) noexcept;
  GcRank(const GcRank&) = delete;
  GcRank& operator=(const GcRank&) = delete;

  // Parses `fasta_path` (plain or gzip; "-" reads stdin) and builds both
  // bitvectors and their rank indices in memory. Parallel: one producer
  // thread, `threads` packers, a hash helper, and the calling thread
  // assembles. The result is byte-identical for any thread count. The source
  // is hashed while it is parsed (one pass, no second read). Throws
  // InputError for unreadable/malformed input.
  void build(const std::string& fasta_path, const BuildOptions& options = {});

  // Writes the index to `<out_dir>/<name>.gcrank/`. The index is written to a
  // uniquely named temporary directory, verified against its checksums, synced
  // to disk, and then renamed into place, so an interrupted save never leaves a
  // half-written index behind. Refuses to replace an existing index unless
  // `overwrite` is true; with overwrite, the existing index is moved aside only
  // after the new one is verified and is put back if the swap fails.
  void save(const std::string& out_dir, const std::string& name,
            bool overwrite = false) const;

  // Loads an index directory (the `<name>.gcrank` directory itself). Warnings
  // (for example the mmap fallback) are reported through `options.log`.
  void load(const std::string& dir, const OpenOptions& options = {});

  // ---- queries: half-open [start, end), 0-based --------------------------
  //
  // The `start, end` overloads use global coordinates over the concatenated
  // contigs (0 is the first base of the first contig). The `contig` overloads
  // resolve a contig by name and take coordinates relative to its start, like
  // a BED interval; a window that leaves the contig, or an unknown name,
  // throws QueryError.

  // Number of unambiguous G/C bases in the window. U counts as T and is never
  // GC; N and the IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V) are
  // never counted either, so this is a lower bound for windows that contain
  // them.
  [[nodiscard]] uint64_t gc_count(uint64_t start, uint64_t end) const;
  [[nodiscard]] uint64_t gc_count(std::string_view contig, uint64_t start, uint64_t end) const;
  // Number of unambiguous A/C/G/T bases in the window; U counts as T, and N
  // and the IUPAC ambiguity codes are excluded.
  [[nodiscard]] uint64_t valid_count(uint64_t start, uint64_t end) const;
  [[nodiscard]] uint64_t valid_count(std::string_view contig, uint64_t start,
                                     uint64_t end) const;
  // GC fraction of the window, computed exactly like `bedtools nuc`: the
  // integer counts are converted to `float` and divided in single precision,
  // and the result is returned as a double holding that single-precision
  // value. Printing it with "%.6f" is therefore identical to bedtools' pct_gc
  // column, down to the last digit (above 2^24 bases per window the float
  // conversion rounds the counts themselves, the same way bedtools does).
  //
  // With remove_ambiguous = false the denominator is the window length, i.e.
  // everything the window spans (bedtools' `seq_len`); with
  // remove_ambiguous = true it is valid_count, so every ambiguous base is
  // removed (bedtools has no such column). The window must then hold at least
  // one A/C/G/T base, otherwise QueryError is thrown, as it is for a
  // zero-length window either way.
  [[nodiscard]] double gc_fraction(uint64_t start, uint64_t end,
                                   bool remove_ambiguous = false) const;
  [[nodiscard]] double gc_fraction(std::string_view contig, uint64_t start, uint64_t end,
                                   bool remove_ambiguous = false) const;

  // All per-window values in one call. This is the efficient choice for a
  // bedtools-style row (pct_gc, num_gc, num_valid, seq_len): it resolves the
  // contig once and reads each rank index once, instead of once per column.
  // Unlike gc_fraction(), an undefined fraction is reported as NaN rather than
  // thrown, so callers that print "nan" can use the result directly.
  [[nodiscard]] WindowStats stats(uint64_t start, uint64_t end,
                                  bool remove_ambiguous = false) const;
  [[nodiscard]] WindowStats stats(std::string_view contig, uint64_t start, uint64_t end,
                                  bool remove_ambiguous = false) const;

  // ---- metadata ----------------------------------------------------------

  [[nodiscard]] uint64_t total_bases() const;
  [[nodiscard]] const std::vector<ContigInfo>& contigs() const;
  [[nodiscard]] const Manifest& manifest() const;
  // Lookup by contig name (a hash lookup); returns nullptr when absent.
  [[nodiscard]] const ContigInfo* contig(std::string_view name) const;

private:
  // Out-of-line implementation state (bitvectors, contigs, manifest), so this
  // header does not pull the internal headers in.
  struct Impl;
  std::unique_ptr<Impl> impl_;

  void check_range(uint64_t start, uint64_t end) const;
  // Shared by both gc_fraction overloads; rejects a zero denominator.
  [[nodiscard]] double fraction(uint64_t start, uint64_t end, bool remove_ambiguous) const;
  // Shared by both stats overloads; the range is already validated.
  [[nodiscard]] WindowStats stats_global(uint64_t start, uint64_t end,
                                         bool remove_ambiguous) const;
  // Resolves `name`, checks the contig-relative window, and returns the
  // equivalent global window. Throws QueryError when the name is unknown, the
  // index is not built, or the window leaves the contig.
  [[nodiscard]] std::pair<uint64_t, uint64_t> resolve_window(std::string_view contig,
                                                             uint64_t start, uint64_t end) const;
  void rebuild_contig_index();
};

}  // namespace gcrank
