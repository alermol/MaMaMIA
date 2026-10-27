#include "gcrank/gc_rank.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <iterator>
#include <limits>
#include <mutex>
#include <span>
#include <string_view>
#include <thread>
#include <utility>

#if defined(_WIN32)
#include <process.h>
#else
#include <fcntl.h>
#include <unistd.h>
#endif

#include "gcrank/detail/bitvector.h"
#include "gcrank/detail/crc32.h"
#include "gcrank/error.h"
#include "gcrank/detail/fasta_reader.h"
#include "gcrank/format.h"
#include "gcrank/detail/input_stream.h"
#include "gcrank/detail/thread_pool.h"

namespace gcrank {

namespace fs = std::filesystem;

unsigned default_threads() {
  const unsigned hardware = std::thread::hardware_concurrency();
  const unsigned available = hardware == 0 ? 4u : hardware;
  return std::min(8u, std::max(1u, available));
}

namespace {

// ---- parallel build pipeline --------------------------------------------
//
//   hash helper (1 thread)  hashes the read buffers, in read order
//   producer (1 thread)     parses FASTA -> bounded ring of slots (chunk_bits bases)
//   packers (N threads)     pack A/C/G/T/N bases into GC and valid words
//   assembler (caller)      appends blocks in order, records contig boundaries
//
// Slots are recycled through the ring; the assembler consumes in submission
// order, so the output does not depend on the thread count. The hash helper
// lives in HashingInputStream, so hashing overlaps the parse instead of
// extending the serial critical path.

struct Slot {
  enum class Kind { ContigBegin, ContigEnd, Block };

  Kind kind = Kind::Block;
  std::string name;            // ContigBegin
  uint64_t nbases = 0;         // Block
  std::vector<char> bases;     // Block, producer filled
  std::vector<uint64_t> gc_words;     // Block, worker filled
  std::vector<uint64_t> valid_words;  // Block, worker filled

  bool free = true;   // producer may fill it
  bool ready = false; // assembler may consume it
};

// Per-byte packing code: bit 0 is the GC bit, bit 1 the valid bit. A table
// lookup keeps the inner loop free of the data-dependent branches that a
// base_is_gc()/base_is_valid() pair turns into on random sequence.
struct PackTable {
  std::array<uint8_t, 256> code{};
  constexpr PackTable() {
    for (int c = 0; c < 256; ++c) {
      const BaseClass cls = base_class(static_cast<char>(c));
      code[static_cast<std::size_t>(c)] =
          static_cast<uint8_t>((base_is_gc(cls) ? 1u : 0u) | (base_is_valid(cls) ? 2u : 0u));
    }
  }
};
inline constexpr PackTable kPackTable{};

// Packs bases into MSB-first 64-bit words (bit 0 of the block = MSB of word 0).
// Padding bits after `n` bases in the last word are zero.
//
// Eight bases are packed per iteration: each one produces one byte for the GC
// word and one for the valid word, which is then shifted into place. This is
// several times faster than shifting a 64-bit accumulator once per base on
// random sequence, and produces identical words.
void pack_block(std::span<const char> bases, std::vector<uint64_t>& gc_words,
                std::vector<uint64_t>& valid_words) {
  const std::size_t n = bases.size();
  const std::size_t nwords = (n + 63) / 64;
  gc_words.assign(nwords, 0);
  valid_words.assign(nwords, 0);

  const unsigned char* bytes = reinterpret_cast<const unsigned char*>(bases.data());
  std::size_t i = 0;
  for (; i + 8 <= n; i += 8) {
    unsigned gc_byte = 0;
    unsigned valid_byte = 0;
    for (unsigned j = 0; j < 8; ++j) {
      const unsigned code = kPackTable.code[bytes[i + j]];
      gc_byte = (gc_byte << 1) | (code & 1u);
      valid_byte = (valid_byte << 1) | (code >> 1);
    }
    const std::size_t word = i >> 6;
    const unsigned shift = 56u - static_cast<unsigned>(i & 63u);
    gc_words[word] |= static_cast<uint64_t>(gc_byte) << shift;
    valid_words[word] |= static_cast<uint64_t>(valid_byte) << shift;
  }

  if (i < n) {
    uint64_t gc_acc = 0;
    uint64_t valid_acc = 0;
    const std::size_t start = i;
    for (; i < n; ++i) {
      const unsigned code = kPackTable.code[bytes[i]];
      gc_acc = (gc_acc << 1) | (code & 1u);
      valid_acc = (valid_acc << 1) | (code >> 1);
    }
    const std::size_t count = n - start;
    const unsigned shift = 64u - static_cast<unsigned>(start & 63u) -
                           static_cast<unsigned>(count);
    gc_words[start >> 6] |= gc_acc << shift;
    valid_words[start >> 6] |= valid_acc << shift;
  }
}

// Thrown inside the producer thread when the assembler is unwinding.
struct PipelineAbort {};

class BuildPipeline {
public:
  BuildPipeline(InputStream& input, uint64_t chunk_bases, unsigned workers, bool verbose,
                const LogSink& log)
      : input_(input),
        chunk_(chunk_bases),
        verbose_(verbose),
        log_(log),
        capacity_(static_cast<uint64_t>(workers) + 2),
        slots_(static_cast<std::size_t>(capacity_)),
        pool_(workers) {}

  void run(BitVector& genome, BitVector& valid, std::vector<ContigInfo>& contigs) {
    genome_ = &genome;
    valid_ = &valid;
    contigs_ = &contigs;
    producer_ = std::jthread([this](std::stop_token stop) { produce(stop); });
    try {
      assemble();
    } catch (...) {
      shutdown();
      throw;
    }
    shutdown();
  }

private:
  void log(LogLevel level, std::string message) {
    if (log_) log_(level, message);
  }

  // ---- producer ----------------------------------------------------------

  void produce(std::stop_token stop) {
    try {
      FastaReader reader(input_, chunk_);
      uint64_t index = 0;
      uint64_t bases = 0;
      auto next_progress = std::chrono::steady_clock::now() + std::chrono::seconds(4);
      for (;;) {
        const FastaReader::Event event = reader.next();
        if (event == FastaReader::Event::End) break;

        Slot& slot = slots_[static_cast<std::size_t>(index % capacity_)];
        {
          std::unique_lock<std::mutex> lock(mutex_);
          cv_free_.wait(lock, [&] {
            return slot.free || failed_.load() || stop.stop_requested();
          });
          if (failed_.load() || stop.stop_requested()) throw PipelineAbort{};
          slot.free = false;
        }

        bool is_block = false;
        switch (event) {
          case FastaReader::Event::ContigBegin:
            slot.kind = Slot::Kind::ContigBegin;
            slot.name = reader.contig_name();
            break;
          case FastaReader::Event::Block:
            slot.kind = Slot::Kind::Block;
            // Take the reader's buffer instead of copying it; the slot returns
            // its old buffer, so both sides keep their capacity.
            reader.take_block(slot.bases);
            slot.nbases = static_cast<uint64_t>(slot.bases.size());
            bases += slot.nbases;
            is_block = true;
            break;
          case FastaReader::Event::ContigEnd:
            slot.kind = Slot::Kind::ContigEnd;
            slot.name.clear();
            break;
          default:
            break;
        }
        publish(index, is_block);
        ++index;

        if (verbose_ && is_block) {
          const auto now = std::chrono::steady_clock::now();
          if (now >= next_progress) {
            next_progress = now + std::chrono::seconds(4);
            log(LogLevel::Info, std::format("parsed {} bases", bases));
          }
        }
      }
      if (verbose_) log(LogLevel::Info, std::format("parsed {} bases", bases));
    } catch (const PipelineAbort&) {
      // The assembler is failing; unwind quietly.
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!producer_error_ && !failed_.load()) producer_error_ = std::current_exception();
    }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      producer_done_ = true;
    }
    cv_ready_.notify_all();
  }

  void publish(uint64_t index, bool is_block) {
    if (is_block) {
      // One pool task per block: the pool's queue is the work queue, so this
      // class keeps no second one.
      {
        std::lock_guard<std::mutex> lock(mutex_);
        produce_index_ = index + 1;
      }
      pool_.submit([this, index] { pack_slot(index); });
    } else {
      {
        std::lock_guard<std::mutex> lock(mutex_);
        produce_index_ = index + 1;
        slots_[static_cast<std::size_t>(index % capacity_)].ready = true;
      }
      cv_ready_.notify_one();
    }
  }

  // ---- workers -----------------------------------------------------------

  void pack_slot(uint64_t index) {
    Slot& slot = slots_[static_cast<std::size_t>(index % capacity_)];
    try {
      pack_block(slot.bases, slot.gc_words, slot.valid_words);
      // Keep the capacity for reuse; this is the buffer the peak-memory
      // budget is about.
      slot.bases.clear();
    } catch (...) {
      std::lock_guard<std::mutex> lock(mutex_);
      if (!worker_error_) worker_error_ = std::current_exception();
      failed_ = true;
      cv_ready_.notify_all();
      return;
    }
    {
      std::lock_guard<std::mutex> lock(mutex_);
      slot.ready = true;
    }
    cv_ready_.notify_one();
  }

  // ---- assembler (calling thread) ----------------------------------------

  void assemble() {
    ContigInfo current;
    for (;;) {
      std::exception_ptr error;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_ready_.wait(lock, [&] {
          if (failed_.load() || producer_error_) return true;
          if (consume_index_ < produce_index_ &&
              slots_[static_cast<std::size_t>(consume_index_ % capacity_)].ready)
            return true;
          return producer_done_ && consume_index_ == produce_index_;
        });
        if (worker_error_) {
          error = worker_error_;
        } else if (producer_error_) {
          error = producer_error_;
        } else if (failed_.load()) {
          error = std::make_exception_ptr(Error("gcrank: internal pipeline failure"));
        } else if (producer_done_ && consume_index_ == produce_index_) {
          break;  // all items consumed
        }
      }
      if (error) std::rethrow_exception(error);

      Slot& slot = slots_[static_cast<std::size_t>(consume_index_ % capacity_)];
      switch (slot.kind) {
        case Slot::Kind::ContigBegin:
          current.name = slot.name;
          current.length = 0;
          current.offset = genome_->size();
          break;
        case Slot::Kind::Block:
          genome_->append_bits(slot.gc_words, slot.nbases);
          valid_->append_bits(slot.valid_words, slot.nbases);
          current.length += slot.nbases;
          break;
        case Slot::Kind::ContigEnd:
          contigs_->push_back(current);
          break;
      }
      {
        std::lock_guard<std::mutex> lock(mutex_);
        slot.ready = false;
        slot.free = true;
        ++consume_index_;
      }
      cv_free_.notify_one();
    }
  }

  void shutdown() {
    producer_.request_stop();
    cv_free_.notify_all();
    cv_ready_.notify_all();
    if (producer_.joinable()) producer_.join();
    try {
      pool_.wait();
      pool_.stop();
    } catch (...) {
      // Worker exceptions are already captured in worker_error_.
    }
  }

  InputStream& input_;
  uint64_t chunk_;
  bool verbose_;
  const LogSink& log_;
  uint64_t capacity_;
  std::vector<Slot> slots_;

  std::mutex mutex_;
  std::condition_variable cv_free_;
  std::condition_variable cv_ready_;
  uint64_t produce_index_ = 0;
  uint64_t consume_index_ = 0;
  bool producer_done_ = false;
  std::exception_ptr producer_error_;
  std::exception_ptr worker_error_;
  std::atomic<bool> failed_{false};

  std::jthread producer_;
  ThreadPool pool_;

  BitVector* genome_ = nullptr;
  BitVector* valid_ = nullptr;
  std::vector<ContigInfo>* contigs_ = nullptr;
};

// ---- index files ---------------------------------------------------------

FileInfo write_contigs(const fs::path& path, const std::vector<ContigInfo>& contigs) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path.string() + "' for writing");
  CRC32 crc;
  uint64_t bytes = 0;
  std::string line;
  for (const ContigInfo& contig : contigs) {
    if (contig.name.find('\n') != std::string::npos)
      throw InputError("contig name contains a newline; refusing to write '" + path.string() +
                       "'");
    line.clear();
    std::format_to(std::back_inserter(line), "{}\t{}\t{}\n", contig.name, contig.length,
                   contig.offset);
    out.write(line.data(), static_cast<std::streamsize>(line.size()));
    if (!out) throw OutputError("write error on '" + path.string() + "'");
    crc.update(line.data(), line.size());
    bytes += line.size();
  }
  out.flush();
  if (!out) throw OutputError("write error on '" + path.string() + "'");
  return FileInfo{bytes, crc.value()};
}

uint64_t parse_u64_field(std::string_view text, const std::string& path, std::string_view what) {
  if (text.empty()) throw FormatError("'" + path + "': empty " + std::string(what));
  uint64_t value = 0;
  const std::from_chars_result parsed =
      std::from_chars(text.data(), text.data() + text.size(), value);
  if (parsed.ec == std::errc::result_out_of_range)
    throw FormatError("'" + path + "': " + std::string(what) + " overflows: '" +
                      std::string(text) + "'");
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
    throw FormatError("'" + path + "': " + std::string(what) + " is not a number: '" +
                      std::string(text) + "'");
  return value;
}

std::vector<ContigInfo> parse_contigs(const fs::path& path, uint64_t expected_contigs,
                                      uint64_t expected_bases) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw IoError("cannot open '" + path.string() + "' for reading");
  const std::string name = path.string();
  std::vector<ContigInfo> contigs;
  std::string line;
  uint64_t next_offset = 0;
  while (std::getline(in, line)) {
    std::string_view view(line);
    if (!view.empty() && view.back() == '\r') view.remove_suffix(1);
    if (view.empty()) continue;
    const std::size_t tab_offset = view.rfind('\t');
    if (tab_offset == std::string_view::npos)
      throw FormatError("'" + name + "': expected '<name>\\t<length>\\t<offset>'");
    const std::size_t tab_length = view.rfind('\t', tab_offset - 1);
    if (tab_length == std::string_view::npos)
      throw FormatError("'" + name + "': expected '<name>\\t<length>\\t<offset>'");
    ContigInfo contig;
    contig.name = view.substr(0, tab_length);
    contig.length = parse_u64_field(view.substr(tab_length + 1, tab_offset - tab_length - 1),
                                    name, "contig length");
    contig.offset = parse_u64_field(view.substr(tab_offset + 1), name, "contig offset");
    if (contig.offset != next_offset)
      throw FormatError("'" + name + "': contig '" + contig.name + "' starts at " +
                        std::to_string(contig.offset) + " but " +
                        std::to_string(next_offset) + " was expected");
    next_offset += contig.length;
    contigs.push_back(std::move(contig));
  }
  if (in.bad()) throw IoError("read error on '" + name + "'");
  if (contigs.size() != expected_contigs)
    throw FormatError("'" + name + "': holds " + std::to_string(contigs.size()) +
                      " contigs but the manifest declares " +
                      std::to_string(expected_contigs));
  if (next_offset != expected_bases)
    throw FormatError("'" + name + "': contigs cover " + std::to_string(next_offset) +
                      " bases but the manifest declares " + std::to_string(expected_bases));
  return contigs;
}

void write_sidecar(const fs::path& path, uint32_t crc) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw OutputError("cannot open '" + path.string() + "' for writing");
  const std::string text = crc32_hex(crc) + "\n";
  out.write(text.data(), static_cast<std::streamsize>(text.size()));
  out.flush();
  if (!out) throw OutputError("write error on '" + path.string() + "'");
}

// Re-reads the data files and compares them against the manifest. This is what
// makes the temporary-directory -> final rename safe. Returns an empty string
// when everything matches.
std::string verify_index_files(const fs::path& dir, const Manifest& manifest) {
  for (const auto& entry : manifest.files) {
    const fs::path path = dir / entry.first;
    std::error_code ec;
    const uint64_t size = fs::file_size(path, ec);
    if (ec)
      return "cannot stat '" + path.string() + "': " + ec.message();
    if (size != entry.second.size)
      return "'" + path.string() + "' holds " + std::to_string(size) + " bytes but " +
             std::to_string(entry.second.size) + " were written";
    uint32_t crc = 0;
    try {
      crc = crc32_file(path.string());
    } catch (const IoError& error) {
      return error.what();
    }
    if (crc != entry.second.crc32)
      return "'" + path.string() + "' failed its checksum validation";
  }
  return std::string();
}

// A token that is unique among concurrent saves: the process id plus a
// per-process counter, so two saves of the same index cannot share a
// temporary directory.
std::string unique_save_token() {
  static std::atomic<uint64_t> counter{0};
#if defined(_WIN32)
  const uint64_t pid = static_cast<uint64_t>(::_getpid());
#else
  const uint64_t pid = static_cast<uint64_t>(::getpid());
#endif
  return std::to_string(pid) + "." + std::to_string(counter.fetch_add(1));
}

// Flushes a file (or, on a best-effort basis, a directory) to stable storage.
// The temp-directory rename is only crash-safe once the bytes it publishes have
// reached the disk, because verify_index_files() re-reads through the page
// cache and proves nothing about durability.
void sync_to_disk(const fs::path& path, bool directory) {
#if defined(_WIN32)
  (void)path;
  (void)directory;  // no portable equivalent; the Windows path is untested
#else
  const int flags = O_RDONLY | (directory ? O_DIRECTORY : 0);
  const int fd = ::open(path.c_str(), flags);
  if (fd < 0)
    throw OutputError("cannot open '" + path.string() + "' to flush it: " + std::strerror(errno));
  const int flushed = ::fsync(fd);
  const int err = errno;
  ::close(fd);
  // Some filesystems reject fsync on a directory; that is not fatal here.
  if (flushed != 0 && !(directory && (err == EINVAL || err == ENOTSUP)))
    throw OutputError("cannot flush '" + path.string() + "': " + std::strerror(err));
#endif
}

// Reads a `*.crc32` sidecar ("<8 lowercase hex digits>\n") and returns the
// checksum; throws FormatError when the file is missing or malformed.
uint32_t read_sidecar(const fs::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw FormatError("cannot open '" + path.string() + "' for reading");
  std::string text;
  in >> text;
  if (text.size() != 8 || in.bad())
    throw FormatError("'" + path.string() + "': expected 8 hex digits");
  uint32_t value = 0;
  const std::from_chars_result parsed =
      std::from_chars(text.data(), text.data() + text.size(), value, 16);
  if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
    throw FormatError("'" + path.string() + "': expected 8 hex digits");
  return value;
}

// The sidecars are small and cheap, so they are cross-checked against the
// manifest on every load, even when the expensive data-file CRCs are skipped.
void verify_sidecars(const fs::path& dir, const Manifest& manifest) {
  for (const auto& [name, info] : manifest.files) {
    const uint32_t sidecar = read_sidecar(dir / (name + kCrcSuffix));
    if (sidecar != info.crc32)
      throw FormatError("'" + (dir / (name + kCrcSuffix)).string() +
                        "' disagrees with the manifest checksum of '" + name + "'");
  }
}

void load_bitvector(BitVector& bitvector, const fs::path& bv_path, const fs::path& rb_path,
                    uint64_t nbits, bool stored_little_endian, bool use_mmap,
                    const LogSink& log) {
  if (use_mmap) {
    try {
      MappedFile bv_file = MappedFile::map(bv_path.string());
      MappedFile rb_file = MappedFile::map(rb_path.string());
      bitvector.map_bv(std::move(bv_file), nbits);
      bitvector.map_rank(std::move(rb_file), nbits, stored_little_endian);
      return;
    } catch (const IoError& error) {
      // A library must not write to stderr; report through the caller's sink.
      if (log) log(LogLevel::Warning, std::string(error.what()) + "; falling back to reading");
    }
  }
  bitvector.load_bv(bv_path.string(), nbits);
  bitvector.load_rank(rb_path.string(), nbits, stored_little_endian);
}

}  // namespace

// ---- GcRank ---------------------------------------------------------------

// The out-of-line state. Keeping it here means the public header shows only
// GcRank, the option structs, ContigInfo, Manifest and the error types.
struct GcRank::Impl {
  // Transparent hash so contig() and the contig-relative queries can look up a
  // string_view without constructing a std::string.
  struct NameHash {
    using is_transparent = void;
    std::size_t operator()(std::string_view name) const {
      return std::hash<std::string_view>{}(name);
    }
  };

  BitVector genome_;
  BitVector valid_;
  std::vector<ContigInfo> contigs_;
  std::unordered_map<std::string, std::size_t, NameHash, std::equal_to<>> contig_index_;
  Manifest manifest_;
  bool built_ = false;
};

GcRank::GcRank() : impl_(std::make_unique<Impl>()) {}
GcRank::~GcRank() = default;
GcRank::GcRank(GcRank&&) noexcept = default;
GcRank& GcRank::operator=(GcRank&&) noexcept = default;

void GcRank::build(const std::string& fasta_path, const BuildOptions& options) {
  if (fasta_path.empty()) throw UsageError("build: empty FASTA path");
  if (options.chunk_bits == 0) throw UsageError("build: chunk_bits must be > 0");
  const unsigned workers = options.threads == 0 ? default_threads() : options.threads;

  // Validate SOURCE_DATE_EPOCH before parsing: it is only needed for the
  // manifest, but a bad value should not surface after an hours-long build.
  const std::string created = created_from_env();

  // Open the source and hash it while parsing: one pass over the input, and no
  // time-of-check/time-of-use window between the hash and the parse. Non-
  // seekable sources ("-" for stdin, pipes, process substitution) work too,
  // because nothing needs the file size up front.
  HashingInputStream input(open_input(fasta_path));

  const auto started = std::chrono::steady_clock::now();
  auto log = [&options](LogLevel level, std::string message) {
    if (options.log) options.log(level, message);
  };
  if (options.verbose)
    log(LogLevel::Info,
        std::format("building from {} with {} worker thread(s), chunk {} bases",
                    input.name(), workers, options.chunk_bits));

  // Build into locals and commit only on success, so a failed build leaves any
  // loaded index usable (the strong exception guarantee).
  BitVector genome;
  BitVector valid;
  std::vector<ContigInfo> contigs;

  BuildPipeline pipeline(input, options.chunk_bits, workers, options.verbose, options.log);
  pipeline.run(genome, valid, contigs);

  genome.finalize();
  valid.finalize();

  Manifest manifest;
  manifest.format_version = kFormatVersion;
  manifest.builder_version = kBuilderVersion;
  manifest.created = created;
  manifest.source.path = fasta_path;
  manifest.source.sha256 = input.sha256_hex();
  manifest.source.file_bytes = input.bytes();
  manifest.source.total_bases = genome.size();
  manifest.source.num_contigs = static_cast<uint64_t>(contigs.size());
  manifest.bit_order = kBitOrderMsb0;
  manifest.endianness = kEndiannessLittle;  // .rb scalars are always little-endian
  manifest.superblock_bits = kSuperblockBits;
  manifest.basic_block_bits = kBasicBlockBits;
  manifest.remove_ambiguous_supported = true;
  manifest.total_gc = genome.rank1(genome.size());
  manifest.total_valid = valid.rank1(valid.size());
  manifest.has_totals = true;
  manifest.semantics = default_semantics();

  const uint64_t total_gc = manifest.total_gc;
  impl_->genome_ = std::move(genome);
  impl_->valid_ = std::move(valid);
  impl_->contigs_ = std::move(contigs);
  rebuild_contig_index();
  impl_->manifest_ = std::move(manifest);
  impl_->built_ = true;

  if (options.verbose) {
    const double seconds =
        std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    log(LogLevel::Info,
        std::format("built {} bases, {} contigs, {} GC in {:.2f}s",
                    impl_->manifest_.source.total_bases, impl_->manifest_.source.num_contigs, total_gc, seconds));
  }
}

void GcRank::save(const std::string& out_dir, const std::string& name, bool overwrite) const {
  if (!impl_->built_) throw Error("save: no index in memory (call build() or load() first)");
  if (name.empty()) throw UsageError("save: index name must not be empty");
  if (name == "." || name == "..") throw UsageError("save: invalid index name '" + name + "'");
  if (name.find('/') != std::string::npos || name.find('\\') != std::string::npos)
    throw UsageError("save: index name must not contain path separators: '" + name + "'");
  if (!impl_->genome_.has_rank() || !impl_->valid_.has_rank())
    throw Error("save: rank indices are not built");

  const fs::path out(out_dir);
  const fs::path final_dir = out / (name + ".gcrank");

  std::error_code ec;
  fs::create_directories(out, ec);
  if (ec)
    throw OutputError("cannot create output directory '" + out.string() + "': " + ec.message());
  const bool target_exists = fs::exists(final_dir, ec);
  if (ec)
    throw OutputError("cannot stat '" + final_dir.string() + "': " + ec.message());
  if (target_exists && !overwrite)
    throw OutputError("'" + final_dir.string() +
                      "' already exists (pass overwrite/--overwrite to replace it)");

  // Unique per call, so two concurrent saves of the same index cannot delete
  // each other's work; the creation is retried in the unlikely case that the
  // name is already taken by a leftover from a crashed process with a reused
  // pid, because that directory must not be reused.
  fs::path tmp_dir;
  bool created = false;
  for (int attempt = 0; attempt < 8 && !created; ++attempt) {
    tmp_dir = out / (name + ".gcrank.tmp." + unique_save_token());
    created = fs::create_directories(tmp_dir, ec);
    if (ec) throw OutputError("cannot create '" + tmp_dir.string() + "': " + ec.message());
  }
  if (!created)
    throw OutputError("cannot create a temporary directory for '" + final_dir.string() + "'");
  fs::path old_dir;  // named freshly on every publish attempt below

  // The existing index is not touched until the new one has been written and
  // verified: it is moved aside only just before the new directory is renamed
  // into place, and moved back if anything fails in between.
  bool old_moved = false;
  bool published = false;
  try {
    Manifest manifest = impl_->manifest_;
    if (manifest.created.empty()) manifest.created = created_from_env();
    // .rb scalars are always little-endian, and the population totals are
    // recomputed from the in-memory rank indices.
    manifest.endianness = kEndiannessLittle;
    manifest.total_gc = impl_->genome_.rank1(impl_->genome_.size());
    manifest.total_valid = impl_->valid_.rank1(impl_->valid_.size());
    manifest.has_totals = true;
    manifest.files.clear();
    manifest.files[kContigsFileName] = write_contigs(tmp_dir / kContigsFileName, impl_->contigs_);
    manifest.files[kGenomeBvFileName] = impl_->genome_.save_bv((tmp_dir / kGenomeBvFileName).string());
    manifest.files[kGenomeRbFileName] = impl_->genome_.save_rank((tmp_dir / kGenomeRbFileName).string());
    manifest.files[kValidBvFileName] = impl_->valid_.save_bv((tmp_dir / kValidBvFileName).string());
    manifest.files[kValidRbFileName] = impl_->valid_.save_rank((tmp_dir / kValidRbFileName).string());
    for (const char* data_file : kDataFiles)
      write_sidecar(tmp_dir / (std::string(data_file) + kCrcSuffix),
                    manifest.files[data_file].crc32);
    manifest.save((tmp_dir / kManifestFileName).string());
    if (const std::string problem = verify_index_files(tmp_dir, manifest); !problem.empty())
      throw OutputError(problem + " after writing; the index was not published");
    for (const auto& entry : fs::directory_iterator(tmp_dir))
      if (entry.is_regular_file()) sync_to_disk(entry.path(), false);
    sync_to_disk(tmp_dir, true);  // the new entries themselves, too

    // Publish the verified directory. With overwrite, the existing index is
    // moved aside first and restored if the rename fails; because another
    // process may be saving the same index concurrently, the move/publish pair
    // is retried a few times rather than assumed to be the only writer.
    std::error_code publish_error;
    for (int attempt = 0; attempt < 5 && !published; ++attempt) {
      const bool present = fs::exists(final_dir, ec);
      if (ec)
        throw OutputError("cannot stat '" + final_dir.string() + "': " + ec.message());
      if (!present) {
        fs::rename(tmp_dir, final_dir, publish_error);
        if (!publish_error) published = true;
        continue;  // retry when another saver published in the meantime
      }
      if (!overwrite)
        throw OutputError("'" + final_dir.string() +
                          "' already exists (pass overwrite/--overwrite to replace it)");
      old_dir = out / (name + ".gcrank.old." + unique_save_token());
      fs::rename(final_dir, old_dir, publish_error);
      if (publish_error) continue;  // lost the race for the existing index; retry
      old_moved = true;
      fs::rename(tmp_dir, final_dir, publish_error);
      if (!publish_error) {
        published = true;
        break;
      }
      // Undo the move and try again; keep the failure from the last attempt.
      std::error_code ignored;
      fs::rename(old_dir, final_dir, ignored);
      if (!ignored) old_moved = false;
    }
    if (!published)
      throw OutputError("cannot publish '" + final_dir.string() + "': " +
                        (publish_error ? publish_error.message()
                                       : "another process is writing the same index"));
    sync_to_disk(out, true);  // make the rename itself durable
    if (old_moved) {
      std::error_code ignored;
      fs::remove_all(old_dir, ignored);  // best effort; a leftover is harmless
    }
  } catch (...) {
    std::error_code ignored;
    if (old_moved && !published) {
      // Restore the previous index. Should even this fail, the old index still
      // exists under old_dir instead of being gone.
      fs::rename(old_dir, final_dir, ignored);
    }
    fs::remove_all(tmp_dir, ignored);
    throw;
  }
}

void GcRank::load(const std::string& dir, const OpenOptions& options) {
  const fs::path root(dir);
  Manifest manifest = Manifest::load((root / kManifestFileName).string());

  // Checksums first: a corrupt index should report a checksum or sidecar
  // mismatch, not a confusing parse error from inside a corrupted data file.
  if (options.validate_crc)
    if (const std::string problem = verify_index_files(root, manifest); !problem.empty())
      throw FormatError(problem);
  verify_sidecars(root, manifest);

  std::vector<ContigInfo> contigs =
      parse_contigs(root / kContigsFileName, manifest.source.num_contigs,
                    manifest.source.total_bases);

  const bool stored_little_endian = manifest.endianness == kEndiannessLittle;
  BitVector genome;
  BitVector valid;
  load_bitvector(genome, root / kGenomeBvFileName, root / kGenomeRbFileName,
                 manifest.source.total_bases, stored_little_endian, options.mmap, options.log);
  load_bitvector(valid, root / kValidBvFileName, root / kValidRbFileName,
                 manifest.source.total_bases, stored_little_endian, options.mmap, options.log);

  // Cheap structural cross-check: a corrupt rank index changes the end-of-vector
  // rank, which must equal the totals stored in the manifest. This catches it
  // even when validate_crc is false.
  if (manifest.has_totals) {
    if (genome.rank1(genome.size()) != manifest.total_gc ||
        valid.rank1(valid.size()) != manifest.total_valid)
      throw FormatError(root.string() +
                        ": rank indices disagree with the manifest population counts");
  }

  impl_->genome_ = std::move(genome);
  impl_->valid_ = std::move(valid);
  impl_->contigs_ = std::move(contigs);
  rebuild_contig_index();
  impl_->manifest_ = std::move(manifest);
  impl_->built_ = true;
}

void GcRank::rebuild_contig_index() {
  impl_->contig_index_.clear();
  impl_->contig_index_.reserve(impl_->contigs_.size());
  for (std::size_t i = 0; i < impl_->contigs_.size(); ++i)
    impl_->contig_index_.emplace(impl_->contigs_[i].name, i);
}

void GcRank::check_range(uint64_t start, uint64_t end) const {
  if (!impl_->built_) throw UsageError("no index in memory (call build() or load() first)");
  if (start > end)
    throw QueryError("invalid window [" + std::to_string(start) + ", " + std::to_string(end) +
                     "): start > end");
  if (end > impl_->genome_.size())
    throw QueryError("invalid window [" + std::to_string(start) + ", " + std::to_string(end) +
                     "): end exceeds the genome length " + std::to_string(impl_->genome_.size()));
}

uint64_t GcRank::gc_count(uint64_t start, uint64_t end) const {
  check_range(start, end);
  return impl_->genome_.rank1(end) - impl_->genome_.rank1(start);
}

uint64_t GcRank::gc_count(std::string_view contig_name, uint64_t start, uint64_t end) const {
  const auto [begin, finish] = resolve_window(contig_name, start, end);
  return impl_->genome_.rank1(finish) - impl_->genome_.rank1(begin);
}

uint64_t GcRank::valid_count(uint64_t start, uint64_t end) const {
  check_range(start, end);
  return impl_->valid_.rank1(end) - impl_->valid_.rank1(start);
}

uint64_t GcRank::valid_count(std::string_view contig_name, uint64_t start, uint64_t end) const {
  const auto [begin, finish] = resolve_window(contig_name, start, end);
  return impl_->valid_.rank1(finish) - impl_->valid_.rank1(begin);
}

// Shared by both gc_fraction overloads; the window is already validated (but
// may still be degenerate), so this only rejects a zero denominator.
double GcRank::fraction(uint64_t start, uint64_t end, bool remove_ambiguous) const {
  const uint64_t gc = impl_->genome_.rank1(end) - impl_->genome_.rank1(start);
  const uint64_t length = end - start;
  uint64_t denominator = length;
  if (remove_ambiguous) {
    denominator = impl_->valid_.rank1(end) - impl_->valid_.rank1(start);
    if (denominator == 0)
      throw QueryError("gc_fraction(remove_ambiguous = true): window [" +
                       std::to_string(start) + ", " + std::to_string(end) +
                       ") holds no unambiguous A/C/G/T bases");
  } else if (length == 0) {
    throw QueryError("gc_fraction: empty window [" + std::to_string(start) + ", " +
                     std::to_string(end) + ")");
  }
  // `bedtools nuc` converts the counts to float and divides in single
  // precision; doing the same makes "%.6f" of the result identical to its
  // pct_gc column.
  const float quotient = static_cast<float>(gc) / static_cast<float>(denominator);
  return static_cast<double>(quotient);
}

double GcRank::gc_fraction(uint64_t start, uint64_t end, bool remove_ambiguous) const {
  check_range(start, end);
  return fraction(start, end, remove_ambiguous);
}

double GcRank::gc_fraction(std::string_view contig_name, uint64_t start, uint64_t end,
                           bool remove_ambiguous) const {
  const auto [begin, finish] = resolve_window(contig_name, start, end);
  return fraction(begin, finish, remove_ambiguous);
}

WindowStats GcRank::stats_global(uint64_t start, uint64_t end, bool remove_ambiguous) const {
  WindowStats result;
  result.length = end - start;
  result.gc = impl_->genome_.rank1(end) - impl_->genome_.rank1(start);
  result.valid = impl_->valid_.rank1(end) - impl_->valid_.rank1(start);
  const uint64_t denominator = remove_ambiguous ? result.valid : result.length;
  if (denominator == 0) {
    result.gc_fraction = std::numeric_limits<double>::quiet_NaN();
    return result;
  }
  // The same single-precision division as gc_fraction(), so the value is
  // bit-identical to bedtools nuc's pct_gc.
  const float quotient =
      static_cast<float>(result.gc) / static_cast<float>(denominator);
  result.gc_fraction = static_cast<double>(quotient);
  return result;
}

WindowStats GcRank::stats(uint64_t start, uint64_t end, bool remove_ambiguous) const {
  check_range(start, end);
  return stats_global(start, end, remove_ambiguous);
}

WindowStats GcRank::stats(std::string_view contig_name, uint64_t start, uint64_t end,
                          bool remove_ambiguous) const {
  const auto [begin, finish] = resolve_window(contig_name, start, end);
  return stats_global(begin, finish, remove_ambiguous);
}

uint64_t GcRank::total_bases() const { return impl_->genome_.size(); }

const std::vector<ContigInfo>& GcRank::contigs() const { return impl_->contigs_; }

const Manifest& GcRank::manifest() const { return impl_->manifest_; }

const ContigInfo* GcRank::contig(std::string_view name) const {
  const auto found = impl_->contig_index_.find(name);
  return found == impl_->contig_index_.end() ? nullptr : &impl_->contigs_[found->second];
}

std::pair<uint64_t, uint64_t> GcRank::resolve_window(std::string_view name, uint64_t start,
                                                     uint64_t end) const {
  if (!impl_->built_) throw UsageError("no index in memory (call build() or load() first)");
  const ContigInfo* info = contig(name);
  if (info == nullptr) throw QueryError("unknown contig '" + std::string(name) + "'");
  if (start > end)
    throw QueryError("invalid window on contig '" + std::string(name) + "' [" +
                     std::to_string(start) + ", " + std::to_string(end) + "): start > end");
  if (end > info->length)
    throw QueryError("invalid window on contig '" + std::string(name) + "' [" +
                     std::to_string(start) + ", " + std::to_string(end) +
                     "): end exceeds the contig length " + std::to_string(info->length));
  return {info->offset + start, info->offset + end};
}

}  // namespace gcrank
