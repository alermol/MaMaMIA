#include "gcrank/detail/fasta_reader.h"

#include <algorithm>
#include <array>
#include <format>

#include "gcrank/error.h"
#include "gcrank/format.h"

namespace gcrank {

namespace {

// Class and normalized form of one input byte in a single table load. The
// parser looks at every sequence byte, so folding the two lookups together
// removes one dependent load from the hot loop.
struct SequenceInfo {
  BaseClass cls = BaseClass::Invalid;
  char normalized = '\0';
};

constexpr std::array<SequenceInfo, 256> make_sequence_info() {
  std::array<SequenceInfo, 256> table{};
  for (int c = 0; c < 256; ++c) {
    const char byte = static_cast<char>(c);
    table[static_cast<std::size_t>(c)] = SequenceInfo{base_class(byte), normalized_base(byte)};
  }
  return table;
}

inline constexpr std::array<SequenceInfo, 256> kSequenceInfo = make_sequence_info();

// Printable rendering of an offending byte for error messages.
std::string describe_byte(char c) {
  const unsigned char u = static_cast<unsigned char>(c);
  if (u >= 0x20 && u < 0x7f) return std::format("'{}'", c);
  return std::format("0x{:02x}", u);
}

// Extra explanation for the bytes that most often show up in FASTA files that
// are not nucleotide sequences.
std::string invalid_byte_hint(char c) {
  if (c == '-' || c == '.' || c == '*')
    return " (alignment gaps are not nucleotide codes; strip them before indexing)";
  return std::string();
}

}  // namespace

FastaReader::FastaReader(std::unique_ptr<InputStream> in, uint64_t block_bases)
    : owned_(std::move(in)), block_bases_(block_bases) {
  if (owned_ == nullptr) throw UsageError("FastaReader: null input stream");
  in_ = owned_.get();
  init();
}

FastaReader::FastaReader(InputStream& in, uint64_t block_bases)
    : in_(&in), block_bases_(block_bases) {
  init();
}

void FastaReader::init() {
  if (block_bases_ == 0) throw UsageError("FastaReader: block size must be > 0");
  lines_ = std::make_unique<LineReader>(*in_);
  // Reserve a sane amount up front; the buffer grows as needed otherwise.
  block_.reserve(static_cast<std::size_t>(std::min<uint64_t>(block_bases_, uint64_t(1) << 20)));
}

FastaReader::~FastaReader() = default;

std::unique_ptr<FastaReader> FastaReader::from_path(const std::string& path,
                                                    uint64_t block_bases) {
  return std::make_unique<FastaReader>(open_input(path), block_bases);
}

void FastaReader::parse_header() {
  // The name is the first whitespace-delimited token (samtools faidx style).
  const std::size_t begin = line_.find_first_not_of(kWhitespaceBytes, 1);
  if (begin == std::string_view::npos) {
    next_name_.clear();
    return;
  }
  const std::size_t end = line_.find_first_of(kWhitespaceBytes, begin);
  next_name_.assign(line_, begin,
                    end == std::string_view::npos ? std::string_view::npos : end - begin);
}

bool FastaReader::all_whitespace(std::string_view piece, std::size_t begin) {
  for (std::size_t i = begin; i < piece.size(); ++i)
    if (kSequenceInfo[static_cast<unsigned char>(piece[i])].cls != BaseClass::Whitespace)
      return false;
  return true;
}

FastaReader::AppendResult FastaReader::append_sequence(std::string_view piece, std::size_t begin,
                                                       std::size_t column_base) {
  AppendResult result;
  std::size_t i = begin;
  std::size_t written = block_.size();
  const std::size_t capacity = static_cast<std::size_t>(block_bases_);
  uint64_t contig_bases = contig_bases_;
  uint64_t total_bases = total_bases_;

  while (i < piece.size() && written < capacity) {
    // Size the buffer for the bytes this pass may consume, then compact
    // whitespace in place. Resizing once per line avoids a capacity check and
    // a push_back call per base.
    const std::size_t end = std::min(piece.size(), i + (capacity - written));
    const std::size_t old_size = block_.size();
    block_.resize(old_size + (end - i));
    char* const dst = block_.data() + old_size;
    std::size_t out = 0;
    for (; i < end; ++i) {
      const char c = piece[i];
      const SequenceInfo& info = kSequenceInfo[static_cast<unsigned char>(c)];
      if (info.cls == BaseClass::Invalid)
        throw InputError(in_->name() + ": invalid character " + describe_byte(c) +
                         " at line " + std::to_string(lines_->line_number()) + ", column " +
                         std::to_string(column_base + i + 1) + invalid_byte_hint(c));
      if (info.cls == BaseClass::Whitespace) continue;
      dst[out++] = info.normalized;
      ++contig_bases;
      ++total_bases;
    }
    block_.resize(old_size + out);
    written = old_size + out;
  }

  contig_bases_ = contig_bases;
  total_bases_ = total_bases;
  result.consumed = i;
  result.full = written == capacity;
  return result;
}

bool FastaReader::fetch_piece() {
  const char* data = nullptr;
  std::size_t len = 0;
  const LineReader::Status status = lines_->next(&data, &len);
  if (status == LineReader::Status::End) {
    piece_ = {};
    have_piece_ = false;
    return false;
  }
  next_header_line_ = lines_->line_number();
  piece_ = data != nullptr ? std::string_view(data, len) : std::string_view{};
  piece_pos_ = 0;
  piece_ends_line_ = status == LineReader::Status::Line;
  piece_starts_line_ = at_line_start_;
  have_piece_ = true;
  at_line_start_ = false;
  return true;
}

void FastaReader::end_piece() {
  have_piece_ = false;
  piece_pos_ = 0;
  if (piece_ends_line_) {
    at_line_start_ = true;
    line_bytes_ = 0;
  }
  piece_ = {};
}

void FastaReader::read_header_line() {
  // The current piece starts with '>' at the start of a line. Header lines are
  // the only lines that are materialized; sequence data is never buffered.
  header_.assign(piece_.data(), piece_.size());
  have_piece_ = false;
  bool done = piece_ends_line_;
  while (!done) {
    const char* data = nullptr;
    std::size_t len = 0;
    const LineReader::Status status = lines_->next(&data, &len);
    if (status == LineReader::Status::End) break;  // unterminated last line
    header_.append(data != nullptr ? data : "", len);
    done = status == LineReader::Status::Line;
  }
  line_ = header_;
  parse_header();
  at_line_start_ = true;
  line_bytes_ = 0;
  header_pending_ = true;
  if (phase_ == Phase::InSequence) pending_contig_end_ = true;
}

FastaReader::Event FastaReader::next() {
  if (block_valid_) {
    // The caller has consumed the previous block.
    block_.clear();
    block_valid_ = false;
  }

  for (;;) {
    // A finished contig is always reported before anything else.
    if (pending_contig_end_) {
      if (!block_.empty()) {
        block_valid_ = true;
        return Event::Block;  // flush the tail of the contig first
      }
      pending_contig_end_ = false;
      if (at_eof_) phase_ = Phase::Done;
      return Event::ContigEnd;
    }

    // A parsed header waiting to be announced.
    if (header_pending_) {
      if (next_name_.empty())
        throw InputError(in_->name() + ": empty contig name at line " +
                         std::to_string(next_header_line_));
      if (!seen_names_.insert(next_name_).second)
        throw InputError(in_->name() + ": duplicate contig name '" + next_name_ +
                         "' at line " + std::to_string(next_header_line_));
      header_pending_ = false;
      contig_name_ = next_name_;
      contig_bases_ = 0;
      phase_ = Phase::InSequence;
      ++num_contigs_;
      return Event::ContigBegin;
    }

    // Consume the piece in hand.
    if (have_piece_) {
      if (phase_ == Phase::InSequence) {
        const AppendResult result = append_sequence(piece_, piece_pos_, line_bytes_);
        line_bytes_ += result.consumed - piece_pos_;
        piece_pos_ = result.consumed;
        if (result.full) {
          block_valid_ = true;
          return Event::Block;
        }
        end_piece();
        continue;
      }
      // Before the first header only blank lines are allowed; the offending
      // line is reported as soon as a non-whitespace byte appears.
      if (!all_whitespace(piece_, piece_pos_))
        throw InputError(in_->name() + ": sequence data at line " +
                         std::to_string(next_header_line_) +
                         " before the first '>' header");
      end_piece();
      continue;
    }

    // Fetch the next piece of input.
    if (!fetch_piece()) {
      if (phase_ == Phase::InSequence) {
        pending_contig_end_ = true;
        at_eof_ = true;
        continue;
      }
      phase_ = Phase::Done;
      return Event::End;
    }

    // A '>' at the start of a line opens a header, which is buffered whole.
    if (piece_starts_line_ && !piece_.empty() && piece_.front() == '>') {
      read_header_line();
      continue;
    }
  }
}

}  // namespace gcrank
