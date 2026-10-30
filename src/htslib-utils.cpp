// htslib helpers used by the BAM-first workflow.
//
// The package talks to BAM files through Rhtslib directly. This file reads
// headers, builds CSI indexes and measures reads. bam_metrics_rg_cpp() serves
// the pipeline; bam_count_ranges_cpp() serves the public countReads() helper
// and arbitrary GRanges.

#include <Rcpp.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct HtsFileCloser {
  void operator()(htsFile* fp) const {
    if (fp != nullptr) hts_close(fp);
  }
};

struct HeaderDeleter {
  void operator()(sam_hdr_t* hdr) const {
    if (hdr != nullptr) sam_hdr_destroy(hdr);
  }
};

struct BamDeleter {
  void operator()(bam1_t* b) const {
    if (b != nullptr) bam_destroy1(b);
  }
};

struct IndexDeleter {
  void operator()(hts_idx_t* idx) const {
    if (idx != nullptr) hts_idx_destroy(idx);
  }
};

struct ItrDeleter {
  void operator()(hts_itr_t* itr) const {
    if (itr != nullptr) hts_itr_destroy(itr);
  }
};

using HtsFilePtr = std::unique_ptr<htsFile, HtsFileCloser>;
using HeaderPtr = std::unique_ptr<sam_hdr_t, HeaderDeleter>;
using BamPtr = std::unique_ptr<bam1_t, BamDeleter>;
using IndexPtr = std::unique_ptr<hts_idx_t, IndexDeleter>;
using ItrPtr = std::unique_ptr<hts_itr_t, ItrDeleter>;

const char* format_name(enum htsExactFormat format) {
  switch (format) {
  case bam:
    return "BAM";
  case sam:
    return "SAM";
  case cram:
    return "CRAM";
  case fasta_format:
    return "FASTA";
  case fastq_format:
    return "FASTQ";
  case vcf:
    return "VCF";
  case bcf:
    return "BCF";
  default:
    return "unknown";
  }
}

}  // namespace

// [[Rcpp::export]]
std::string htslib_version_cpp() {
  return std::string(hts_version());
}

// [[Rcpp::export]]
Rcpp::List bam_header_info_cpp(const std::string& path) {
  HtsFilePtr fp(hts_open(path.c_str(), "r"));
  if (!fp) {
    Rcpp::stop("Could not open BAM file '%s': %s", path,
               errno != 0 ? std::strerror(errno) : "unknown htslib error");
  }

  const htsFormat* format = hts_get_format(fp.get());
  if (format == nullptr) {
    Rcpp::stop("Could not determine the format of '%s'", path);
  }
  if (format->format == unknown_format || format->format == binary_format ||
      format->format == text_format || format->format == empty_format) {
    Rcpp::stop(
        "Could not determine the format of '%s'; the file may be truncated, "
        "corrupt, or not a BAM file",
        path);
  }
  if (format->format != bam) {
    Rcpp::stop(
        "Unsupported input format '%s' for '%s'; only coordinate-sorted BAM "
        "files are accepted (convert SAM/CRAM to BAM first)",
        format_name(format->format), path);
  }

  HeaderPtr hdr(sam_hdr_read(fp.get()));
  if (!hdr) {
    Rcpp::stop(
        "Could not read the BAM header of '%s'; the file may be truncated, "
        "corrupt, or not a BAM file",
        path);
  }

  const int n_targets = sam_hdr_nref(hdr.get());
  if (n_targets <= 0) {
    Rcpp::stop("BAM file '%s' contains no target sequences in its header", path);
  }

  Rcpp::CharacterVector target_names(n_targets);
  Rcpp::NumericVector target_lengths(n_targets);
  for (int i = 0; i < n_targets; ++i) {
    const char* name = sam_hdr_tid2name(hdr.get(), i);
    if (name == nullptr) {
      Rcpp::stop("BAM header of '%s' has an unnamed target at index %d", path, i);
    }
    target_names[i] = name;
    target_lengths[i] = static_cast<double>(sam_hdr_tid2len(hdr.get(), i));
  }

  Rcpp::DataFrame targets = Rcpp::DataFrame::create(
      Rcpp::Named("name") = target_names,
      Rcpp::Named("length") = target_lengths);

  const char* header_text = sam_hdr_str(hdr.get());
  return Rcpp::List::create(
      Rcpp::Named("path") = path,
      Rcpp::Named("format") = "BAM",
      Rcpp::Named("targets") = targets,
      Rcpp::Named("header_text") = header_text == nullptr ? "" : header_text,
      Rcpp::Named("n_targets") = n_targets);
}

// [[Rcpp::export]]
bool bam_index_build_cpp(const std::string& path,
                         const std::string& index_path,
                         int min_shift = 14,
                         int threads = 1) {
  if (min_shift < 1) {
    Rcpp::stop(
        "min_shift must be at least 1; BAI indexes (min_shift = 0) are not "
        "supported, use the default CSI min_shift = 14");
  }
  if (threads < 1) {
    Rcpp::stop("threads must be at least 1");
  }
  if (path.empty()) {
    Rcpp::stop("A BAM path is required");
  }
  if (index_path.empty()) {
    Rcpp::stop("An explicit CSI index path is required");
  }

  errno = 0;
  const int ret = sam_index_build3(path.c_str(), index_path.c_str(), min_shift,
                                   threads);
  const int saved_errno = errno;
  if (ret != 0) {
    Rcpp::stop(
        "Could not build a CSI index for '%s' (htslib error code %d): %s",
        path, ret,
        saved_errno != 0
            ? std::strerror(saved_errno)
            : "the BAM must be coordinate-sorted and readable, and the "
              "output directory must be writable");
  }
  return true;
}

// Count alignments and covered bases per sample in exact-N windows.
//
// This is the breadth backend. For every window it builds one bitset per
// sample and marks every reference base covered by an aligned block (M, =, X
// and, unless `count_deletions` is false, D). N, I, S, H and P do not cover a
// reference base, so a spliced read leaves a gap. A window's `count` is the
// number of alignments with at least one aligned base in the window, which is
// the semantics of bedtools coverage -counts.
//
// `read_length <= 0` infers the aligned query length from the first passing
// read. One sample counts every read and ignores RG; several samples require an
// RG tag naming a declared read group.
// [[Rcpp::export]]
Rcpp::List bam_metrics_rg_cpp(
    const std::string& path,
    Rcpp::CharacterVector contig_names,
    Rcpp::NumericVector contig_lengths,
    int n_per_chromosome,
    Rcpp::CharacterVector read_group_ids,
    Rcpp::IntegerVector sample_index,
    Rcpp::CharacterVector sample_names,
    int mapq,
    int exclude_flags,
    int read_length,
    bool count_deletions) {
  const int n_contigs = contig_names.size();
  const int n_samples = sample_names.size();
  const bool single_sample = n_samples == 1;
  if (n_contigs < 1) {
    Rcpp::stop("At least one contig is required for read counting");
  }
  if (n_samples < 1) {
    Rcpp::stop("At least one sample is required for read counting");
  }
  if (n_per_chromosome < 1) {
    Rcpp::stop("n_per_chromosome must be at least 1");
  }
  if (contig_lengths.size() != n_contigs) {
    Rcpp::stop("contig_names and contig_lengths must have the same length");
  }
  if (read_group_ids.size() != sample_index.size()) {
    Rcpp::stop("read_group_ids and sample_index must have the same length");
  }

  std::unordered_map<std::string, int> rg_to_sample;
  for (int i = 0; i < read_group_ids.size(); ++i) {
    const std::string id = Rcpp::as<std::string>(read_group_ids[i]);
    const int sample = sample_index[i];
    if (id.empty()) {
      Rcpp::stop("Read group IDs must not be empty");
    }
    if (sample < 1 || sample > n_samples) {
      Rcpp::stop("sample_index must point into sample_names");
    }
    if (!rg_to_sample.emplace(id, sample - 1).second) {
      Rcpp::stop("Duplicate read group ID '%s'", id);
    }
  }

  HtsFilePtr fp(hts_open(path.c_str(), "r"));
  if (!fp) {
    Rcpp::stop("Could not open BAM file '%s': %s", path,
               errno != 0 ? std::strerror(errno) : "unknown htslib error");
  }
  HeaderPtr hdr(sam_hdr_read(fp.get()));
  if (!hdr) {
    Rcpp::stop("Could not read the BAM header of '%s'", path);
  }
  IndexPtr idx(sam_index_load(fp.get(), path.c_str()));
  if (!idx) {
    Rcpp::stop(
        "Could not load an index for BAM '%s'; build a CSI index first", path);
  }

  const hts_pos_t n_windows =
      static_cast<hts_pos_t>(n_contigs) * n_per_chromosome;
  if (n_windows > static_cast<hts_pos_t>(INT32_MAX)) {
    Rcpp::stop("The analysis window set is too large for per-sample counting");
  }
  Rcpp::IntegerMatrix counts(n_samples, static_cast<int>(n_windows));
  Rcpp::IntegerMatrix covered(n_samples, static_cast<int>(n_windows));
  Rcpp::IntegerMatrix mapq_high(n_samples, static_cast<int>(n_windows));
  std::fill(counts.begin(), counts.end(), 0);
  std::fill(covered.begin(), covered.end(), 0);
  std::fill(mapq_high.begin(), mapq_high.end(), 0);

  BamPtr record(bam_init1());
  if (!record) {
    Rcpp::stop("Could not allocate a BAM record");
  }

  std::vector<std::vector<uint64_t>> bits(n_samples);
  int inferred_read_length = read_length > 0 ? read_length : 0;

  for (int contig = 0; contig < n_contigs; ++contig) {
    const std::string name = Rcpp::as<std::string>(contig_names[contig]);
    const double expected_length = contig_lengths[contig];
    const int tid = sam_hdr_name2tid(hdr.get(), name.c_str());
    if (tid < 0) {
      Rcpp::stop("Contig '%s' is not present in the BAM header of '%s'", name,
                 path);
    }
    const hts_pos_t header_length = sam_hdr_tid2len(hdr.get(), tid);
    if (header_length != static_cast<hts_pos_t>(expected_length)) {
      Rcpp::stop(
          "Contig '%s' is %lld bp in BAM '%s' but %lld bp in the analysis "
          "windows",
          name, static_cast<long long>(header_length), path,
          static_cast<long long>(expected_length));
    }
    if (header_length <= 0) {
      Rcpp::stop("Contig '%s' in BAM '%s' has no length", name, path);
    }

    for (int window = 0; window < n_per_chromosome; ++window) {
      // Match makeAnalysisWindows() exactly: floor(i * L / n) with integer
      // arithmetic, not a floating-point division.
      const hts_pos_t start =
          (static_cast<hts_pos_t>(window) * header_length) /
          n_per_chromosome;
      const hts_pos_t stop =
          (static_cast<hts_pos_t>(window + 1) * header_length) /
          n_per_chromosome;
      if (stop <= start) {
        continue;
      }
      const size_t n_words =
          static_cast<size_t>((stop - start + 63) / 64);
      for (int sample = 0; sample < n_samples; ++sample) {
        bits[sample].assign(n_words, 0);
      }

      ItrPtr itr(sam_itr_queryi(idx.get(), tid, start, stop));
      if (!itr) {
        Rcpp::stop("Could not iterate contig '%s' in BAM '%s'", name, path);
      }
      int ret = 0;
      while ((ret = sam_itr_next(fp.get(), itr.get(), record.get())) >= 0) {
        const bam1_t* b = record.get();
        if ((b->core.flag & exclude_flags) != 0) {
          continue;
        }
        if (b->core.qual < mapq) {
          continue;
        }
        if (b->core.tid < 0) {
          continue;
        }

        int sample = 0;
        if (!single_sample) {
          const uint8_t* rg_aux = bam_aux_get(b, "RG");
          if (rg_aux == nullptr) {
            Rcpp::stop(
                "Read '%s' in BAM '%s' has no RG tag; multi-sample BAMs need "
                "an RG tag on every read",
                bam_get_qname(b), path);
          }
          const char* rg = bam_aux2Z(rg_aux);
          if (rg == nullptr) {
            Rcpp::stop(
                "Read '%s' in BAM '%s' has a malformed RG tag; RG must be a "
                "string naming an @RG ID",
                bam_get_qname(b), path);
          }
          const auto found = rg_to_sample.find(std::string(rg));
          if (found == rg_to_sample.end()) {
            Rcpp::stop(
                "Read '%s' in BAM '%s' references read group '%s', which is "
                "not declared in the @RG header lines",
                bam_get_qname(b), path, rg);
          }
          sample = found->second;
        }

        const uint32_t* cigar = bam_get_cigar(b);
        const int n_cigar = b->core.n_cigar;
        hts_pos_t ref_pos = b->core.pos;
        bool overlaps = false;
        int aligned_query = 0;
        for (int k = 0; k < n_cigar; ++k) {
          const int op = bam_cigar_op(cigar[k]);
          const int len = bam_cigar_oplen(cigar[k]);
          const bool covers = op == BAM_CMATCH || op == BAM_CEQUAL ||
                              op == BAM_CDIFF ||
                              (count_deletions && op == BAM_CDEL);
          if (covers) {
            const hts_pos_t block_start = ref_pos;
            const hts_pos_t block_stop = ref_pos + len;
            if (block_stop > start && block_start < stop) {
              overlaps = true;
              const hts_pos_t clipped_start = std::max(block_start, start);
              const hts_pos_t clipped_stop = std::min(block_stop, stop);
              const size_t offset =
                  static_cast<size_t>(clipped_start - start);
              size_t first = offset / 64;
              const size_t last =
                  static_cast<size_t>(clipped_stop - start - 1) / 64;
              const uint64_t head =
                  ~uint64_t{0} << (offset % 64);
              const size_t tail_shift =
                  63 - static_cast<size_t>(clipped_stop - start - 1) % 64;
              const uint64_t tail = ~uint64_t{0} >> tail_shift;
              if (first == last) {
                bits[sample][first] |=
                    head & tail;
              } else {
                bits[sample][first] |= head;
                for (size_t w = first + 1; w < last; ++w) {
                  bits[sample][w] = ~uint64_t{0};
                }
                bits[sample][last] |= tail;
              }
            }
          }
          if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            aligned_query += len;
          }
          if ((bam_cigar_type(op) & 2) != 0) {
            ref_pos += len;
          }
        }
        if (!overlaps) {
          continue;
        }
        if (inferred_read_length <= 0 && aligned_query > 0) {
          inferred_read_length = aligned_query;
        }
        const hts_pos_t column =
            static_cast<hts_pos_t>(contig) * n_per_chromosome + window;
        int& value = counts(sample, static_cast<int>(column));
        if (value == INT32_MAX) {
          Rcpp::stop("Read counts exceed R's integer range for sample '%s'",
                     Rcpp::as<std::string>(sample_names[sample]));
        }
        ++value;
        if (b->core.qual >= 30) {
          int& high = mapq_high(sample, static_cast<int>(column));
          if (high == INT32_MAX) {
            Rcpp::stop(
                "High-MAPQ read counts exceed R's integer range for sample "
                "'%s'",
                Rcpp::as<std::string>(sample_names[sample]));
          }
          ++high;
        }
      }
      if (ret < -1) {
        Rcpp::stop("Could not read alignments from contig '%s' in BAM '%s'",
                   name, path);
      }

      for (int sample = 0; sample < n_samples; ++sample) {
        uint64_t observed = 0;
        for (size_t w = 0; w < bits[sample].size(); ++w) {
          observed += static_cast<uint64_t>(
              __builtin_popcountll(bits[sample][w]));
        }
        const hts_pos_t column =
            static_cast<hts_pos_t>(contig) * n_per_chromosome + window;
        if (observed > INT32_MAX) {
          Rcpp::stop("Covered bases exceed R's integer range for sample '%s'",
                     Rcpp::as<std::string>(sample_names[sample]));
        }
        covered(sample, static_cast<int>(column)) =
            static_cast<int>(observed);
      }
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("counts") = counts,
      Rcpp::Named("covered") = covered,
      Rcpp::Named("mapq_high") = mapq_high,
      Rcpp::Named("read_length") = inferred_read_length > 0
                                       ? static_cast<double>(inferred_read_length)
                                       : NA_REAL);
}

// Count reads in arbitrary half-open ranges.
//
// The ranges arrive in GRanges order. `range_contig` is a 1-based index into
// `contig_names`; starts are 0-based and ends are exclusive. Each range gets
// its own index iterator, so overlapping ranges are counted independently.
// Reads are assigned by their 5' end: a forward read uses its leftmost
// coordinate and a reverse read uses the last aligned base. `exclude_flags` follows samtools `view -F`: a read is dropped
// when any mask bit is set.
// [[Rcpp::export]]
Rcpp::IntegerVector bam_count_ranges_cpp(
    const std::string& path,
    Rcpp::CharacterVector contig_names,
    Rcpp::NumericVector contig_lengths,
    Rcpp::IntegerVector range_contig,
    Rcpp::NumericVector range_start,
    Rcpp::NumericVector range_end,
    int mapq,
    int exclude_flags) {
  const int n_contigs = contig_names.size();
  const R_xlen_t n_ranges = range_contig.size();
  if (n_contigs < 1) {
    Rcpp::stop("At least one contig is required for read counting");
  }
  if (contig_lengths.size() != n_contigs) {
    Rcpp::stop("contig_names and contig_lengths must have the same length");
  }
  if (range_start.size() != n_ranges || range_end.size() != n_ranges) {
    Rcpp::stop("range_contig, range_start and range_end must have the same length");
  }
  if (n_ranges == 0) {
    return Rcpp::IntegerVector(0);
  }

  HtsFilePtr fp(hts_open(path.c_str(), "r"));
  if (!fp) {
    Rcpp::stop("Could not open BAM file '%s': %s", path,
               errno != 0 ? std::strerror(errno) : "unknown htslib error");
  }
  HeaderPtr hdr(sam_hdr_read(fp.get()));
  if (!hdr) {
    Rcpp::stop("Could not read the BAM header of '%s'", path);
  }
  IndexPtr idx(sam_index_load(fp.get(), path.c_str()));
  if (!idx) {
    Rcpp::stop(
        "Could not load an index for BAM '%s'; build a CSI index first", path);
  }

  std::vector<int> tids(n_contigs, -1);
  std::vector<hts_pos_t> lengths(n_contigs, 0);
  for (int contig = 0; contig < n_contigs; ++contig) {
    const std::string name = Rcpp::as<std::string>(contig_names[contig]);
    const double expected_length = contig_lengths[contig];
    if (!std::isfinite(expected_length) || expected_length <= 0) {
      Rcpp::stop("Contig '%s' has an invalid length", name);
    }
    const int tid = sam_hdr_name2tid(hdr.get(), name.c_str());
    if (tid < 0) {
      Rcpp::stop("Contig '%s' is not present in the BAM header of '%s'", name,
                 path);
    }
    const hts_pos_t header_length = sam_hdr_tid2len(hdr.get(), tid);
    if (header_length != static_cast<hts_pos_t>(expected_length)) {
      Rcpp::stop(
          "Contig '%s' is %lld bp in BAM '%s' but %lld bp in the count ranges",
          name, static_cast<long long>(header_length), path,
          static_cast<long long>(expected_length));
    }
    if (header_length <= 0) {
      Rcpp::stop("Contig '%s' in BAM '%s' has no length", name, path);
    }
    tids[contig] = tid;
    lengths[contig] = header_length;
  }

  Rcpp::IntegerVector counts(n_ranges);
  std::fill(counts.begin(), counts.end(), 0);
  BamPtr record(bam_init1());
  if (!record) {
    Rcpp::stop("Could not allocate a BAM record");
  }

  for (R_xlen_t i = 0; i < n_ranges; ++i) {
    const int contig = range_contig[i] - 1;
    if (contig < 0 || contig >= n_contigs) {
      Rcpp::stop("range_contig must point into contig_names");
    }
    const double start_value = range_start[i];
    const double end_value = range_end[i];
    if (!std::isfinite(start_value) || !std::isfinite(end_value) ||
        start_value != std::floor(start_value) ||
        end_value != std::floor(end_value)) {
      Rcpp::stop("Range %lld has non-integer coordinates",
                 static_cast<long long>(i + 1));
    }
    const hts_pos_t start = static_cast<hts_pos_t>(start_value);
    const hts_pos_t end = static_cast<hts_pos_t>(end_value);
    if (start < 0 || end <= start || end > lengths[contig]) {
      Rcpp::stop(
          "Range %lld is outside contig '%s' in BAM '%s'",
          static_cast<long long>(i + 1),
          Rcpp::as<std::string>(contig_names[contig]), path);
    }

    ItrPtr itr(sam_itr_queryi(idx.get(), tids[contig], start, end));
    if (!itr) {
      Rcpp::stop(
          "Could not iterate BAM '%s' over range %lld on contig '%s'", path,
          static_cast<long long>(i + 1),
          Rcpp::as<std::string>(contig_names[contig]));
    }
    int ret = 0;
    while ((ret = sam_itr_next(fp.get(), itr.get(), record.get())) >= 0) {
      const bam1_t* b = record.get();
      if ((b->core.flag & exclude_flags) != 0) {
        continue;
      }
      if (b->core.qual < mapq) {
        continue;
      }
      if (b->core.tid < 0) {
        continue;
      }
      const hts_pos_t pos0 = (b->core.flag & BAM_FREVERSE) != 0
                                 ? bam_endpos(b) - 1
                                 : b->core.pos;
      if (pos0 < start || pos0 >= end) {
        continue;
      }
      if (counts[i] == INT32_MAX) {
        Rcpp::stop("Read count exceeds R's integer range for range %lld",
                   static_cast<long long>(i + 1));
      }
      ++counts[i];
    }
    if (ret < -1) {
      Rcpp::stop(
          "Could not read alignments from BAM '%s' over range %lld", path,
          static_cast<long long>(i + 1));
    }
  }
  return counts;
}
