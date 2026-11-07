// htslib helpers used by the BAM analysis workflow.
//
// The application talks to BAM files through Rhtslib directly. This file reads
// headers, builds CSI indexes and measures reads. bam_metrics_rg_cpp() serves
// the pipeline; bam_count_ranges_cpp() serves the public countReads() helper
// and arbitrary GRanges.

#include <Rcpp.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
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
        "files are accepted (convert SAM/CRAM to BAM before uploading)",
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

// Median polish for the cohort score matrix.
//
// Alternating row and column medians remove the shared window effect and the
// per-sample offset. This is the C++ equivalent of the previous R loop and is
// used by medianPolish().
// [[Rcpp::export]]
Rcpp::List median_polish_cpp(Rcpp::NumericMatrix x, int max_iter, double tol,
                             bool na_rm) {
  const int n_rows = x.nrow();
  const int n_cols = x.ncol();
  if (max_iter < 1) {
    Rcpp::stop("max_iter must be a positive integer");
  }
  std::vector<double> residuals(x.begin(), x.end());
  std::vector<double> row_effects(n_rows, 0.0);
  std::vector<double> column_effects(n_cols, 0.0);
  std::vector<double> row_medians(n_rows, 0.0);
  std::vector<double> column_medians(n_cols, 0.0);
  std::vector<double> values;
  values.reserve(static_cast<size_t>(std::max(n_rows, n_cols)));

  const auto median_of = [&](std::vector<double>& input) {
    if (input.empty()) {
      return 0.0;
    }
    std::sort(input.begin(), input.end());
    const size_t n = input.size();
    if ((n & 1U) != 0U) {
      return input[n / 2U];
    }
    return (input[n / 2U - 1U] + input[n / 2U]) / 2.0;
  };

  double mu = 0.0;
  bool converged = false;
  int iteration = 0;

  for (int step = 1; step <= max_iter; ++step) {
    iteration = step;
    const double old = [&]() {
      double total = std::abs(mu);
      for (double value : row_effects) {
        total += std::abs(value);
      }
      for (double value : column_effects) {
        total += std::abs(value);
      }
      return total;
    }();

    for (int row = 0; row < n_rows; ++row) {
      values.clear();
      bool missing = false;
      for (int column = 0; column < n_cols; ++column) {
        const double value = residuals[static_cast<size_t>(row) +
                                       static_cast<size_t>(n_rows) * column];
        if (!std::isnan(value)) {
          values.push_back(value);
        } else if (!na_rm) {
          missing = true;
        }
      }
      double median = missing ? 0.0 : median_of(values);
      if (!std::isfinite(median)) {
        median = 0.0;
      }
      row_medians[row] = median;
      row_effects[row] += median;
      for (int column = 0; column < n_cols; ++column) {
        const size_t index = static_cast<size_t>(row) +
                             static_cast<size_t>(n_rows) * column;
        if (std::isfinite(residuals[index])) {
          residuals[index] -= median;
        }
      }
    }

    for (int column = 0; column < n_cols; ++column) {
      values.clear();
      bool missing = false;
      for (int row = 0; row < n_rows; ++row) {
        const size_t index = static_cast<size_t>(row) +
                             static_cast<size_t>(n_rows) * column;
        const double value = residuals[index];
        if (!std::isnan(value)) {
          values.push_back(value);
        } else if (!na_rm) {
          missing = true;
        }
      }
      double median = missing ? 0.0 : median_of(values);
      if (!std::isfinite(median)) {
        median = 0.0;
      }
      column_medians[column] = median;
      column_effects[column] += median;
      for (int row = 0; row < n_rows; ++row) {
        const size_t index = static_cast<size_t>(row) +
                             static_cast<size_t>(n_rows) * column;
        if (std::isfinite(residuals[index])) {
          residuals[index] -= median;
        }
      }
    }

    values.clear();
    for (double value : row_medians) {
      values.push_back(value);
    }
    for (double value : column_medians) {
      values.push_back(value);
    }
    mu += median_of(values);

    double current = std::abs(mu);
    for (double value : row_effects) {
      current += std::abs(value);
    }
    for (double value : column_effects) {
      current += std::abs(value);
    }
    if (std::isfinite(current) && std::abs(old - current) < tol) {
      converged = true;
      break;
    }
  }

  Rcpp::NumericMatrix out(n_rows, n_cols);
  std::copy(residuals.begin(), residuals.end(), out.begin());
  return Rcpp::List::create(
      Rcpp::Named("residuals") = out,
      Rcpp::Named("row") = row_effects,
      Rcpp::Named("column") = column_effects,
      Rcpp::Named("mu") = mu,
      Rcpp::Named("iterations") = iteration,
      Rcpp::Named("converged") = converged);
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
// `threads` is the number of counting workers, not the htslib decompression
// pool. Each worker opens its own BAM handle and index and walks a contiguous
// block of windows. Windows are split across workers, and each worker keeps
// only a window-sized bitmap per sample, so memory scales with window length
// instead of contig length.
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
    bool count_deletions,
    int threads) {
  const int n_contigs = contig_names.size();
  const int n_samples = sample_names.size();
  const bool map_read_groups = read_group_ids.size() > 0;
  if (threads < 1) {
    Rcpp::stop("threads must be at least 1");
  }
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

  // Extract every R object before any worker starts: workers must not call
  // into the R API.
  std::vector<std::string> names(n_contigs);
  std::vector<double> expected_lengths(n_contigs);
  for (int contig = 0; contig < n_contigs; ++contig) {
    names[contig] = Rcpp::as<std::string>(contig_names[contig]);
    expected_lengths[contig] = contig_lengths[contig];
  }

  std::unordered_map<std::string, int> rg_to_sample;
  for (int i = 0; i < read_group_ids.size(); ++i) {
    const std::string id = Rcpp::as<std::string>(read_group_ids[i]);
    const int sample = sample_index[i];
    if (id.empty()) {
      Rcpp::stop("Read group IDs must not be empty");
    }
    // 0 is the sentinel for an @RG line that belongs to a sample the caller
    // did not request; its reads are skipped.
    if (sample < 0 || sample > n_samples) {
      Rcpp::stop("sample_index must point into sample_names");
    }
    if (!rg_to_sample.emplace(id, sample - 1).second) {
      Rcpp::stop("Duplicate read group ID '%s'", id);
    }
  }

  // Validate the header and resolve target ids once. Each worker then uses the
  // same tids with its own header and index.
  std::vector<int> tids(n_contigs, -1);
  std::vector<hts_pos_t> contig_sizes(n_contigs, 0);
  {
    HtsFilePtr fp(hts_open(path.c_str(), "r"));
    if (!fp) {
      Rcpp::stop("Could not open BAM file '%s': %s", path,
                 errno != 0 ? std::strerror(errno) : "unknown htslib error");
    }
    HeaderPtr hdr(sam_hdr_read(fp.get()));
    if (!hdr) {
      Rcpp::stop("Could not read the BAM header of '%s'", path);
    }
    for (int contig = 0; contig < n_contigs; ++contig) {
      const double expected_length = expected_lengths[contig];
      if (!std::isfinite(expected_length) || expected_length <= 0) {
        Rcpp::stop("Contig '%s' has an invalid length", names[contig]);
      }
      const int tid = sam_hdr_name2tid(hdr.get(), names[contig].c_str());
      if (tid < 0) {
        Rcpp::stop("Contig '%s' is not present in the BAM header of '%s'",
                   names[contig], path);
      }
      const hts_pos_t header_length = sam_hdr_tid2len(hdr.get(), tid);
      if (header_length != static_cast<hts_pos_t>(expected_length)) {
        Rcpp::stop(
            "Contig '%s' is %lld bp in BAM '%s' but %lld bp in the analysis "
            "windows",
            names[contig], static_cast<long long>(header_length), path,
            static_cast<long long>(expected_length));
      }
      if (header_length <= 0) {
        Rcpp::stop("Contig '%s' in BAM '%s' has no length", names[contig],
                   path);
      }
      tids[contig] = tid;
      contig_sizes[contig] = header_length;
    }
    IndexPtr idx(sam_index_load(fp.get(), path.c_str()));
    if (!idx) {
      Rcpp::stop(
          "Could not load an index for BAM '%s'; build a CSI index first",
          path);
    }
  }

  const hts_pos_t n_windows =
      static_cast<hts_pos_t>(n_contigs) * n_per_chromosome;
  if (n_windows > static_cast<hts_pos_t>(INT32_MAX)) {
    Rcpp::stop("The analysis window set is too large for per-sample counting");
  }
  Rcpp::IntegerMatrix counts(n_samples, static_cast<int>(n_windows));
  Rcpp::IntegerMatrix covered(n_samples, static_cast<int>(n_windows));
  Rcpp::IntegerMatrix mapq_high(n_samples, static_cast<int>(n_windows));
  // R's allocMatrix() returns zero-filled memory; filling again here would
  // touch n_samples * n_windows cells three extra times.
  int* counts_ptr = counts.begin();
  int* covered_ptr = covered.begin();
  int* mapq_ptr = mapq_high.begin();

  struct MetricsTask {
    int contig;
    int window_begin;
    int window_end;
  };
  const int window_block = std::max(1, std::min(32, n_per_chromosome));
  std::vector<MetricsTask> tasks;
  tasks.reserve(
      static_cast<size_t>(n_contigs) *
      static_cast<size_t>((n_per_chromosome + window_block - 1) /
                          window_block));
  for (int contig = 0; contig < n_contigs; ++contig) {
    for (int window = 0; window < n_per_chromosome; window += window_block) {
      tasks.push_back(MetricsTask{
          contig, window,
          std::min(window + window_block, n_per_chromosome)});
    }
  }

  std::vector<int> task_inferred_read_length(tasks.size(), 0);
  std::atomic<size_t> next_task{0};
  std::atomic<bool> failed{false};
  std::mutex error_mutex;
  std::string error_message;
  const auto set_error = [&](const std::string& message) {
    bool expected = false;
    if (failed.compare_exchange_strong(expected, true)) {
      std::lock_guard<std::mutex> lock(error_mutex);
      error_message = message;
    }
  };

  const auto run_worker = [&]() {
    HtsFilePtr worker_file(hts_open(path.c_str(), "r"));
    if (!worker_file) {
      set_error("Could not open BAM file '" + path + "': " +
                (errno != 0 ? std::strerror(errno)
                            : "unknown htslib error"));
      return;
    }
    HeaderPtr worker_header(sam_hdr_read(worker_file.get()));
    if (!worker_header) {
      set_error("Could not read the BAM header of '" + path + "'");
      return;
    }
    IndexPtr worker_index(
        sam_index_load(worker_file.get(), path.c_str()));
    if (!worker_index) {
      set_error("Could not load an index for BAM '" + path +
                "'; build a CSI index first");
      return;
    }
    BamPtr record(bam_init1());
    if (!record) {
      set_error("Could not allocate a BAM record");
      return;
    }

    std::vector<std::vector<uint64_t>> bits(n_samples);
    std::vector<size_t> first_word(n_samples, 0);
    std::vector<size_t> last_word(n_samples, 0);
    std::vector<int> window_count(n_samples, 0);
    std::vector<int> window_mapq(n_samples, 0);
    int inferred_read_length = read_length > 0 ? read_length : 0;

    while (!failed.load(std::memory_order_relaxed)) {
      const size_t task_index =
          next_task.fetch_add(1, std::memory_order_relaxed);
      if (task_index >= tasks.size()) {
        break;
      }
      const MetricsTask& task = tasks[task_index];
      const int contig = task.contig;
      const hts_pos_t contig_length = contig_sizes[contig];

      for (int window = task.window_begin; window < task.window_end;
           ++window) {
        if (failed.load(std::memory_order_relaxed)) {
          return;
        }
        const hts_pos_t start =
            (static_cast<hts_pos_t>(window) * contig_length) /
            n_per_chromosome;
        const hts_pos_t stop =
            (static_cast<hts_pos_t>(window + 1) * contig_length) /
            n_per_chromosome;
        if (stop <= start) {
          continue;
        }
        const size_t window_words =
            static_cast<size_t>((stop - start + 63) / 64);
        for (int sample = 0; sample < n_samples; ++sample) {
          if (bits[sample].size() < window_words) {
            bits[sample].resize(window_words, 0);
          }
          first_word[sample] = window_words;
          last_word[sample] = 0;
          window_count[sample] = 0;
          window_mapq[sample] = 0;
        }

        ItrPtr itr(sam_itr_queryi(
            worker_index.get(), tids[contig], start, stop));
        if (!itr) {
          set_error("Could not iterate contig '" + names[contig] +
                    "' in BAM '" + path + "'");
          return;
        }
        int ret = 0;
        long long seen = 0;
        while ((ret = sam_itr_next(
                    worker_file.get(), itr.get(), record.get())) >= 0) {
          if ((++seen & 0xFFFFF) == 0 &&
              failed.load(std::memory_order_relaxed)) {
            return;
          }
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
          if (map_read_groups) {
            const uint8_t* rg_aux = bam_aux_get(b, "RG");
            if (rg_aux == nullptr) {
              set_error(
                  "Read '" + std::string(bam_get_qname(b)) +
                  "' in BAM '" + path +
                  "' has no RG tag; multi-sample BAMs need an RG tag on "
                  "every read");
              return;
            }
            const char* rg = bam_aux2Z(rg_aux);
            if (rg == nullptr) {
              set_error(
                  "Read '" + std::string(bam_get_qname(b)) +
                  "' in BAM '" + path +
                  "' has a malformed RG tag; RG must be a string naming an "
                  "@RG ID");
              return;
            }
            const auto found = rg_to_sample.find(std::string(rg));
            if (found == rg_to_sample.end()) {
              set_error(
                  "Read '" + std::string(bam_get_qname(b)) +
                  "' in BAM '" + path + "' references read group '" + rg +
                  "', which is not declared in the @RG header lines");
              return;
            }
            sample = found->second;
            if (sample < 0) {
              // This read group belongs to a sample the caller did not
              // request.
              continue;
            }
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
                const hts_pos_t clipped_start =
                    std::max(block_start, start);
                const hts_pos_t clipped_stop =
                    std::min(block_stop, stop);
                const size_t offset =
                    static_cast<size_t>(clipped_start - start);
                const size_t first = offset / 64;
                const size_t last =
                    static_cast<size_t>(clipped_stop - start - 1) / 64;
                const uint64_t head = ~uint64_t{0} << (offset % 64);
                const size_t tail_shift =
                    63 -
                    static_cast<size_t>(clipped_stop - start - 1) % 64;
                const uint64_t tail = ~uint64_t{0} >> tail_shift;
                if (first < first_word[sample]) {
                  first_word[sample] = first;
                }
                if (last > last_word[sample]) {
                  last_word[sample] = last;
                }
                if (first == last) {
                  bits[sample][first] |= head & tail;
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
          ++window_count[sample];
          if (window_count[sample] == INT32_MAX) {
            set_error(
                "Read counts exceed R's integer range for a window in BAM '" +
                path + "'");
            return;
          }
          if (b->core.qual >= 30) {
            ++window_mapq[sample];
            if (window_mapq[sample] == INT32_MAX) {
              set_error(
                  "High-MAPQ read counts exceed R's integer range for a "
                  "window in BAM '" +
                  path + "'");
              return;
            }
          }
        }
        if (ret < -1) {
          set_error("Could not read alignments from contig '" +
                    names[contig] + "' in BAM '" + path + "'");
          return;
        }

        const hts_pos_t column =
            static_cast<hts_pos_t>(contig) * n_per_chromosome + window;
        for (int sample = 0; sample < n_samples; ++sample) {
          uint64_t observed = 0;
          const size_t first = first_word[sample];
          const size_t last = last_word[sample];
          if (first <= last) {
            for (size_t w = first; w <= last; ++w) {
              observed += static_cast<uint64_t>(
                  __builtin_popcountll(bits[sample][w]));
            }
            std::fill(
                bits[sample].begin() + first,
                bits[sample].begin() + last + 1, 0);
          }
          if (observed > INT32_MAX) {
            set_error(
                "Covered bases exceed R's integer range for a window in BAM '" +
                path + "'");
            return;
          }
          const size_t index =
              static_cast<size_t>(sample) +
              static_cast<size_t>(n_samples) *
                  static_cast<size_t>(column);
          counts_ptr[index] = window_count[sample];
          covered_ptr[index] = static_cast<int>(observed);
          mapq_ptr[index] = window_mapq[sample];
        }
      }
      task_inferred_read_length[task_index] = inferred_read_length;
    }
  };

  unsigned hardware = std::thread::hardware_concurrency();
  if (hardware == 0) {
    hardware = static_cast<unsigned>(threads);
  }
  const size_t worker_count = std::max<size_t>(
      1, std::min({static_cast<size_t>(threads), tasks.size(),
                   static_cast<size_t>(hardware)}));

  if (worker_count <= 1) {
    run_worker();
  } else {
    std::vector<std::thread> pool;
    pool.reserve(worker_count);
    try {
      for (size_t worker = 0; worker < worker_count; ++worker) {
        pool.emplace_back(run_worker);
      }
    } catch (const std::exception& e) {
      set_error(std::string("Could not start counting threads: ") + e.what());
    }
    for (std::thread& worker : pool) {
      worker.join();
    }
  }

  if (failed.load()) {
    Rcpp::stop("%s", error_message);
  }
  Rcpp::checkUserInterrupt();

  int inferred_read_length = read_length > 0 ? read_length : 0;
  if (inferred_read_length <= 0) {
    for (size_t task = 0; task < task_inferred_read_length.size(); ++task) {
      if (task_inferred_read_length[task] > 0) {
        inferred_read_length = task_inferred_read_length[task];
        break;
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
    int exclude_flags,
    int threads) {
  const int n_contigs = contig_names.size();
  const R_xlen_t n_ranges = range_contig.size();
  if (threads < 1) {
    Rcpp::stop("threads must be at least 1");
  }
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
  if (threads > 1 && hts_set_threads(fp.get(), threads) != 0) {
    Rcpp::stop("Could not enable %d htslib worker thread(s) for '%s'", threads,
               path);
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
  // Rcpp vectors are zero-filled by allocVector(); no extra fill needed.
  BamPtr record(bam_init1());
  if (!record) {
    Rcpp::stop("Could not allocate a BAM record");
  }

  for (R_xlen_t i = 0; i < n_ranges; ++i) {
    if ((i & 0x3F) == 0) {
      Rcpp::checkUserInterrupt();
    }
    if (range_contig[i] == NA_INTEGER) {
      Rcpp::stop("range_contig must not be NA");
    }
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
    long long seen = 0;
    while ((ret = sam_itr_next(fp.get(), itr.get(), record.get())) >= 0) {
      if ((++seen & 0xFFFFF) == 0) {
        Rcpp::checkUserInterrupt();
      }
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
