// htslib helpers used by the BAM-first workflow.
//
// The package talks to BAM files through Rhtslib directly (D8): no
// Rsamtools dependency. This file only reads headers and builds CSI indexes;
// read counting stays in bamsignals.

#include <Rcpp.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cerrno>
#include <climits>
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

// Minimal compile/load spike kept from T0.2. It proves that Rcpp, Rhtslib and
// the C++23 build flags produce a loadable shared library.
// [[Rcpp::export]]
int build_stub_cpp() {
  return 1;
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

// Count reads per sample in one BAM.
//
// Reads are assigned to the window containing their 5' end, matching
// bamsignals::bamCount(): a forward read uses its leftmost coordinate and a
// reverse read uses the last aligned base. `exclude_flags` follows samtools
// `view -F`: a read is dropped when any mask bit is set. Reads are assigned to
// a sample through their `RG` tag, which must name an `@RG ID` from the header.
// `sample_index` is 1-based and aligned to `read_group_ids`; the returned
// matrix has one row per `sample_names` entry and one column per window.
// [[Rcpp::export]]
Rcpp::IntegerMatrix bam_count_rg_cpp(
    const std::string& path,
    Rcpp::CharacterVector contig_names,
    Rcpp::NumericVector contig_lengths,
    int n_per_chromosome,
    Rcpp::CharacterVector read_group_ids,
    Rcpp::IntegerVector sample_index,
    Rcpp::CharacterVector sample_names,
    int mapq,
    int exclude_flags,
    bool allow_missing_rg) {
  const int n_contigs = contig_names.size();
  const int n_samples = sample_names.size();
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
    Rcpp::stop(
        "read_group_ids and sample_index must have the same length");
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
  std::fill(counts.begin(), counts.end(), 0);
  BamPtr record(bam_init1());
  if (!record) {
    Rcpp::stop("Could not allocate a BAM record");
  }

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

    ItrPtr itr(sam_itr_queryi(idx.get(), tid, 0, header_length));
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

      int sample = -1;
      const uint8_t* rg_aux = bam_aux_get(b, "RG");
      if (rg_aux != nullptr) {
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
      } else if (allow_missing_rg && n_samples == 1) {
        sample = 0;
      } else {
        Rcpp::stop(
            "Read '%s' in BAM '%s' has no RG tag; multi-sample BAMs need an "
            "RG tag on every read",
            bam_get_qname(b), path);
      }

      const hts_pos_t pos0 = (b->core.flag & BAM_FREVERSE) != 0
                                 ? bam_endpos(b) - 1
                                 : b->core.pos;
      if (pos0 < 0 || pos0 >= header_length) {
        continue;
      }
      const hts_pos_t window =
          (pos0 * n_per_chromosome) / header_length;
      if (window < 0 || window >= n_per_chromosome) {
        continue;
      }
      const hts_pos_t column =
          static_cast<hts_pos_t>(contig) * n_per_chromosome + window;
      int& value = counts(sample, static_cast<int>(column));
      if (value == INT32_MAX) {
        Rcpp::stop("Read counts exceed R's integer range for sample '%s'",
                   Rcpp::as<std::string>(sample_names[sample]));
      }
      ++value;
    }
    if (ret < -1) {
      Rcpp::stop("Could not read alignments from contig '%s' in BAM '%s'",
                 name, path);
    }
  }
  return counts;
}
