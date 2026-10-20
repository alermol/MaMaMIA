// htslib helpers used by the BAM-first workflow.
//
// The package talks to BAM files through Rhtslib directly (D8): no
// Rsamtools dependency. This file only reads headers and builds CSI indexes;
// read counting stays in bamsignals.

#include <Rcpp.h>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <cerrno>
#include <cstring>
#include <memory>
#include <string>
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

using HtsFilePtr = std::unique_ptr<htsFile, HtsFileCloser>;
using HeaderPtr = std::unique_ptr<sam_hdr_t, HeaderDeleter>;

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
  if (min_shift < 0) {
    Rcpp::stop("min_shift must be non-negative");
  }
  if (threads < 1) {
    Rcpp::stop("threads must be at least 1");
  }
  if (index_path.empty()) {
    Rcpp::stop("An explicit CSI index path is required");
  }

  const int ret = sam_index_build3(path.c_str(), index_path.c_str(), min_shift,
                                   threads);
  if (ret != 0) {
    Rcpp::stop(
        "Could not build a CSI index for '%s' (htslib error code %d); the BAM "
        "must be coordinate-sorted and readable, and the output directory "
        "must be writable",
        path, ret);
  }
  return true;
}
