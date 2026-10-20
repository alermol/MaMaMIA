// Rcpp wrapper around the vendored gcrank library.
//
// An index is opened once per session and kept alive through an external
// pointer; the R side owns an LRU cache keyed by the index directory. Queries
// use contig-relative, half-open, 0-based coordinates, exactly like BED.

#include <Rcpp.h>

#include "gcrank/error.h"
#include "gcrank/gc_rank.h"

#include <algorithm>
#include <string>
#include <vector>

namespace {

gcrank::GcRank* checked_index(SEXP xptr) {
  Rcpp::XPtr<gcrank::GcRank> index(xptr);
  if (index.get() == nullptr) {
    Rcpp::stop("The gcrank index handle is no longer valid");
  }
  return index.get();
}

}  // namespace

// [[Rcpp::export]]
SEXP gcrank_open_cpp(const std::string& index_dir, bool mmap, bool validate_crc) {
  gcrank::GcRank* index = nullptr;
  try {
    index = new gcrank::GcRank();
    gcrank::OpenOptions options;
    options.mmap = mmap;
    options.validate_crc = validate_crc;
    index->load(index_dir, options);
    return Rcpp::XPtr<gcrank::GcRank>(index, true);
  } catch (const std::exception& e) {
    delete index;
    Rcpp::stop("Could not open gcrank index '%s': %s", index_dir, e.what());
  }
  return R_NilValue;
}

// [[Rcpp::export]]
void gcrank_close_cpp(SEXP xptr) {
  gcrank::GcRank* index = nullptr;
  {
    Rcpp::XPtr<gcrank::GcRank> handle(xptr);
    index = handle.get();
  }
  if (index != nullptr) {
    R_ClearExternalPtr(xptr);
    delete index;
  }
  return;
}

// [[Rcpp::export]]
Rcpp::DataFrame gcrank_contigs_cpp(SEXP xptr) {
  gcrank::GcRank* index = checked_index(xptr);
  const std::vector<gcrank::ContigInfo>& contigs = index->contigs();

  Rcpp::CharacterVector names(contigs.size());
  Rcpp::NumericVector lengths(contigs.size());
  Rcpp::NumericVector offsets(contigs.size());
  for (std::size_t i = 0; i < contigs.size(); ++i) {
    names[i] = contigs[i].name;
    lengths[i] = static_cast<double>(contigs[i].length);
    offsets[i] = static_cast<double>(contigs[i].offset);
  }
  return Rcpp::DataFrame::create(
      Rcpp::Named("name") = names,
      Rcpp::Named("length") = lengths,
      Rcpp::Named("offset") = offsets);
}

// [[Rcpp::export]]
Rcpp::NumericVector gcrank_gc_cpp(SEXP xptr,
                                  Rcpp::CharacterVector contig,
                                  Rcpp::NumericVector start,
                                  Rcpp::NumericVector end,
                                  bool remove_ambiguous = false) {
  gcrank::GcRank* index = checked_index(xptr);

  const R_xlen_t n = std::max(
      contig.size(), std::max(start.size(), end.size()));
  if (n == 0) {
    return Rcpp::NumericVector(0);
  }
  if (contig.size() != 1 && contig.size() != n) {
    Rcpp::stop("`contig` must have length 1 or %lld", static_cast<long long>(n));
  }
  if (start.size() != 1 && start.size() != n) {
    Rcpp::stop("`start` must have length 1 or %lld", static_cast<long long>(n));
  }
  if (end.size() != 1 && end.size() != n) {
    Rcpp::stop("`end` must have length 1 or %lld", static_cast<long long>(n));
  }

  Rcpp::NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    const R_xlen_t ci = contig.size() == 1 ? 0 : i;
    const R_xlen_t si = start.size() == 1 ? 0 : i;
    const R_xlen_t ei = end.size() == 1 ? 0 : i;

    if (Rcpp::CharacterVector::is_na(contig[ci])) {
      Rcpp::stop("Contig name at position %lld is NA", static_cast<long long>(i + 1));
    }
    if (!R_finite(start[si]) || !R_finite(end[ei])) {
      Rcpp::stop("Window coordinates must be finite (position %lld)",
                 static_cast<long long>(i + 1));
    }
    if (start[si] < 0 || end[ei] < 0 ||
        start[si] > 9007199254740991.0 || end[ei] > 9007199254740991.0) {
      Rcpp::stop("Window coordinates are out of range (position %lld)",
                 static_cast<long long>(i + 1));
    }
    if (end[ei] <= start[si]) {
      Rcpp::stop("Window end must be greater than start (position %lld)",
                 static_cast<long long>(i + 1));
    }

    const std::string name = Rcpp::as<std::string>(contig[ci]);
    const std::uint64_t s = static_cast<std::uint64_t>(start[si]);
    const std::uint64_t e = static_cast<std::uint64_t>(end[ei]);

    try {
      out[i] = index->gc_fraction(name, s, e, remove_ambiguous);
    } catch (const std::exception& error) {
      Rcpp::stop("gcrank query failed for %s:%lld-%lld: %s", name,
                 static_cast<long long>(s), static_cast<long long>(e),
                 error.what());
    }
  }
  return out;
}
