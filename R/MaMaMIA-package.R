#' MaMaMIA: breadth-based introgression analysis
#'
#' MaMaMIA turns one or more coordinate-sorted BAM files into donor/recipient
#' introgression calls. It parses BAM headers with Rhtslib, matches them against
#' a YAML assembly registry, builds CSI indexes, tiles contigs into exact
#' windows, reads GC content through the vendored gcrank library, measures
#' breadth and counts in one htslib pass, scores donor/recipient homoeologs,
#' applies window masks, segments with CBS and calibrates thresholds by segment
#' length.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib MaMaMIA, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom Rhtslib pkgconfig
#' @importFrom rlang .data
NULL
