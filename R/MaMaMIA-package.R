#' MaMaMIA: BAM-first introgression analysis
#'
#' MaMaMIA turns one or more coordinate-sorted BAM files into pairwise
#' donor-recipient introgression calls. The BAM-first workflow parses BAM
#' headers with Rhtslib, matches them against a YAML assembly registry builds
#' CSI indexes, tiles contigs into exact windows, reads GC content through the
#' vendored gcrank library, counts reads with bamsignals, corrects GC bias,
#' optionally normalizes several samples to dosage, and finally segments the
#' pairwise difference profiles with CBS.
#'
#' @keywords internal
"_PACKAGE"

#' @useDynLib MaMaMIA, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom rlang .data
NULL
