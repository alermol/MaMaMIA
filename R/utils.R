## Small internal utilities shared by the pipeline.

#' Stable hash of an R object
#'
#' Serialises an object to a temporary file and returns its MD5. Used for
#' provenance and mask/calibration versions.
#' @noRd
hash_object <- function(x) {
    path <- tempfile()
    on.exit(unlink(path), add = TRUE)
    saveRDS(x, path, version = 2)
    unname(tools::md5sum(path))
}

#' Null-coalescing helper
#' @noRd
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}
