## Small internal utilities shared by the pipeline.

#' Stable hash of an R object
#'
#' Serialises an object to a temporary file and returns its MD5. Used for
#' provenance and mask/calibration versions.
#' @noRd
hash_object <- function(x) {
    rlang::hash(x)
}

#' Null-coalescing helper
#' @noRd
`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

#' Apply over independent tasks with the configured worker budget
#'
#' Forks on Unix when more than one worker is allowed. Falls back to a serial
#' `lapply()` on Windows and when the caller keeps the default of one core.
#' Errors from forked workers come back as `try-error` values so the caller can
#' preserve the original condition message.
#' @noRd
.mamamia_lapply <- function(X, FUN, ..., workers = NULL) {
    if (is.null(workers)) {
        workers <- getOption("MaMaMIA.cores", 1L)
    }
    workers <- suppressWarnings(as.integer(workers))
    if (length(workers) != 1L || is.na(workers) || workers < 1L) {
        workers <- 1L
    }
    if (.Platform$OS.type == "unix" && workers > 1L && length(X) > 1L) {
        return(parallel::mclapply(
            X, FUN, ...,
            mc.cores = min(workers, length(X)),
            mc.preschedule = TRUE
        ))
    }
    lapply(X, FUN, ...)
}
