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

#' Cores available to this process
#'
#' `parallel::detectCores()` reports host CPUs even when a container has a
#' lower CPU quota. Read the cgroup limit so the default worker budget matches
#' the CPUs the process can actually use.
#' @noRd
mamamia_available_cores <- function() {
    detected <- parallel::detectCores()
    if (length(detected) != 1L || !is.numeric(detected) ||
        is.na(detected) || !is.finite(detected) || detected < 1) {
        detected <- 1
    }
    limit <- Inf
    cpu_max <- "/sys/fs/cgroup/cpu.max"
    if (file.exists(cpu_max)) {
        fields <- strsplit(
            readLines(cpu_max, n = 1L, warn = FALSE), "\\s+"
        )[[1L]]
        if (length(fields) >= 2L && !identical(fields[[1L]], "max")) {
            quota <- suppressWarnings(as.numeric(fields[[1L]]))
            period <- suppressWarnings(as.numeric(fields[[2L]]))
            if (is.finite(quota) && is.finite(period) && period > 0) {
                limit <- quota / period
            }
        }
    } else {
        quota_file <- "/sys/fs/cgroup/cpu/cpu.cfs_quota_us"
        period_file <- "/sys/fs/cgroup/cpu/cpu.cfs_period_us"
        if (file.exists(quota_file) && file.exists(period_file)) {
            quota <- suppressWarnings(as.numeric(
                readLines(quota_file, n = 1L, warn = FALSE)
            ))
            period <- suppressWarnings(as.numeric(
                readLines(period_file, n = 1L, warn = FALSE)
            ))
            if (is.finite(quota) && is.finite(period) && quota > 0 &&
                period > 0) {
                limit <- quota / period
            }
        }
    }
    if (is.finite(limit)) {
        detected <- min(detected, max(1, floor(limit)))
    }
    as.integer(max(1L, detected))
}

#' Default worker budget for BAM calculations
#'
#' Leaves one core for the operating system and the Shiny event loop. The
#' result is capped at the pipeline's 256-core limit. Callers can still pass an
#' explicit `cores` value or set `options(MaMaMIA.cores = n)` to override it.
#' @noRd
mamamia_default_cores <- function(available = mamamia_available_cores()) {
    if (length(available) != 1L || !is.numeric(available) ||
        is.na(available) || !is.finite(available) || available < 1) {
        return(1L)
    }
    as.integer(min(256L, max(1L, available - 1L)))
}

#' Apply over independent tasks with the configured worker budget
#'
#' Forks on Unix when more than one worker is allowed. Falls back to a serial
#' `lapply()` on Windows and when only one worker is available.
#' Errors from forked workers come back as `try-error` values so the caller can
#' preserve the original condition message.
#' @noRd
.mamamia_lapply <- function(X, FUN, ..., workers = NULL) {
    if (is.null(workers)) {
        workers <- getOption("MaMaMIA.cores", mamamia_default_cores())
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
