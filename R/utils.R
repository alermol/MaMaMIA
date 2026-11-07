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

#' Compare key for natural ordering
#'
#' Pads every run of digits to a fixed width so that a lexicographic comparison
#' puts `sample2` before `sample10`. Non-numeric runs are lower-cased, so the
#' order does not depend on capitalisation.
#' @noRd
mamamia_natural_key <- function(x) {
    vapply(x, function(value) {
        if (is.na(value)) {
            return("")
        }
        parts <- strsplit(
            value,
            "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)",
            perl = TRUE
        )[[1L]]
        if (length(parts) == 0L) {
            return("")
        }
        digits <- grepl("^[0-9]+$", parts)
        if (any(digits)) {
            parts[digits] <- paste0(
                strrep("0", pmax(0L, 20L - nchar(parts[digits]))),
                parts[digits]
            )
        }
        paste(tolower(parts), collapse = "\u0001")
    }, character(1), USE.NAMES = FALSE)
}

#' Sort labels in natural order
#'
#' Sorts names containing numbers the way a person reads them: the digit runs
#' are compared as numbers and the text is compared case-insensitively.
#' Missing values stay at the end and ties keep their original order.
#' @noRd
mamamia_natural_sort <- function(x) {
    x <- as.character(x)
    if (length(x) < 2L) {
        return(x)
    }
    key <- mamamia_natural_key(x)
    key[is.na(x)] <- "\uffff"
    x[order(key, tolower(x), x, seq_along(x), na.last = TRUE)]
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

#' Memory limit available to this process
#'
#' `memory.max` (cgroup v2) and `memory.limit_in_bytes` (cgroup v1) cap a
#' container below the host's RAM. Bundle export forks one process per worker,
#' so the limit decides how many workers can run without swapping.
#' @noRd
mamamia_memory_limit <- function() {
    v2 <- "/sys/fs/cgroup/memory.max"
    if (file.exists(v2)) {
        value <- readLines(v2, n = 1L, warn = FALSE)
        if (length(value) == 1L && !identical(trimws(value), "max")) {
            limit <- suppressWarnings(as.numeric(value))
            if (is.finite(limit) && limit > 0) {
                return(limit)
            }
        }
        return(Inf)
    }
    v1 <- "/sys/fs/cgroup/memory/memory.limit_in_bytes"
    if (file.exists(v1)) {
        limit <- suppressWarnings(as.numeric(
            readLines(v1, n = 1L, warn = FALSE)
        ))
        ## cgroup v1 often reports a sentinel close to the 64-bit maximum.
        if (is.finite(limit) && limit > 0 && limit < 2^60) {
            return(limit)
        }
    }
    Inf
}

#' Memory-based cap on forked bundle workers
#'
#' Each bundle worker carries a copy-on-write view of the Shiny process heap.
#' Forking more workers than fit in the container memory limit pushes the run
#' into swap, which is far slower than using fewer workers. Reserve memory for
#' the app and the operating system, then allow roughly `per_worker` bytes per
#' worker.
#' @noRd
mamamia_memory_worker_cap <- function(
    memory_limit = mamamia_memory_limit(),
    reserve = 4 * 1024^3,
    per_worker = 2 * 1024^3
) {
    if (length(memory_limit) != 1L || !is.numeric(memory_limit) ||
        is.na(memory_limit) || memory_limit <= 0) {
        return(1L)
    }
    if (!is.finite(memory_limit)) {
        return(Inf)
    }
    usable <- memory_limit - reserve
    if (!is.finite(usable) || usable < per_worker) {
        return(1L)
    }
    as.integer(max(1, floor(usable / per_worker)))
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
