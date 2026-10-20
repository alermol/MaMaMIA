## GC content from a gcrank index through the vendored C++ library (T2.3).
##
## An index is opened once and kept in an LRU cache keyed by its normalized
## directory. With mmap = TRUE (the default) the index files are memory-mapped:
## resident memory grows only as pages are touched. Queries use contig-relative,
## half-open, 0-based coordinates, exactly like BED.

.MaMaMIAgcrankCache <- new.env(parent = emptyenv())
.MaMaMIAgcrankState <- new.env(parent = emptyenv())
.MaMaMIAgcrankState$counter <- 0L
.MaMaMIAgcrankCacheMax <- 2L

.gcrank_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAGcrankError", fmt, ...)
}

.gcrank_tick <- function() {
    .MaMaMIAgcrankState$counter <- .MaMaMIAgcrankState$counter + 1L
    .MaMaMIAgcrankState$counter
}

.gcrank_logical <- function(x, what) {
    if (length(x) != 1L || is.na(x) || !is.logical(x)) {
        .gcrank_error("`%s` must be true or false", what)
    }
    x
}

.validate_gcrank_index <- function(index_dir) {
    if (!is.character(index_dir) || length(index_dir) != 1L ||
        is.na(index_dir) || !nzchar(index_dir)) {
        .gcrank_error("`index_dir` must be a single non-empty directory path")
    }
    index_dir <- path.expand(index_dir)
    if (!dir.exists(index_dir)) {
        .gcrank_error("gcrank index does not exist: %s", index_dir)
    }
    if (!file.exists(file.path(index_dir, "contigs.tsv")) ||
        !file.exists(file.path(index_dir, "manifest.json"))) {
        .gcrank_error(
            "Directory '%s' is not a gcrank index: expected contigs.tsv and manifest.json",
            index_dir
        )
    }
    normalizePath(index_dir, winslash = "/", mustWork = TRUE)
}

.gcrank_evict <- function() {
    keys <- ls(envir = .MaMaMIAgcrankCache, all.names = TRUE)
    while (length(keys) > .MaMaMIAgcrankCacheMax) {
        last_used <- vapply(keys, function(key) {
            get(key, envir = .MaMaMIAgcrankCache, inherits = FALSE)$last_used
        }, numeric(1))
        victim <- keys[which.min(last_used)]
        entry <- get(victim, envir = .MaMaMIAgcrankCache, inherits = FALSE)
        try(gcrank_close_cpp(entry$xptr), silent = TRUE)
        rm(list = victim, envir = .MaMaMIAgcrankCache)
        keys <- setdiff(keys, victim)
    }
    invisible(NULL)
}

.gcrank_handle <- function(index_dir, mmap = TRUE, validate_crc = FALSE) {
    index_dir <- .validate_gcrank_index(index_dir)
    mmap <- .gcrank_logical(mmap, "mmap")
    validate_crc <- .gcrank_logical(validate_crc, "validate_crc")

    if (exists(index_dir, envir = .MaMaMIAgcrankCache, inherits = FALSE)) {
        entry <- get(index_dir, envir = .MaMaMIAgcrankCache, inherits = FALSE)
        entry$last_used <- .gcrank_tick()
        assign(index_dir, entry, envir = .MaMaMIAgcrankCache)
        return(entry$xptr)
    }

    xptr <- tryCatch(
        gcrank_open_cpp(index_dir, mmap, validate_crc),
        error = function(e) {
            .gcrank_error("%s", conditionMessage(e))
        }
    )
    assign(
        index_dir,
        list(
            xptr = xptr,
            mmap = mmap,
            validate_crc = validate_crc,
            last_used = .gcrank_tick()
        ),
        envir = .MaMaMIAgcrankCache
    )
    .gcrank_evict()
    xptr
}

#' Read contig metadata from a gcrank index
#'
#' @param index_dir Path to the `*.gcrank` index directory.
#'
#' @return A data frame with `name`, `length` and `offset`.
#' @export
gcrankContigs <- function(index_dir) {
    xptr <- .gcrank_handle(index_dir, mmap = TRUE, validate_crc = FALSE)
    gcrank_contigs_cpp(xptr)
}

#' Query GC content from a gcrank index
#'
#' Returns the GC fraction for half-open, 0-based, contig-relative intervals.
#' The query is vectorized over `contig`, `start` and `end`; each argument may
#' have length 1 or the common length. With `remove_ambiguous = FALSE` (the
#' default) the denominator is the window length, matching `bedtools nuc`;
#' with `TRUE`, ambiguous bases are removed from the denominator.
#'
#' With `mmap = TRUE` (the default) the index is memory-mapped and resident
#' memory grows as pages are touched. The index handle is cached per normalized
#' directory, with at most two indexes kept open.
#'
#' @param index_dir Path to the `*.gcrank` index directory.
#' @param contig Contig name(s).
#' @param start Window start(s), 0-based.
#' @param end Window end(s), exclusive.
#' @param remove_ambiguous Remove ambiguous bases from the denominator.
#' @param mmap Memory-map the index instead of reading it.
#' @param validate_crc Re-read and verify every data file against the manifest.
#'
#' @return A numeric vector of GC fractions.
#' @export
gcrankGc <- function(index_dir, contig, start, end,
                     remove_ambiguous = FALSE, mmap = TRUE,
                     validate_crc = FALSE) {
    index_dir <- .validate_gcrank_index(index_dir)
    remove_ambiguous <- .gcrank_logical(remove_ambiguous, "remove_ambiguous")
    mmap <- .gcrank_logical(mmap, "mmap")
    validate_crc <- .gcrank_logical(validate_crc, "validate_crc")

    if (is.factor(contig)) {
        contig <- as.character(contig)
    }
    if (!is.atomic(contig) || is.list(contig) || anyNA(contig)) {
        .gcrank_error("`contig` must be an atomic vector without NA values")
    }
    contig <- as.character(contig)

    if (!is.numeric(start) || !is.numeric(end)) {
        .gcrank_error("`start` and `end` must be numeric vectors")
    }
    start <- as.numeric(start)
    end <- as.numeric(end)

    n <- max(length(contig), length(start), length(end))
    if (n == 0L) {
        return(numeric(0))
    }
    for (value in list(contig = contig, start = start, end = end)) {
        if (length(value) != 1L && length(value) != n) {
            .gcrank_error(
                "`contig`, `start` and `end` must have length 1 or %d",
                n
            )
        }
    }
    if (length(contig) == 1L) contig <- rep(contig, n)
    if (length(start) == 1L) start <- rep(start, n)
    if (length(end) == 1L) end <- rep(end, n)

    if (any(!nzchar(contig))) {
        .gcrank_error("`contig` must not contain empty names")
    }
    if (any(!is.finite(start)) || any(!is.finite(end))) {
        .gcrank_error("Window coordinates must be finite")
    }
    if (any(start < 0) || any(end < 0) ||
        any(start > 9007199254740991) || any(end > 9007199254740991)) {
        .gcrank_error("Window coordinates are out of range")
    }
    if (any(start != floor(start)) || any(end != floor(end))) {
        .gcrank_error("Window coordinates must be whole numbers")
    }
    if (any(end <= start)) {
        bad <- which(end <= start)[[1L]]
        .gcrank_error(
            "Window end must be greater than start (position %d: %s:%s-%s)",
            bad, contig[[bad]],
            format(start[[bad]], scientific = FALSE),
            format(end[[bad]], scientific = FALSE)
        )
    }

    xptr <- .gcrank_handle(index_dir, mmap = mmap, validate_crc = validate_crc)
    tryCatch(
        gcrank_gc_cpp(xptr, contig, start, end, remove_ambiguous),
        error = function(e) {
            .gcrank_error("%s", conditionMessage(e))
        }
    )
}

#' Clear the gcrank index cache
#'
#' Closes every cached index handle. Useful in tests and before handing a
#' long-lived process to another user.
#'
#' @return `NULL`, invisibly.
#' @export
gcrankClearCache <- function() {
    keys <- ls(envir = .MaMaMIAgcrankCache, all.names = TRUE)
    for (key in keys) {
        entry <- get(key, envir = .MaMaMIAgcrankCache, inherits = FALSE)
        try(gcrank_close_cpp(entry$xptr), silent = TRUE)
    }
    rm(list = keys, envir = .MaMaMIAgcrankCache)
    .MaMaMIAgcrankState$counter <- 0L
    invisible(NULL)
}
