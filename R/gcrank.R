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

#' Extract GC content for analysis windows
#'
#' Routes every target contig to the gcrank index configured for it, calls
#' [gcrankGc()] once per index with all of that index's windows, and assembles
#' the result in the same order as the input window table. Windows are 0-based
#' half-open, contig-relative intervals.
#'
#' @param windows An `AnalysisWindows` list with a `table` element, or a data
#'   frame with `chr_id`, `start` and `stop`.
#' @param assembly An `AssemblyConfig` list whose contigs carry `index_path`.
#' @param mmap,validate_crc Passed to [gcrankGc()].
#'
#' @return A data frame with `chr_id`, `start`, `stop` and `gc`, one row per
#'   input window in input order. GC values are validated to be non-missing and
#'   inside `[0, 1]`.
#' @noRd
extractWindowGc <- function(windows, assembly, mmap = TRUE,
                            validate_crc = FALSE) {
    if (is.list(windows) && !is.data.frame(windows) &&
        is.data.frame(windows$table)) {
        table <- windows$table
    } else if (is.data.frame(windows)) {
        table <- windows
    } else {
        .window_error(
            "`windows` must be an AnalysisWindows list or a data frame with chr_id, start and stop"
        )
    }
    required <- c("chr_id", "start", "stop")
    missing <- setdiff(required, names(table))
    if (length(missing) > 0L) {
        .window_error(
            "The window table must have columns: %s",
            paste(required, collapse = ", ")
        )
    }
    if (nrow(table) == 0L) {
        return(data.frame(
            chr_id = character(0),
            start = integer(0),
            stop = integer(0),
            gc = numeric(0),
            stringsAsFactors = FALSE
        ))
    }

    targets <- .assembly_target_contigs(assembly)
    key <- if (is.null(assembly$key)) "<unknown>" else assembly$key
    position <- match(as.character(table$chr_id), targets$chr_id)
    if (anyNA(position)) {
        unknown <- unique(as.character(table$chr_id)[is.na(position)])[[1L]]
        .window_error(
            "Assembly '%s' has no target contig '%s' for the supplied windows",
            key, unknown
        )
    }

    start <- suppressWarnings(as.numeric(table$start))
    stop <- suppressWarnings(as.numeric(table$stop))
    if (anyNA(start) || anyNA(stop) || any(!is.finite(start)) ||
        any(!is.finite(stop))) {
        .window_error("Window coordinates must be finite numbers")
    }
    if (any(start != floor(start)) || any(stop != floor(stop))) {
        .window_error("Window coordinates must be whole numbers")
    }
    if (any(start < 0) || any(stop <= start)) {
        bad <- which(start < 0 | stop <= start)[[1L]]
        .window_error(
            "Invalid window %s:%s-%s: start must be non-negative and smaller than stop",
            table$chr_id[[bad]],
            format(start[[bad]], scientific = FALSE),
            format(stop[[bad]], scientific = FALSE)
        )
    }
    contig_length <- targets$length[position]
    if (any(stop > contig_length)) {
        bad <- which(stop > contig_length)[[1L]]
        .window_error(
            "Window %s:%s-%s exceeds the contig length (%s bp)",
            table$chr_id[[bad]],
            format(start[[bad]], scientific = FALSE),
            format(stop[[bad]], scientific = FALSE),
            format(contig_length[[bad]], scientific = FALSE)
        )
    }

    index_path <- targets$index_path[position]
    if (anyNA(index_path) || any(!nzchar(index_path))) {
        bad <- which(is.na(index_path) | !nzchar(index_path))[[1L]]
        .gcrank_error(
            "Assembly '%s' has no gc_index for contig '%s'",
            key, table$chr_id[[bad]]
        )
    }

    gc <- rep(NA_real_, nrow(table))
    groups <- split(seq_len(nrow(table)), index_path)
    for (index_dir in names(groups)) {
        idx <- groups[[index_dir]]
        values <- tryCatch(
            gcrankGc(
                index_dir,
                as.character(table$chr_id[idx]),
                start[idx],
                stop[idx],
                remove_ambiguous = FALSE,
                mmap = mmap,
                validate_crc = validate_crc
            ),
            error = function(e) {
                .gcrank_error(
                    "GC extraction failed for index '%s': %s",
                    index_dir, conditionMessage(e)
                )
            }
        )
        if (!is.numeric(values) || length(values) != length(idx)) {
            .gcrank_error(
                "GC extraction for index '%s' returned %d values for %d windows",
                index_dir, length(values), length(idx)
            )
        }
        if (anyNA(values) || any(!is.finite(values))) {
            bad <- which(is.na(values) | !is.finite(values))[[1L]]
            .gcrank_error(
                "GC content is missing for %s:%s-%s in index '%s'",
                table$chr_id[idx[[bad]]],
                format(start[idx[[bad]]], scientific = FALSE),
                format(stop[idx[[bad]]], scientific = FALSE),
                index_dir
            )
        }
        if (any(values < 0 | values > 1)) {
            bad <- which(values < 0 | values > 1)[[1L]]
            .gcrank_error(
                "GC content %.6f is outside [0, 1] for %s:%s-%s in index '%s'",
                values[[bad]],
                table$chr_id[idx[[bad]]],
                format(start[idx[[bad]]], scientific = FALSE),
                format(stop[idx[[bad]]], scientific = FALSE),
                index_dir
            )
        }
        gc[idx] <- values
    }

    if (anyNA(gc)) {
        .gcrank_error("GC content is missing for %d window(s)", sum(is.na(gc)))
    }

    data.frame(
        chr_id = as.character(table$chr_id),
        start = as.integer(start),
        stop = as.integer(stop),
        gc = gc,
        stringsAsFactors = FALSE
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
