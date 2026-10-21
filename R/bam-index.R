## CSI index management (D3, D8, D12).
##
## Every BAM is read through a CSI index with min_shift = 14. A BAI is never
## created or accepted. When the BAM directory is read-only the caller stages
## the BAM in session temp and indexes the staged copy.

.directory_writable <- function(path) {
    if (!dir.exists(path)) {
        return(FALSE)
    }
    if (.Platform$OS.type == "unix") {
        ## Check permission bits rather than file.access(): root can write to a
        ## 0555 directory, but the application must still treat it as
        ## read-only and stage the BAM instead.
        mode <- suppressWarnings(as.integer(file.info(path)$mode))
        if (length(mode) == 1L && !is.na(mode)) {
            return(bitwAnd(mode, 146L) != 0L) # 0222
        }
    }
    file.access(path, mode = 2L) == 0L
}

.bam_index_magic <- function(index_path) {
    if (!file.exists(index_path)) {
        return(FALSE)
    }
    info <- file.info(index_path)
    if (is.na(info$size) || info$size < 4) {
        return(FALSE)
    }
    ## CSI files are BGZF-compressed, so the "CSI\1" magic is only visible
    ## after decompression.
    magic <- suppressWarnings(tryCatch({
        connection <- gzfile(index_path, "rb")
        on.exit(close(connection), add = TRUE)
        readBin(connection, what = "raw", n = 4L)
    }, error = function(e) raw(0)))
    identical(magic, as.raw(c(0x43, 0x53, 0x49, 0x01))) # "CSI\1"
}

.index_integer <- function(x, what, minimum = 1L,
                           maximum = .Machine$integer.max) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
        !is.finite(x) || x != floor(x) || x < minimum || x > maximum) {
        .mamamia_error(
            "MaMaMIAIndexError",
            "`%s` must be a single integer between %d and %d",
            what, as.integer(minimum),
            as.integer(base::min(maximum, .Machine$integer.max))
        )
    }
    as.integer(x)
}

#' Ensure a valid CSI index exists for a BAM
#'
#' Returns an existing current `.csi` index or builds one atomically next to
#' the BAM with `min_shift = 14`. BAI indexes are never created or returned.
#' If the BAM directory is not writable, the function stops with guidance to
#' stage the BAM in session temp first (see the internal `stageBam()` helper).
#'
#' @param path Path to a BAM file.
#' @param min_shift CSI minimum shift; must be a positive integer. The default
#'   `14L` supports chromosomes longer than 512 Mb.
#' @param threads Number of htslib indexing threads.
#'
#' @return The path to the CSI index.
#'
#' @examples
#' demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
#' bam <- file.path(tempdir(), "demo.bam")
#' file.copy(file.path(demo, "demo.bam"), bam, overwrite = TRUE)
#' index <- ensureBamIndex(bam)
#' file.exists(index)
#'
#' @export
ensureBamIndex <- function(path, min_shift = 14L, threads = 1L) {
    path <- .validate_bam_path(path, class = "MaMaMIAIndexError")
    min_shift <- .index_integer(
        min_shift, "min_shift", minimum = 1L, maximum = 31L
    )
    threads <- .index_integer(threads, "threads", minimum = 1L)

    index_path <- paste0(path, ".csi")
    bam_mtime <- file.info(path)$mtime
    writable <- .directory_writable(dirname(path))

    if (file.exists(index_path)) {
        index_info <- file.info(index_path)
        current <- !is.na(index_info$mtime) && !is.na(bam_mtime) &&
            index_info$mtime >= bam_mtime
        valid <- .bam_index_magic(index_path)
        if (current && valid) {
            return(index_path)
        }
        if (!writable) {
            .mamamia_error(
                "MaMaMIAIndexError",
                "BAM '%s' has a corrupt, stale or unreadable CSI index: '%s'; the directory is not writable, so stage the BAM in session temp and index the staged copy",
                path, index_path
            )
        }
    }

    if (!writable) {
        .mamamia_error(
            "MaMaMIAIndexError",
            "Could not build a CSI index for '%s': the directory '%s' is not writable; stage the BAM in session temp and index the staged copy",
            path, dirname(path)
        )
    }

    temporary <- tempfile(
        pattern = paste0(basename(path), ".csi."),
        tmpdir = dirname(path)
    )
    on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)

    tryCatch(
        bam_index_build_cpp(path, temporary, min_shift, threads),
        error = function(e) {
            .mamamia_error(
                "MaMaMIAIndexError",
                "Could not build a CSI index for '%s': %s",
                path, conditionMessage(e)
            )
        }
    )

    if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
        file.info(temporary)$size <= 0) {
        .mamamia_error(
            "MaMaMIAIndexError",
            "Could not build a CSI index for '%s': htslib produced no index",
            path
        )
    }
    if (!.bam_index_magic(temporary)) {
        .mamamia_error(
            "MaMaMIAIndexError",
            "Could not build a CSI index for '%s': the output is not a CSI file",
            path
        )
    }

    ## R's file.rename() does not replace an existing file on Windows.
    if (.Platform$OS.type == "windows" && file.exists(index_path)) {
        unlink(index_path)
    }
    if (!file.rename(temporary, index_path)) {
        .mamamia_error(
            "MaMaMIAIndexError",
            "Could not build a CSI index for '%s': could not move '%s' to '%s'",
            path, temporary, index_path
        )
    }
    index_path
}

#' Stage a BAM in session temp for indexing
#'
#' Creates a private staging directory, symlinks the BAM into it when the
#' platform supports symlinks and copies it otherwise, and returns the staged
#' path. The returned path is deliberately **not** passed through
#' `normalizePath()`: htslib must create the CSI next to the symlink in session
#' temp rather than next to the read-only source BAM. The staging directory is
#' new and empty, so it cannot contain a BAI sibling that htslib might prefer.
#'
#' @param path Path to a BAM file.
#' @param session_dir Session temporary directory.
#'
#' @return The staged BAM path.
#' @noRd
stageBam <- function(path, session_dir = tempdir()) {
    path <- .validate_bam_path(path, class = "MaMaMIAIndexError")
    session_dir <- path.expand(session_dir)
    if (!dir.exists(session_dir)) {
        if (!dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)) {
            .mamamia_error(
                "MaMaMIAIndexError",
                "Could not create the staging directory '%s'", session_dir
            )
        }
    }
    session_dir <- normalizePath(session_dir, winslash = "/", mustWork = TRUE)
    stage_dir <- tempfile(pattern = "mamamia-bam-", tmpdir = session_dir)
    if (!dir.create(stage_dir, showWarnings = FALSE)) {
        .mamamia_error(
            "MaMaMIAIndexError",
            "Could not create the staging directory '%s'", stage_dir
        )
    }
    staged <- file.path(stage_dir, basename(path))

    linked <- FALSE
    if (.Platform$OS.type == "unix") {
        linked <- isTRUE(file.symlink(path, staged))
    }
    if (!linked) {
        copied <- file.copy(path, staged, overwrite = FALSE, copy.mode = FALSE)
        if (!isTRUE(copied)) {
            unlink(stage_dir, recursive = TRUE, force = TRUE)
            .mamamia_error(
                "MaMaMIAIndexError",
                "Could not stage BAM '%s' in '%s'", path, stage_dir
            )
        }
    }
    staged
}
