## Multi-BAM preflight validation (T1.4).
##
## Preflight is deliberately cheap: it reads file metadata and headers only.
## No index is built and no alignment is read, so every structural problem
## (wrong assembly, duplicate sample, limits, path escapes) fails before the
## expensive stages start. Index status is reported for the UI; the pipeline
## itself calls ensureBamIndex() later.

.preflight_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAPreflightError", fmt, ...)
}

.strip_trailing_slash <- function(path) {
    if (nchar(path) > 1L && grepl("[/\\\\]$", path)) {
        sub("[/\\\\]+$", "", path)
    } else {
        path
    }
}

.path_under <- function(path, root) {
    path <- .strip_trailing_slash(path)
    root <- .strip_trailing_slash(root)
    if (!nzchar(root)) {
        root <- "/"
    }
    if (.Platform$OS.type == "windows") {
        path <- tolower(path)
        root <- tolower(root)
    }
    if (identical(path, root)) {
        return(TRUE)
    }
    prefix <- if (identical(root, "/") || grepl("[/\\\\]$", root)) {
        root
    } else {
        paste0(root, "/")
    }
    startsWith(path, prefix)
}

.file_readable <- function(path) {
    if (.Platform$OS.type == "unix") {
        ## Permission bits, not file.access(): root can read a 0000 file, but
        ## the application must still treat it as unreadable.
        mode <- suppressWarnings(as.integer(file.info(path)$mode))
        if (length(mode) == 1L && !is.na(mode)) {
            return(bitwAnd(mode, 292L) != 0L) # 0444
        }
    }
    file.access(path, mode = 4L) == 0L
}

.session_temp_dir <- function() {
    getOption("MaMaMIA.session_dir", tempdir())
}

.check_bam_containment <- function(path, bam_root, session_dir = .session_temp_dir()) {
    if (is.null(bam_root)) {
        return(invisible(path))
    }
    if (!is.character(bam_root) || length(bam_root) != 1L ||
        is.na(bam_root) || !nzchar(bam_root)) {
        .mamamia_error(
            "MaMaMIAPathError", "`server.bam_root` must be a single directory path"
        )
    }
    bam_root <- path.expand(bam_root)
    if (!dir.exists(bam_root)) {
        .mamamia_error(
            "MaMaMIAPathError",
            "server.bam_root '%s' does not exist or is not a directory", bam_root
        )
    }

    resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
    root <- normalizePath(bam_root, winslash = "/", mustWork = TRUE)
    if (.path_under(resolved, root)) {
        return(invisible(path))
    }

    ## Uploads and staged copies live in session temp. A symlink there that
    ## points outside both the session directory and bam_root is still
    ## rejected because normalizePath() resolves it first.
    session <- normalizePath(session_dir, winslash = "/", mustWork = FALSE)
    if (.path_under(resolved, session)) {
        return(invisible(path))
    }

    .mamamia_error(
        "MaMaMIAPathError",
        "BAM path '%s' is outside server.bam_root '%s'", path, bam_root
    )
}

.csi_status <- function(path) {
    index_path <- paste0(path, ".csi")
    writable <- .directory_writable(dirname(path))
    status <- "missing"

    if (file.exists(index_path)) {
        if (!.bam_index_magic(index_path)) {
            status <- "corrupt"
        } else {
            index_mtime <- file.info(index_path)$mtime
            bam_mtime <- file.info(path)$mtime
            current <- !is.na(index_mtime) && !is.na(bam_mtime) &&
                as.numeric(index_mtime) >= as.numeric(bam_mtime)
            status <- if (current) "valid" else "stale"
        }
    }

    action <- if (identical(status, "valid")) {
        "reuse"
    } else if (writable) {
        if (identical(status, "missing")) "create" else "rebuild"
    } else {
        "stage"
    }

    list(index_path = index_path, status = status, action = action)
}

.preflight_warnings <- function(header, targets, index, path, assembly_key) {
    warnings <- character(0)

    if (identical(header$sample_source, "basename")) {
        warnings <- c(
            warnings,
            "No @RG SM value; the sample name is derived from the file name."
        )
    }

    extra <- setdiff(header$targets$name, targets$chr_id)
    if (length(extra) > 0L) {
        warnings <- c(
            warnings,
            sprintf(
                "The BAM header contains %d target(s) not in assembly '%s'; they will be ignored.",
                length(extra), assembly_key
            )
        )
    }

    if (file.exists(paste0(path, ".bai"))) {
        warnings <- c(
            warnings,
            "A BAI index is present but ignored; MaMaMIA requires a CSI index."
        )
    }

    index_warning <- switch(index$action,
        create = "A CSI index will be created next to the BAM.",
        rebuild = "The existing CSI index is stale or corrupt and will be rebuilt.",
        stage = "The BAM directory is not writable; the BAM will be staged in session temp and indexed there.",
        NULL
    )
    if (!is.null(index_warning)) {
        warnings <- c(warnings, index_warning)
    }

    paste(warnings, collapse = " ")
}

#' Validate BAM inputs before an analysis
#'
#' Reads only file metadata and BAM headers. Every path must exist, be a
#' non-empty regular file, be readable, and (when `server.bam_root` is set) live
#' under that root or in session temporary storage. File count and total size
#' limits from `config$server` are enforced. Headers are parsed once, sample
#' names are resolved per D2, duplicate samples are rejected, and every BAM must
#' match the same configured assembly.
#'
#' No index is built and no alignment is read. The returned table reports the
#' CSI status and the action the pipeline will take: `reuse`, `create`,
#' `rebuild`, or `stage` (for read-only BAM directories).
#'
#' @param paths Character vector of BAM paths.
#' @param config A validated configuration list from [readMaMaMIAConfig()], or a
#'   path to a YAML configuration.
#' @param assembly Optional assembly key to validate instead of auto-selecting.
#'
#' @return A data frame with one row per BAM and columns `path`, `file_size`,
#'   `sample`, `sample_source`, `assembly`, `n_targets`,
#'   `n_assembly_targets`, `index_path`, `index_status`, `index_action`,
#'   `warnings` and `valid`.
#' @export
preflightBams <- function(paths, config, assembly = NULL) {
    if (is.character(config)) {
        config <- readMaMaMIAConfig(config)
    }
    if (!is.list(config) || !is.list(config$assemblies) ||
        length(config$assemblies) == 0L) {
        .config_error(
            "`config` must be a validated MaMaMIA configuration with at least one assembly"
        )
    }
    if (!is.character(paths) || length(paths) == 0L || anyNA(paths) ||
        any(!nzchar(paths))) {
        .preflight_error("`paths` must be a non-empty character vector of BAM paths")
    }

    server <- if (is.null(config$server)) {
        .config_server_defaults()
    } else {
        config$server
    }
    defaults <- .config_server_defaults()
    if (is.null(server$max_bams)) {
        server$max_bams <- defaults$max_bams
    }
    if (is.null(server$max_total_bytes)) {
        server$max_total_bytes <- defaults$max_total_bytes
    }
    if (is.null(server$bam_root)) {
        server$bam_root <- defaults$bam_root
    }

    n <- length(paths)
    if (n > server$max_bams) {
        .mamamia_error(
            "MaMaMIALimitError",
            "Too many BAM files: %d (limit %d)", n, server$max_bams
        )
    }

    checked <- character(n)
    sizes <- numeric(n)
    for (i in seq_len(n)) {
        path <- .validate_bam_path(paths[[i]])
        if (!isTRUE(utils::file_test("-f", path))) {
            .bam_header_error("BAM path is not a regular file: %s", path)
        }
        if (!.file_readable(path)) {
            .bam_header_error("BAM file is not readable: %s", path)
        }
        .check_bam_containment(path, server$bam_root)
        size <- file.info(path)$size
        if (length(size) != 1L || is.na(size) || !is.finite(size)) {
            .bam_header_error("Could not determine the size of BAM '%s'", path)
        }
        checked[[i]] <- path
        sizes[[i]] <- as.numeric(size)
    }

    total_size <- sum(sizes)
    if (total_size > server$max_total_bytes) {
        .mamamia_error(
            "MaMaMIALimitError",
            "Total BAM size %s exceeds %s bytes",
            format(total_size, scientific = FALSE),
            format(server$max_total_bytes, scientific = FALSE)
        )
    }

    headers <- lapply(checked, bamHeaderInfo)
    samples <- vapply(headers, `[[`, character(1), "sample")
    duplicated_samples <- unique(samples[duplicated(samples)])
    if (length(duplicated_samples) > 0L) {
        sample <- duplicated_samples[[1L]]
        sample_paths <- checked[samples == sample]
        .mamamia_error(
            "MaMaMIASampleError",
            "Duplicate sample name '%s' in input BAMs: %s",
            sample, paste(sprintf("'%s'", sample_paths), collapse = ", ")
        )
    }

    matches <- matchAssemblyCohort(headers, config, assembly)
    keys <- vapply(matches, `[[`, character(1), "assembly")
    assembly_key <- keys[[1L]]
    targets <- .assembly_targets_for_match(
        matches[[1L]]$assembly_config, assembly_key
    )

    indexes <- lapply(checked, .csi_status)
    warnings <- vapply(seq_len(n), function(i) {
        .preflight_warnings(
            headers[[i]], targets, indexes[[i]], checked[[i]], assembly_key
        )
    }, character(1))

    out <- data.frame(
        path = checked,
        file_size = sizes,
        sample = samples,
        sample_source = vapply(headers, `[[`, character(1), "sample_source"),
        assembly = rep(assembly_key, n),
        n_targets = vapply(headers, `[[`, integer(1), "n_targets"),
        n_assembly_targets = rep(as.integer(nrow(targets)), n),
        index_path = vapply(indexes, `[[`, character(1), "index_path"),
        index_status = vapply(indexes, `[[`, character(1), "status"),
        index_action = vapply(indexes, `[[`, character(1), "action"),
        warnings = warnings,
        valid = rep(TRUE, n),
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
}

#' List server-side BAM files under a configured root
#'
#' Returns regular `.bam` files below `root`, bounded by depth and file count.
#' Symlinks are resolved; anything that escapes `root` is dropped. This is the
#' listing used by the Shiny server-side picker.
#'
#' @param root Directory to scan.
#' @param max_depth Maximum directory depth (1 = files directly in `root`).
#' @param max_files Maximum number of files to return.
#'
#' @return A sorted character vector of absolute BAM paths. The
#'   `"truncated"` attribute is `TRUE` when more files were found than
#'   `max_files`.
#' @noRd
.list_bam_files <- function(root, max_depth = 3L, max_files = 5000L) {
    if (!is.character(root) || length(root) != 1L || is.na(root) ||
        !nzchar(root)) {
        .mamamia_error("MaMaMIAPathError", "`root` must be a single directory path")
    }
    root <- path.expand(root)
    if (!dir.exists(root)) {
        .mamamia_error(
            "MaMaMIAPathError",
            "BAM root '%s' does not exist or is not a directory", root
        )
    }
    max_depth <- .index_integer(max_depth, "max_depth", minimum = 1L)
    max_files <- .index_integer(max_files, "max_files", minimum = 1L)
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)

    files <- list.files(
        root,
        pattern = "(?i)\\.bam$",
        recursive = TRUE,
        full.names = TRUE,
        all.files = FALSE,
        include.dirs = FALSE,
        no.. = TRUE
    )
    if (length(files) == 0L) {
        out <- character(0)
        attr(out, "truncated") <- FALSE
        return(out)
    }

    relative <- substring(files, nchar(root) + 2L)
    depth <- lengths(strsplit(relative, "/", fixed = TRUE))
    files <- files[depth <= max_depth]
    files <- files[utils::file_test("-f", files)]
    if (length(files) == 0L) {
        out <- character(0)
        attr(out, "truncated") <- FALSE
        return(out)
    }

    resolved <- normalizePath(files, winslash = "/", mustWork = TRUE)
    keep <- vapply(resolved, .path_under, logical(1), root = root)
    resolved <- sort(unique(resolved[keep]))

    truncated <- length(resolved) > max_files
    if (truncated) {
        resolved <- utils::head(resolved, max_files)
    }
    attr(resolved, "truncated") <- truncated
    resolved
}
