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
#' non-empty regular file and be readable. File count and total size limits from
#' `config$server` are enforced. Headers are parsed once, sample names are
#' resolved per D2, duplicate samples are rejected, and every BAM must match the
#' same configured assembly.
#'
#' BAMs are supplied only through Shiny uploads; there is no server-side BAM
#' directory. No index is built and no alignment is read. The returned table
#' reports the CSI status and the action the pipeline will take: `reuse`,
#' `create`,
#' `rebuild`, or `stage` (when the containing directory is not writable).
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
