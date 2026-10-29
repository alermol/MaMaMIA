## BAM-first Shiny application (T6.1-T6.4).
##
## The app is deliberately thin: configuration, preflight, run orchestration
## and cohort export live here, while measurement, scoring, masking,
## segmentation and calibration stay in the pipeline. The uploaded BAM is
## staged by symlink into a per-session directory before preflight, so the
## pipeline never writes an index into the upload directory and the session
## can remove the staged file and index when it ends. The dashboard accepts one
## BAM per run; a multi-sample BAM is read through its @RG SM tags. A bundled
## pre-computed example is available without a BAM upload.

## Error taxonomy ------------------------------------------------------------

.bam_error_suggestions <- c(
    config = paste(
        "Check MAMAMIA_CONFIG and the YAML species/assembly registry;",
        "validate it with readMaMaMIAConfig()."
    ),
    BAM = paste(
        "Check that the file is a coordinate-sorted BAM;",
        "re-upload it or choose another file."
    ),
    assembly = paste(
        "Choose a species pair or assembly override that matches the BAM header,",
        "or use a BAM that contains every configured target."
    ),
    index = paste(
        "Check that the session temporary directory is writable and has",
        "enough free space, then re-upload the BAM."
    ),
    windows = paste(
        "Choose a lower resolution;",
        "a contig is shorter than the requested window count."
    ),
    GC = paste(
        "Check that every gcrank index in the configuration exists",
        "and has a supported format_version."
    ),
    metrics = paste(
        "Check MAPQ and the exclude flags;",
        "verify that the BAM is intact and indexed."
    ),
    score = paste(
        "Check that donor and recipient windows share chromosome groups",
        "and that every group has at least one donor and one recipient."
    ),
    mask = paste(
        "Rebuild the panel mask with the same window grid,",
        "or check the blacklist BED coordinates."
    ),
    segmentation = paste(
        "Check that both the donor and the recipient side",
        "have enough finite windows."
    ),
    calibration = paste(
        "Check the calibration length bins and the control segments;",
        "a very small cohort needs the configured single-sample floors."
    ),
    limits = paste(
        "Reduce the number of BAMs or their total size; the limits come from",
        "server.max_bams and server.max_total_bytes."
    ),
    analysis = paste(
        "Check the input mode, the registry defaults (MAPQ, exclude flags,",
        "fit limits) and the server log."
    )
)

#' @noRd
bam_error_stage <- function(e) {
    if (inherits(e, "MaMaMIAConfigError")) {
        return("config")
    }
    if (inherits(e, "MaMaMIAAssemblyMatchError") ||
        inherits(e, "MaMaMIAHeaderMismatchError")) {
        return("assembly")
    }
    if (inherits(e, "MaMaMIAIndexError")) {
        return("index")
    }
    if (inherits(e, "MaMaMIAWindowError")) {
        return("windows")
    }
    if (inherits(e, "MaMaMIAGcrankError")) {
        return("GC")
    }
    if (inherits(e, "MaMaMIACountError") ||
        inherits(e, "MaMaMIAWindowMetricsError")) {
        return("metrics")
    }
    if (inherits(e, "MaMaMIAScoringError")) {
        return("score")
    }
    if (inherits(e, "MaMaMIAMaskError")) {
        return("mask")
    }
    if (inherits(e, "MaMaMIACalibrationError")) {
        return("calibration")
    }
    if (inherits(e, "MaMaMIASegmentationError")) {
        return("segmentation")
    }
    if (inherits(e, "MaMaMIALimitError")) {
        return("limits")
    }
    if (inherits(e, "MaMaMIABamHeaderError") ||
        inherits(e, "MaMaMIAUnsupportedFormatError") ||
        inherits(e, "MaMaMIAPreflightError") ||
        inherits(e, "MaMaMIAPathError")) {
        return("BAM")
    }

    message <- conditionMessage(e)
    if (inherits(e, "MaMaMIASampleError")) {
        stage <- sub(
            "^Sample '[^']*' \\(([^)]*)\\): .*$", "\\1", message,
            perl = TRUE
        )
        if (!identical(stage, message)) {
            return(switch(stage,
                metrics = "metrics",
                score = "score",
                mask = "mask",
                segmentation = "segmentation",
                calibration = "calibration",
                stage
            ))
        }
        return("analysis")
    }

    "analysis"
}

#' Resolve the input file(s) responsible for an error
#'
#' Prefers a `path` in the condition payload, then the sample name embedded in
#' a `MaMaMIASampleError` message, then the caller's fallback list. File output
#' always uses the original source path when a preflight table is available,
#' never the staged session copy.
#'
#' @param e A condition.
#' @param files Fallback file paths.
#' @param table Optional preflight table with `sample` and `source` columns.
#'
#' @return A unique character vector of file paths (possibly length zero).
#' @noRd
bam_error_files <- function(e, files = NULL, table = NULL) {
    data <- if (inherits(e, "condition")) e$data else NULL
    if (is.list(data) && !is.null(data$path)) {
        paths <- as.character(data$path)
        paths <- paths[!is.na(paths) & nzchar(paths)]
        if (length(paths) > 0L) {
            return(unique(paths))
        }
    }

    message <- conditionMessage(e)
    sample <- sub("^Sample '([^']+)' .*$", "\\1", message, perl = TRUE)
    if (!identical(sample, message) && is.data.frame(table) &&
        all(c("sample", "source") %in% names(table))) {
        rows <- table$sample == sample
        if (any(rows)) {
            return(unique(table$source[rows]))
        }
    }

    files <- if (is.null(files)) character() else as.character(files)
    files <- files[!is.na(files) & nzchar(files)]
    unique(files)
}

#' Render a compact `File:`/`Files:` line
#' @noRd
bam_error_file_label <- function(files, max_files = 5L) {
    files <- as.character(files)
    files <- files[!is.na(files) & nzchar(files)]
    if (length(files) == 0L) {
        return(NULL)
    }
    hidden <- length(files) - max_files
    if (hidden > 0L) {
        files <- c(
            utils::head(files, max_files),
            sprintf("(+%d more)", hidden)
        )
    }
    sprintf(
        "%s: %s",
        if (length(files) == 1L) "File" else "Files",
        paste(files, collapse = ", ")
    )
}

#' Describe a user-facing error with a stage, file and a suggested fix
#'
#' @param e A condition.
#' @param files Fallback file paths when the condition does not name one.
#' @param table Optional preflight table used to map a sample to its source
#'   BAM.
#'
#' @return A list with `stage`, `class`, `message`, `files` and `suggestion`.
#' @noRd
bam_error_info <- function(e, files = NULL, table = NULL) {
    stage <- bam_error_stage(e)
    suggestion <- .bam_error_suggestions[[stage]]
    if (is.null(suggestion) || is.na(suggestion)) {
        suggestion <- .bam_error_suggestions[["analysis"]]
    }
    list(
        stage = stage,
        class = class(e)[[1L]],
        message = conditionMessage(e),
        files = bam_error_files(e, files = files, table = table),
        suggestion = suggestion
    )
}

#' @noRd
bam_error_text <- function(info) {
    lines <- sprintf("[%s] %s", info$stage, info$message)
    file_label <- bam_error_file_label(info$files)
    if (!is.null(file_label)) {
        lines <- c(lines, file_label)
    }
    paste(c(lines, paste("Suggested fix:", info$suggestion)), collapse = "\n")
}

## Configuration -------------------------------------------------------------

#' Load and normalize the application configuration
#'
#' `config` may be a path to a YAML registry, a validated configuration list,
#' or `NULL`. In the latter case `MAMAMIA_CONFIG` is used; when that is empty
#' too, the returned state is not `ok` and the app renders a configuration
#' error card instead of crashing.
#'
#' @param config Configuration path, list or `NULL`.
#'
#' @return A list with `ok`, `config`, `path`, `source` and `error`.
#' @noRd
bam_launcher_state <- function(config = NULL) {
    source <- "argument"
    if (is.null(config)) {
        path <- Sys.getenv("MAMAMIA_CONFIG", unset = "")
        if (!nzchar(path)) {
            return(list(
                ok = FALSE,
                config = NULL,
                path = NULL,
                source = "missing",
                error = list(
                    stage = "config",
                    class = "MaMaMIAConfigError",
                    message = paste(
                        "No MaMaMIA configuration was supplied: set",
                        "MAMAMIA_CONFIG or pass `config` to runMaMaMIA()."
                    ),
                    suggestion = .bam_error_suggestions[["config"]]
                )
            ))
        }
        config <- path
        source <- "environment"
    }

    parsed <- tryCatch(
        if (is.character(config)) {
            readMaMaMIAConfig(config)
        } else if (is.list(config) && !is.null(config$defaults) &&
            !is.null(config$server) && !is.null(config$assemblies) &&
            !is.null(config$species)) {
            config
        } else {
            validateMaMaMIAConfig(config)
        },
        error = function(e) e
    )
    if (inherits(parsed, "condition")) {
        return(list(
            ok = FALSE,
            config = NULL,
            path = if (is.character(config)) config else NULL,
            source = source,
            error = bam_error_info(
                parsed,
                files = if (is.character(config)) config else NULL
            )
        ))
    }
    list(
        ok = TRUE,
        config = parsed,
        path = attr(parsed, "config_path"),
        source = source,
        error = NULL
    )
}

#' @noRd
bam_species_choices <- function(config) {
    species <- config$species
    if (is.null(species) || length(species) == 0L) {
        return(NULL)
    }
    keys <- names(species)
    stats::setNames(
        keys,
        vapply(keys, function(key) {
            sprintf("%s (%s)", species[[key]]$label, key)
        }, character(1))
    )
}

#' Every component declared for a species
#'
#' The dashboard never preselects a component: it analyses all of them in one
#' assembly. The order follows the registry so labels read A, B, D or S as the
#' config declares them.
#' @noRd
bam_species_components <- function(species) {
    if (is.null(species) || is.null(species$contigs)) {
        return(character())
    }
    unique(as.character(species$contigs$component))
}

#' Resolve a species pair over all of its components
#' @noRd
bam_resolve_species_pair <- function(config, recipient, donor) {
    recipient_species <- config$species[[recipient]]
    donor_species <- config$species[[donor]]
    resolveSpeciesAssembly(
        config,
        donor = donor,
        recipient = recipient,
        donor_components = bam_species_components(donor_species),
        recipient_components = bam_species_components(recipient_species)
    )
}

## Sidebar choices -----------------------------------------------------------

#' @noRd
bam_resolution_choices <- function(defaults = .config_defaults()) {
    allowed <- sort(unique(as.integer(defaults$allowed_windows)))
    if (length(allowed) == 0L) {
        allowed <- 1000L
    }
    labels <- vapply(allowed, function(n) {
        switch(as.character(n),
            `100` = "Low (100 per chromosome)",
            `1000` = "Medium (1,000 per chromosome)",
            `10000` = "High (10,000 per chromosome)",
            sprintf(
                "%s per chromosome",
                format(n, big.mark = ",", scientific = FALSE)
            )
        )
    }, character(1))
    stats::setNames(as.character(allowed), labels)
}

#' @noRd
bam_default_resolution <- function(defaults = .config_defaults()) {
    allowed <- sort(unique(as.integer(defaults$allowed_windows)))
    if (length(allowed) == 0L) {
        allowed <- 1000L
    }
    value <- suppressWarnings(as.integer(defaults$windows_per_chromosome))
    if (is.na(value) || !value %in% allowed) {
        value <- if (1000L %in% allowed) 1000L else allowed[[1L]]
    }
    as.character(value)
}

#' @noRd
bam_resolution_value <- function(value, defaults = .config_defaults()) {
    allowed <- sort(unique(as.integer(defaults$allowed_windows)))
    if (length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
        .mamamia_error("MaMaMIAWindowError", "Choose an analysis resolution")
    }
    resolution <- suppressWarnings(as.integer(value))
    if (is.na(resolution) || !resolution %in% allowed) {
        .mamamia_error(
            "MaMaMIAWindowError",
            "Resolution must be one of: %s",
            paste(allowed, collapse = ", ")
        )
    }
    resolution
}

#' @noRd
bam_integer <- function(value, label, minimum, maximum) {
    if (length(value) != 1L || is.na(value) || !is.numeric(value) ||
        !is.finite(value) || value != floor(value) ||
        value < minimum || value > maximum) {
        .mamamia_error(
            "MaMaMIAAnalysisError",
            "%s must be a whole number between %s and %s",
            label,
            format(minimum, scientific = FALSE),
            format(maximum, scientific = FALSE)
        )
    }
    as.integer(value)
}

#' @noRd
bam_max_cores <- function() {
    cores <- parallel::detectCores()
    if (is.na(cores) || cores < 1L) {
        return(1L)
    }
    as.integer(cores)
}

#' Resolve the single uploaded BAM
#'
#' The dashboard analyses one BAM per run. A multi-sample BAM is still one
#' file, so several rows from `input$bams` always mean several files and are
#' rejected rather than silently truncated. The client already trims a
#' multi-file drag-and-drop; this is the server-side contract.
#'
#' @param uploads The data frame Shiny stores in `input$bams`.
#'
#' @return A character vector with zero or one BAM path.
#' @noRd
bam_single_upload_path <- function(uploads) {
    if (!is.data.frame(uploads) || !"datapath" %in% names(uploads)) {
        return(character())
    }
    paths <- as.character(uploads$datapath)
    paths <- paths[!is.na(paths) & nzchar(paths)]
    if (length(paths) > 1L) {
        .mamamia_error(
            "MaMaMIAUploadError",
            paste(
                "The dashboard analyses one BAM per run;",
                "%d files were supplied. A multi-sample BAM is one file."
            ),
            length(paths)
        )
    }
    paths
}

## BAM staging ---------------------------------------------------------------

## Staged files below this size are copied rather than symlinked.
.MaMaMIACopyThreshold <- 1024^2

#' Stage BAM files in a session directory
#'
#' Symlinks each BAM (and its CSI index when present) into `session_dir` on
#' platforms that support symlinks, and copies it otherwise. Small files (below
#' 1 MiB) are copied even when symlinks are available: the copy is cheap, and
#' temporary-directory cleanup that follows and chmods symlinks then cannot
#' touch the source file. The staged directory is writable either way, so
#' indexing and analysis never modify the source directory.
#'
#' @param paths Uploaded BAM paths from Shiny.
#' @param session_dir Private per-session directory.
#'
#' @return A data frame with `path` (staged), `source` and `dir`.
#' @noRd
bam_stage_paths <- function(paths, session_dir) {
    if (!is.character(paths) || length(paths) == 0L || anyNA(paths) ||
        any(!nzchar(paths))) {
        .mamamia_error(
            "MaMaMIAPathError",
            "`paths` must be a non-empty character vector of BAM paths"
        )
    }
    session_dir <- path.expand(session_dir)
    if (!dir.exists(session_dir)) {
        if (!dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)) {
            .mamamia_error(
                "MaMaMIAIndexError",
                "Could not create the session directory '%s'", session_dir
            )
        }
    }
    stage_dir <- tempfile(pattern = "bams-", tmpdir = session_dir)
    if (!dir.create(stage_dir, showWarnings = FALSE)) {
        .mamamia_error(
            "MaMaMIAIndexError",
            "Could not create the staging directory '%s'", stage_dir
        )
    }

    stage_one <- function(source, target) {
        size <- file.info(source)$size
        if (length(size) == 1L && !is.na(size) &&
            size < .MaMaMIACopyThreshold) {
            copied <- isTRUE(file.copy(
                source, target,
                overwrite = FALSE, copy.mode = FALSE
            ))
            if (copied) {
                return(TRUE)
            }
        }
        linked <- FALSE
        if (.Platform$OS.type == "unix") {
            linked <- isTRUE(suppressWarnings(file.symlink(source, target)))
        }
        if (!linked) {
            linked <- isTRUE(file.copy(
                source, target,
                overwrite = FALSE, copy.mode = FALSE
            ))
        }
        linked
    }

    staged <- character(length(paths))
    sources <- character(length(paths))
    for (i in seq_along(paths)) {
        source <- path.expand(paths[[i]])
        if (!file.exists(source) || dir.exists(source)) {
            unlink(stage_dir, recursive = TRUE, force = TRUE)
            .mamamia_error(
                "MaMaMIAPathError",
                "BAM file does not exist or is not a file: %s", paths[[i]]
            )
        }
        ## A relative symlink would be resolved against the staging directory,
        ## not the working directory, so stage absolute paths only.
        source <- normalizePath(source, winslash = "/", mustWork = TRUE)
        name <- sprintf(
            "%03d-%s", i,
            gsub("[^A-Za-z0-9._-]+", "_", basename(source))
        )
        target <- file.path(stage_dir, name)
        if (!stage_one(source, target)) {
            unlink(stage_dir, recursive = TRUE, force = TRUE)
            .mamamia_error(
                "MaMaMIAIndexError",
                paste(
                    "Could not stage BAM '%s' in '%s';",
                    "upload the file instead"
                ),
                source, stage_dir
            )
        }
        index <- paste0(source, ".csi")
        if (file.exists(index)) {
            stage_one(index, paste0(target, ".csi"))
        }
        staged[[i]] <- target
        sources[[i]] <- source
    }

    data.frame(
        path = staged,
        source = sources,
        dir = stage_dir,
        stringsAsFactors = FALSE
    )
}

#' Format a preflight state for the dashboard calculation log
#'
#' The footer log is the only place that reports preflight status, detected
#' assemblies and problems, so the messages are deliberately explicit.
#' @noRd
bam_preflight_log <- function(pf) {
    if (is.null(pf)) {
        return("Preflight: no state available.")
    }
    if (isTRUE(pf$empty)) {
        return("Preflight: waiting for a BAM upload.")
    }
    if (!isTRUE(pf$ok)) {
        info <- pf$error
        if (is.null(info)) {
            return("Preflight failed.")
        }
        lines <- sprintf(
            "Preflight failed at stage '%s': %s", info$stage, info$message
        )
        file_label <- bam_error_file_label(info$files)
        if (!is.null(file_label)) {
            lines <- c(lines, file_label)
        }
        lines <- c(lines, sprintf("Suggested fix: %s", info$suggestion))
        candidates <- pf$candidates
        if (is.data.frame(candidates) && nrow(candidates) > 0L) {
            labels <- candidates$label
            if (is.null(labels) || length(labels) != nrow(candidates)) {
                labels <- candidates$key
            }
            labels[is.na(labels) | !nzchar(labels)] <-
                candidates$key[is.na(labels) | !nzchar(labels)]
            lines <- c(lines, paste0(
                "Assembly candidates: ",
                paste(sprintf(
                    "%s (%d/%d targets, %.1f%%)",
                    labels, candidates$n_matched,
                    candidates$n_targets, 100 * candidates$score
                ), collapse = "; ")
            ))
        }
        return(lines)
    }

    table <- pf$table
    n_bams <- if (is.null(table) || nrow(table) == 0L) {
        0L
    } else {
        length(unique(table$path))
    }
    n_samples <- if (is.null(table)) 0L else nrow(table)
    lines <- character()
    lines <- c(lines, sprintf(
        "Preflight OK: %d BAM file%s, %d sample%s.",
        n_bams, if (n_bams == 1L) "" else "s",
        n_samples, if (n_samples == 1L) "" else "s"
    ))
    if (!is.null(table) && nrow(table) > 0L) {
        warnings <- unique(table$warnings[nzchar(table$warnings)])
        if (length(warnings) > 0L) {
            lines <- c(lines, paste0("Preflight note: ", warnings))
        }
    }

    assembly <- pf$assembly_config
    if (!is.null(assembly)) {
        label <- assembly$label
        if (is.null(label) || length(label) != 1L || is.na(label) ||
            !nzchar(label)) {
            label <- assembly$key
        }
        lines <- c(lines, sprintf("Detected assembly: %s.", label))
    }
    lines
}

#' Format one timestamped calculation-log entry
#'
#' Pure helper used both when the log is seeded during server creation and by
#' the reactive `log_line()` writer. Returning `NULL` for empty text lets `c()`
#' drop the entry.
#' @noRd
bam_log_entry <- function(message, time = Sys.time()) {
    if (length(message) != 1L || is.na(message) || !nzchar(message)) {
        return(NULL)
    }
    sprintf("[%s] %s", format(time, "%H:%M:%S"), message)
}

#' Initial calculation-log lines for a new session
#'
#' This must not read or write a reactive value: it runs during server
#' creation, before Shiny has entered a reactive consumer context.
#' @noRd
bam_initial_log <- function(time = Sys.time()) {
    c(
        bam_log_entry("MaMaMIA session started.", time),
        bam_log_entry("Configuration validated.", time),
        bam_log_entry("Samples are detected from the BAM header.", time)
    )
}

## Preflight presentation ----------------------------------------------------

#' @noRd
bam_preflight_display <- function(table) {
    if (is.null(table) || nrow(table) == 0L) {
        return(data.frame(
            file = character(), sample = character(), assembly = character(),
            targets = character(), index = character(),
            warnings = character(), stringsAsFactors = FALSE
        ))
    }
    source <- if ("source" %in% names(table)) table$source else table$path
    data.frame(
        file = basename(source),
        sample = table$sample,
        assembly = table$assembly,
        targets = sprintf("%d/%d", table$n_targets, table$n_assembly_targets),
        index = table$index_action,
        warnings = table$warnings,
        stringsAsFactors = FALSE
    )
}

#' Score every configured assembly against one BAM header
#'
#' Uses the same candidate scoring as `matchAssembly()`, without raising the
#' matching error, so the UI can show the closest assemblies.
#'
#' @noRd
bam_assembly_candidates <- function(path, config) {
    header <- tryCatch(bamHeaderInfo(path), error = function(e) NULL)
    if (is.null(header)) {
        return(NULL)
    }
    targets <- data.frame(
        name = as.character(header$targets$name),
        length = as.numeric(header$targets$length),
        stringsAsFactors = FALSE
    )
    candidates <- lapply(seq_along(config$assemblies), function(i) {
        candidate <- .match_candidate(
            targets, config$assemblies[[i]],
            names(config$assemblies)[[i]], header$path
        )
        candidate$header_target_names <- targets$name
        candidate
    })
    names(candidates) <- names(config$assemblies)
    .candidate_table(candidates)
}

#' Extract assembly candidates from an error, or compute them from the headers
#' @noRd
bam_error_candidates <- function(e, paths, config) {
    candidates <- e$data$candidates
    if (!is.null(candidates)) {
        return(candidates)
    }
    if (!inherits(e, "MaMaMIAAssemblyMatchError") &&
        !inherits(e, "MaMaMIAHeaderMismatchError")) {
        return(NULL)
    }
    for (path in paths) {
        candidates <- bam_assembly_candidates(path, config)
        if (!is.null(candidates)) {
            return(candidates)
        }
    }
    NULL
}

#' @noRd
bam_preflight_failure <- function(e, paths, config, files = NULL) {
    list(
        ok = FALSE,
        empty = FALSE,
        table = NULL,
        staged = NULL,
        assembly = NULL,
        assembly_config = NULL,
        config = config,
        key = NULL,
        unpaired = NULL,
        error = bam_error_info(
            e,
            files = if (is.null(files)) paths else files
        ),
        candidates = bam_error_candidates(e, paths, config)
    )
}

#' @noRd
bam_assembly_windows <- function(assembly, resolutions) {
    targets <- .assembly_target_contigs(assembly)
    lengths <- targets$length
    vapply(resolutions, function(n) {
        if (any(lengths < n)) {
            NA_integer_
        } else {
            as.integer(length(lengths) * n)
        }
    }, integer(1L))
}

#' Explain why an assembly cannot use a resolution
#'
#' Mirrors the check in [makeAnalysisWindows()] so the dashboard can disable
#' the Run button before an analysis starts. Returns `NULL` when every target
#' contig can be tiled into the requested number of windows.
#'
#' @param assembly An `AssemblyConfig` list.
#' @param resolution Windows per chromosome.
#'
#' @return `NULL` or a one-element character message.
#' @noRd
bam_resolution_problem <- function(assembly, resolution) {
    if (is.null(assembly) || length(resolution) != 1L || is.na(resolution)) {
        return(NULL)
    }
    resolution <- suppressWarnings(as.numeric(resolution))
    if (!is.finite(resolution) || resolution < 1) {
        return(NULL)
    }
    targets <- tryCatch(
        .assembly_target_contigs(assembly),
        error = function(e) NULL
    )
    if (is.null(targets) || !all(c("chr_id", "length") %in% names(targets))) {
        return(NULL)
    }
    short <- targets$length < resolution
    if (!any(short)) {
        return(NULL)
    }
    contigs <- sprintf(
        "'%s' (%s bp)",
        targets$chr_id[short],
        format(targets$length[short], big.mark = ",", scientific = FALSE)
    )
    sprintf(
        "Contig%s %s %s shorter than the requested %s windows; choose a lower resolution.",
        if (sum(short) == 1L) "" else "s",
        paste(contigs, collapse = ", "),
        if (sum(short) == 1L) "is" else "are",
        format(resolution, big.mark = ",", scientific = FALSE)
    )
}

#' @noRd
bam_assembly_info <- function(assembly, resolutions = NULL) {
    if (is.null(resolutions)) {
        resolutions <- sort(unique(as.integer(
            .config_defaults()$allowed_windows
        )))
    }
    resolutions <- as.integer(resolutions)
    donor <- assembly$contigs$donor
    recipient <- assembly$contigs$recipient
    pairs <- assembly$pairs
    list(
        key = assembly$key,
        label = assembly$label,
        description = assembly$description,
        donor = sprintf("%s (%s)", donor$id, donor$name),
        recipient = sprintf("%s (%s)", recipient$id, recipient$name),
        pairs = if (length(pairs) == 0L) {
            character()
        } else {
            sprintf("%s -> %s", names(pairs), unname(pairs))
        },
        windows = data.frame(
            resolution = resolutions,
            n_windows = bam_assembly_windows(assembly, resolutions),
            stringsAsFactors = FALSE
        )
    )
}

#' @noRd
bam_assembly_summary_ui <- function(assembly, resolutions = NULL) {
    info <- bam_assembly_info(assembly, resolutions)
    windows <- if (all(is.na(info$windows$n_windows))) {
        "unavailable"
    } else {
        paste(
            vapply(seq_len(nrow(info$windows)), function(i) {
                n <- info$windows$n_windows[[i]]
                if (is.na(n)) {
                    sprintf(
                        "%s: a contig is too short",
                        format(info$windows$resolution[[i]],
                            big.mark = ",", scientific = FALSE
                        )
                    )
                } else {
                    sprintf(
                        "%s: %s windows",
                        format(info$windows$resolution[[i]],
                            big.mark = ",", scientific = FALSE
                        ),
                        format(n, big.mark = ",", scientific = FALSE)
                    )
                }
            }, character(1)),
            collapse = "; "
        )
    }
    shiny::tags$dl(
        shiny::tags$dt("Assembly"),
        shiny::tags$dd(sprintf("%s (%s)", info$label, info$key)),
        if (nzchar(info$description)) {
            shiny::tags$dd(class = "text-muted", info$description)
        },
        shiny::tags$dt("Donor contigs"),
        shiny::tags$dd(paste(info$donor, collapse = ", ")),
        shiny::tags$dt("Recipient contigs"),
        shiny::tags$dd(paste(info$recipient, collapse = ", ")),
        shiny::tags$dt("Pairs"),
        shiny::tags$dd(paste(info$pairs, collapse = "; ")),
        shiny::tags$dt("Analysis windows"),
        shiny::tags$dd(windows)
    )
}

## Run parameters and progress ------------------------------------------------

#' Validate and resolve the run controls
#'
#' @param input The Shiny input list (or a list with the same fields).
#' @param config A validated configuration.
#' @param samples Sample names from preflight.
#'
#' @return A list with the resolved `analyzeBamInput()` arguments.
#' @noRd
bam_run_parameters <- function(input, config, samples) {
    if (!is.character(samples) || length(samples) == 0L) {
        .mamamia_error("MaMaMIAAnalysisError", "At least one sample is required")
    }
    defaults <- config$defaults
    resolution <- bam_resolution_value(input$resolution, defaults)
    mapq <- bam_integer(defaults$mapq, "MAPQ", 0L, 255L)
    exclude_flags <- bam_integer(
        defaults$exclude_flags, "Exclude flags", 0L, 4095L
    )
    cores <- bam_integer(
        min(30L, bam_max_cores()), "Cores", 1L, bam_max_cores()
    )
    donor_aggregation <- input$donor_aggregation
    if (is.null(donor_aggregation) || length(donor_aggregation) != 1L ||
        is.na(donor_aggregation) || !nzchar(donor_aggregation)) {
        donor_aggregation <- defaults$donor_aggregation
    }
    donor_aggregation <- as.character(donor_aggregation)
    if (!donor_aggregation %in% c("mean", "max", "sum")) {
        .mamamia_error(
            "MaMaMIAAnalysisError",
            "`donor_aggregation` must be one of: mean, max, sum"
        )
    }
    control_samples <- input$control_samples
    if (is.null(control_samples) || length(control_samples) == 0L ||
        length(samples) < 2L) {
        control_samples <- character(0)
    } else {
        control_samples <- unique(as.character(control_samples))
        control_samples <- control_samples[nzchar(control_samples)]
        unknown <- setdiff(control_samples, samples)
        if (length(unknown) > 0L) {
            .mamamia_error(
                "MaMaMIAAnalysisError",
                "Unknown control sample(s): %s", paste(unknown, collapse = ", ")
            )
        }
    }
    calibration <- if (length(control_samples) == 0L) {
        NULL
    } else {
        list(controls = control_samples)
    }
    list(
        resolution = resolution,
        mapq = mapq,
        exclude_flags = exclude_flags,
        cores = cores,
        metric = defaults$metric,
        score = defaults$score,
        donor_aggregation = donor_aggregation,
        control_samples = control_samples,
        calibration = calibration
    )
}

#' @noRd
bam_progress_stages <- function(n_samples, ...) {
    n_samples <- max(1L, as.integer(n_samples))
    c(
        "preflight", "index", "windows", "gc",
        rep("metrics", n_samples),
        "score", "mask", "segmentation", "calibration"
    )
}

.bam_progress_labels <- c(
    preflight = "Checking BAM headers",
    index = "Ensuring CSI indexes",
    windows = "Generating analysis windows",
    gc = "Reading GC content",
    metrics = "Measuring breadth and counts",
    score = "Scoring donor/recipient windows",
    mask = "Masking recurrent windows",
    segmentation = "Segmenting score profiles",
    calibration = "Calibrating call thresholds"
)

#' @noRd
bam_progress_label <- function(stage, detail = list()) {
    if (!is.list(detail)) {
        detail <- list()
    }
    ## `[[` on a named character vector throws "subscript out of bounds" for
    ## an unknown name. Progress stages are an open set (a pipeline stage, a
    ## sample count, a future event), so an unknown one must not abort the
    ## analysis.
    label <- unname(.bam_progress_labels[stage])
    if (length(label) != 1L || is.na(label)) {
        label <- as.character(stage)[[1L]]
    }

    count <- function(x) {
        if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x)) {
            return(NULL)
        }
        format(x, big.mark = ",", scientific = FALSE)
    }
    sample <- detail$sample
    sample_suffix <- if (!is.null(sample) && length(sample) == 1L &&
        !is.na(sample) && nzchar(sample)) {
        paste0(" for ", sample)
    } else {
        ""
    }

    switch(stage,
        preflight = {
            n <- count(detail$n_bams)
            if (is.null(n)) {
                return(label)
            }
            n_samples <- count(detail$n_samples)
            sprintf(
                "%s (%s BAM file%s%s)",
                label, n,
                if (isTRUE(detail$n_bams == 1)) "" else "s",
                if (!is.null(n_samples)) {
                    sprintf(", %s sample%s", n_samples,
                        if (isTRUE(detail$n_samples == 1)) "" else "s"
                    )
                } else {
                    ""
                }
            )
        },
        index = {
            n <- count(detail$n_bams)
            staged <- count(detail$staged)
            if (is.null(n) || is.null(staged)) {
                return(label)
            }
            sprintf("%s (%s of %s staged)", label, staged, n)
        },
        windows = {
            n <- count(detail$n_windows)
            per <- count(detail$n_per_chromosome)
            if (is.null(n)) {
                return(label)
            }
            sprintf(
                "%s (%s window%s%s)",
                label, n,
                if (isTRUE(detail$n_windows == 1)) "" else "s",
                if (!is.null(per)) sprintf(", %s per chromosome", per) else ""
            )
        },
        gc = {
            n <- count(detail$n_windows)
            if (is.null(n)) {
                return(label)
            }
            sprintf("%s for %s windows", label, n)
        },
        metrics = {
            total <- count(detail$total)
            if (is.null(total)) {
                return(paste0(label, sample_suffix))
            }
            sprintf("%s%s (%s reads)", label, sample_suffix, total)
        },
        score = {
            mode <- detail$mode
            metric <- detail$metric
            if (is.null(mode) || is.null(metric) || length(mode) != 1L ||
                is.na(mode) || !nzchar(mode)) {
                return(label)
            }
            sprintf("%s (%s, %s)", label, mode, metric)
        },
        mask = {
            fraction <- detail$masked_fraction
            if (is.null(fraction) || length(fraction) != 1L ||
                !is.finite(fraction)) {
                return(label)
            }
            sprintf("%s (%s%% of windows)", label, round(100 * fraction, 1))
        },
        calibration = {
            n <- count(detail$n_calls)
            if (is.null(n)) {
                return(label)
            }
            sprintf("%s (%s calibrated call%s)", label, n,
                if (isTRUE(detail$n_calls == 1)) "" else "s"
            )
        },
        segmentation = {
            n <- count(detail$n_segments)
            if (is.null(n)) {
                return(paste0(label, sample_suffix))
            }
            sprintf("%s%s (%s segment%s)", label, sample_suffix, n,
                if (isTRUE(detail$n_segments == 1)) "" else "s"
            )
        },
        paste0(label, sample_suffix)
    )
}

#' @noRd
bam_progress_callback <- function(n_samples, total_steps = NULL) {
    stages <- length(bam_progress_stages(n_samples))
    total <- suppressWarnings(as.integer(total_steps))
    if (length(total) != 1L || is.na(total) || total < 1L) {
        total <- stages
    }
    state <- new.env(parent = emptyenv())
    state$done <- 0L
    function(stage, detail = list()) {
        state$done <- state$done + 1L
        shiny::incProgress(
            1 / total,
            detail = bam_progress_label(stage, detail)
        )
        invisible(NULL)
    }
}

## Cohort summary and selection ---------------------------------------------

#' @noRd
bam_sample_choices <- function(result) {
    if (is.null(result) || is.null(result$samples) ||
        !"sample" %in% names(result$samples)) {
        return(character())
    }
    unique(as.character(result$samples$sample))
}

#' One row per sample with measured reads, calls and mask burden
#' @noRd
bam_cohort_summary <- function(result) {
    samples <- bam_sample_choices(result)
    if (length(samples) == 0L) {
        return(data.frame(
            sample = character(), reads_donor = numeric(),
            reads_recipient = numeric(), reads_total = numeric(),
            n_windows = integer(), n_calls = integer(),
            n_high_confidence = integer(), masked_fraction = numeric(),
            mode = character(), stringsAsFactors = FALSE
        ))
    }
    metrics <- result$metrics
    reads <- function(role) {
        vapply(samples, function(sample) {
            sum(
                as.numeric(metrics$count[
                    metrics$sample == sample & metrics$role == role
                ]),
                na.rm = TRUE
            )
        }, numeric(1))
    }
    donor <- reads("donor")
    recipient <- reads("recipient")
    meta <- result$samples
    position <- match(samples, meta$sample)
    data.frame(
        sample = samples,
        reads_donor = donor,
        reads_recipient = recipient,
        reads_total = donor + recipient,
        n_windows = as.integer(meta$n_windows[position]),
        n_calls = as.integer(meta$n_calls[position]),
        n_high_confidence = as.integer(meta$n_high_confidence[position]),
        masked_fraction = as.numeric(meta$masked_fraction[position]),
        mode = as.character(meta$mode[position]),
        stringsAsFactors = FALSE
    )
}

#' Coverage warnings for zero or near-zero samples
#'
#' A side is near zero when it has fewer than `near_zero` reads per analysis
#' window on average.
#'
#' @param result A `BamAnalysisResult`.
#' @param near_zero Minimum average reads per window.
#'
#' @return A character vector of actionable warnings.
#' @noRd
bam_coverage_warnings <- function(result, near_zero = 1) {
    samples <- bam_sample_choices(result)
    warnings <- character()
    if (length(samples) == 0L) {
        return(warnings)
    }
    for (sample in samples) {
        for (side in c("donor", "recipient")) {
            values <- as.numeric(result$metrics$count[
                result$metrics$sample == sample &
                    result$metrics$role == side
            ])
            n_windows <- length(values)
            total <- sum(values, na.rm = TRUE)
            if (n_windows == 0L) {
                next
            }
            if (total <= 0) {
                warnings <- c(warnings, sprintf(
                    "Sample '%s' has no reads on the %s targets.",
                    sample, side
                ))
            } else if (total / n_windows < near_zero) {
                warnings <- c(warnings, sprintf(
                    "Sample '%s' has near-zero coverage on the %s side (%s reads over %d windows).",
                    sample, side, format(total, scientific = FALSE), n_windows
                ))
            }
        }
    }
    warnings
}

#' Cache completed analyses for one Shiny session
#'
#' Returning to an input that was already analysed should reuse that result
#' instead of repeating the whole pipeline. `get()` returns `NULL` on a miss;
#' `set()` evicts the oldest entry once the current limit is reached;
#' `reserve()` raises the budget before a run adds its entry.
#'
#' @noRd
bam_result_cache <- function(max_entries = 3L) {
    max_entries <- max(1L, as.integer(max_entries))
    store <- new.env(parent = emptyenv())
    order <- character()

    valid_key <- function(key) {
        length(key) == 1L && !is.na(key) && nzchar(key)
    }

    list(
        get = function(key) {
            if (!valid_key(key) ||
                !exists(key, envir = store, inherits = FALSE)) {
                return(NULL)
            }
            get(key, envir = store, inherits = FALSE)
        },
        set = function(key, value) {
            if (!valid_key(key) || is.null(value)) {
                return(invisible(NULL))
            }
            assign(key, value, envir = store)
            order <<- c(order[order != key], key)
            while (length(order) > max_entries) {
                oldest <- order[[1L]]
                if (exists(oldest, envir = store, inherits = FALSE)) {
                    rm(list = oldest, envir = store)
                }
                order <<- order[-1L]
            }
            invisible(NULL)
        },
        ## A run can need more entries than the default session budget; grow
        ## the budget before the result arrives.
        reserve = function(n) {
            n <- suppressWarnings(as.integer(n))
            if (length(n) == 1L && !is.na(n) && n > max_entries) {
                max_entries <<- n
            }
            invisible(NULL)
        },
        size = function() length(order),
        limit = function() max_entries
    )
}

#' Key for a completed analysis
#'
#' Use the resolved assembly key rather than the raw species selectors: an
#' invalid widget value that falls back to the species default must not
#' invalidate a result that is still valid. The key also covers the staged BAM
#' paths, sample names and every registry default or Advanced control that
#' reaches [analyzeBamInput()] and changes the result.
#'
#' @noRd
bam_analysis_key <- function(preflight, resolution, donor_aggregation = NULL,
                             controls = NULL) {
    if (!isTRUE(preflight$ok) || is.null(preflight$table) ||
        !all(c("path", "sample") %in% names(preflight$table))) {
        return(NULL)
    }
    defaults <- preflight$config$defaults
    if (is.null(defaults)) {
        defaults <- list()
    }
    if (is.null(donor_aggregation) || length(donor_aggregation) != 1L ||
        is.na(donor_aggregation) || !nzchar(donor_aggregation)) {
        donor_aggregation <- defaults$donor_aggregation
    }
    samples <- unique(as.character(preflight$table$sample))
    controls <- unique(as.character(controls))
    controls <- controls[nzchar(controls)]
    controls <- sort(intersect(controls, samples))
    if (length(controls) > 0L && setequal(controls, samples)) {
        controls <- character(0)
    }
    resolution_key <- tryCatch(
        as.character(bam_resolution_value(resolution, defaults)),
        error = function(e) as.character(resolution)
    )
    hash_object(list(
        example = isTRUE(preflight$example),
        paths = unique(as.character(preflight$table$path)),
        samples = samples,
        assembly = as.character(preflight$assembly),
        resolution = resolution_key,
        mapq = as.integer(defaults$mapq),
        exclude_flags = as.integer(defaults$exclude_flags),
        metric = as.character(defaults$metric),
        score = as.character(defaults$score),
        donor_aggregation = as.character(donor_aggregation),
        controls = controls
    ))
}

#' @noRd
bam_bedpe_filename <- function(sample) {
    if (is.null(sample) || length(sample) != 1L || is.na(sample) ||
        !nzchar(sample)) {
        return("mamamia-selected-introgressions.bedpe")
    }
    safe <- gsub("[^A-Za-z0-9._-]+", "-", sample)
    safe <- gsub("^[.-]+|[.-]+$", "", safe)
    if (!nzchar(safe)) {
        safe <- "sample"
    }
    paste0("mamamia-", safe, "-selected-introgressions.bedpe")
}

## UI ------------------------------------------------------------------------

#' @noRd
## Dashboard layout ----------------------------------------------------------

## The page is a three-column, three-row grid. The 2rem first column is a
## gutter; the header and footer span the sidebar and main columns. The grid
## keeps the log visible while the result tabs scroll independently.
#' @noRd
bam_dashboard_css <- function() {
    "
    .container-fluid {
        height: 100vh;
        padding: 0 !important;
    }
    html, body { height: 100%; }
    body { overflow: hidden; }

    .mamamia-header {
        display: flex;
        align-items: baseline;
        gap: 1rem;
        min-width: 0;
        height: 100%;
        padding: 0 .9rem;
        border-bottom: 1px solid var(--bs-border-color);
        background: var(--bs-tertiary-bg);
    }

    .mamamia-title {
        margin: 0;
        font-size: 1.95rem;
        font-weight: 800;
        line-height: 1.05;
        letter-spacing: -0.02em;
    }

    .mamamia-subtitle {
        margin: 0;
        color: var(--bs-secondary-color);
        font-size: .82rem;
    }

    .mamamia-sidebar {
        min-width: 0;
        min-height: 0;
        overflow-x: hidden;
        overflow-y: auto;
        padding: .65rem .75rem;
        border-right: 1px solid var(--bs-border-color);
        background: var(--bs-body-bg);
        font-size: .88rem;
        line-height: 1.35;
    }

    .mamamia-main {
        display: flex;
        flex-direction: column;
        min-width: 0;
        min-height: 0;
        overflow: hidden;
        padding: .35rem .5rem;
    }

    .mamamia-main > .tabbable,
    .mamamia-main > .card {
        flex: 1 1 auto;
        min-height: 0;
        margin-bottom: 0;
    }
    .mamamia-main .card-header {
        flex: 0 0 auto;
        padding: .2rem .45rem;
    }
    .mamamia-main .card-body {
        min-height: 0;
        overflow: auto;
    }
    .mamamia-main .tab-content {
        min-height: 0;
        overflow: hidden;
    }
    .mamamia-main .tab-pane {
        height: 100%;
        overflow: auto;
    }
    .mamamia-main .nav-link {
        padding: .35rem .6rem;
        font-size: .86rem;
    }

    .mamamia-footer {
        display: flex;
        flex-direction: column;
        gap: .15rem;
        min-width: 0;
        min-height: 0;
        overflow: hidden;
        padding: .25rem .5rem .3rem;
        border-top: 1px solid var(--bs-border-color);
        background: var(--bs-tertiary-bg);
    }

    .mamamia-footer .shiny-text-output {
        flex: 1 1 auto;
        margin: 0;
        padding: .15rem .4rem;
        overflow-y: auto;
        font-size: .7rem;
        line-height: 1.28;
        white-space: pre-wrap;
        background: var(--bs-body-bg);
        border: 1px solid var(--bs-border-color);
        border-radius: .25rem;
    }

    .mamamia-footer .footer-status {
        display: flex;
        flex: 0 0 auto;
        align-items: center;
        gap: .4rem;
        font-size: .74rem;
        line-height: 1.2;
        min-height: 1.1rem;
    }

    .mamamia-footer .mamamia-footer-status {
        flex: 1 1 auto;
        min-width: 0;
    }

    .mamamia-footer .mamamia-footer-toggle {
        flex: 0 0 auto;
        padding: .1rem .45rem;
        font-size: .7rem;
        line-height: 1.2;
    }

    /* Collapsing the log gives its height back to the sidebar and main area. */
    #mamamia-grid.mamamia-footer-collapsed {
        grid-template-rows: 75px 1fr 2.25rem;
    }
    #mamamia-grid.mamamia-footer-collapsed .mamamia-footer .shiny-text-output {
        display: none;
    }
    #mamamia-grid.mamamia-footer-collapsed #explorer-score_plot {
        height: calc(100vh - 341px) !important;
    }
    @media (max-width: 500px) {
        #mamamia-grid.mamamia-footer-collapsed {
            grid-template-rows: 85px 350px 2.25rem 350px;
        }
    }

    .mamamia-sidebar .shiny-input-container {
        width: 100% !important;
        margin-bottom: .6rem;
    }
    .mamamia-sidebar .form-label,
    .mamamia-sidebar .control-label,
    .mamamia-sidebar label {
        margin-bottom: .15rem;
        font-size: .8rem;
        font-weight: 600;
        line-height: 1.25;
    }
    .mamamia-sidebar .form-control,
    .mamamia-sidebar .form-select {
        min-height: 1.9rem;
        height: auto;
        padding: .22rem .4rem;
        font-size: .85rem;
        line-height: 1.3;
    }
    .mamamia-sidebar .btn {
        padding: .25rem .5rem;
        font-size: .82rem;
        line-height: 1.25;
    }
    .mamamia-sidebar .accordion-button {
        padding: .45rem .55rem;
        font-size: .84rem;
        font-weight: 600;
    }
    .mamamia-sidebar .accordion-body {
        padding: .55rem .6rem;
    }
    .mamamia-run-panel {
        margin-top: .75rem;
        padding-top: .75rem;
        border-top: 1px solid var(--bs-border-color);
    }
    .mamamia-run-panel .btn {
        width: 100%;
    }
    .mamamia-main > .mamamia-filter-bar {
        flex: 0 0 auto;
    }
    .mamamia-filter-bar {
        display: flex;
        align-items: flex-end;
        gap: .7rem;
        padding: .45rem .55rem;
        border: 1px solid var(--bs-border-color);
        border-radius: .375rem;
        background: var(--bs-tertiary-bg);
    }
    .mamamia-sample-picker {
        flex: 0 0 auto;
    }
    .mamamia-sample-picker .shiny-input-container,
    .mamamia-sample-picker .form-group {
        width: 12rem;
        margin-bottom: 0;
    }
    .mamamia-filter-heading {
        flex: 0 0 auto;
        padding-bottom: .35rem;
        font-size: .8rem;
        font-weight: 600;
        color: var(--bs-secondary-color);
        white-space: nowrap;
    }
    .mamamia-filter-controls {
        display: flex;
        flex: 1 1 auto;
        flex-wrap: wrap;
        align-items: flex-end;
        gap: .7rem;
        min-width: 0;
    }
    .mamamia-filter-controls .shiny-input-container {
        width: auto;
        min-width: 12rem;
        margin-bottom: 0;
    }
    .mamamia-filter-controls .form-check {
        margin-bottom: .35rem;
        white-space: nowrap;
    }
    .mamamia-segment-controls {
        display: flex;
        flex-wrap: wrap;
        align-items: flex-end;
        gap: .7rem;
        margin-bottom: .55rem;
    }
    .mamamia-segment-controls .shiny-input-container {
        margin-bottom: 0;
    }
    .mamamia-segment-controls .form-group {
        margin-bottom: 0;
    }
    .mamamia-export-panel {
        max-width: 48rem;
    }
    .mamamia-export-panel h5 {
        margin-top: .4rem;
        font-size: .95rem;
    }
    .mamamia-karyotype {
        height: 100%;
        min-height: 0;
    }
    .mamamia-karyotype .shiny-input-container {
        width: auto;
        margin-bottom: 0;
    }
    .mamamia-karyotype .control-label {
        margin-bottom: 0;
    }
    .mamamia-karyotype .radio-inline {
        margin-right: .75rem;
    }
    .mamamia-karyotype .html-widget {
        min-height: 0;
    }
    .mamamia-sidebar .form-check {
        margin-bottom: .25rem;
        min-height: 1.25rem;
    }
    .mamamia-sidebar .form-check-label,
    .mamamia-sidebar .radio label,
    .mamamia-sidebar .checkbox label {
        font-size: .8rem;
        font-weight: 400;
        line-height: 1.3;
    }
    .mamamia-sidebar .radio-inline + .radio-inline {
        margin-left: .6rem;
    }
    .mamamia-sidebar .table {
        margin-bottom: .25rem;
        font-size: .7rem;
    }
    .mamamia-sidebar .text-muted.small,
    .mamamia-sidebar small {
        font-size: .75rem;
        line-height: 1.35;
    }
    "
}

#' Single-file BAM chooser
#'
#' The dashboard runs one BAM per analysis, so the control never offers a
#' batch. `multiple = FALSE` removes multi-select from the native file picker;
#' the launcher's drop guard trims a multi-file drag-and-drop to its first file
#' and the server rejects a multi-file payload. A multi-sample BAM is a single
#' file and remains supported.
#'
#' @return A `shiny.tag` containing the labelled file input.
#' @noRd
bam_single_file_input <- function() {
    shiny::div(
        class = "mamamia-upload",
        shiny::fileInput(
            "bams",
            "BAM file",
            multiple = FALSE,
            accept = ".bam",
            width = "100%",
            buttonLabel = "Choose BAM file",
            placeholder = "No BAM selected"
        )
    )
}

#' @noRd
bam_sidebar_species_ui <- function(config, selection = NULL) {
    species_choices <- bam_species_choices(config)
    species_keys <- names(config$species)
    if (length(species_keys) < 2L) {
        return(shiny::p(
            class = "text-muted small",
            paste(
                "This registry has no species section, so the assembly is",
                "detected from the BAM header."
            )
        ))
    }

    recipient <- if (!is.null(selection$recipient)) {
        selection$recipient
    } else {
        species_keys[[1L]]
    }
    donor <- if (!is.null(selection$donor)) {
        selection$donor
    } else {
        species_keys[[2L]]
    }
    if (!recipient %in% species_keys) {
        recipient <- species_keys[[1L]]
    }
    if (!donor %in% species_keys) {
        donor <- species_keys[[2L]]
    }

    shiny::tagList(
        shiny::selectInput(
            "recipient_species", "Recipient species",
            choices = species_choices,
            selected = recipient
        ),
        shiny::selectInput(
            "donor_species", "Donor species",
            choices = species_choices,
            selected = donor
        ),
        shiny::uiOutput("component_summary_ui"),
        shiny::actionButton(
            "swap_species", "Swap donor and recipient",
            class = "btn-sm btn-outline-secondary",
            style = "width: 100%"
        )
    )
}

#' @noRd
bam_sidebar_details_ui <- function(config) {
    defaults <- config$defaults
    shiny::tagList(
        bslib::accordion(
            id = "bam_sidebar_accordion",
            multiple = FALSE,
            open = "Analysis",
            bslib::accordion_panel(
                "Analysis",
                shiny::radioButtons(
                    "resolution", "Resolution",
                    choices = bam_resolution_choices(defaults),
                    selected = bam_default_resolution(defaults)
                ),
            ),
            bslib::accordion_panel(
                "Assembly",
                shiny::uiOutput("species_ui")
            ),
            bslib::accordion_panel(
                "Advanced",
                shiny::selectInput(
                    "donor_aggregation", "Donor aggregation",
                    choices = c(
                        "Max" = "max",
                        "Mean" = "mean",
                        "Sum" = "sum"
                    ),
                    selected = defaults$donor_aggregation
                ),
                shiny::p(
                    class = "text-muted small mb-0",
                    sprintf(
                        paste(
                            "How donor components are combined per window.",
                            "Registry default: %s."
                        ),
                        defaults$donor_aggregation
                    )
                ),
                shiny::uiOutput("control_samples_ui")
            )
        ),
        shiny::div(
            class = "mamamia-run-panel",
            shiny::uiOutput("run_ui")
        )
    )
}

## The BAM-first result area. The view filters stay visible above all tabs;
## call filters belong to the segment table and exports have their own tab.
#' @noRd
bam_results_ui <- function() {
    ns <- shiny::NS("explorer")
    shiny::tagList(
        shiny::div(
            class = "mamamia-filter-bar",
            shiny::uiOutput(
                "sample_selector_ui", class = "mamamia-sample-picker"
            ),
            shiny::div(class = "mamamia-filter-heading", "Filters"),
            shiny::div(
                class = "mamamia-filter-controls",
                explorer_view_filters_ui("explorer")
            )
        ),
        bslib::navset_card_tab(
            id = ns("tabs"),
            bslib::nav_panel(
                "Karyotype",
                explorer_karyotype_ui(
                    "explorer", height = "calc(100vh - 390px)"
                )
            ),
            bslib::nav_panel(
                "Score profile",
                shiny::plotOutput(
                    ns("score_plot"), height = "calc(100vh - 335px)"
                )
            ),
            bslib::nav_panel(
                "Calls",
                shiny::div(
                    class = "mamamia-segment-controls",
                    explorer_call_filters_ui("explorer"),
                    explorer_selection_ui("explorer")
                ),
                DT::DTOutput(ns("calls"))
            ),
            bslib::nav_panel(
                "Export",
                shiny::div(
                    class = "mamamia-export-panel",
                    shiny::h5("Selected calls"),
                    explorer_export_ui("explorer"),
                    shiny::h5("Cohort"),
                    shiny::uiOutput("cohort_ui")
                )
            )
        )
    )
}

#' Grid layout template for the dashboard
#'
#' The 2rem first column is the gap size understood by `gridlayout`, so the
#' page has a 200 px sidebar, a fluid main area, a 150 px header and a 120 px
#' footer with a uniform 2rem gap.
#' @noRd
bam_dashboard_layout <- function() {
    paste(
        "|0.5rem|240px   |1fr      |",
        "|75px  |header  |header   |",
        "|1fr   |sidebar |main     |",
        "|100px |footer  |footer   |",
        sep = "\n"
    )
}

bam_launcher_ui <- function(state) {
    if (!isTRUE(state$ok)) {
        file_label <- bam_error_file_label(state$error$files)
        return(bslib::page_sidebar(
            title = "MaMaMIA",
            bslib::card(
                bslib::card_header("Configuration error"),
                shiny::p(state$error$message),
                if (!is.null(file_label)) {
                    shiny::p(class = "text-muted", file_label)
                },
                shiny::p(
                    class = "text-muted",
                    paste("Suggested fix:", state$error$suggestion)
                ),
                shiny::p(
                    class = "text-muted",
                    paste(
                        "The production registry is mounted read-only and",
                        "passed through the MAMAMIA_CONFIG environment",
                        "variable."
                    )
                )
            )
        ))
    }

    example_available <- bam_example_available()

    bslib::page_fluid(
        title = "MaMaMIA",
        theme = bslib::bs_theme(version = 5),
        shiny::tags$head(
            shiny::tags$style(shiny::HTML(bam_dashboard_css())),
            shiny::tags$script(shiny::HTML(
                paste(
                    "$(document).on('shiny:value', function(e) {",
                    "  if (e.target.id !== 'calculation_log') return;",
                    "  var el = document.getElementById('calculation_log');",
                    "  if (el) { el.scrollTop = el.scrollHeight; }",
                    "});",
                    "Shiny.addCustomMessageHandler('mamamia-log-scroll',",
                    "function(x) {",
                    "  setTimeout(function() {",
                    "    var el = document.getElementById('calculation_log');",
                    "    if (el) { el.scrollTop = el.scrollHeight; }",
                    "  }, 0);",
                    "});"
                )
            )),
            shiny::tags$script(shiny::HTML(
                paste(
                    "$(document).on('click', '#footer_toggle', function() {",
                    "  var grid = document.getElementById('mamamia-grid');",
                    "  if (!grid) return;",
                    "  var collapsed =",
                    "    grid.classList.toggle('mamamia-footer-collapsed');",
                    "  this.setAttribute(",
                    "    'aria-expanded', collapsed ? 'false' : 'true');",
                    "  this.textContent = collapsed ? 'Show log' : 'Hide log';",
                    "  this.setAttribute(",
                    "    'title',",
                    "    collapsed ?",
                    "      'Expand the calculation log' :",
                    "      'Collapse the calculation log');",
                    "  var log = document.getElementById('calculation_log');",
                    "  if (!collapsed && log) { log.scrollTop = log.scrollHeight; }",
                    "  setTimeout(function() {",
                    "    window.dispatchEvent(new Event('resize'));",
                    "  }, 0);",
                    "});"
                )
            )),
            shiny::tags$script(shiny::HTML(
                paste(
                    "(function() {",
                    "  var id = 'bams';",
                    "  function onlyFirst(files) {",
                    "    if (!files || files.length <= 1) return files;",
                    "    var transfer = new DataTransfer();",
                    "    transfer.items.add(files[0]);",
                    "    return transfer.files;",
                    "  }",
                    "  function uploadInput(node) {",
                    "    var zone = node && node.closest('.mamamia-upload');",
                    "    if (!zone) return null;",
                    "    var input = zone.querySelector('input[type=\"file\"]');",
                    "    return input && input.id === id ? input : null;",
                    "  }",
                    "  document.addEventListener('change', function(event) {",
                    "    var input = event.target;",
                    "    if (!input || input.id !== id || !input.files) return;",
                    "    if (input.files.length > 1) {",
                    "      input.files = onlyFirst(input.files);",
                    "    }",
                    "  }, true);",
                    "  document.addEventListener('drop', function(event) {",
                    "    var input = uploadInput(event.target);",
                    "    var files = event.dataTransfer && event.dataTransfer.files;",
                    "    if (!input || !files || files.length <= 1) return;",
                    "    event.preventDefault();",
                    "    event.stopPropagation();",
                    "    input.files = onlyFirst(files);",
                    "    input.dispatchEvent(new Event('change', { bubbles: true }));",
                    "  }, true);",
                    "})();"
                )
            ))
        ),
        gridlayout::grid_container(
            id = "mamamia-grid",
            layout = bam_dashboard_layout(),
            gridlayout::grid_place(
                "header",
                shiny::tags$header(
                    class = "mamamia-header",
                    shiny::div(
                        shiny::h1("MaMaMIA", class = "mamamia-title"),
                        shiny::p(
                            class = "mamamia-subtitle",
                            paste(
                                "Max Mapping for Meticulous Introgression",
                                "Analysis"
                            )
                        )
                    )
                )
            ),
            gridlayout::grid_place(
                "sidebar",
                shiny::tags$div(
                    class = "mamamia-sidebar",
                    shiny::uiOutput("bam_mode_note"),
                    shiny::uiOutput("bam_upload_ui"),
                    if (example_available) {
                        shiny::actionButton(
                            "example_data", "Example data",
                            class = "btn-outline-secondary btn-sm",
                            style = "width: 100%"
                        )
                    },
                    shiny::uiOutput("sidebar_recovery"),
                    shiny::uiOutput("sidebar_body")
                )
            ),
            gridlayout::grid_place(
                "main",
                shiny::tags$main(
                    class = "mamamia-main",
                    bam_results_ui()
                )
            ),
            gridlayout::grid_place(
                "footer",
                shiny::tags$footer(
                    class = "mamamia-footer",
                    shiny::div(
                        class = "footer-status",
                        shiny::div(
                            class = "mamamia-footer-status",
                            shiny::uiOutput("footer_status")
                        ),
                        shiny::tags$button(
                            id = "footer_toggle",
                            type = "button",
                            class = paste(
                                "btn btn-sm btn-outline-secondary",
                                "mamamia-footer-toggle"
                            ),
                            `aria-expanded` = "false",
                            `aria-controls` = "calculation_log",
                            title = "Expand the calculation log",
                            "Show log"
                        )
                    ),
                    shiny::verbatimTextOutput("calculation_log")
                )
            )
        ),
        shiny::tags$script(shiny::HTML(
            "document.getElementById('mamamia-grid')",
            ".classList.add('mamamia-footer-collapsed');"
        ))
    )
}

## Server --------------------------------------------------------------------

#' @noRd
bam_launcher_server <- function(state) {
    if (!isTRUE(state$ok)) {
        return(function(input, output, session) {})
    }
    function(input, output, session) {
        config <- state$config
        session_dir <- tempfile(
            pattern = "mamamia-session-", tmpdir = tempdir()
        )
        dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
        session$onSessionEnded(function() {
            unlink(session_dir, recursive = TRUE, force = TRUE)
        })

        example_available <- bam_example_available()
        example_mode <- shiny::reactiveVal(FALSE)
        sidebar_unlocked <- shiny::reactiveVal(FALSE)
        upload_epoch <- shiny::reactiveVal(0L)

        ## The bundled example has 20 kb contigs, so its 2-8 kb calls would
        ## disappear under the 5 Mb production filter. Use a toy-scale filter
        ## while example mode is active and restore the default afterwards.
        shiny::observeEvent(example_mode(), {
            shiny::updateNumericInput(
                session, "explorer-min_width",
                value = if (isTRUE(example_mode())) 0.001 else 5
            )
        }, ignoreInit = TRUE)

        result <- shiny::reactiveVal(NULL)
        last_error <- shiny::reactiveVal(NULL)
        result_cache <- bam_result_cache()
        active_key <- shiny::reactiveVal(NULL)

        log_entries <- shiny::reactiveVal(bam_initial_log())
        log_line <- function(...) {
            entry <- bam_log_entry(paste0(...))
            if (is.null(entry)) {
                return(invisible(NULL))
            }
            log_entries(c(log_entries(), entry))
            try(
                session$sendCustomMessage("mamamia-log-scroll", NULL),
                silent = TRUE
            )
            invisible(NULL)
        }

        output$bam_mode_note <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                return(shiny::p(
                    class = "text-muted small mb-2",
                    "Example data: one sample."
                ))
            }
            pf <- preflight_state()
            if (isTRUE(pf$empty)) {
                return(shiny::p(
                    class = "text-muted small mb-2",
                    paste(
                        "One BAM per run. Samples are detected from the @RG SM",
                        "tags in its header."
                    )
                ))
            }
            if (!isTRUE(pf$ok) || is.null(pf$table)) {
                return(shiny::p(
                    class = "text-danger small mb-2",
                    "Sample detection failed; see the calculation log."
                ))
            }
            n_samples <- nrow(pf$table)
            shiny::p(
                class = "text-muted small mb-2",
                sprintf(
                    "%d sample%s detected from @RG SM tags%s",
                    n_samples,
                    if (n_samples == 1L) "" else "s",
                    if (n_samples > 1L) {
                        "; cohort median-polish correction is used."
                    } else {
                        "; single-sample GC correction is used."
                    }
                )
            )
        })

        output$bam_upload_ui <- shiny::renderUI({
            upload_epoch()
            bam_single_file_input()
        })

        shiny::observeEvent(input$example_data, {
            if (!example_available) {
                return()
            }
            ## The action button starts at zero. Ignoring that initial value
            ## (instead of using `ignoreInit`) keeps the observer working when
            ## the button is the first input a test session receives.
            if (is.null(input$example_data) ||
                length(input$example_data) != 1L ||
                is.na(input$example_data) ||
                input$example_data < 1) {
                return()
            }
            if (isTRUE(example_mode())) {
                example_mode(FALSE)
                sidebar_unlocked(FALSE)
                result(NULL)
                last_error(NULL)
                upload_epoch(upload_epoch() + 1L)
                shiny::updateActionButton(
                    session, "example_data", label = "Example data"
                )
                shiny::updateRadioButtons(
                    session, "resolution",
                    choices = bam_resolution_choices(config$defaults),
                    selected = bam_default_resolution(config$defaults)
                )
                log_line(
                    "Example data cleared. Upload a BAM to start a new run."
                )
                return()
            }
            example_mode(TRUE)
            sidebar_unlocked(TRUE)
            result(NULL)
            last_error(NULL)
            upload_epoch(upload_epoch() + 1L)
            shiny::updateActionButton(
                session, "example_data", label = "Clear example data"
            )
            shiny::updateRadioButtons(
                session, "resolution",
                choices = bam_example_resolution_choices(),
                selected = bam_example_default_resolution()
            )
            log_line(bam_example_note())
            log_line("Example data selected.")
        }, ignoreInit = FALSE)

        shiny::observeEvent(input$unlock_setup, {
            sidebar_unlocked(TRUE)
            log_line("Setup panels opened after a preflight failure.")
        }, ignoreInit = TRUE)

        uploaded_paths <- shiny::reactive({
            tryCatch(
                bam_single_upload_path(input$bams),
                MaMaMIAUploadError = function(e) {
                    log_line(conditionMessage(e))
                    shiny::showNotification(
                        conditionMessage(e),
                        type = "error",
                        duration = 10
                    )
                    character()
                }
            )
        })

        shiny::observeEvent(uploaded_paths(), {
            if (length(uploaded_paths()) > 0L && isTRUE(example_mode())) {
                example_mode(FALSE)
                sidebar_unlocked(FALSE)
                shiny::updateActionButton(
                    session, "example_data", label = "Example data"
                )
                shiny::updateRadioButtons(
                    session, "resolution",
                    choices = bam_resolution_choices(config$defaults),
                    selected = bam_default_resolution(config$defaults)
                )
                log_line(
                    "Uploaded BAM selected; example data is no longer active."
                )
            }
            if (length(uploaded_paths()) == 0L && !isTRUE(example_mode())) {
                sidebar_unlocked(FALSE)
            }
        }, ignoreInit = TRUE, ignoreNULL = FALSE)

        selected_paths <- shiny::reactive({
            paths <- uploaded_paths()
            paths <- paths[!is.na(paths) & nzchar(paths)]
            utils::head(unique(paths), 1L)
        })

        staging <- local({
            cache <- NULL
            function(paths) {
                key <- paste(paths, collapse = "\n")
                if (!is.null(cache) && identical(cache$key, key)) {
                    return(cache$value)
                }
                if (!is.null(cache)) {
                    unlink(cache$value$dir, recursive = TRUE, force = TRUE)
                }
                value <- bam_stage_paths(paths, session_dir)
                cache <<- list(key = key, value = value)
                value
            }
        })

        staged_bams <- shiny::reactive({
            paths <- selected_paths()
            if (length(paths) == 0L) {
                return(NULL)
            }
            staging(paths)
        })

        species_keys <- names(config$species)
        species_available <- length(species_keys) >= 2L
        default_recipient <- if (length(species_keys) > 0L) {
            species_keys[[1L]]
        } else {
            NULL
        }
        default_donor <- if (length(species_keys) > 1L) {
            species_keys[[2L]]
        } else {
            default_recipient
        }

        recipient_species <- shiny::reactive({
            if (isTRUE(example_mode())) {
                return(bam_example_species_selection()$recipient)
            }
            value <- input$recipient_species
            if (is.null(value) || length(value) != 1L || is.na(value) ||
                !value %in% species_keys) {
                return(default_recipient)
            }
            as.character(value)
        })
        donor_species <- shiny::reactive({
            if (isTRUE(example_mode())) {
                return(bam_example_species_selection()$donor)
            }
            value <- input$donor_species
            if (is.null(value) || length(value) != 1L || is.na(value) ||
                !value %in% species_keys) {
                return(default_donor)
            }
            as.character(value)
        })

        shiny::observeEvent(input$swap_species, {
            if (!species_available) {
                return()
            }
            recipient_key <- recipient_species()
            donor_key <- donor_species()
            if (is.null(recipient_key) || is.null(donor_key)) {
                return()
            }
            shiny::updateSelectInput(
                session, "recipient_species", selected = donor_key
            )
            shiny::updateSelectInput(
                session, "donor_species", selected = recipient_key
            )
        })

        output$component_summary_ui <- shiny::renderUI({
            species_config <- if (isTRUE(example_mode())) {
                bam_example_config()
            } else {
                config
            }
            if (isTRUE(example_mode())) {
                selection <- bam_example_species_selection()
                recipient_key <- selection$recipient
                donor_key <- selection$donor
            } else {
                recipient_key <- recipient_species()
                donor_key <- donor_species()
            }
            recipient <- species_config$species[[recipient_key]]
            donor <- species_config$species[[donor_key]]
            if (is.null(recipient) || is.null(donor)) {
                return(NULL)
            }
            recipient_components <- bam_species_components(recipient)
            donor_components <- bam_species_components(donor)
            shiny::p(
                class = "text-muted small mb-2",
                sprintf(
                    "All components: recipient %s; donor %s.",
                    paste(recipient_components, collapse = "+"),
                    paste(donor_components, collapse = "+")
                )
            )
        })

        assembly_selection <- shiny::reactive({
            if (!species_available) {
                ## Header-based auto-detection is the only mode for a registry
                ## without species.
                return(list(
                    ok = TRUE, error = NULL, config = config, key = NULL,
                    assembly = NULL, unpaired = NULL
                ))
            }

            empty <- list(
                ok = FALSE, error = NULL, config = config,
                key = NULL, assembly = NULL, unpaired = NULL
            )
            recipient_key <- recipient_species()
            donor_key <- donor_species()
            if (is.null(recipient_key) || is.null(donor_key)) {
                empty$error <- tryCatch(
                    .config_error(
                        "Choose a recipient species and a donor species"
                    ),
                    error = function(e) e
                )
                return(empty)
            }
            recipient <- config$species[[recipient_key]]
            donor <- config$species[[donor_key]]
            if (is.null(recipient) || is.null(donor)) {
                empty$error <- tryCatch(
                    .config_error(
                        "Unknown species selection ('%s' recipient, '%s' donor)",
                        recipient_key, donor_key
                    ),
                    error = function(e) e
                )
                return(empty)
            }
            resolved <- tryCatch(
                bam_resolve_species_pair(config, recipient_key, donor_key),
                error = function(e) e
            )
            if (inherits(resolved, "condition")) {
                empty$error <- resolved
                return(empty)
            }
            list(
                ok = TRUE, error = NULL, config = resolved$config,
                key = resolved$key, assembly = resolved$assembly,
                unpaired = resolved$unpaired
            )
        })

        preflight_state <- shiny::reactive({
            if (example_mode()) {
                return(tryCatch(
                    bam_example_preflight_state(),
                    error = function(e) {
                        bam_preflight_failure(e, character(), config)
                    }
                ))
            }
            selection <- assembly_selection()
            if (!isTRUE(selection$ok)) {
                failure <- bam_preflight_failure(
                    selection$error, character(), selection$config
                )
                failure$config <- selection$config
                failure$key <- selection$key
                failure$unpaired <- selection$unpaired
                return(failure)
            }
            paths <- selected_paths()
            if (length(paths) == 0L) {
                return(list(
                    ok = FALSE, empty = TRUE, table = NULL, staged = NULL,
                    assembly = NULL, assembly_config = NULL,
                    config = selection$config, key = selection$key,
                    unpaired = selection$unpaired,
                    error = NULL, candidates = NULL
                ))
            }
            staged <- tryCatch(staged_bams(), error = function(e) e)
            if (inherits(staged, "condition")) {
                return(bam_preflight_failure(
                    staged, paths, selection$config, files = paths
                ))
            }
            table <- tryCatch(
                preflightBams(
                    staged$path, selection$config,
                    assembly = selection$key
                ),
                error = function(e) e
            )
            if (inherits(table, "condition")) {
                return(bam_preflight_failure(
                    table, staged$path, selection$config,
                    files = staged$source
                ))
            }
            table$source <- staged$source
            key <- table$assembly[[1L]]
            list(
                ok = TRUE, empty = FALSE, table = table,
                staged = staged$path, assembly = key,
                assembly_config = selection$config$assemblies[[key]],
                config = selection$config, key = key,
                unpaired = selection$unpaired,
                error = NULL, candidates = NULL
            )
        })

        ## The setup controls stay hidden until the first successful
        ## preflight. A failure can unlock them explicitly through the
        ## recovery panel, so a wrong species choice can be corrected.
        shiny::observe({
            pf <- preflight_state()
            if (isTRUE(pf$ok) && !isTRUE(sidebar_unlocked())) {
                sidebar_unlocked(TRUE)
            }
        })

        preflight_log_signature <- shiny::reactive({
            pf <- preflight_state()
            list(
                ok = isTRUE(pf$ok),
                empty = isTRUE(pf$empty),
                example = isTRUE(pf$example),
                assembly = pf$assembly,
                samples = if (is.null(pf$table)) NULL else pf$table$sample,
                warnings = if (is.null(pf$table)) NULL else pf$table$warnings,
                error = if (is.null(pf$error)) {
                    NULL
                } else {
                    paste(
                        pf$error$stage, pf$error$message, pf$error$suggestion
                    )
                },
                candidates = if (is.null(pf$candidates)) {
                    NULL
                } else {
                    pf$candidates$key
                }
            )
        })
        shiny::observeEvent(preflight_log_signature(), {
            for (line in bam_preflight_log(preflight_state())) {
                log_line(line)
            }
        }, ignoreInit = FALSE)

        output$sidebar_body <- shiny::renderUI({
            if (!isTRUE(sidebar_unlocked())) {
                return(NULL)
            }
            bam_sidebar_details_ui(config)
        })

        output$control_samples_ui <- shiny::renderUI({
            pf <- preflight_state()
            if (!isTRUE(pf$ok) || is.null(pf$table) ||
                !"sample" %in% names(pf$table)) {
                return(NULL)
            }
            samples <- unique(as.character(pf$table$sample))
            samples <- samples[!is.na(samples) & nzchar(samples)]
            if (length(samples) < 2L) {
                return(shiny::p(
                    class = "text-muted small mb-0",
                    "Control samples apply only to multi-sample cohorts."
                ))
            }
            selected <- shiny::isolate(input$control_samples)
            selected <- intersect(as.character(selected), samples)
            shiny::tagList(
                shiny::selectizeInput(
                    "control_samples", "Calibration controls",
                    choices = samples,
                    selected = selected,
                    multiple = TRUE,
                    options = list(
                        placeholder = "All samples (default)",
                        plugins = list("remove_button")
                    )
                ),
                shiny::p(
                    class = "text-muted small mb-0",
                    paste(
                        "Used only to calibrate call thresholds.",
                        "Leave empty to pool all samples."
                    )
                )
            )
        })

        output$sidebar_recovery <- shiny::renderUI({
            if (isTRUE(sidebar_unlocked())) {
                return(NULL)
            }
            pf <- preflight_state()
            if (isTRUE(pf$empty) || isTRUE(pf$ok)) {
                return(NULL)
            }
            info <- pf$error
            if (is.null(info)) {
                return(NULL)
            }
            shiny::tagList(
                shiny::p(
                    class = "text-danger small mb-1",
                    sprintf(
                        "Preflight failed (stage: %s): %s",
                        info$stage, info$message
                    )
                ),
                shiny::actionButton(
                    "unlock_setup", "Show setup",
                    class = "btn-sm btn-outline-secondary"
                )
            )
        })

        output$footer_status <- shiny::renderUI({
            pf <- preflight_state()
            if (isTRUE(pf$ok)) {
                table <- pf$table
                n_bams <- length(unique(table$path))
                label <- "Preflight OK"
                badge <- "text-bg-success"
                detail <- sprintf(
                    "%d BAM file%s",
                    n_bams,
                    if (n_bams == 1L) "" else "s"
                )
            } else if (isTRUE(pf$empty)) {
                label <- "Waiting"
                badge <- "text-bg-secondary"
                detail <- "Upload a BAM or use the example data."
            } else {
                label <- "Preflight failed"
                badge <- "text-bg-danger"
                detail <- if (is.null(pf$error)) {
                    "See the calculation log."
                } else {
                    sprintf(
                        "Stage '%s': %s",
                        pf$error$stage, pf$error$message
                    )
                }
            }
            if (isTRUE(example_mode())) {
                detail <- paste("Example data |", detail)
            }
            shiny::div(
                class = "d-flex align-items-center gap-2",
                shiny::span(class = paste("badge", badge), label),
                shiny::span(class = "text-truncate", title = detail, detail)
            )
        })

        output$calculation_log <- shiny::renderText({
            paste(log_entries(), collapse = "\n")
        })

        output$species_ui <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                return(bam_sidebar_species_ui(
                    bam_example_config(),
                    selection = bam_example_species_selection()
                ))
            }
            bam_sidebar_species_ui(config)
        })

        output$preflight_table <- shiny::renderTable({
            pf <- preflight_state()
            if (isTRUE(pf$empty)) {
                return(data.frame(
                    status = paste(
                        "Upload a BAM or select the bundled example",
                        "to see the preflight."
                    ),
                    stringsAsFactors = FALSE
                ))
            }
            if (!isTRUE(pf$ok)) {
                return(data.frame(
                    status = "Preflight failed; see the calculation log.",
                    stringsAsFactors = FALSE
                ))
            }
            bam_preflight_display(pf$table)
        },
        striped = TRUE, bordered = TRUE, hover = TRUE, spacing = "xs",
        rownames = FALSE, width = "100%"
        )

        output$assembly_summary <- shiny::renderUI({
            pf <- preflight_state()
            if (isTRUE(pf$empty)) {
                return(shiny::p(
                    class = "text-muted",
                    paste(
                        "Upload a BAM or select the bundled example",
                        "to detect an assembly."
                    )
                ))
            }
            if (!isTRUE(pf$ok)) {
                return(shiny::p(
                    class = "text-muted",
                    "No assembly was detected. See the calculation log."
                ))
            }
            assembly <- pf$assembly_config
            if (is.null(assembly)) {
                return(shiny::p(
                    class = "text-muted",
                    "No assembly is available for the current selection."
                ))
            }
            shiny::tagList(
                bam_assembly_summary_ui(
                    assembly,
                    if (!is.null(pf$config)) {
                        pf$config$defaults$allowed_windows
                    } else {
                        config$defaults$allowed_windows
                    }
                ),
                if (!is.null(pf$unpaired) && nrow(pf$unpaired) > 0L) {
                    shiny::p(
                        class = "text-muted",
                        sprintf(
                            "Not paired (%d): %s",
                            nrow(pf$unpaired),
                            paste(
                                sprintf(
                                    "%s %s (%s)",
                                    pf$unpaired$role, pf$unpaired$contig,
                                    pf$unpaired$reason
                                ),
                                collapse = "; "
                            )
                        )
                    )
                }
            )
        })

        example_resolution_problem <- function(resolution) {
            if (!example_mode()) {
                return(NULL)
            }
            resolutions <- bam_example_resolutions()
            if (length(resolution) != 1L || is.na(resolution) ||
                !resolution %in% resolutions) {
                return(sprintf(
                    "The bundled example supports resolution(s): %s.",
                    paste(resolutions, collapse = ", ")
                ))
            }
            NULL
        }

        output$run_ui <- shiny::renderUI({
            pf <- preflight_state()
            ready <- isTRUE(pf$ok)
            problem <- NULL
            if (ready) {
                resolution <- suppressWarnings(as.integer(input$resolution))
                problem <- bam_resolution_problem(
                    pf$assembly_config, resolution
                )
                if (is.null(problem)) {
                    problem <- example_resolution_problem(resolution)
                }
                ready <- is.null(problem)
            }
            shiny::tagList(
                shiny::tags$button(
                    id = "run", type = "button",
                    class = paste(
                        "btn btn-primary action-button",
                        if (!ready) "disabled"
                    ),
                    disabled = if (ready) NULL else "",
                    shiny::span(class = "action-label", "Run analysis")
                ),
                if (!is.null(problem)) {
                    shiny::p(class = "text-danger", problem)
                }
            )
        })

        ## A completed analysis is cached under its resolved assembly and run
        ## parameters. Returning to a combination that was already analysed
        ## restores that result instead of clearing it and recomputing every
        ## sample. A combination that has never been run still needs one run.
        analysis_key <- shiny::reactive({
            bam_analysis_key(
                preflight_state(),
                input$resolution,
                input$donor_aggregation,
                input$control_samples
            )
        })

        shiny::observeEvent(analysis_key(), {
            key <- analysis_key()
            cached <- result_cache$get(key)
            if (!is.null(cached)) {
                active_key(key)
                result(cached)
                last_error(NULL)
                log_line(
                    "Restored a cached analysis result for these inputs; ",
                    "no recalculation was needed."
                )
                return()
            }
            had_result <- !is.null(result())
            active_key(NULL)
            result(NULL)
            last_error(NULL)
            if (had_result) {
                log_line(
                    "Analysis inputs changed; the previous result was cleared. ",
                    "Run the analysis for the new inputs."
                )
            }
        }, ignoreInit = TRUE)

        ## Samples are known as soon as preflight succeeds, so the selector can
        ## appear before the run and keep the chosen sample through the run.
        result_samples <- shiny::reactive({
            res <- result()
            if (!is.null(res)) {
                samples <- bam_sample_choices(res)
                if (length(samples) > 0L) {
                    return(samples)
                }
            }
            pf <- preflight_state()
            if (!isTRUE(pf$ok) || is.null(pf$table) ||
                !"sample" %in% names(pf$table)) {
                return(character())
            }
            unique(as.character(pf$table$sample))
        })

        shiny::observeEvent(result_samples(), {
            samples <- result_samples()
            if (length(samples) == 0L) {
                return()
            }
            current <- shiny::isolate(input$result_sample)
            selected <- if (!is.null(current) && length(current) == 1L &&
                !is.na(current) && current %in% samples) {
                current
            } else {
                samples[[1L]]
            }
            shiny::updateSelectInput(
                session, "result_sample",
                choices = samples, selected = selected
            )
        }, ignoreNULL = FALSE)

        result_sample <- shiny::reactive({
            samples <- result_samples()
            if (length(samples) == 0L) {
                return(NULL)
            }
            value <- input$result_sample
            if (is.null(value) || length(value) != 1L || is.na(value) ||
                !value %in% samples) {
                return(samples[[1L]])
            }
            value
        })

        result_view <- shiny::reactive({
            result()
        })

        shiny::observeEvent(input$run, {
            pf <- preflight_state()
            if (!isTRUE(pf$ok)) {
                log_line(paste(
                    "Run refused: preflight is not ready.",
                    "Fix the problems reported in this log first."
                ))
                shiny::showNotification(
                    "Fix the preflight problems before running the analysis.",
                    type = "warning", duration = 8
                )
                return()
            }
            samples <- unique(pf$table$sample)
            params <- tryCatch(
                bam_run_parameters(input, pf$config, samples),
                error = function(e) e
            )
            if (inherits(params, "condition")) {
                info <- bam_error_info(
                    params, files = pf$table$source, table = pf$table
                )
                last_error(info)
                log_line(sprintf(
                    "Run parameters rejected at stage '%s': %s",
                    info$stage, info$message
                ))
                log_line(sprintf("Suggested fix: %s", info$suggestion))
                shiny::showNotification(
                    bam_error_text(info),
                    type = "error", duration = 12
                )
                return()
            }

            resolution_problem <- bam_resolution_problem(
                pf$assembly_config, params$resolution
            )
            if (is.null(resolution_problem)) {
                resolution_problem <- example_resolution_problem(
                    params$resolution
                )
            }
            if (!is.null(resolution_problem)) {
                info <- bam_error_info(
                    .mamamia_error(
                        "MaMaMIAWindowError", "%s", resolution_problem
                    ),
                    files = pf$table$source, table = pf$table
                )
                last_error(info)
                log_line(sprintf(
                    "Run refused at stage '%s': %s",
                    info$stage, info$message
                ))
                shiny::showNotification(
                    bam_error_text(info),
                    type = "error", duration = 12
                )
                return()
            }

            previous <- result()
            key <- analysis_key()
            cached <- result_cache$get(key)
            if (!is.null(cached)) {
                if (!identical(active_key(), key)) {
                    log_line(
                        "Using the cached analysis result; no recalculation ",
                        "was needed."
                    )
                    shiny::showNotification(
                        "Using the cached analysis result for these inputs.",
                        type = "message", duration = 5
                    )
                    active_key(key)
                    result(cached)
                }
                last_error(NULL)
                return()
            }
            log_line(sprintf(
                paste(
                    "Starting analysis: %d sample%s, %s windows/chromosome,",
                    "metric '%s', score '%s', donor aggregation '%s'."
                ),
                length(samples),
                if (length(samples) == 1L) "" else "s",
                format(params$resolution, big.mark = ",", scientific = FALSE),
                params$metric, params$score, params$donor_aggregation
            ))
            if (length(params$control_samples) > 0L) {
                log_line(paste(
                    "Calibration controls:",
                    paste(params$control_samples, collapse = ", ")
                ))
            }

            result_cache$reserve(result_cache$size() + 1L)

            stages_per_run <- length(bam_progress_stages(length(samples)))
            total_steps <- stages_per_run
            run_analysis <- function(assembly_config, assembly_key, progress) {
                if (isTRUE(pf$example)) {
                    example_result <- bam_example_result(params$resolution)
                    bam_example_progress(
                        example_result, progress, pause = 0.1
                    )
                    return(example_result)
                }
                analyzeBamInput(
                    pf$staged,
                    assembly_config,
                    assembly = assembly_key,
                    windows_per_chromosome = params$resolution,
                    metric = params$metric,
                    score = params$score,
                    donor_aggregation = params$donor_aggregation,
                    calibration = params$calibration,
                    mapq = params$mapq,
                    exclude_flags = params$exclude_flags,
                    cores = params$cores,
                    progress = progress
                )
            }
            progress_for <- function(prefix = "") {
                raw_progress <- bam_progress_callback(
                    length(samples),
                    total_steps = total_steps
                )
                function(stage, detail = list()) {
                    log_line(paste0(
                        prefix, bam_progress_label(stage, detail)
                    ))
                    raw_progress(stage, detail)
                }
            }

            new_result <- tryCatch(
                shiny::withProgress(
                    message = "MaMaMIA analysis", value = 0,
                    {
                        primary <- run_analysis(
                            pf$config, pf$assembly, progress_for()
                        )
                        result_cache$set(key, primary)
                        active_key(key)
                        primary
                    }
                ),
                error = function(e) {
                    info <- bam_error_info(
                        e, files = pf$table$source, table = pf$table
                    )
                    last_error(info)
                    log_line(sprintf(
                        "Run failed at stage '%s': %s",
                        info$stage, info$message
                    ))
                    log_line(sprintf("Suggested fix: %s", info$suggestion))
                    shiny::showNotification(
                        bam_error_text(info),
                        type = "error", duration = 20
                    )
                    NULL
                }
            )
            if (is.null(new_result)) {
                result(previous)
                return()
            }
            last_error(NULL)
            result_cache$set(key, new_result)
            active_key(key)
            result(new_result)
            log_line(sprintf(
                "Analysis finished: %d sample(s), %d calibrated call(s).",
                length(bam_sample_choices(new_result)),
                sum(new_result$samples$n_calls, na.rm = TRUE)
            ))
            run_summary <- sprintf(
                "%d sample(s), %d calibrated call(s)",
                length(bam_sample_choices(new_result)),
                sum(new_result$samples$n_calls, na.rm = TRUE)
            )
            shiny::showNotification(
                run_summary,
                type = "message",
                duration = 6
            )
            coverage <- bam_coverage_warnings(new_result)
            for (warning in coverage) {
                log_line(paste("Coverage warning:", warning))
            }
        })

        ## The sample selector belongs to the main result area: it controls the
        ## karyotype, score profile, call table and exports. It is only useful
        ## when the input BAM contains several @RG SM samples.
        output$sample_selector_ui <- shiny::renderUI({
            samples <- result_samples()
            if (length(samples) <= 1L) {
                return(NULL)
            }
            shiny::selectInput(
                "result_sample", "Sample",
                choices = samples,
                selected = shiny::isolate(result_sample())
            )
        })

        output$cohort_ui <- shiny::renderUI({
            res <- result()
            if (is.null(res)) {
                return(shiny::p(
                    class = "text-muted small",
                    "Run the analysis to enable the exports."
                ))
            }
            samples <- bam_sample_choices(res)
            warnings <- bam_coverage_warnings(res)
            multiple_samples <- length(samples) > 1L
            shiny::tagList(
                shiny::div(
                    class = "table-responsive",
                    shiny::tableOutput("cohort_table")
                ),
                if (length(warnings) > 0L) {
                    shiny::div(
                        class = "alert alert-warning small",
                        shiny::tags$ul(lapply(warnings, shiny::tags$li))
                    )
                },
                shiny::div(
                    class = "d-flex flex-column gap-1",
                    shiny::downloadButton(
                        "download_sample",
                        "Download BEDPE (current sample)",
                        class = "btn-sm"
                    ),
                    if (multiple_samples) {
                        shiny::downloadButton(
                            "download_cohort",
                            "Download selected BEDPE (all samples)",
                            class = "btn-sm"
                        )
                    }
                ),
                shiny::p(
                    class = "text-muted small",
                    paste(
                        "The per-sample export uses the rows selected for the",
                        "current sample, or all calibrated calls when none are",
                        "selected. The cohort export includes only selected",
                        "rows and adds a `sample` column."
                    )
                )
            )
        })

        output$cohort_table <- shiny::renderTable({
            res <- result()
            shiny::validate(
                shiny::need(
                    !is.null(res),
                    "Run the analysis to see the cohort summary"
                )
            )
            bam_cohort_summary(res)
        },
        striped = TRUE, bordered = TRUE, spacing = "xs",
        rownames = FALSE, digits = 0
        )

        output$download_sample <- shiny::downloadHandler(
            filename = function() {
                bam_bedpe_filename(result_sample())
            },
            content = function(file) {
                res <- result()
                sample <- result_sample()
                if (is.null(res) || is.null(sample)) {
                    stop(
                        "Run the analysis before downloading",
                        call. = FALSE
                    )
                }
                writeCallsToBedpe(
                    res, file,
                    samples = sample,
                    pass_only = TRUE
                )
            }
        )

        output$download_cohort <- shiny::downloadHandler(
            filename = function() {
                "mamamia-cohort-selected-introgressions.bedpe"
            },
            content = function(file) {
                res <- result()
                if (is.null(res)) {
                    stop(
                        "Run the analysis before downloading",
                        call. = FALSE
                    )
                }
                writeCallsToBedpe(
                    res, file,
                    pass_only = TRUE
                )
            }
        )

        explorer_server(
            "explorer",
            result = result_view,
            sample = result_sample
        )
    }
}

#' @noRd
bam_launcher_app <- function(config = NULL) {
    require_dashboard_packages()
    state <- bam_launcher_state(config)
    ui <- bam_launcher_ui(state)
    server <- bam_launcher_server(state)
    shiny::shinyApp(ui, server)
}
