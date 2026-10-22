## BAM-first Shiny application (T6.1-T6.4).
##
## The app is deliberately thin: configuration, preflight, run orchestration
## and cohort export live here, while counting, correction, normalization and
## segmentation stay in the pipeline. Uploaded BAMs are staged by symlink into
## a per-session directory before preflight, so the pipeline never writes an
## index into the upload directory and the session can remove every staged
## file and index when it ends. A bundled pre-computed example is available
## without a BAM upload.

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
    counts = paste(
        "Check MAPQ and the exclude flags;",
        "verify that the BAM is intact and indexed."
    ),
    correction = paste(
        "Too few valid or ideal windows for the GC fit: lower the resolution,",
        "raise the max fit windows cap, or inspect coverage and GC outliers."
    ),
    normalization = paste(
        "Pick a reference sample from the input BAMs,",
        "or use median normalization."
    ),
    segmentation = paste(
        "Check that both the donor and the recipient side",
        "have enough ideal windows."
    ),
    limits = paste(
        "Reduce the number of BAMs or their total size; the limits come from",
        "server.max_bams and server.max_total_bytes."
    ),
    analysis = paste(
        "Check the run controls (MAPQ, exclude flags, cores, max fit windows)",
        "and the server log."
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
    if (inherits(e, "MaMaMIACountError")) {
        return("counts")
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
                counts = "counts",
                RCA = "counts",
                "GC correction" = "correction",
                segmentation = "segmentation",
                stage
            ))
        }
        return("analysis")
    }
    if (inherits(e, "MaMaMIAAnalysisError") &&
        grepl("^(Normalization failed|Reference sample)", message)) {
        return("normalization")
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
bam_config_labels <- function(labels, max_items = 4L) {
    if (length(labels) <= max_items) {
        return(paste(labels, collapse = "; "))
    }
    paste(
        c(
            utils::head(labels, max_items),
            sprintf("+%d more", length(labels) - max_items)
        ),
        collapse = "; "
    )
}

#' @noRd
bam_config_summary <- function(config) {
    species <- listMaMaMIASpecies(config)
    if (nrow(species) > 0L) {
        labels <- bam_config_labels(
            sprintf("%s (%s)", species$label, species$key)
        )
        return(sprintf("%d species: %s", nrow(species), labels))
    }

    assemblies <- listMaMaMIAAssemblies(config)
    if (nrow(assemblies) > 0L) {
        labels <- bam_config_labels(
            sprintf("%s (%s)", assemblies$label, assemblies$key)
        )
        return(sprintf(
            "%d assembl%s, detected from the BAM header: %s",
            nrow(assemblies),
            if (nrow(assemblies) == 1L) "y" else "ies",
            labels
        ))
    }
    "No species or assemblies configured"
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

#' Resolve the subgenome used for a species selector
#'
#' Falls back to the species default (or its first subgenome) when the current
#' input is missing or belongs to a previously selected species.
#' @noRd
bam_effective_subgenome <- function(species, requested) {
    available <- unique(as.character(species$contigs$subgenome))
    if (length(available) == 0L) {
        return(NULL)
    }
    if (!is.null(requested) && length(requested) == 1L &&
        !is.na(requested) && requested %in% available) {
        return(as.character(requested))
    }
    if (!is.null(species$default_subgenome) &&
        species$default_subgenome %in% available) {
        return(species$default_subgenome)
    }
    available[[1L]]
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

#' @noRd


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
            lines <- c(lines, paste0(
                "Assembly candidates: ",
                paste(sprintf(
                    "%s (%d/%d targets, %.1f%%)",
                    candidates$key, candidates$n_matched,
                    candidates$n_targets, 100 * candidates$score
                ), collapse = "; ")
            ))
        }
        return(lines)
    }

    table <- pf$table
    n <- if (is.null(table)) 0L else nrow(table)
    lines <- character()
    if (isTRUE(pf$example)) {
        lines <- c(
            lines,
            "Example data loaded: pre-computed tables, no BAM is read."
        )
    }
    lines <- c(lines, sprintf(
        "Preflight OK: %d BAM%s, assembly '%s'.",
        n, if (n == 1L) "" else "s", pf$assembly
    ))
    if (!is.null(table) && n > 0L) {
        lines <- c(lines, sprintf(
            "Samples: %s", paste(unique(table$sample), collapse = ", ")
        ))
        warnings <- unique(table$warnings[nzchar(table$warnings)])
        if (length(warnings) > 0L) {
            lines <- c(lines, paste0("Preflight note: ", warnings))
        }
    }

    assembly <- pf$assembly_config
    if (!is.null(assembly)) {
        lines <- c(lines, sprintf(
            "Detected assembly: %s (%s); %d contig pair%s.",
            assembly$label, assembly$key, length(assembly$pairs),
            if (length(assembly$pairs) == 1L) "" else "s"
        ))
        if (length(assembly$pairs) > 0L) {
            lines <- c(lines, paste0(
                "Pairs: ",
                paste(
                    sprintf(
                        "%s -> %s", names(assembly$pairs), unname(assembly$pairs)
                    ),
                    collapse = "; "
                )
            ))
        }
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
bam_initial_log <- function(config, time = Sys.time()) {
    c(
        bam_log_entry("MaMaMIA session started.", time),
        bam_log_entry(bam_config_summary(config), time),
        bam_log_entry("Input mode: single BAM.", time)
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

#' @noRd
bam_candidate_display <- function(candidates) {
    if (is.null(candidates) || nrow(candidates) == 0L) {
        return(data.frame(
            assembly = character(), label = character(),
            targets = character(), score = character(),
            stringsAsFactors = FALSE
        ))
    }
    data.frame(
        assembly = candidates$key,
        label = candidates$label,
        targets = sprintf("%d/%d", candidates$n_matched, candidates$n_targets),
        score = sprintf("%.1f%%", 100 * candidates$score),
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
        .mamamia_error("MaMaMIAAnalysisError", "At least one BAM is required")
    }
    resolution <- bam_resolution_value(input$resolution, config$defaults)
    mapq <- bam_integer(input$mapq, "MAPQ", 0L, 255L)
    exclude_flags <- bam_integer(input$exclude_flags, "Exclude flags", 0L, 4095L)
    cores <- bam_integer(input$cores, "Cores", 1L, bam_max_cores())
    max_fit_windows <- input$max_fit_windows
    if (length(max_fit_windows) != 1L || is.na(max_fit_windows) ||
        !is.numeric(max_fit_windows) || max_fit_windows <= 0 ||
        (!is.infinite(max_fit_windows) &&
            max_fit_windows != floor(max_fit_windows))) {
        .mamamia_error(
            "MaMaMIAAnalysisError",
            "Max fit windows must be a positive whole number"
        )
    }
    if (is.finite(max_fit_windows) &&
        max_fit_windows > .Machine$integer.max) {
        max_fit_windows <- .Machine$integer.max
    }

    reference <- input$reference_sample
    if (length(samples) == 1L) {
        normalization <- "none"
        reference <- NULL
    } else if (length(reference) == 1L && !is.na(reference) &&
        nzchar(reference)) {
        if (!reference %in% samples) {
            .mamamia_error(
                "MaMaMIAAnalysisError",
                "Reference sample '%s' is not one of: %s",
                reference, paste(samples, collapse = ", ")
            )
        }
        normalization <- "reference"
    } else {
        normalization <- "median"
        reference <- NULL
    }

    list(
        resolution = resolution,
        mapq = mapq,
        exclude_flags = exclude_flags,
        cores = cores,
        max_fit_windows = max_fit_windows,
        normalization = normalization,
        reference_sample = reference
    )
}

#' @noRd
bam_progress_stages <- function(n_samples, normalization = "none") {
    n_samples <- max(1L, as.integer(n_samples))
    c(
        "preflight", "index", "windows", "gc",
        rep("counts", n_samples),
        rep("correction", n_samples),
        if (!identical(normalization, "none")) "normalization",
        rep("segmentation", n_samples)
    )
}

.bam_progress_labels <- c(
    preflight = "Checking BAM headers",
    index = "Ensuring CSI indexes",
    windows = "Generating analysis windows",
    gc = "Reading GC content",
    counts = "Counting reads",
    correction = "Correcting GC bias",
    normalization = "Normalizing samples",
    segmentation = "Segmenting chromosome pairs"
)

#' @noRd
bam_progress_label <- function(stage, detail = list()) {
    if (!is.list(detail)) {
        detail <- list()
    }
    label <- .bam_progress_labels[[stage]]
    if (is.null(label) || is.na(label)) {
        label <- stage
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
            assembly <- detail$assembly
            sprintf(
                "%s (%s file%s%s)",
                label, n,
                if (isTRUE(detail$n_bams == 1)) "" else "s",
                if (!is.null(assembly) && length(assembly) == 1L &&
                    !is.na(assembly) && nzchar(assembly)) {
                    sprintf(", assembly '%s'", assembly)
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
        counts = {
            total <- count(detail$total)
            if (is.null(total)) {
                return(paste0(label, sample_suffix))
            }
            sprintf("%s%s (%s reads)", label, sample_suffix, total)
        },
        correction = {
            ideal <- count(detail$n_ideal)
            if (is.null(ideal)) {
                return(paste0(label, sample_suffix))
            }
            sprintf("%s%s (%s ideal windows)", label, sample_suffix, ideal)
        },
        normalization = {
            n <- count(detail$n_samples)
            mode <- detail$mode
            if (is.null(n) || is.null(mode) || length(mode) != 1L ||
                is.na(mode) || !nzchar(mode)) {
                return(label)
            }
            reference <- detail$reference
            detail_text <- if (identical(mode, "reference") &&
                !is.null(reference) && length(reference) == 1L &&
                !is.na(reference) && nzchar(reference)) {
                sprintf("reference %s", reference)
            } else {
                mode
            }
            sprintf("%s (%s sample%s, %s)", label, n,
                if (isTRUE(detail$n_samples == 1)) "" else "s",
                detail_text
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
bam_progress_callback <- function(n_samples, normalization = "none") {
    total <- length(bam_progress_stages(n_samples, normalization))
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
    if (is.null(result) || is.null(result$isa)) {
        return(character())
    }
    names(result$isa)
}

#' One row per sample with reads on targets, calls and normalization
#' @noRd
bam_cohort_summary <- function(result) {
    samples <- bam_sample_choices(result)
    if (length(samples) == 0L) {
        return(data.frame(
            sample = character(),
            reads_donor = numeric(),
            reads_recipient = numeric(),
            reads_total = numeric(),
            n_segments = integer(),
            n_putative_introgressions = integer(),
            normalization = character(),
            reference = character(),
            stringsAsFactors = FALSE
        ))
    }

    side_reads <- function(counts) {
        if (is.null(counts) || is.null(dim(counts))) {
            return(rep(NA_real_, length(samples)))
        }
        present <- samples %in% colnames(counts)
        out <- rep(NA_real_, length(samples))
        out[present] <- as.numeric(colSums(counts[, samples[present], drop = FALSE]))
        out
    }
    reads_donor <- side_reads(result$counts$donor)
    reads_recipient <- side_reads(result$counts$recipient)

    n_segments <- vapply(
        result$isa[samples],
        function(x) as.integer(nrow(x$segments)),
        integer(1L)
    )
    n_putative <- vapply(
        result$isa[samples],
        function(x) {
            as.integer(sum(x$segments$putative_introgression, na.rm = TRUE))
        },
        integer(1L)
    )

    meta <- result$samples
    mode <- result$normalization$mode
    if (is.null(mode) || length(mode) != 1L || is.na(mode)) {
        mode <- "none"
    }
    normalization <- rep(as.character(mode), length(samples))
    reference <- rep(FALSE, length(samples))
    if (is.data.frame(meta) && all(c("sample", "normalization") %in% names(meta))) {
        position <- match(samples, meta$sample)
        known <- !is.na(position)
        normalization[known] <- as.character(meta$normalization[position[known]])
        if ("reference" %in% names(meta)) {
            reference[known] <- as.logical(meta$reference[position[known]])
        }
    }

    data.frame(
        sample = samples,
        reads_donor = reads_donor,
        reads_recipient = reads_recipient,
        reads_total = reads_donor + reads_recipient,
        n_segments = n_segments,
        n_putative_introgressions = n_putative,
        normalization = normalization,
        reference = ifelse(reference, "yes", ""),
        stringsAsFactors = FALSE
    )
}

#' Coverage warnings for zero or near-zero samples
#'
#' A side is near zero when it has fewer than `near_zero` reads per analysis
#' window on average. Reference samples with all-zero windows are reported too,
#' because the dose there follows the pseudocount rule.
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
    donor <- result$counts$donor
    recipient <- result$counts$recipient
    reference <- rep(FALSE, length(samples))
    if (is.data.frame(result$samples) &&
        all(c("sample", "reference") %in% names(result$samples))) {
        position <- match(samples, result$samples$sample)
        known <- !is.na(position)
        reference[known] <- as.logical(result$samples$reference[position[known]])
    }
    reference[is.na(reference)] <- FALSE

    side_values <- function(counts, sample) {
        if (is.null(counts) || is.null(dim(counts)) ||
            !sample %in% colnames(counts)) {
            return(numeric())
        }
        as.numeric(counts[, sample])
    }

    for (i in seq_along(samples)) {
        sample <- samples[[i]]
        for (side in c("donor", "recipient")) {
            values <- if (identical(side, "donor")) {
                side_values(donor, sample)
            } else {
                side_values(recipient, sample)
            }
            n_windows <- length(values)
            total <- sum(values, na.rm = TRUE)
            if (n_windows == 0L) {
                next
            }
            if (total <= 0) {
                warnings <- c(
                    warnings,
                    sprintf(
                        "Sample '%s' has no reads on the %s targets.",
                        sample, side
                    )
                )
            } else if (total / n_windows < near_zero) {
                warnings <- c(
                    warnings,
                    sprintf(
                        "Sample '%s' has near-zero coverage on the %s side (%s reads over %d windows).",
                        sample, side,
                        format(total, scientific = FALSE), n_windows
                    )
                )
            }
        }
        if (isTRUE(reference[[i]])) {
            zero_donor <- sum(side_values(donor, sample) == 0)
            zero_recipient <- sum(side_values(recipient, sample) == 0)
            if (zero_donor > 0L || zero_recipient > 0L) {
                warnings <- c(
                    warnings,
                    sprintf(
                        paste0(
                            "Reference sample '%s' has %d donor and %d recipient",
                            " all-zero window(s); the dose there follows the",
                            " pseudocount rule."
                        ),
                        sample, zero_donor, zero_recipient
                    )
                )
            }
        }
    }
    warnings
}

#' @noRd
bam_selection_map <- function(result) {
    samples <- bam_sample_choices(result)
    stats::setNames(rep(list(integer(0)), length(samples)), samples)
}

#' @noRd
bam_selection_update <- function(map, sample, rows) {
    if (is.null(sample) || length(sample) != 1L || is.na(sample) ||
        !nzchar(sample)) {
        return(map)
    }
    if (is.null(map) || !is.list(map)) {
        map <- list()
    }
    rows <- if (is.null(rows)) integer(0) else suppressWarnings(as.integer(rows))
    rows <- rows[!is.na(rows)]
    map[[as.character(sample)]] <- sort(unique(rows))
    map
}

## A per-sample export uses the rows selected for that sample. When nothing
## is selected yet, `NULL` asks writeCohortToBedpe() for the sample's putative
## calls, which is more useful than an empty header.
#' @noRd
bam_sample_selection <- function(map, sample) {
    if (is.null(map) || !is.list(map) || is.null(sample) ||
        length(sample) != 1L || is.na(sample)) {
        return(NULL)
    }
    rows <- map[[as.character(sample)]]
    rows <- rows[!is.na(rows)]
    if (length(rows) == 0L) {
        return(NULL)
    }
    stats::setNames(
        list(sort(unique(as.integer(rows)))),
        as.character(sample)
    )
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
        align-items: center;
        min-width: 0;
        height: 100%;
        padding: 0 1.5rem;
        border-bottom: 1px solid var(--bs-border-color);
        background: var(--bs-tertiary-bg);
    }

    .mamamia-title {
        margin: 0;
        font-size: clamp(2rem, 4vw, 3.25rem);
        font-weight: 800;
        line-height: 1;
        letter-spacing: -0.02em;
    }

    .mamamia-subtitle {
        margin: .35rem 0 0;
        color: var(--bs-secondary-color);
        font-size: .95rem;
    }

    .mamamia-sidebar {
        min-width: 0;
        min-height: 0;
        overflow-x: hidden;
        overflow-y: auto;
        padding: .55rem .6rem;
        border-right: 1px solid var(--bs-border-color);
        background: var(--bs-body-bg);
    }

    .mamamia-main {
        min-width: 0;
        min-height: 0;
        overflow: auto;
        padding: .75rem 1rem;
    }

    .mamamia-footer {
        display: flex;
        flex-direction: column;
        gap: .25rem;
        min-width: 0;
        min-height: 0;
        overflow: hidden;
        padding: .35rem .75rem .45rem;
        border-top: 1px solid var(--bs-border-color);
        background: var(--bs-tertiary-bg);
    }

    .mamamia-footer .shiny-text-output {
        flex: 1 1 auto;
        margin: 0;
        padding: .15rem .35rem;
        overflow-y: auto;
        font-size: .72rem;
        line-height: 1.3;
        white-space: pre-wrap;
        background: var(--bs-body-bg);
        border: 1px solid var(--bs-border-color);
        border-radius: .25rem;
    }

    .mamamia-footer .footer-status { font-size: .76rem; min-height: 1.25rem; }

    .mamamia-sidebar .shiny-input-container {
        width: 100% !important;
        margin-bottom: .45rem;
    }
    .mamamia-sidebar .form-label,
    .mamamia-sidebar .control-label,
    .mamamia-sidebar label {
        margin-bottom: .1rem;
        font-size: .72rem;
    }
    .mamamia-sidebar .form-control,
    .mamamia-sidebar .form-select {
        height: auto;
        padding: .18rem .35rem;
        font-size: .78rem;
    }
    .mamamia-sidebar .btn {
        padding: .18rem .4rem;
        font-size: .74rem;
    }
    .mamamia-sidebar .accordion-button {
        padding: .35rem .45rem;
        font-size: .78rem;
    }
    .mamamia-sidebar .accordion-body { padding: .45rem; }
    .mamamia-sidebar .radio label,
    .mamamia-sidebar .checkbox label { font-size: .72rem; }
    .mamamia-sidebar .table { margin-bottom: .25rem; font-size: .65rem; }

    .mamamia-main .nav-link { font-size: .9rem; }
    .mamamia-main .tab-content { padding-top: .75rem; }
    "
}

#' @noRd
bam_sidebar_species_ui <- function(config) {
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

    species_controls <- shiny::tagList(
        shiny::selectInput(
            "recipient_species", "Recipient species",
            choices = species_choices,
            selected = species_keys[[1L]]
        ),
        shiny::uiOutput("recipient_subgenome_ui"),
        shiny::selectInput(
            "donor_species", "Donor species",
            choices = species_choices,
            selected = species_keys[[2L]]
        ),
        shiny::uiOutput("donor_subgenome_ui"),
        shiny::actionButton(
            "swap_species", "Swap donor and recipient", class = "btn-sm"
        ),
        shiny::p(
            class = "text-muted small",
            paste(
                "Contigs are paired by chromosome group within the selected",
                "subgenomes; only paired contigs enter the analysis."
            )
        )
    )

    if (length(config$assemblies) > 0L) {
        shiny::tagList(
            shiny::checkboxInput(
                "auto_detect",
                "Detect the assembly from the BAM header",
                value = FALSE
            ),
            shiny::conditionalPanel("!input.auto_detect", species_controls)
        )
    } else {
        species_controls
    }
}

#' @noRd
bam_sidebar_details_ui <- function(config) {
    defaults <- config$defaults
    bslib::accordion(
        id = "bam_sidebar_accordion",
        multiple = FALSE,
        open = FALSE,
        bslib::accordion_panel(
            "Analysis",
            shiny::radioButtons(
                "resolution", "Resolution",
                choices = bam_resolution_choices(defaults),
                selected = bam_default_resolution(defaults)
            ),
            shiny::uiOutput("reference_ui"),
            shiny::uiOutput("run_ui"),
            shiny::p(
                class = "text-muted small",
                paste(
                    "The analysis runs in this session and blocks until it",
                    "finishes; a disconnect is handled after the current",
                    "stage returns."
                )
            )
        ),
        bslib::accordion_panel(
            "Assembly",
            shiny::uiOutput("assembly_summary"),
            shiny::uiOutput("species_ui"),
            shiny::div(
                class = "table-responsive",
                shiny::tableOutput("preflight_table")
            )
        ),
        bslib::accordion_panel(
            "Advanced",
            shiny::numericInput(
                "mapq", "Minimum MAPQ",
                value = defaults$mapq, min = 0, max = 255, step = 1
            ),
            shiny::numericInput(
                "exclude_flags", "Exclude flags (SAM mask)",
                value = defaults$exclude_flags,
                min = 0, max = 4095, step = 1
            ),
            shiny::numericInput(
                "cores", "Cores for the GC-bias fit",
                value = default_cores(), min = 1L, max = bam_max_cores()
            ),
            shiny::numericInput(
                "max_fit_windows", "Max fit windows",
                value = defaults$max_fit_windows,
                min = 1, step = 1000
            )
        ),
        bslib::accordion_panel(
            "Filters",
            explorer_controls_ui("explorer")
        ),
        bslib::accordion_panel(
            "Export",
            shiny::uiOutput("cohort_ui")
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
        "|      |        |         |",
        "|------|--------|---------|",
        "|2rem  |200px   |1fr      |",
        "|150px |header  |header   |",
        "|1fr   |sidebar |main     |",
        "|120px |footer  |footer   |",
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

    config <- state$config
    example_available <- bam_example_available()

    bslib::page_fluid(
        title = "MaMaMIA",
        theme = bslib::bs_theme(version = 5),
        shiny::tags$head(
            shiny::tags$style(shiny::HTML(bam_dashboard_css())),
            shiny::tags$script(shiny::HTML(
                paste(
                    "Shiny.addCustomMessageHandler('mamamia-log-scroll',",
                    "function(x) {",
                    "  var el = document.getElementById('calculation_log');",
                    "  if (el) { el.scrollTop = el.scrollHeight; }",
                    "});"
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
                    shiny::selectInput(
                        "bam_mode", "Input mode",
                        choices = c(
                            "Single BAM" = "single",
                            "Multiple BAMs" = "multi"
                        ),
                        selected = "single",
                        width = "100%"
                    ),
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
                    shiny::uiOutput("main_status"),
                    explorer_results_ui("explorer")
                )
            ),
            gridlayout::grid_place(
                "footer",
                shiny::tags$footer(
                    class = "mamamia-footer",
                    shiny::div(
                        class = "footer-status",
                        shiny::uiOutput("footer_status")
                    ),
                    shiny::verbatimTextOutput("calculation_log")
                )
            )
        )
    )
}

## Server --------------------------------------------------------------------

#' @noRd
bam_launcher_server <- function(state, cache = NULL) {
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

        result <- shiny::reactiveVal(NULL)
        last_error <- shiny::reactiveVal(NULL)
        selection_map <- shiny::reactiveVal(list())

        log_entries <- shiny::reactiveVal(bam_initial_log(config))
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

        bam_mode <- shiny::reactive({
            value <- input$bam_mode
            if (length(value) == 1L && !is.na(value) &&
                value %in% c("single", "multi")) {
                return(as.character(value))
            }
            "single"
        })

        output$bam_mode_note <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                return(shiny::p(
                    class = "text-muted small mb-2",
                    paste(
                        "Example data: one synthetic sample; normalization",
                        "is off."
                    )
                ))
            }
            if (identical(bam_mode(), "multi")) {
                return(shiny::p(
                    class = "text-muted small mb-2",
                    paste(
                        "Cohort run: choose a reference sample or use median",
                        "normalization."
                    )
                ))
            }
            shiny::p(
                class = "text-muted small mb-2",
                paste(
                    "Single-sample run: normalization is off and one BAM",
                    "is analysed."
                )
            )
        })

        output$bam_upload_ui <- shiny::renderUI({
            upload_epoch()
            multi <- identical(bam_mode(), "multi")
            shiny::fileInput(
                "bams",
                if (multi) "BAM files" else "BAM file",
                multiple = multi,
                accept = ".bam",
                width = "100%"
            )
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
            log_line(paste(
                "Example data selected: the normal pipeline stages are",
                "replayed from pre-computed tables."
            ))
        }, ignoreInit = FALSE)

        shiny::observeEvent(input$bam_mode, {
            sidebar_unlocked(FALSE)
            upload_epoch(upload_epoch() + 1L)
            log_line(sprintf(
                "Input mode: %s.",
                if (identical(bam_mode(), "multi")) {
                    "multiple BAMs"
                } else {
                    "single BAM"
                }
            ))
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$unlock_setup, {
            sidebar_unlocked(TRUE)
            log_line("Setup panels opened after a preflight failure.")
        }, ignoreInit = TRUE)

        uploaded_paths <- shiny::reactive({
            uploads <- input$bams
            if (is.data.frame(uploads) && "datapath" %in% names(uploads)) {
                paths <- as.character(uploads$datapath)
                paths <- paths[!is.na(paths) & nzchar(paths)]
                return(paths)
            }
            character()
        })

        shiny::observeEvent(uploaded_paths(), {
            if (length(uploaded_paths()) > 0L && isTRUE(example_mode())) {
                example_mode(FALSE)
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
            paths <- unique(paths)
            if (!identical(bam_mode(), "multi")) {
                paths <- utils::head(paths, 1L)
            }
            paths
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
            swapped <- pending_subgenomes()
            if (!is.null(swapped)) {
                return(swapped$recipient_species)
            }
            value <- input$recipient_species
            if (is.null(value) || length(value) != 1L || is.na(value) ||
                !value %in% species_keys) {
                return(default_recipient)
            }
            as.character(value)
        })
        donor_species <- shiny::reactive({
            swapped <- pending_subgenomes()
            if (!is.null(swapped)) {
                return(swapped$donor_species)
            }
            value <- input$donor_species
            if (is.null(value) || length(value) != 1L || is.na(value) ||
                !value %in% species_keys) {
                return(default_donor)
            }
            as.character(value)
        })

        auto_detect <- shiny::reactive({
            if (!species_available) {
                return(TRUE)
            }
            isTRUE(input$auto_detect)
        })

        ## Subgenome choices to apply after swapping the two species. They are
        ## keyed by species so they never leak into an unrelated selection.
        pending_subgenomes <- shiny::reactiveVal(NULL)
        shiny::observeEvent(input$swap_species, {
            if (!species_available) {
                return()
            }
            recipient_key <- recipient_species()
            donor_key <- donor_species()
            recipient <- config$species[[recipient_key]]
            donor <- config$species[[donor_key]]
            if (is.null(recipient) || is.null(donor)) {
                return()
            }
            pending_subgenomes(list(
                recipient_species = donor_key,
                recipient_subgenome = bam_effective_subgenome(
                    donor, input$donor_subgenome
                ),
                donor_species = recipient_key,
                donor_subgenome = bam_effective_subgenome(
                    recipient, input$recipient_subgenome
                )
            ))
            shiny::updateSelectInput(
                session, "recipient_species", selected = donor_key
            )
            shiny::updateSelectInput(
                session, "donor_species", selected = recipient_key
            )
        })

        ## A manual species change ends the swapped state; the selects updated
        ## by the swap itself match the pending values and keep it.
        shiny::observeEvent(
            list(input$recipient_species, input$donor_species),
            {
                swapped <- pending_subgenomes()
                if (is.null(swapped)) {
                    return()
                }
                if (!identical(input$recipient_species, swapped$recipient_species) ||
                    !identical(input$donor_species, swapped$donor_species)) {
                    pending_subgenomes(NULL)
                }
            },
            ignoreInit = TRUE
        )
        shiny::observeEvent(input$recipient_subgenome, {
            swapped <- pending_subgenomes()
            if (!is.null(swapped) &&
                identical(recipient_species(), swapped$recipient_species) &&
                !identical(input$recipient_subgenome,
                    swapped$recipient_subgenome)) {
                pending_subgenomes(NULL)
            }
        }, ignoreInit = TRUE)
        shiny::observeEvent(input$donor_subgenome, {
            swapped <- pending_subgenomes()
            if (!is.null(swapped) &&
                identical(donor_species(), swapped$donor_species) &&
                !identical(input$donor_subgenome, swapped$donor_subgenome)) {
                pending_subgenomes(NULL)
            }
        }, ignoreInit = TRUE)

        effective_subgenome <- function(role) {
            swapped <- pending_subgenomes()
            if (!is.null(swapped) &&
                identical(
                    swapped[[paste0(role, "_species")]],
                    if (identical(role, "recipient")) {
                        recipient_species()
                    } else {
                        donor_species()
                    }
                )) {
                return(swapped[[paste0(role, "_subgenome")]])
            }
            input[[paste0(role, "_subgenome")]]
        }

        species_subgenome_ui <- function(species_key, input_id, label, role) {
            species <- config$species[[species_key]]
            if (is.null(species)) {
                return(NULL)
            }
            available <- unique(as.character(species$contigs$subgenome))
            if (length(available) <= 1L) {
                return(NULL)
            }
            pending <- pending_subgenomes()
            current <- NULL
            if (!is.null(pending) &&
                identical(pending[[paste0(role, "_species")]], species_key)) {
                current <- pending[[paste0(role, "_subgenome")]]
            }
            if (is.null(current) || length(current) != 1L || is.na(current) ||
                !current %in% available) {
                current <- shiny::isolate(input[[input_id]])
            }
            if (is.null(current) || length(current) != 1L || is.na(current) ||
                !current %in% available) {
                current <- bam_effective_subgenome(species, NULL)
            }
            shiny::selectInput(
                input_id, label, choices = available, selected = current
            )
        }
        output$recipient_subgenome_ui <- shiny::renderUI({
            species_subgenome_ui(
                recipient_species(), "recipient_subgenome",
                "Recipient subgenome", "recipient"
            )
        })
        output$donor_subgenome_ui <- shiny::renderUI({
            species_subgenome_ui(
                donor_species(), "donor_subgenome", "Donor subgenome", "donor"
            )
        })

        assembly_selection <- shiny::reactive({
            if (auto_detect()) {
                ## Header-based auto-detection. It is the only mode for a
                ## registry without species, and an explicit option when the
                ## registry also declares predefined assemblies.
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
                resolveSpeciesAssembly(
                    config,
                    donor = donor_key,
                    recipient = recipient_key,
                    donor_subgenome = bam_effective_subgenome(
                        donor, effective_subgenome("donor")
                    ),
                    recipient_subgenome = bam_effective_subgenome(
                        recipient, effective_subgenome("recipient")
                    )
                ),
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

        output$main_status <- shiny::renderUI({
            res <- result()
            if (is.null(res)) {
                return(shiny::p(
                    class = "text-muted mb-2",
                    paste(
                        "The circos, coverage and segment tabs populate after",
                        "a successful analysis."
                    )
                ))
            }
            samples <- bam_sample_choices(res)
            shiny::p(
                class = "text-muted mb-2",
                sprintf(
                    "%d sample%s loaded: %s",
                    length(samples),
                    if (length(samples) == 1L) "" else "s",
                    paste(samples, collapse = ", ")
                )
            )
        })

        output$footer_status <- shiny::renderUI({
            pf <- preflight_state()
            if (isTRUE(pf$ok)) {
                label <- "Preflight OK"
                badge <- "text-bg-success"
                detail <- sprintf(
                    "%d BAM%s; assembly '%s'; sample%s: %s",
                    nrow(pf$table),
                    if (nrow(pf$table) == 1L) "" else "s",
                    pf$assembly,
                    if (nrow(pf$table) == 1L) "" else "s",
                    paste(unique(pf$table$sample), collapse = ", ")
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
                return(shiny::p(
                    class = "text-muted small", bam_example_note()
                ))
            }
            bam_sidebar_species_ui(config)
        })

        output$preflight_table <- shiny::renderTable({
            pf <- preflight_state()
            if (isTRUE(pf$empty)) {
                return(data.frame(
                    status = paste(
                        "Upload BAM files or select the bundled example",
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
                        "Upload BAM files or select the bundled example",
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

        output$reference_ui <- shiny::renderUI({
            if (!identical(bam_mode(), "multi")) {
                return(NULL)
            }
            pf <- preflight_state()
            if (isTRUE(pf$example) || !isTRUE(pf$ok) || nrow(pf$table) < 2L) {
                return(NULL)
            }
            samples <- unique(pf$table$sample)
            choices <- c(
                "None (median normalization)" = "",
                stats::setNames(samples, samples)
            )
            current <- shiny::isolate(input$reference_sample)
            if (is.null(current) || !current %in% unname(choices)) {
                current <- ""
            }
            shiny::selectInput(
                "reference_sample", "Reference sample",
                choices = choices, selected = current
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

        ## Changing any input that affects the result invalidates the previous
        ## one, so a stale BEDPE can never be downloaded for a different input.
        run_key <- shiny::reactive({
            list(
                mode = bam_mode(),
                example = example_mode(),
                paths = selected_paths(),
                resolution = input$resolution,
                auto_detect = auto_detect(),
                recipient_species = recipient_species(),
                donor_species = donor_species(),
                recipient_subgenome = effective_subgenome("recipient"),
                donor_subgenome = effective_subgenome("donor"),
                reference = input$reference_sample,
                mapq = input$mapq,
                exclude_flags = input$exclude_flags,
                max_fit_windows = input$max_fit_windows
            )
        })
        shiny::observeEvent(run_key(), {
            result(NULL)
            last_error(NULL)
        }, ignoreInit = TRUE)

        shiny::observeEvent(result(), {
            selection_map(bam_selection_map(result()))
            samples <- bam_sample_choices(result())
            if (length(samples) > 0L) {
                shiny::updateSelectInput(
                    session, "result_sample",
                    choices = samples, selected = samples[[1L]]
                )
            }
        }, ignoreNULL = FALSE)

        result_sample <- shiny::reactive({
            samples <- bam_sample_choices(result())
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

        result_isa <- shiny::reactive({
            res <- result()
            sample <- result_sample()
            if (is.null(res) || is.null(sample)) {
                return(NULL)
            }
            res$isa[[sample]]
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
            log_line(sprintf(
                "Starting analysis: %d sample%s, %s windows/chromosome, normalization '%s'%s.",
                length(samples),
                if (length(samples) == 1L) "" else "s",
                format(params$resolution, big.mark = ",", scientific = FALSE),
                params$normalization,
                if (!is.null(params$reference_sample)) {
                    sprintf(", reference '%s'", params$reference_sample)
                } else {
                    ""
                }
            ))
            raw_progress <- bam_progress_callback(
                length(samples), params$normalization
            )
            progress <- function(stage, detail = list()) {
                log_line(bam_progress_label(stage, detail))
                raw_progress(stage, detail)
            }
            new_result <- tryCatch(
                shiny::withProgress(
                    message = "MaMaMIA analysis", value = 0,
                    {
                        if (isTRUE(pf$example)) {
                            example_result <- bam_example_result(
                                params$resolution
                            )
                            bam_example_progress(
                                example_result, progress, pause = 0.1
                            )
                            example_result
                        } else {
                            analyzeBamInput(
                                pf$staged,
                                pf$config,
                                assembly = pf$assembly,
                                windows_per_chromosome = params$resolution,
                                reference_sample = params$reference_sample,
                                normalization = params$normalization,
                                mapq = params$mapq,
                                exclude_flags = params$exclude_flags,
                                cores = params$cores,
                                max_fit_windows = params$max_fit_windows,
                                cache = cache,
                                progress = progress
                            )
                        }
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
            result(new_result)
            log_line(sprintf(
                "Analysis finished: %d sample(s), %d segment(s), %d putative introgression(s).",
                length(new_result$isa),
                sum(new_result$samples$n_segments),
                sum(new_result$samples$n_putative_introgressions)
            ))
            shiny::showNotification(
                sprintf(
                    "%d sample(s), %d segment(s), %d putative introgression(s)",
                    length(new_result$isa),
                    sum(new_result$samples$n_segments),
                    sum(new_result$samples$n_putative_introgressions)
                ),
                type = "message", duration = 6
            )
            coverage <- bam_coverage_warnings(new_result)
            for (warning in coverage) {
                log_line(paste("Coverage warning:", warning))
            }
            if (length(coverage) > 0L) {
                shiny::showNotification(
                    paste(coverage, collapse = "\n"),
                    type = "warning", duration = 15
                )
            }
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
            multi <- identical(bam_mode(), "multi")
            shiny::tagList(
                if (multi && length(samples) > 1L) {
                    shiny::selectInput(
                        "result_sample", "Sample",
                        choices = samples,
                        selected = shiny::isolate(result_sample())
                    )
                },
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
                    if (multi && length(samples) > 1L) {
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
                        "current sample, or all putative calls when none are",
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
                writeCohortToBedpe(
                    res, file,
                    samples = sample,
                    selected = bam_sample_selection(selection_map(), sample)
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
                writeCohortToBedpe(
                    res, file,
                    selected = selection_map()
                )
            }
        )

        explorer_server(
            "explorer",
            ISA = result_isa,
            sample = result_sample,
            selection = shiny::reactive({
                rows <- selection_map()[[result_sample()]]
                if (is.null(rows)) integer() else as.integer(rows)
            }),
            on_selection = function(rows, sample) {
                selection_map(
                    bam_selection_update(selection_map(), sample, rows)
                )
            }
        )
    }
}

#' @noRd
bam_launcher_app <- function(config = NULL, cache = NULL) {
    require_dashboard_packages()
    state <- bam_launcher_state(config)
    ui <- bam_launcher_ui(state)
    server <- bam_launcher_server(state, cache)
    shiny::shinyApp(ui, server)
}
