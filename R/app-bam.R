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
    labels <- vapply(keys, function(key) {
        label <- species[[key]]$label
        if (is.null(label) || length(label) != 1L || is.na(label) ||
            !nzchar(label)) {
            key
        } else {
            label
        }
    }, character(1))
    if (anyDuplicated(labels)) {
        repeated <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
        labels[repeated] <- sprintf("%s (%s)", labels[repeated], keys[repeated])
    }
    stats::setNames(keys, labels)
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
            `5000` = "High (5,000 per chromosome)",
            `10000` = "Extra high (10,000 per chromosome)",
            sprintf(
                "%s per chromosome",
                format(n, big.mark = ",", scientific = FALSE)
            )
        )
    }, character(1))
    stats::setNames(as.character(allowed), labels)
}

#' @noRd
bam_resolution_range <- function(defaults = .config_defaults(),
                                 step = 10L) {
    allowed <- sort(unique(as.integer(defaults$allowed_windows)))
    allowed <- allowed[!is.na(allowed) & allowed > 0L]
    if (length(allowed) == 0L) {
        allowed <- 1000L
    }
    list(
        min = allowed[[1L]],
        max = allowed[[length(allowed)]],
        ticks = allowed,
        step = as.integer(step)
    )
}

#' @noRd
bam_resolution_slider_ui <- function(defaults = .config_defaults(),
                                     locked = FALSE) {
    range <- bam_resolution_range(defaults)
    positions <- if (range$max > range$min) {
        100 * (range$ticks - range$min) / (range$max - range$min)
    } else {
        rep(50, length(range$ticks))
    }
    labels <- lapply(seq_along(range$ticks), function(i) {
        shiny::tags$span(
            class = if (i == 1L) {
                "first"
            } else if (i == length(range$ticks)) {
                "last"
            } else {
                NULL
            },
            style = sprintf("left:%.3f%%;", positions[[i]]),
            format(range$ticks[[i]], big.mark = ",", scientific = FALSE)
        )
    })
    exact <- shiny::numericInput(
        "resolution_number", NULL,
        value = bam_default_resolution(defaults),
        min = range$min,
        max = range$max,
        step = range$step,
        width = "6rem"
    )
    exact_query <- htmltools::tagQuery(exact)$find("input")
    if (isTRUE(locked)) {
        exact_query <- exact_query$addAttrs(
            disabled = "disabled", `aria-disabled` = "true"
        )
    }
    exact <- exact_query$
        addAttrs(
            `aria-label` = "Exact resolution in windows per chromosome",
            title = "Exact resolution in windows per chromosome"
        )$
        allTags()
    shiny::div(
        class = "mamamia-resolution-slider",
        shiny::div(
            class = "mamamia-resolution-head",
            shiny::tags$label(
                class = "control-label", `for` = "resolution", "Resolution"
            ),
            exact
        ),
        shiny::sliderInput(
            "resolution", NULL,
            min = range$min,
            max = range$max,
            value = bam_default_resolution(defaults),
            step = range$step,
            ticks = FALSE,
            round = TRUE
        ),
        shiny::div(class = "mamamia-resolution-ticks", labels)
    )
}

#' @noRd
bam_example_resolution_slider_ui <- function(defaults = .config_defaults()) {
    resolution <- as.integer(bam_example_default_resolution())
    display <- defaults
    display$allowed_windows <- sort(unique(c(
        as.integer(defaults$allowed_windows), resolution
    )))
    display$windows_per_chromosome <- resolution
    shiny::tagList(
        shiny::div(
            class = "mamamia-resolution-locked",
            bam_resolution_slider_ui(display, locked = TRUE),
            shiny::actionButton(
                "resolution_locked", "",
                class = "mamamia-resolution-lock",
                title = "Resolution is fixed for the bundled example",
                `aria-label` = "Resolution is fixed for the bundled example"
            )
        ),
        shiny::p(
            class = "mamamia-resolution-note",
            sprintf(
                paste(
                    "Fixed at %s windows per chromosome.",
                    "Upload a BAM to change it."
                ),
                format(
                    as.integer(resolution),
                    big.mark = ",", scientific = FALSE
                )
            )
        )
    )
}

#' @noRd
bam_default_resolution <- function(defaults = .config_defaults()) {
    range <- bam_resolution_range(defaults)
    valid <- function(value) {
        !is.na(value) && value >= range$min && value <= range$max &&
            (value - range$min) %% range$step == 0L
    }
    value <- suppressWarnings(as.integer(defaults$windows_per_chromosome))
    if (!valid(value)) {
        value <- 1000L
    }
    if (!valid(value)) {
        value <- range$min
    }
    value
}

#' @noRd
bam_resolution_value <- function(value, defaults = .config_defaults()) {
    range <- bam_resolution_range(defaults)
    if (length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
        .mamamia_error("MaMaMIAWindowError", "Choose an analysis resolution")
    }
    resolution <- suppressWarnings(as.integer(value))
    if (is.na(resolution) || resolution < range$min ||
        resolution > range$max ||
        (resolution - range$min) %% range$step != 0L) {
        .mamamia_error(
            "MaMaMIAWindowError",
            "Resolution must be between %s and %s in steps of %s",
            format(range$min, big.mark = ",", scientific = FALSE),
            format(range$max, big.mark = ",", scientific = FALSE),
            range$step
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
    lines <- character()
    lines <- c(lines, sprintf(
        "Preflight OK: %d BAM file%s.",
        n_bams, if (n_bams == 1L) "" else "s"
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
        bam_log_entry("BAM header parsed.", time)
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
#' @param resolution Optional resolved resolution; overrides `input$resolution`.
#'
#' @return A list with the resolved `analyzeBamInput()` arguments.
#' @noRd
bam_run_parameters <- function(input, config, samples, resolution = NULL) {
    if (!is.character(samples) || length(samples) == 0L) {
        .mamamia_error("MaMaMIAAnalysisError", "At least one sample is required")
    }
    defaults <- config$defaults
    if (is.null(resolution)) {
        resolution <- input$resolution
    }
    resolution <- bam_resolution_value(resolution, defaults)
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
    preflight = "Checking BAM",
    index = "Indexing BAM",
    windows = "Building windows",
    gc = "Reading GC content",
    metrics = "Counting reads",
    score = "Scoring windows",
    mask = "Masking windows",
    segmentation = "Segmenting profiles",
    calibration = "Calibrating calls"
)

#' @noRd
bam_progress_label <- function(stage, detail = list()) {
    ## Progress text stays short: the counts, sample names and other detail
    ## fields remain available to the pipeline but are not shown in the toast.
    ## `[` on a named character vector returns NA for an unknown name, so a
    ## future stage cannot abort the analysis.
    label <- unname(.bam_progress_labels[stage])
    if (length(label) != 1L || is.na(label)) {
        label <- as.character(stage)[[1L]]
    }
    label
}

#' @noRd
bam_progress_callback <- function(n_samples, total_steps = NULL,
                                  handler = NULL) {
    stages <- length(bam_progress_stages(n_samples))
    total <- suppressWarnings(as.integer(total_steps))
    if (length(total) != 1L || is.na(total) || total < 1L) {
        total <- stages
    }
    state <- new.env(parent = emptyenv())
    state$done <- 0L
    function(stage, detail = list()) {
        state$done <- state$done + 1L
        if (is.function(handler)) {
            handler(stage, detail, min(1, state$done / total))
        } else {
            shiny::incProgress(
                1 / total,
                detail = bam_progress_label(stage, detail)
            )
        }
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
    if (isTRUE(preflight$example)) {
        donor_aggregation <- defaults$donor_aggregation
        controls <- character(0)
    }
    resolution_key <- if (isTRUE(preflight$example)) {
        as.character(bam_example_default_resolution())
    } else {
        tryCatch(
            as.character(bam_resolution_value(resolution, defaults)),
            error = function(e) as.character(resolution)
        )
    }
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

## The bundled Latin subset keeps the interface independent of the host's
## installed fonts. Inter is distributed under the SIL Open Font License; the
## licence text ships next to the font.
#' @noRd
mamamia_font_dependency <- function() {
    path <- system.file("www", package = "MaMaMIA")
    if (!nzchar(path) || !dir.exists(path)) {
        return(NULL)
    }
    htmltools::htmlDependency(
        name = "mamamia-font",
        version = "1.0.0",
        src = c(file = path),
        stylesheet = "mamamia-font.css",
        all_files = TRUE
    )
}

#' @noRd
mamamia_package_version <- function() {
    tryCatch(
        as.character(utils::packageVersion("MaMaMIA")),
        error = function(e) ""
    )
}

#' @noRd
mamamia_asset_url <- function(file) {
    dep <- mamamia_font_dependency()
    if (is.null(dep)) {
        return(file)
    }
    paste0(dep$name, "-", dep$version, "/", file)
}

#' @noRd
mamamia_favicon <- function() {
    shiny::tags$link(
        rel = "icon",
        type = "image/svg+xml",
        href = mamamia_asset_url("mamamia-logo.svg")
    )
}

## A warm, paper-like surface with a single deep-teal accent keeps the
## scientific content in front. The chart palette stays blue-red so scores and
## calls remain distinguishable from the interface accent.
#' @noRd
bam_dashboard_theme <- function() {
    bslib::bs_theme(
        version = 5,
        primary = "#0f5c53",
        secondary = "#6d7168",
        success = "#2f7d5d",
        info = "#2c6e8f",
        warning = "#b7791f",
        danger = "#b23a3a",
        bg = "#f6f5f0",
        fg = "#1b1e1a",
        base_font = paste(
            "'Inter'", "system-ui", "-apple-system", "'Segoe UI'", "Roboto",
            "'Helvetica Neue'", "Arial", "sans-serif",
            sep = ", "
        ),
        code_font = paste(
            "ui-monospace", "SFMono-Regular", "Menlo", "Consolas",
            "'Liberation Mono'", "monospace",
            sep = ", "
        ),
        "border-color" = "#e4e1d7",
        "border-radius" = ".5rem",
        "font-size-base" = ".95rem",
        "headings-font-weight" = "700"
    )
}

## The page is a three-column, three-row grid. The 2rem first column is a
## gutter; the header and footer span the sidebar and main columns. The grid
## keeps the log visible while the result tabs scroll independently.
#' @noRd
bam_dashboard_css <- function() {
    "
    :root {
        --mm-paper: #f6f5f0;
        --mm-surface: #ffffff;
        --mm-sidebar: #fbfaf6;
        --mm-ink: #1b1e1a;
        --mm-muted: #6d7168;
        --mm-line: #e4e1d7;
        --mm-line-strong: #d2cec0;
        --mm-accent: #0f5c53;
        --mm-accent-rgb: 15, 92, 83;
        --mm-accent-soft: rgba(15, 92, 83, .07);
        --mm-warm: #b7791f;
        --mm-shadow: 0 1px 2px rgba(27, 30, 26, .04),
            0 10px 30px rgba(27, 30, 26, .05);
    }
    .container-fluid {
        height: 100vh;
        padding: 0 !important;
    }
    html, body { height: 100%; }
    body {
        overflow: hidden;
        background-color: var(--mm-paper);
        background-image: radial-gradient(
            rgba(27, 30, 26, .028) .6px, transparent .6px
        );
        background-size: 18px 18px;
        color: var(--mm-ink);
        -webkit-font-smoothing: antialiased;
        text-rendering: optimizeLegibility;
    }

    .mamamia-header {
        position: relative;
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 1rem;
        min-width: 0;
        height: 100%;
        padding: 0 .9rem;
        border-bottom: 1px solid var(--mm-line);
        background: var(--mm-surface);
    }
    .mamamia-header::after {
        content: '';
        position: absolute;
        right: 0;
        bottom: -1px;
        left: 0;
        height: 1px;
        background: linear-gradient(
            90deg,
            var(--mm-accent),
            rgba(183, 121, 31, .5) 32%,
            transparent 72%
        );
    }
    .mamamia-brand {
        display: flex;
        align-items: center;
        gap: .65rem;
        min-width: 0;
    }
    .mamamia-mark {
        display: grid;
        flex: 0 0 auto;
        place-items: center;
        width: 2.05rem;
        height: 2.05rem;
        border-radius: .62rem;
        background: var(--mm-ink);
        box-shadow: 0 2px 6px rgba(27, 30, 26, .18);
        color: var(--mm-paper);
    }
    .mamamia-mark svg {
        display: block;
        width: 1.3rem;
        height: 1.3rem;
    }
    .mamamia-brand-text { min-width: 0; }
    .mamamia-title-row {
        display: flex;
        align-items: baseline;
        gap: .45rem;
    }
    .mamamia-title {
        margin: 0;
        font-size: 1.55rem;
        font-weight: 800;
        line-height: 1;
        letter-spacing: -.035em;
    }
    .mamamia-version {
        padding: .05rem .35rem;
        border: 1px solid var(--mm-line);
        border-radius: 999px;
        color: var(--mm-muted);
        font-size: .58rem;
        font-weight: 650;
        letter-spacing: .04em;
        white-space: nowrap;
    }
    .mamamia-subtitle {
        margin: .2rem 0 0;
        color: var(--mm-muted);
        font-size: .7rem;
        letter-spacing: .012em;
    }
    .mamamia-header-note {
        display: flex;
        align-items: center;
        gap: .45rem;
        color: var(--mm-muted);
        font-size: .64rem;
        font-weight: 600;
        letter-spacing: .1em;
        text-transform: uppercase;
        white-space: nowrap;
    }
    .mamamia-header-note::before {
        content: '';
        width: .38rem;
        height: .38rem;
        border-radius: 50%;
        background: var(--mm-accent);
        box-shadow: 0 0 0 3px var(--mm-accent-soft);
    }
    @media (max-width: 980px) {
        .mamamia-header-note { display: none; }
    }

    .mamamia-sidebar {
        min-width: 0;
        min-height: 0;
        overflow-x: hidden;
        overflow-y: auto;
        padding: .75rem .85rem 1rem;
        border-right: 1px solid var(--mm-line);
        background: var(--mm-sidebar);
        font-size: .85rem;
        line-height: 1.35;
        scrollbar-width: thin;
        scrollbar-color: var(--mm-line-strong) transparent;
    }
    .mamamia-sidebar::-webkit-scrollbar { width: 7px; }
    .mamamia-sidebar::-webkit-scrollbar-thumb {
        border: 2px solid var(--mm-sidebar);
        border-radius: 999px;
        background: var(--mm-line-strong);
    }
    .mamamia-sidebar-section { margin-top: .15rem; }
    .mamamia-sidebar .mamamia-section-label {
        margin: .1rem 0 .45rem;
        color: var(--mm-muted);
        font-size: .6rem;
        font-weight: 750;
        letter-spacing: .13em;
        text-transform: uppercase;
    }
    .mamamia-sidebar .mamamia-sidebar-intro {
        margin-bottom: .55rem;
        color: var(--mm-muted);
        font-size: .7rem;
        line-height: 1.4;
    }
    .mamamia-sidebar .mamamia-sidebar-intro p {
        margin: 0;
    }
    .mamamia-sidebar .mamamia-recovery {
        margin: .1rem 0 .8rem;
        padding: .72rem .78rem .78rem;
        border: 1px solid rgba(183, 121, 31, .26);
        border-left: 3px solid var(--mm-warm);
        border-radius: .55rem;
        background: linear-gradient(
            180deg,
            rgba(183, 121, 31, .095),
            rgba(183, 121, 31, .035)
        );
        box-shadow: 0 1px 2px rgba(27, 30, 26, .03);
    }
    .mamamia-recovery-head {
        display: flex;
        align-items: center;
        gap: .42rem;
        margin-bottom: .35rem;
    }
    .mamamia-recovery-icon {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        width: 1.25rem;
        height: 1.25rem;
        border-radius: 50%;
        background: rgba(183, 121, 31, .16);
        color: var(--mm-warm);
        font-size: .62rem;
    }
    .mamamia-recovery-title {
        color: var(--mm-ink);
        font-size: .74rem;
        font-weight: 700;
        line-height: 1.25;
    }
    .mamamia-recovery-message {
        margin: 0 0 .4rem;
        color: #675027;
        font-size: .67rem;
        line-height: 1.42;
        overflow-wrap: anywhere;
    }
    .mamamia-recovery-note {
        margin: 0 0 .45rem;
        color: var(--mm-muted);
        font-size: .64rem;
        line-height: 1.4;
    }
    .mamamia-recovery-suggestion {
        margin: 0 0 .55rem;
        color: var(--mm-muted);
        font-size: .65rem;
        line-height: 1.4;
    }
    .mamamia-recovery-suggestion strong {
        color: #675027;
        font-weight: 700;
    }
    .mamamia-recovery-action {
        width: 100%;
        padding: .32rem .5rem;
        border: 1px solid rgba(183, 121, 31, .38);
        border-radius: .4rem;
        background: #fffaf0;
        color: #6b4f18;
        font-size: .68rem;
        font-weight: 650;
        line-height: 1.25;
        transition: background-color .15s ease, border-color .15s ease;
    }
    .mamamia-recovery-action:hover {
        border-color: rgba(183, 121, 31, .6);
        background: #fdf2dc;
        color: #563d10;
    }

    .mamamia-main {
        display: flex;
        flex-direction: column;
        min-width: 0;
        min-height: 0;
        overflow: hidden;
        padding: .6rem .7rem .65rem;
    }

    .mamamia-main > .tabbable,
    .mamamia-main > .card {
        flex: 1 1 auto;
        min-height: 0;
        margin-bottom: 0;
    }
    .mamamia-main > .card {
        overflow: hidden;
        border: 1px solid var(--mm-line);
        border-radius: .7rem;
        background: var(--mm-surface);
        box-shadow: var(--mm-shadow);
    }
    .mamamia-main .card-header {
        flex: 0 0 auto;
        padding: .3rem .6rem 0;
        border-bottom: 1px solid var(--mm-line);
        background: transparent;
    }
    .mamamia-main .card-body {
        min-height: 0;
        padding: 0;
        overflow: auto;
    }
    .mamamia-main .tab-content {
        min-height: 0;
        overflow: hidden;
    }
    /* bslib writes `padding:0` inline on fillable tab panes. Keep the tab
       content aligned with the tab labels instead of the card border. */
    .mamamia-main .tab-pane {
        height: 100%;
        padding: .65rem .7rem !important;
        overflow: auto;
    }
    .mamamia-main .nav-tabs {
        gap: .1rem;
        border-bottom: 0;
    }
    .mamamia-main .nav-link {
        padding: .48rem .72rem;
        border: 0;
        border-bottom: 2px solid transparent;
        border-radius: 0;
        color: var(--mm-muted);
        font-size: .8rem;
        font-weight: 600;
        line-height: 1.2;
    }
    .mamamia-main .nav-link:hover {
        border-bottom-color: var(--mm-line-strong);
        color: var(--mm-ink);
    }
    .mamamia-main .nav-link.active {
        border-bottom-color: var(--mm-accent);
        background: transparent;
        color: var(--mm-ink);
    }
    .mamamia-main .nav-link:focus-visible {
        outline: 2px solid var(--mm-accent);
        outline-offset: -2px;
    }

    .mamamia-footer {
        display: flex;
        flex-direction: column;
        gap: .2rem;
        min-width: 0;
        min-height: 0;
        overflow: hidden;
        padding: .25rem .85rem .35rem;
        border-top: 1px solid var(--mm-line);
        background: var(--mm-sidebar);
    }

    .mamamia-footer .shiny-text-output {
        flex: 1 1 auto;
        margin: 0;
        padding: .2rem .45rem;
        overflow-y: auto;
        border: 1px solid var(--mm-line);
        border-radius: .35rem;
        background: var(--mm-paper);
        color: #474a42;
        font-family: ui-monospace, SFMono-Regular, Menlo, Consolas,
            'Liberation Mono', monospace;
        font-size: .66rem;
        line-height: 1.35;
        white-space: pre-wrap;
    }

    .mamamia-footer .footer-status {
        display: flex;
        flex: 0 0 auto;
        align-items: center;
        gap: .5rem;
        min-height: 1.35rem;
        font-size: .72rem;
        line-height: 1.2;
    }

    .mamamia-footer .mamamia-footer-status {
        flex: 1 1 auto;
        min-width: 0;
    }
    .mamamia-footer .badge {
        padding: .32em .6em;
        border-radius: 999px;
        font-size: .58rem;
        font-weight: 750;
        letter-spacing: .07em;
        text-transform: uppercase;
        vertical-align: .08em;
    }

    .mamamia-footer .mamamia-footer-toggle {
        flex: 0 0 auto;
        padding: .15rem .45rem;
        border: 0;
        border-radius: .35rem;
        background: transparent;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.2;
    }
    .mamamia-footer .mamamia-footer-toggle:hover,
    .mamamia-footer .mamamia-footer-toggle:focus {
        background: var(--mm-paper);
        color: var(--mm-ink);
    }

    /* Collapsing the log gives its height back to the sidebar and main area. */
    #mamamia-grid.mamamia-footer-collapsed {
        grid-template-rows: 75px 1fr 2.25rem;
    }
    #mamamia-grid.mamamia-footer-collapsed .mamamia-footer .shiny-text-output {
        display: none;
    }
    @media (max-width: 500px) {
        #mamamia-grid.mamamia-footer-collapsed {
            grid-template-rows: 85px 350px 2.25rem 350px;
        }
    }

    .mamamia-sidebar .shiny-input-container {
        width: 100% !important;
        margin-bottom: .55rem;
    }
    .mamamia-sidebar .form-label,
    .mamamia-sidebar .control-label,
    .mamamia-sidebar .shiny-input-container > label {
        margin-bottom: .18rem;
        color: var(--mm-muted);
        font-size: .67rem;
        font-weight: 650;
        letter-spacing: .03em;
        line-height: 1.25;
    }
    .mamamia-sidebar .form-control,
    .mamamia-sidebar .form-select {
        min-height: 2rem;
        height: auto;
        padding: .28rem .5rem;
        border-color: var(--mm-line);
        border-radius: .42rem;
        background-color: var(--mm-surface);
        font-size: .79rem;
        line-height: 1.3;
    }
    .mamamia-sidebar .form-control:focus,
    .mamamia-sidebar .form-select:focus {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-sidebar .btn {
        padding: .25rem .5rem;
        font-size: .78rem;
        line-height: 1.25;
    }
    .mamamia-sidebar .mamamia-ghost-btn {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: .35rem;
        width: 100%;
        height: 1.95rem;
        margin-top: .45rem;
        border: 1px solid var(--mm-line);
        border-radius: .42rem;
        background: transparent;
        color: var(--mm-ink);
        font-size: .73rem;
        font-weight: 600;
    }
    .mamamia-sidebar .mamamia-ghost-btn:hover,
    .mamamia-sidebar .mamamia-ghost-btn:focus-visible {
        border-color: var(--mm-line-strong);
        background: var(--mm-paper);
        color: var(--mm-ink);
    }

    .mamamia-upload .form-group { margin-bottom: 0; }
    .mamamia-upload .control-label {
        position: absolute !important;
        width: 1px !important;
        height: 1px !important;
        padding: 0 !important;
        margin: -1px !important;
        overflow: hidden !important;
        clip: rect(0, 0, 0, 0) !important;
        white-space: nowrap !important;
        border: 0 !important;
    }
    .mamamia-upload-head {
        display: flex;
        align-items: center;
        gap: .5rem;
        margin-bottom: .45rem;
        cursor: pointer;
    }
    .mamamia-upload-icon {
        display: grid;
        flex: 0 0 auto;
        place-items: center;
        width: 1.65rem;
        height: 1.65rem;
        border-radius: .45rem;
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
        font-size: .74rem;
    }
    .mamamia-upload-title {
        display: block;
        font-size: .76rem;
        font-weight: 650;
        line-height: 1.2;
    }
    .mamamia-upload-hint {
        display: block;
        color: var(--mm-muted);
        font-size: .64rem;
        line-height: 1.3;
    }
    .mamamia-upload .input-group {
        gap: .4rem;
        padding: .5rem;
        border: 1px dashed var(--mm-line-strong);
        border-radius: .55rem;
        background: var(--mm-surface);
        cursor: pointer;
        transition: border-color .15s ease, background-color .15s ease;
    }
    .mamamia-upload .input-group:hover,
    .mamamia-upload.is-dragover .input-group {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
    }
    .mamamia-upload .input-group-btn {
        display: block;
        width: 100%;
        margin: 0;
    }
    .mamamia-upload .btn-file {
        display: flex;
        align-items: center;
        justify-content: center;
        width: 100%;
        height: 1.85rem;
        padding: 0 .55rem;
        border: 0;
        border-radius: .38rem;
        background: var(--mm-ink);
        color: var(--mm-paper);
        font-size: .7rem;
        font-weight: 650;
        letter-spacing: .01em;
        box-shadow: none;
    }
    .mamamia-upload .btn-file:hover,
    .mamamia-upload .btn-file:focus {
        background: var(--mm-accent);
        color: #fff;
    }
    .mamamia-upload input.form-control[readonly] {
        flex: 1 1 100%;
        min-height: 0;
        height: auto;
        padding: 0 .1rem;
        border: 0;
        border-radius: 0;
        background: transparent;
        color: var(--mm-muted);
        font-size: .64rem;
        line-height: 1.3;
        text-align: center;
        box-shadow: none;
    }
    .mamamia-upload .progress {
        height: 3px;
        margin: .35rem 0 0;
        border-radius: 999px;
        background: var(--mm-line);
    }
    /* Shiny writes the file name (and then 'Upload complete') into the
       progress bar. Keep the bar as a pure fill indicator. */
    .mamamia-upload .progress-bar {
        overflow: hidden;
        background-color: var(--mm-accent);
        color: transparent;
        font-size: 0;
        line-height: 0;
        text-indent: -9999px;
        white-space: nowrap;
    }

    .mamamia-sidebar .accordion {
        --bs-accordion-bg: transparent;
        --bs-accordion-btn-bg: transparent;
        --bs-accordion-active-bg: transparent;
        --bs-accordion-border-color: var(--mm-line);
        --bs-accordion-btn-focus-box-shadow: none;
        margin-top: .6rem;
    }
    .mamamia-sidebar .accordion-item {
        border: 0;
        border-top: 1px solid var(--mm-line);
        background: transparent;
    }
    .mamamia-sidebar .accordion-button {
        padding: .62rem .1rem;
        background: transparent;
        color: var(--mm-muted);
        font-size: .61rem;
        font-weight: 750;
        letter-spacing: .12em;
        line-height: 1.2;
        text-transform: uppercase;
        box-shadow: none;
    }
    .mamamia-sidebar .accordion-button:not(.collapsed) {
        color: var(--mm-ink);
    }
    .mamamia-sidebar .accordion-button::after {
        width: .62rem;
        height: .62rem;
        background-size: .62rem;
        opacity: .55;
    }
    .mamamia-sidebar .accordion-body {
        padding: .05rem .1rem .8rem;
    }
    .mamamia-sidebar .shiny-options-group {
        display: flex;
        flex-direction: column;
        gap: .2rem;
    }
    .mamamia-sidebar .radio-inline {
        display: flex;
        align-items: center;
        gap: .4rem;
        margin: 0;
        font-size: .76rem;
    }
    .mamamia-sidebar .radio-inline input,
    .mamamia-sidebar .form-check-input {
        accent-color: var(--mm-accent);
    }
    .mamamia-sidebar .form-check-label,
    .mamamia-sidebar .radio label,
    .mamamia-sidebar .checkbox label {
        font-size: .76rem;
        font-weight: 450;
        line-height: 1.35;
    }
    .mamamia-sidebar .text-muted.small,
    .mamamia-sidebar small {
        color: var(--mm-muted) !important;
        font-size: .69rem;
        line-height: 1.4;
    }
    .mamamia-sidebar .selectize-input {
        min-height: 2rem;
        padding: .25rem .4rem;
        border-color: var(--mm-line);
        border-radius: .42rem;
        background: var(--mm-surface);
        font-size: .74rem;
    }
    .mamamia-sidebar .selectize-control.multi .selectize-input > div {
        border: 0;
        border-radius: .3rem;
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
        font-size: .68rem;
    }
    .mamamia-run-panel {
        margin-top: .35rem;
        padding-top: .7rem;
        border-top: 1px solid var(--mm-line);
    }
    .mamamia-run-panel .btn {
        width: 100%;
    }
    .mamamia-run-panel #run {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: .45rem;
        height: 2.3rem;
        border: 0;
        border-radius: .45rem;
        background: var(--mm-accent);
        color: #ffffff;
        font-size: .76rem;
        font-weight: 650;
        letter-spacing: .015em;
        box-shadow: 0 2px 5px rgba(15, 92, 83, .22);
        transition: background-color .15s ease, transform .1s ease;
    }
    .mamamia-run-panel #run:hover:not(.disabled) {
        background: #0b4a43;
    }
    .mamamia-run-panel #run:active:not(.disabled) {
        transform: translateY(1px);
    }
    .mamamia-run-panel #run.disabled {
        box-shadow: none;
        opacity: .42;
    }

    .mamamia-main > .mamamia-filter-bar {
        flex: 0 0 auto;
    }
    .mamamia-filter-bar {
        display: flex;
        align-items: flex-end;
        gap: 0;
        padding: .1rem 0 .55rem;
        border: 0;
        border-bottom: 1px solid var(--mm-line);
        border-radius: 0;
        background: transparent;
    }
    .mamamia-sample-picker {
        flex: 0 0 auto;
    }
    /* Never hide the empty wrapper: Shiny suspends display:none outputs, so
       the picker could not render once a multi-sample BAM arrives. The
       populated input carries its own spacing so the empty state is flush. */
    .mamamia-sample-picker .shiny-input-container,
    .mamamia-sample-picker .form-group {
        width: 10.5rem;
        margin-bottom: 0;
        margin-right: .55rem;
    }
    .mamamia-sample-picker .control-label {
        margin-bottom: .15rem;
        color: var(--mm-muted);
        font-size: .61rem;
        font-weight: 700;
        letter-spacing: .1em;
        text-transform: uppercase;
    }
    .mamamia-sample-picker .form-select {
        min-height: 1.9rem;
        padding: .2rem .6rem;
        border-color: var(--mm-line);
        border-radius: .45rem;
        background-color: var(--mm-surface);
        font-size: .72rem;
    }
    .mamamia-filter-controls {
        display: flex;
        flex: 1 1 auto;
        flex-wrap: wrap;
        align-items: center;
        gap: .4rem;
        min-width: 0;
    }
    .mamamia-filter-pills {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: .45rem;
        width: 100%;
    }
    .mamamia-filter-pills .shiny-input-container {
        width: auto;
        min-width: 0;
        margin-bottom: 0;
        font-size: 0;
    }
    .mamamia-filter-pills .bootstrap-select {
        flex: 0 0 auto;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle {
        min-height: 1.9rem;
        height: 1.9rem;
        padding: 0 .78rem;
        border: 1px solid var(--mm-line);
        border-radius: .45rem;
        background: var(--mm-surface);
        color: var(--mm-muted);
        font-size: .72rem;
        font-weight: 600;
        box-shadow: none;
        transition: border-color .15s ease, color .15s ease,
            background-color .15s ease;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle .filter-option {
        display: flex;
        align-items: center;
        height: 100%;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle::after {
        width: 5px;
        height: 5px;
        margin: 0 0 0 .42rem;
        border: 0;
        border-right: 1.5px solid currentColor;
        border-bottom: 1.5px solid currentColor;
        transform: rotate(45deg) translateY(-1px);
        transition: transform .15s ease;
    }
    .mamamia-filter-pills .bootstrap-select.show > .dropdown-toggle::after {
        transform: rotate(225deg) translateY(-1px);
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle:hover {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle:focus {
        outline: none !important;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle:focus-visible {
        outline: 2px solid var(--mm-accent) !important;
        outline-offset: 2px !important;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-toggle:not(.bs-placeholder) {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-filter-pills .bootstrap-select.show > .dropdown-toggle {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-filter-pills .bootstrap-select .filter-option-inner-inner {
        font-size: .72rem;
    }
    .mamamia-filter-pills .bootstrap-select > .dropdown-menu {
        width: 10rem;
        min-width: 10rem;
        max-width: 10rem;
        padding: .28rem 0;
        border-color: var(--mm-line);
        border-radius: .45rem;
        box-shadow: 0 .45rem 1.2rem rgba(27, 30, 26, .13);
        font-size: .72rem;
    }
    .mamamia-filter-pills .bootstrap-select .dropdown-item {
        padding: .3rem .8rem;
    }
    .mamamia-filter-pills .bootstrap-select .dropdown-item.active,
    .mamamia-filter-pills .bootstrap-select .dropdown-item:active {
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-filter-pills .mamamia-filter-switch {
        display: flex;
        align-items: center;
        min-height: 1.9rem;
        margin-left: auto;
    }
    .mamamia-filter-pills .form-check.form-switch {
        display: flex;
        align-items: center;
        gap: .42rem;
        min-height: 1.9rem;
        margin: 0;
        padding-left: 2.4em;
        font-size: .72rem;
        white-space: nowrap;
    }
    .mamamia-filter-pills .form-check-input {
        margin-top: 0;
        border-color: var(--mm-line-strong);
        background-color: var(--mm-surface);
    }
    .mamamia-filter-pills .form-check-input:checked {
        border-color: var(--mm-accent);
        background-color: var(--mm-accent);
    }
    .mamamia-filter-pills .form-check-input:focus {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-filter-pills .form-check-label {
        color: var(--mm-muted);
        font-size: .72rem;
        font-weight: 600;
    }
    .mamamia-segment-controls {
        display: flex;
        flex-wrap: wrap;
        align-items: flex-end;
        gap: .5rem;
        margin-bottom: .55rem;
        padding-bottom: .55rem;
        border-bottom: 1px solid var(--mm-line);
    }
    .mamamia-segment-controls .shiny-input-container {
        margin-bottom: 0;
    }
    .mamamia-segment-controls .form-group {
        margin-bottom: 0;
    }
    .mamamia-segment-controls .control-label {
        margin-bottom: .18rem;
        color: var(--mm-muted);
        font-size: .64rem;
        font-weight: 650;
        letter-spacing: .03em;
    }
    .mamamia-segment-controls .form-control {
        height: 2rem;
        padding: .2rem .5rem;
        border-color: var(--mm-line);
        border-radius: .4rem;
        font-size: .76rem;
    }
    .mamamia-segment-controls .btn {
        height: 2rem;
        padding: 0 .65rem;
        border: 1px solid var(--mm-line);
        border-radius: .4rem;
        background: var(--mm-surface);
        color: var(--mm-ink);
        font-size: .72rem;
        font-weight: 600;
    }
    .mamamia-segment-controls .btn:hover {
        border-color: var(--mm-accent);
        background: var(--mm-accent-soft);
        color: var(--mm-accent);
    }
    .mamamia-segment-controls .d-flex { gap: .4rem !important; }
    .mamamia-export-panel {
        max-width: 56rem;
    }
    .mamamia-export-panel h5 {
        margin: 1rem 0 .4rem;
        color: var(--mm-muted);
        font-size: .62rem;
        font-weight: 750;
        letter-spacing: .12em;
        text-transform: uppercase;
    }
    .mamamia-export-panel h5:first-child { margin-top: 0; }
    .mamamia-export-panel .btn {
        border-radius: .4rem;
        font-size: .78rem;
    }
    .mamamia-export-panel .btn.disabled,
    .mamamia-export-panel .btn:disabled {
        opacity: .45;
        pointer-events: none;
    }
    .mamamia-export-panel #explorer-selection_status {
        margin-bottom: .35rem;
        color: var(--mm-muted);
        font-size: .74rem;
    }
    .mamamia-export-panel .table-responsive {
        margin-bottom: .7rem;
        scrollbar-width: thin;
        scrollbar-color: var(--mm-line-strong) transparent;
    }
    .mamamia-export-panel .table-responsive::-webkit-scrollbar {
        height: 7px;
    }
    .mamamia-export-panel .table-responsive::-webkit-scrollbar-thumb {
        border-radius: 999px;
        background: var(--mm-line-strong);
    }
    .mamamia-export-panel .table-responsive::-webkit-scrollbar-track {
        background: transparent;
    }
    .mamamia-export-panel .table {
        margin-bottom: 0;
        font-size: .72rem;
        white-space: nowrap;
    }
    .mamamia-export-panel .table th {
        color: var(--mm-muted);
        font-size: .62rem;
        font-weight: 700;
        letter-spacing: .04em;
    }
    .mamamia-export-panel .table td,
    .mamamia-export-panel .table th {
        padding: .35rem .5rem;
    }
    .mamamia-karyotype {
        display: flex;
        flex-direction: column;
        height: 100%;
        min-height: 0;
    }
    .mamamia-karyotype > div:first-child {
        flex: 0 0 auto;
        justify-content: space-between;
        gap: .5rem !important;
        padding: 0 0 .5rem !important;
    }
    .mamamia-karyotype-plot {
        position: relative;
        display: flex;
        flex: 1 1 auto;
        flex-direction: column;
        min-height: 0;
    }
    .mamamia-karyotype-plot > .html-widget {
        flex: 1 1 auto;
        min-height: 0;
        height: 100% !important;
    }
    .mamamia-analysis-loading {
        position: absolute;
        inset: 0;
        z-index: 4;
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        gap: .85rem;
        padding: 1.5rem;
        text-align: center;
        background: var(--mm-surface);
        opacity: 0;
        visibility: hidden;
        transition: opacity .25s ease, visibility .25s ease;
        pointer-events: none;
    }
    .mamamia-analysis-loading.is-running {
        opacity: 1;
        visibility: visible;
    }
    .mamamia-analysis-loader {
        position: relative;
        width: 92px;
        height: 92px;
    }
    .mamamia-analysis-loader::before {
        content: '';
        position: absolute;
        inset: 12%;
        border-radius: 50%;
        background: var(--mm-accent-soft);
        animation: mamamia-loader-halo 2.4s ease-in-out infinite;
    }
    .mamamia-analysis-loader svg {
        position: relative;
        display: block;
        width: 100%;
        height: 100%;
        overflow: visible;
    }
    .mamamia-analysis-orbit {
        fill: none;
        stroke: var(--mm-line-strong);
        stroke-width: 1.4;
        stroke-dasharray: 4 7;
        transform-box: view-box;
        transform-origin: 40px 40px;
        animation: mamamia-loader-orbit 10s linear infinite;
    }
    .mamamia-analysis-strand {
        fill: none;
        stroke: var(--mm-accent);
        stroke-width: 2.6;
        stroke-linecap: round;
        animation: mamamia-loader-strand 2.4s ease-in-out infinite;
    }
    .mamamia-analysis-strand-warm {
        stroke: var(--mm-warm);
        animation-delay: .3s;
    }
    .mamamia-analysis-chiasma {
        fill: var(--mm-accent);
        transform-box: view-box;
        transform-origin: 40px 40px;
        animation: mamamia-loader-chiasma 1.6s ease-in-out infinite;
    }
    .mamamia-analysis-eyebrow {
        color: var(--mm-muted);
        font-size: .62rem;
        font-weight: 700;
        letter-spacing: .14em;
        text-transform: uppercase;
    }
    .mamamia-analysis-stage {
        color: var(--mm-ink);
        font-size: 1rem;
        font-weight: 650;
        line-height: 1.3;
    }
    .mamamia-analysis-track {
        width: min(240px, 72%);
        height: 3px;
        overflow: hidden;
        border-radius: 999px;
        background: var(--mm-line);
    }
    .mamamia-analysis-fill {
        display: block;
        width: 0;
        height: 100%;
        border-radius: inherit;
        background: linear-gradient(90deg, var(--mm-accent), var(--mm-warm));
        transition: width .25s ease;
    }
    @keyframes mamamia-loader-orbit {
        to { transform: rotate(360deg); }
    }
    @keyframes mamamia-loader-strand {
        0%, 100% { opacity: .55; }
        50% { opacity: 1; }
    }
    @keyframes mamamia-loader-chiasma {
        0%, 100% { transform: scale(.75); opacity: .55; }
        50% { transform: scale(1.25); opacity: 1; }
    }
    @keyframes mamamia-loader-halo {
        0%, 100% { transform: scale(.85); opacity: .45; }
        50% { transform: scale(1.1); opacity: .9; }
    }
    @media (prefers-reduced-motion: reduce) {
        .mamamia-analysis-orbit,
        .mamamia-analysis-strand,
        .mamamia-analysis-chiasma,
        .mamamia-analysis-loader::before {
            animation: none;
        }
        .mamamia-analysis-fill { transition: none; }
    }
    .mamamia-karyotype .shiny-input-container {
        width: auto;
        margin-bottom: 0;
    }
    .mamamia-karyotype .control-label {
        position: absolute !important;
        width: 1px !important;
        height: 1px !important;
        padding: 0 !important;
        margin: -1px !important;
        overflow: hidden !important;
        clip: rect(0, 0, 0, 0) !important;
        white-space: nowrap !important;
        border: 0 !important;
    }
    .mamamia-karyotype .shiny-options-group {
        display: inline-flex;
        gap: 2px;
        padding: 2px;
        border: 1px solid var(--mm-line);
        border-radius: .45rem;
        background: var(--mm-paper);
        font-size: 0;
    }
    .mamamia-karyotype .radio-inline {
        margin: 0;
        padding: 0;
        font-size: 0;
    }
    .mamamia-karyotype .radio-inline input {
        position: absolute;
        width: 1px;
        height: 1px;
        margin: 0;
        opacity: 0;
        pointer-events: none;
    }
    .mamamia-karyotype .radio-inline span {
        display: block;
        padding: .3rem .8rem;
        border-radius: .32rem;
        color: var(--mm-muted);
        font-size: .7rem;
        font-weight: 600;
        line-height: 1.25;
        text-align: center;
        cursor: pointer;
        transition: background-color .15s ease, color .15s ease;
    }
    .mamamia-karyotype .radio-inline input:checked + span {
        background: var(--mm-surface);
        box-shadow: 0 1px 2px rgba(27, 30, 26, .12);
        color: var(--mm-ink);
    }
    .mamamia-karyotype .radio-inline input:focus-visible + span {
        outline: 2px solid var(--mm-accent);
        outline-offset: 1px;
    }
    .mamamia-karyotype .html-widget {
        flex: 1 1 auto;
        min-height: 0;
        height: 100% !important;
    }
    .mamamia-tooltip-modal.show {
        display: flex !important;
        align-items: center;
        justify-content: center;
        height: 100vh;
        height: 100dvh;
        padding: 1rem;
        overflow: hidden;
    }
    .mamamia-tooltip-modal .modal-dialog {
        width: min(800px, calc(100vw - 2rem));
        max-width: calc(100vw - 2rem);
        max-height: calc(100vh - 2rem);
        margin: 0;
    }
    .mamamia-tooltip-modal .modal-content {
        display: flex;
        flex-direction: column;
        max-height: calc(100vh - 2rem);
    }
    .mamamia-tooltip-modal .modal-header,
    .mamamia-tooltip-modal .modal-footer {
        flex: 0 0 auto;
    }
    .mamamia-tooltip-modal .modal-body {
        min-height: 0;
        overflow-y: auto;
        overscroll-behavior: contain;
    }
    .mamamia-tooltip-help {
        padding-right: .35rem;
    }
    .mamamia-help-lead {
        margin: 0 0 .55rem;
        color: var(--mm-ink);
        font-size: .82rem;
        line-height: 1.5;
    }
    .mamamia-help-callout {
        margin: 0 0 1rem;
        padding: .55rem .7rem;
        border-left: 3px solid var(--mm-accent);
        border-radius: .35rem;
        background: var(--mm-accent-soft);
        color: var(--mm-ink);
        font-size: .78rem;
        line-height: 1.45;
    }
    .mamamia-help-section + .mamamia-help-section {
        margin-top: 1rem;
    }
    .mamamia-help-intro {
        margin: 0 0 .4rem;
        color: var(--mm-muted);
        font-size: .78rem;
        line-height: 1.45;
    }
    .mamamia-help-section h6 {
        margin: 0 0 .3rem;
        color: var(--mm-muted);
        font-size: .66rem;
        font-weight: 700;
        letter-spacing: .07em;
        text-transform: uppercase;
    }
    .mamamia-help-list {
        display: grid;
        grid-template-columns: 10rem minmax(0, 1fr);
        align-items: start;
        gap: .28rem .9rem;
        margin: 0;
        text-align: left;
    }
    .mamamia-help-list dt {
        color: var(--mm-ink);
        font-size: .78rem;
        font-weight: 650;
    }
    .mamamia-help-list dd {
        margin: 0;
        color: var(--mm-muted);
        font-size: .78rem;
        line-height: 1.4;
    }
    .mamamia-help-list dd strong {
        color: var(--mm-ink);
        font-weight: 650;
    }
    .mamamia-sidebar .form-check {
        margin-bottom: .25rem;
        min-height: 1.25rem;
    }
    .mamamia-sidebar .radio-inline + .radio-inline {
        margin-left: .6rem;
    }
    .mamamia-resolution-slider {
        position: relative;
    }
    .mamamia-resolution-ticks {
        position: relative;
        height: .9rem;
        margin: -.15rem .7rem .15rem;
        color: var(--mm-muted);
        font-size: .62rem;
        line-height: 1;
    }
    .mamamia-resolution-ticks span {
        position: absolute;
        white-space: nowrap;
    }
    .mamamia-resolution-ticks span.last {
        transform: translateX(-100%);
    }
    .mamamia-resolution-head {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: .45rem;
        margin-bottom: .1rem;
    }
    .mamamia-resolution-head > .control-label {
        flex: 0 0 auto;
        margin: 0;
        color: var(--mm-muted);
        font-size: .61rem;
        font-weight: 700;
        letter-spacing: .1em;
        text-transform: uppercase;
    }
    .mamamia-resolution-head .form-group {
        flex: 0 0 6rem;
        width: 6rem !important;
        margin-bottom: 0;
    }
    .mamamia-resolution-head .form-control {
        min-height: 1.7rem;
        padding: .15rem .45rem;
        border-color: var(--mm-line);
        border-radius: .45rem;
        background-color: var(--mm-surface);
        font-size: .72rem;
    }
    .mamamia-resolution-slider .shiny-label-null {
        display: none;
    }
    .mamamia-resolution-locked {
        position: relative;
    }
    .mamamia-resolution-locked .shiny-input-container {
        margin-bottom: 0;
        opacity: .55;
        filter: grayscale(1);
    }
    .mamamia-resolution-locked .irs {
        pointer-events: none;
    }
    .mamamia-resolution-lock {
        position: absolute;
        top: 0;
        right: 0;
        bottom: 0;
        left: 0;
        width: 100%;
        height: 100%;
        padding: 0;
        border: 0;
        border-radius: .4rem;
        background: transparent !important;
        color: transparent !important;
        cursor: not-allowed !important;
        transition: none !important;
        z-index: 5;
    }
    .mamamia-resolution-lock:hover,
    .mamamia-resolution-lock:focus,
    .mamamia-resolution-lock:active,
    .mamamia-resolution-lock.active {
        background: transparent !important;
        border-color: transparent !important;
        box-shadow: none !important;
        color: transparent !important;
        cursor: not-allowed !important;
    }
    .mamamia-resolution-lock:focus-visible {
        outline: 2px solid rgba(var(--mm-accent-rgb), .45);
        outline-offset: 2px;
    }
    .mamamia-resolution-note {
        margin: .1rem 0 0;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.35;
    }
    .mamamia-fixed-control {
        display: flex;
        align-items: center;
        min-height: 2rem;
        padding: .28rem .5rem;
        border: 1px solid var(--mm-line);
        border-radius: .42rem;
        background: var(--mm-paper);
        color: var(--mm-muted);
        font-size: .79rem;
        line-height: 1.3;
        cursor: not-allowed;
    }
    .mamamia-fixed-note {
        margin: .15rem 0 .55rem;
        color: var(--mm-muted);
        font-size: .68rem;
        line-height: 1.35;
    }
    .mamamia-sidebar .table {
        margin-bottom: .25rem;
        font-size: .68rem;
    }

    /* Call table */
    .mamamia-main .dataTables_wrapper {
        font-size: .75rem;
    }
    .mamamia-main .dataTables_wrapper .row:first-child {
        align-items: center;
        margin-bottom: .45rem;
    }
    .mamamia-main .dataTables_length label,
    .mamamia-main .dataTables_filter label {
        color: var(--mm-muted);
        font-size: .68rem;
        font-weight: 550;
    }
    .mamamia-main .dataTables_length .form-select {
        height: 1.85rem;
        padding: .1rem 1.6rem .1rem .45rem;
        border-color: var(--mm-line);
        border-radius: .38rem;
        font-size: .72rem;
    }
    .mamamia-main .dataTables_wrapper .dataTables_filter input {
        height: 1.85rem;
        padding: .1rem .6rem;
        border: 1px solid var(--mm-line);
        border-radius: 999px;
        background: var(--mm-surface);
        font-size: .72rem;
        box-shadow: none;
    }
    .mamamia-main .dataTables_wrapper .dataTables_filter input:focus {
        border-color: var(--mm-accent);
        box-shadow: 0 0 0 .15rem var(--mm-accent-soft);
    }
    .mamamia-main table.dataTable {
        border-collapse: separate;
        border-spacing: 0;
    }
    .mamamia-main table.dataTable thead th {
        padding: .45rem .55rem;
        border-bottom: 1px solid var(--mm-line-strong);
        background: var(--mm-paper);
        color: var(--mm-muted);
        font-size: .6rem;
        font-weight: 750;
        letter-spacing: .06em;
        text-transform: uppercase;
        white-space: nowrap;
    }
    .mamamia-main table.dataTable tbody td {
        padding: .4rem .55rem;
        border-bottom: 1px solid var(--mm-line);
        font-size: .72rem;
        vertical-align: middle;
        white-space: nowrap;
    }
    .mamamia-main table.dataTable.table-striped > tbody > tr:nth-of-type(odd) > * {
        --bs-table-bg-type: var(--mm-paper);
    }
    .mamamia-main table.dataTable tbody tr:hover > * {
        --bs-table-bg-state: var(--mm-accent-soft);
    }
    .mamamia-main table.dataTable tbody tr.selected > * {
        --bs-table-bg-state: var(--mm-accent-soft);
        box-shadow: inset 0 0 0 9999px transparent;
        color: var(--mm-accent);
    }
    .mamamia-main table.dataTable thead tr.filters th {
        padding: .2rem .3rem .35rem;
        background: var(--mm-paper);
    }
    .mamamia-main table.dataTable thead tr.filters input {
        height: 1.65rem;
        padding: .1rem .35rem;
        border: 1px solid var(--mm-line);
        border-radius: .3rem;
        font-size: .66rem;
    }
    .mamamia-main .dataTables_info,
    .mamamia-main .dataTables_paginate {
        margin-top: .55rem;
        color: var(--mm-muted);
        font-size: .68rem;
    }
    .mamamia-main .page-link {
        padding: .2rem .5rem;
        border-color: var(--mm-line);
        color: var(--mm-muted);
        font-size: .7rem;
    }
    .mamamia-main .active > .page-link,
    .mamamia-main .page-link.active {
        border-color: var(--mm-accent);
        background: var(--mm-accent);
    }
    .mamamia-main .dataTables_empty {
        padding: 2rem 1rem !important;
        color: var(--mm-muted);
        font-size: .8rem;
    }

    .mamamia-main .shiny-output-error-validation {
        display: flex;
        align-items: center;
        justify-content: center;
        height: 100%;
        min-height: 12rem;
        padding: 2rem;
        color: var(--mm-muted);
        font-size: .85rem;
        text-align: center;
    }

    /* Analysis progress */
    .shiny-progress-notification {
        padding: .6rem .8rem;
        border: 1px solid var(--mm-line);
        border-radius: .6rem;
        background: var(--mm-surface);
        box-shadow: 0 .6rem 1.6rem rgba(27, 30, 26, .16);
    }
    .shiny-progress-notification .progress {
        height: 4px;
        margin: .3rem 0 .1rem;
        border-radius: 999px;
        background: var(--mm-line);
    }
    .shiny-progress-notification .progress-bar {
        background-color: var(--mm-accent);
    }
    .shiny-progress-notification .progress-text {
        color: var(--mm-ink);
        font-size: .78rem;
    }

    /* Toasts. bslib doubles the class and id selectors, so match that
       specificity to keep the compact spacing. */
    #shiny-notification-panel#shiny-notification-panel {
        right: 1rem;
        bottom: 1rem;
        width: min(23rem, calc(100vw - 2rem));
        padding: 0;
        z-index: 99999;
    }
    .shiny-notification.shiny-notification {
        padding: .5rem 2rem .5rem .75rem;
        margin: 0 0 .4rem;
        border: 1px solid var(--mm-line);
        border-left: 3px solid var(--mm-accent);
        border-radius: .6rem;
        background: var(--mm-surface);
        box-shadow: 0 .6rem 1.6rem rgba(27, 30, 26, .16);
        color: var(--mm-ink);
        font-size: .78rem;
        line-height: 1.4;
    }
    .shiny-notification.shiny-notification:last-child {
        margin-bottom: 0;
    }
    .shiny-notification.shiny-notification .shiny-progress-notification {
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        box-shadow: none;
    }
    .shiny-notification.shiny-notification .shiny-notification-close {
        top: .25rem;
        right: .45rem;
        width: 1.5rem;
        height: 1.5rem;
        padding: 0;
        color: var(--mm-muted);
        font-size: 1rem;
        font-weight: 400;
    }
    .shiny-notification.shiny-notification .shiny-notification-close:hover {
        color: var(--mm-ink);
        font-weight: 400;
    }

    ::selection {
        background: rgba(var(--mm-accent-rgb), .18);
    }
    .btn:focus-visible,
    .form-control:focus-visible,
    .form-select:focus-visible {
        outline: 2px solid var(--mm-accent);
        outline-offset: 2px;
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
        shiny::div(
            class = "mamamia-upload-head",
            shiny::span(
                class = "mamamia-upload-icon",
                `aria-hidden` = "true",
                shiny::icon("dna")
            ),
            shiny::div(
                class = "mamamia-upload-copy",
                shiny::span(class = "mamamia-upload-title", "BAM file"),
                shiny::span(
                    class = "mamamia-upload-hint",
                    "Drag & drop or browse"
                )
            )
        ),
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
bam_fixed_control_ui <- function(label, value) {
    shiny::div(
        class = "form-group shiny-input-container",
        shiny::tags$label(class = "control-label", label),
        shiny::tags$div(
            class = "form-control mamamia-fixed-control",
            `aria-disabled` = "true",
            value
        )
    )
}

#' @noRd
bam_sidebar_recovery_ui <- function(info, show_action = TRUE) {
    if (!is.list(info)) {
        return(NULL)
    }
    text_value <- function(x) {
        if (is.null(x) || length(x) != 1L || is.na(x)) {
            ""
        } else {
            as.character(x)
        }
    }
    stage <- text_value(info$stage)
    assembly <- identical(stage, "assembly")
    message <- text_value(info$message)
    suggestion <- text_value(info$suggestion)
    recovery_copy <- if (assembly) {
        if (grepl("missing target", message, ignore.case = TRUE)) {
            list(
                body = paste(
                    "This BAM does not contain the contigs required by",
                    "the selected assembly."
                ),
                note = "The calculation log lists the missing targets."
            )
        } else if (grepl(
            "does not match any configured assembly", message,
            ignore.case = TRUE
        )) {
            list(
                body = paste(
                    "This BAM's header does not match any assembly",
                    "in the registry."
                ),
                note = paste(
                    "Choose a matching species pair or add the assembly",
                    "to the registry."
                )
            )
        } else if (grepl("Unknown assembly", message, ignore.case = TRUE)) {
            list(
                body = paste(
                    "The BAM header names an assembly that is not in",
                    "the registry."
                ),
                note = paste(
                    "Add the assembly to the registry or choose a matching",
                    "species pair."
                )
            )
        } else {
            list(
                body = "This BAM does not match a configured assembly.",
                note = "Check the BAM header and the assembly setup."
            )
        }
    } else {
        list(
            body = if (nchar(message) > 220) {
                paste0(substr(message, 1, 217), "...")
            } else {
                message
            },
            note = NULL
        )
    }
    shiny::div(
        class = "mamamia-recovery",
        shiny::div(
            class = "mamamia-recovery-head",
            shiny::span(
                class = "mamamia-recovery-icon",
                shiny::icon(if (assembly) "dna" else "sliders")
            ),
            shiny::span(
                class = "mamamia-recovery-title",
                if (assembly) {
                    "No matching assembly"
                } else {
                    "Setup needs attention"
                }
            )
        ),
        if (nzchar(recovery_copy$body)) {
            shiny::p(class = "mamamia-recovery-message", recovery_copy$body)
        },
        if (!is.null(recovery_copy$note) && nzchar(recovery_copy$note)) {
            shiny::p(class = "mamamia-recovery-note", recovery_copy$note)
        },
        if (nzchar(suggestion)) {
            shiny::p(
                class = "mamamia-recovery-suggestion",
                shiny::strong("Next step: "),
                suggestion
            )
        },
        if (isTRUE(show_action)) {
            shiny::actionButton(
                "unlock_setup",
                if (assembly) "Choose assembly" else "Review setup",
                icon = shiny::icon("arrow-right"),
                class = "mamamia-recovery-action"
            )
        }
    )
}

#' @noRd
bam_sidebar_species_ui <- function(config, selection = NULL, locked = FALSE) {
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

    if (isTRUE(locked)) {
        species_label <- function(key) {
            matched <- species_choices[species_choices == key]
            if (length(matched) == 0L) {
                key
            } else {
                names(matched)[[1L]]
            }
        }
        return(shiny::tagList(
            bam_fixed_control_ui(
                "Recipient species", species_label(recipient)
            ),
            bam_fixed_control_ui(
                "Donor species", species_label(donor)
            ),
            shiny::p(
                class = "mamamia-fixed-note",
                paste(
                    "Fixed for the bundled example.",
                    "Upload a BAM to change the assembly."
                )
            )
        ))
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
        shiny::actionButton(
            "swap_species", "Swap donor and recipient",
            class = "btn-sm btn-outline-secondary",
            style = "width: 100%"
        )
    )
}

#' @noRd
bam_sidebar_details_ui <- function(config) {
    shiny::tagList(
        bslib::accordion(
            id = "bam_sidebar_accordion",
            multiple = FALSE,
            open = "Analysis",
            bslib::accordion_panel(
                "Analysis",
                shiny::uiOutput("resolution_ui"),
            ),
            bslib::accordion_panel(
                "Assembly",
                shiny::uiOutput("species_ui")
            ),
            bslib::accordion_panel(
                "Advanced",
                shiny::uiOutput("donor_aggregation_ui"),
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
            shiny::div(
                class = "mamamia-filter-controls",
                explorer_view_filters_ui("explorer")
            )
        ),
        bslib::navset_card_tab(
            id = ns("tabs"),
            bslib::nav_panel(
                "Karyotype",
                explorer_karyotype_ui("explorer", height = "100%")
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
        "|0.5rem|264px   |1fr      |",
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
        theme = bam_dashboard_theme(),
        shiny::tags$head(
            mamamia_font_dependency(),
            mamamia_favicon(),
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
                    "});",
                    "function mamamiaResultTabs(enabled) {",
                    "  var links = document.querySelectorAll(",
                    "    '#explorer-tabs .nav-link');",
                    "  Array.prototype.forEach.call(links, function(link) {",
                    "    var label = (link.textContent || '').trim();",
                    "    if (label !== 'Calls' && label !== 'Export') return;",
                    "    if (enabled) {",
                    "      link.classList.remove('disabled');",
                    "      link.removeAttribute('aria-disabled');",
                    "      link.removeAttribute('tabindex');",
                    "    } else {",
                    "      link.classList.add('disabled');",
                    "      link.setAttribute('aria-disabled', 'true');",
                    "      link.setAttribute('tabindex', '-1');",
                    "    }",
                    "  });",
                    "  if (!enabled) {",
                    "    var active = document.querySelector(",
                    "      '#explorer-tabs .nav-link.active');",
                    "    if (active &&",
                    "        (active.textContent || '').trim() !== 'Karyotype') {",
                    "      var first = document.querySelector(",
                    "        '#explorer-tabs .nav-link');",
                    "      if (first) first.click();",
                    "    }",
                    "  }",
                    "}",
                    "$(function() { mamamiaResultTabs(false); });",
                    "Shiny.addCustomMessageHandler('mamamia-result-tabs',",
                    "function(x) { mamamiaResultTabs(!!x.enabled); });",
                    "function mamamiaAnalysisState(x) {",
                    "  var loading = document.querySelectorAll(",
                    "    '.mamamia-analysis-loading');",
                    "  Array.prototype.forEach.call(loading, function(node) {",
                    "    if (x.running) {",
                    "      node.style.display = 'flex';",
                    "      node.classList.add('is-running');",
                    "      node.setAttribute('aria-hidden', 'false');",
                    "    } else {",
                    "      node.classList.remove('is-running');",
                    "      node.setAttribute('aria-hidden', 'true');",
                    "      node.style.display = 'none';",
                    "    }",
                    "  });",
                    "  if (x.stage) {",
                    "    var stage = document.querySelectorAll(",
                    "      '.mamamia-analysis-stage');",
                    "    Array.prototype.forEach.call(stage, function(node) {",
                    "      node.textContent = x.stage;",
                    "    });",
                    "  }",
                    "  if (x.progress !== null && x.progress !== undefined) {",
                    "    var fill = document.querySelectorAll(",
                    "      '.mamamia-analysis-fill');",
                    "    var width = Math.round(",
                    "      100 * Math.max(0, Math.min(1, x.progress)));",
                    "    Array.prototype.forEach.call(fill, function(node) {",
                    "      node.style.width = width + '%';",
                    "    });",
                    "  }",
                    "  var run = document.getElementById('run');",
                    "  if (run) {",
                    "    if (x.running) {",
                    "      run.classList.add('disabled');",
                    "      run.setAttribute('disabled', 'disabled');",
                    "      run.setAttribute('aria-busy', 'true');",
                    "    } else {",
                    "      run.classList.remove('disabled');",
                    "      run.removeAttribute('disabled');",
                    "      run.removeAttribute('aria-busy');",
                    "    }",
                    "  }",
                    "}",
                    "Shiny.addCustomMessageHandler('mamamia-analysis-state',",
                    "function(x) { mamamiaAnalysisState(x); });"
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
                    "    var zone = uploadZone(node);",
                    "    if (!zone) return null;",
                    "    var input = zone.querySelector('input[type=\"file\"]');",
                    "    return input && input.id === id ? input : null;",
                    "  }",
                    "  function uploadZone(node) {",
                    "    return node && node.closest ?",
                    "      node.closest('.mamamia-upload') : null;",
                    "  }",
                    "  document.addEventListener('change', function(event) {",
                    "    var input = event.target;",
                    "    if (!input || input.id !== id || !input.files) return;",
                    "    if (input.files.length > 1) {",
                    "      input.files = onlyFirst(input.files);",
                    "    }",
                    "  }, true);",
                    "  document.addEventListener('dragover', function(event) {",
                    "    var zone = uploadZone(event.target);",
                    "    if (!zone) return;",
                    "    event.preventDefault();",
                    "    if (event.dataTransfer) {",
                    "      event.dataTransfer.dropEffect = 'copy';",
                    "    }",
                    "    zone.classList.add('is-dragover');",
                    "  }, true);",
                    "  document.addEventListener('dragleave', function(event) {",
                    "    var zone = uploadZone(event.target);",
                    "    if (!zone) return;",
                    "    if (!zone.contains(event.relatedTarget)) {",
                    "      zone.classList.remove('is-dragover');",
                    "    }",
                    "  }, true);",
                    "  document.addEventListener('click', function(event) {",
                    "    var zone = uploadZone(event.target);",
                    "    if (!zone) return;",
                    "    var own = event.target.closest(",
                    "      'input[type=\"file\"], .btn-file, label');",
                    "    if (own) return;",
                    "    var input = zone.querySelector('input[type=\"file\"]');",
                    "    if (input) input.click();",
                    "  }, true);",
                    "  document.addEventListener('drop', function(event) {",
                    "    var zone = uploadZone(event.target);",
                    "    if (zone) zone.classList.remove('is-dragover');",
                    "    var input = uploadInput(event.target);",
                    "    var files = event.dataTransfer && event.dataTransfer.files;",
                    "    if (!input || !files || files.length < 1) return;",
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
                        class = "mamamia-brand",
                        shiny::span(
                            class = "mamamia-mark",
                            `aria-hidden` = "true",
                            shiny::HTML(paste0(
                                '<svg viewBox="0 0 32 32" fill="none" ',
                                'xmlns="http://www.w3.org/2000/svg">',
                                '<path d="M11 6 C 11 9, 11 12, 16 16',
                                ' C 21 20, 21 23, 21 26" stroke="currentColor" ',
                                'stroke-width="2.2" stroke-linecap="round"/>',
                                '<path d="M21 6 C 21 9, 21 12, 16 16',
                                ' C 11 20, 11 23, 11 26" stroke="#d9a441" ',
                                'stroke-width="2.2" stroke-linecap="round"/>',
                                '</svg>'
                            ))
                        ),
                        shiny::div(
                            class = "mamamia-brand-text",
                            shiny::div(
                                class = "mamamia-title-row",
                                shiny::h1("MaMaMIA", class = "mamamia-title"),
                                shiny::span(
                                    class = "mamamia-version",
                                    mamamia_package_version()
                                )
                            ),
                            shiny::p(
                                class = "mamamia-subtitle",
                                paste(
                                    "Max Mapping for Meticulous Introgression",
                                    "Analysis"
                                )
                            )
                        )
                    ),
                    shiny::span(
                        class = "mamamia-header-note",
                        "Breadth-aware introgression mapping"
                    )
                )
            ),
            gridlayout::grid_place(
                "sidebar",
                shiny::tags$div(
                    class = "mamamia-sidebar",
                    shiny::div(
                        class = "mamamia-sidebar-intro",
                        shiny::uiOutput("bam_mode_note")
                    ),
                    shiny::div(
                        class = "mamamia-sidebar-section",
                        shiny::div(class = "mamamia-section-label", "Input"),
                        shiny::uiOutput("bam_upload_ui"),
                        if (example_available) {
                            shiny::actionButton(
                                "example_data", "Example data",
                                class = "btn mamamia-ghost-btn",
                                title = paste(
                                    "Load the bundled example dataset"
                                )
                            )
                        }
                    ),
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

        ## The bundled example has 200 kb contigs, so its 2 kb calls would
        ## disappear under the 5 Mb production filter. Use a toy-scale filter
        ## while example mode is active and restore the default afterwards.
        shiny::observeEvent(example_mode(), {
            shiny::updateNumericInput(
                session, "explorer-min_width",
                value = if (isTRUE(example_mode())) 0.001 else 5
            )
            if (isTRUE(example_mode())) {
                shiny::updateSliderInput(
                    session, "resolution",
                    value = as.integer(bam_example_default_resolution())
                )
            }
        }, ignoreInit = TRUE)

        result <- shiny::reactiveVal(NULL)
        shiny::observeEvent(result(), {
            session$sendCustomMessage(
                "mamamia-result-tabs",
                list(enabled = !is.null(result()))
            )
        }, ignoreInit = FALSE, ignoreNULL = FALSE)
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

        send_analysis_state <- function(running, stage = NULL,
                                        progress = NULL) {
            try(
                session$sendCustomMessage(
                    "mamamia-analysis-state",
                    list(
                        running = isTRUE(running),
                        stage = stage,
                        progress = progress
                    )
                ),
                silent = TRUE
            )
            invisible(NULL)
        }

        output$bam_mode_note <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                return(shiny::p(
                    class = "text-muted small mb-2",
                    "Example data: events on eight chromosomes."
                ))
            }
            pf <- preflight_state()
            if (isTRUE(pf$empty)) {
                return(shiny::p(
                    class = "text-muted small mb-2",
                    paste(
                        "One BAM per run. The BAM header is matched to the",
                        "configured assembly."
                    )
                ))
            }
            if (!isTRUE(pf$ok) || is.null(pf$table)) {
                return(NULL)
            }
            shiny::p(
                class = "text-muted small mb-2",
                "BAM ready. The assembly was matched from the header."
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
            log_line(bam_example_note())
            log_line("Example data selected.")
        }, ignoreInit = FALSE)

        shiny::observeEvent(input$unlock_setup, {
            pf <- preflight_state()
            panel <- if (!is.null(pf$error) &&
                identical(as.character(pf$error$stage), "assembly")) {
                "Assembly"
            } else {
                "Analysis"
            }
            sidebar_unlocked(TRUE)
            session$onFlushed(function() {
                try(
                    bslib::accordion_panel_open(
                        "bam_sidebar_accordion", panel
                    ),
                    silent = TRUE
                )
            }, once = TRUE)
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

        output$resolution_ui <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                return(bam_example_resolution_slider_ui(config$defaults))
            }
            bam_resolution_slider_ui(config$defaults)
        })

        ## The slider stays the canonical control, but a valid typed value
        ## wins immediately so a run click cannot use a stale slider position.
        resolution_value <- shiny::reactive({
            if (isTRUE(example_mode())) {
                return(as.integer(bam_example_default_resolution()))
            }
            defaults <- config$defaults
            value <- suppressWarnings(as.integer(input$resolution_number))
            if (length(value) == 1L && !is.na(value)) {
                typed <- tryCatch(
                    bam_resolution_value(value, defaults),
                    error = function(e) NULL
                )
                if (!is.null(typed)) {
                    return(typed)
                }
            }
            tryCatch(
                bam_resolution_value(input$resolution, defaults),
                error = function(e) bam_default_resolution(defaults)
            )
        })

        notify_example_resolution <- function() {
            resolution <- bam_example_default_resolution()
            shiny::showNotification(
                sprintf(
                    paste(
                        "The bundled example is fixed at %s windows",
                        "per chromosome. Upload a BAM to choose",
                        "another resolution."
                    ),
                    format(
                        as.integer(resolution),
                        big.mark = ",", scientific = FALSE
                    )
                ),
                type = "warning",
                duration = 7
            )
        }

        shiny::observeEvent(input$resolution_locked, {
            notify_example_resolution()
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$resolution, {
            if (!isTRUE(example_mode())) {
                return()
            }
            expected <- suppressWarnings(
                as.integer(bam_example_default_resolution())
            )
            current <- suppressWarnings(as.integer(input$resolution))
            if (length(expected) == 1L && !is.na(expected) &&
                !identical(current, expected)) {
                shiny::updateSliderInput(
                    session, "resolution", value = expected
                )
                notify_example_resolution()
            }
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$resolution, {
            value <- suppressWarnings(as.integer(input$resolution))
            if (length(value) == 1L && !is.na(value)) {
                shiny::updateNumericInput(
                    session, "resolution_number", value = value
                )
            }
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$resolution_number, {
            value <- suppressWarnings(as.integer(input$resolution_number))
            if (length(value) != 1L || is.na(value)) {
                current <- suppressWarnings(as.integer(input$resolution))
                if (length(current) == 1L && !is.na(current)) {
                    shiny::updateNumericInput(
                        session, "resolution_number", value = current
                    )
                }
                return()
            }
            if (isTRUE(example_mode())) {
                expected <- suppressWarnings(
                    as.integer(bam_example_default_resolution())
                )
                if (length(expected) == 1L && !is.na(expected) &&
                    !identical(value, expected)) {
                    shiny::updateNumericInput(
                        session, "resolution_number", value = expected
                    )
                    notify_example_resolution()
                }
                return()
            }
            valid <- tryCatch(
                bam_resolution_value(value, config$defaults),
                error = function(e) e
            )
            if (inherits(valid, "condition")) {
                current <- suppressWarnings(as.integer(input$resolution))
                if (length(current) == 1L && !is.na(current)) {
                    shiny::updateNumericInput(
                        session, "resolution_number", value = current
                    )
                }
                shiny::showNotification(
                    conditionMessage(valid),
                    type = "warning", duration = 7
                )
                return()
            }
            if (!identical(
                valid,
                suppressWarnings(as.integer(input$resolution))
            )) {
                shiny::updateSliderInput(
                    session, "resolution", value = valid
                )
            }
        }, ignoreInit = TRUE)

        output$donor_aggregation_ui <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                value <- bam_example_config()$defaults$donor_aggregation
                label <- switch(
                    as.character(value),
                    max = "Max",
                    mean = "Mean",
                    sum = "Sum",
                    as.character(value)
                )
                return(shiny::tagList(
                    bam_fixed_control_ui("Donor aggregation", label),
                    shiny::p(
                        class = "mamamia-fixed-note",
                        "Fixed for the bundled example."
                    )
                ))
            }
            defaults <- config$defaults
            shiny::tagList(
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
                )
            )
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
            pf <- preflight_state()
            if (isTRUE(pf$empty) || isTRUE(pf$ok) || is.null(pf$error)) {
                return(NULL)
            }
            bam_sidebar_recovery_ui(
                pf$error,
                show_action = !isTRUE(sidebar_unlocked())
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
                stage <- if (is.null(pf$error)) NULL else pf$error$stage
                assembly_error <- identical(as.character(stage), "assembly")
                label <- if (assembly_error) {
                    "Needs assembly"
                } else {
                    "Preflight failed"
                }
                badge <- if (assembly_error) {
                    "text-bg-warning"
                } else {
                    "text-bg-danger"
                }
                detail <- if (is.null(pf$error)) {
                    "See the calculation log."
                } else if (assembly_error) {
                    "Choose a matching assembly in the sidebar."
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
                    selection = bam_example_species_selection(),
                    locked = TRUE
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
                resolution <- resolution_value()
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
                    shiny::span(class = "action-label", "Run analysis"),
                    shiny::icon("arrow-right", `aria-hidden` = "true")
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
                resolution_value(),
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
                bam_run_parameters(
                    input, pf$config, samples,
                    resolution = resolution_value()
                ),
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
            if (isTRUE(pf$example)) {
                params$resolution <- as.integer(
                    bam_example_default_resolution()
                )
                params$donor_aggregation <-
                    pf$config$defaults$donor_aggregation
                params$control_samples <- character(0)
                params$calibration <- NULL
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
                    "Starting analysis: %s windows/chromosome,",
                    "metric '%s', score '%s', donor aggregation '%s'."
                ),
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
                bam_progress_callback(
                    length(samples),
                    total_steps = total_steps,
                    handler = function(stage, detail, progress) {
                        label <- paste0(
                            prefix, bam_progress_label(stage, detail)
                        )
                        log_line(label)
                        send_analysis_state(TRUE, label, progress)
                    }
                )
            }

            send_analysis_state(TRUE, "Preparing analysis", 0)
            session$onFlushed(
                function() send_analysis_state(FALSE),
                once = TRUE
            )
            new_result <- tryCatch(
                {
                    primary <- run_analysis(
                        pf$config, pf$assembly, progress_for()
                    )
                    result_cache$set(key, primary)
                    active_key(key)
                    primary
                },
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
            n_calls <- sum(new_result$samples$n_calls, na.rm = TRUE)
            call_word <- if (n_calls == 1L) {
                "calibrated call"
            } else {
                "calibrated calls"
            }
            log_line(sprintf(
                "Analysis finished: %d %s.",
                n_calls, call_word
            ))
            run_summary <- sprintf("%d %s", n_calls, call_word)
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
        ## karyotype, call table and exports. It is only useful when the input
        ## BAM contains several @RG SM samples.
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
            multiple_samples <- length(samples) > 1L
            shiny::tagList(
                shiny::div(
                    class = "table-responsive",
                    shiny::tableOutput("cohort_table")
                ),
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
