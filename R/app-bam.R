## BAM-first Shiny application (T6.1-T6.4).
##
## The app is deliberately thin: configuration, preflight, run orchestration
## and cohort export live here, while counting, correction, normalization and
## segmentation stay in the pipeline. Uploads and server-side BAMs are staged
## by symlink into a per-session directory before preflight, so the pipeline
## never writes an index under `server.bam_root` and the session can remove
## every staged file and index when it ends.

## Error taxonomy ------------------------------------------------------------

.bam_error_suggestions <- c(
    config = paste(
        "Check MAMAMIA_CONFIG and the YAML assembly registry;",
        "validate it with readMaMaMIAConfig()."
    ),
    BAM = paste(
        "Check that the file is a coordinate-sorted BAM;",
        "re-upload it or pick another file."
    ),
    assembly = paste(
        "Pick an assembly override that matches the BAM header, or use a BAM",
        "that contains every configured target."
    ),
    index = paste(
        "The BAM directory may be read-only: upload the file instead,",
        "or rebuild the CSI index."
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
        "Fewer than 2 ideal windows: lower the resolution, raise the max",
        "fit windows cap, or inspect coverage and GC outliers."
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
    analysis = "See the message above and the server log."
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
        grepl("^Normalization failed", message)) {
        return("normalization")
    }
    "analysis"
}

#' Describe a user-facing error with a stage and a suggested fix
#'
#' @param e A condition.
#'
#' @return A list with `stage`, `class`, `message` and `suggestion`.
#' @noRd
bam_error_info <- function(e) {
    stage <- bam_error_stage(e)
    suggestion <- .bam_error_suggestions[[stage]]
    if (is.null(suggestion) || is.na(suggestion)) {
        suggestion <- .bam_error_suggestions[["analysis"]]
    }
    list(
        stage = stage,
        class = class(e)[[1L]],
        message = conditionMessage(e),
        suggestion = suggestion
    )
}

#' @noRd
bam_error_text <- function(info) {
    sprintf(
        "[%s] %s\nSuggested fix: %s",
        info$stage, info$message, info$suggestion
    )
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
            !is.null(config$server) && !is.null(config$assemblies)) {
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
            error = bam_error_info(parsed)
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
bam_config_summary <- function(config) {
    assemblies <- listMaMaMIAAssemblies(config)
    labels <- if (nrow(assemblies) == 0L) {
        "none"
    } else {
        paste(
            sprintf("%s (%s)", assemblies$label, assemblies$key),
            collapse = "; "
        )
    }
    sprintf(
        "%d assembl%s: %s",
        nrow(assemblies),
        if (nrow(assemblies) == 1L) "y" else "ies",
        labels
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
        stop("Choose an analysis resolution", call. = FALSE)
    }
    resolution <- suppressWarnings(as.integer(value))
    if (is.na(resolution) || !resolution %in% allowed) {
        stop(
            "Resolution must be one of: ",
            paste(allowed, collapse = ", "),
            call. = FALSE
        )
    }
    resolution
}

#' @noRd
bam_integer <- function(value, label, minimum, maximum) {
    if (length(value) != 1L || is.na(value) || !is.numeric(value) ||
        !is.finite(value) || value != floor(value) ||
        value < minimum || value > maximum) {
        stop(
            sprintf(
                "%s must be a whole number between %s and %s",
                label,
                format(minimum, scientific = FALSE),
                format(maximum, scientific = FALSE)
            ),
            call. = FALSE
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
bam_path_labels <- function(paths, root = NULL) {
    if (is.null(root) || !nzchar(root)) {
        return(paths)
    }
    root <- sub(
        "/+$", "",
        normalizePath(root, winslash = "/", mustWork = FALSE)
    )
    vapply(paths, function(path) {
        resolved <- normalizePath(path, winslash = "/", mustWork = FALSE)
        prefix <- paste0(root, "/")
        if (startsWith(resolved, prefix)) {
            substring(resolved, nchar(prefix) + 1L)
        } else {
            resolved
        }
    }, character(1), USE.NAMES = FALSE)
}

#' List the server-side BAM picker choices for one session
#'
#' The scan is bounded by `max_depth` and `max_files`; the listing is computed
#' once per session and never refreshed automatically while the session lives.
#'
#' @param config A validated configuration.
#' @param initial Optional paths to include and preselect (the `bam_paths`
#'   argument of the launcher).
#'
#' @return A list with `choices` (a named character vector), `selected`,
#'   `truncated` and `error`.
#' @noRd
bam_server_choices <- function(config, initial = NULL, max_depth = 3L,
                               max_files = 2000L) {
    root <- config$server$bam_root
    listing_error <- NULL
    listed <- character()
    if (!is.null(root)) {
        listed <- tryCatch(
            .list_bam_files(root, max_depth = max_depth, max_files = max_files),
            error = function(e) {
                listing_error <<- conditionMessage(e)
                character()
            }
        )
    }
    truncated <- isTRUE(attr(listed, "truncated"))

    initial <- if (is.null(initial)) character() else as.character(initial)
    initial <- initial[!is.na(initial) & nzchar(initial)]
    paths <- unique(c(initial, listed))
    list(
        choices = stats::setNames(paths, bam_path_labels(paths, root)),
        selected = initial[initial %in% paths],
        truncated = truncated,
        error = listing_error
    )
}

#' @noRd
bam_listing_note <- function(listing) {
    if (!is.null(listing$error)) {
        return(paste("Could not list server BAMs:", listing$error))
    }
    if (isTRUE(listing$truncated)) {
        return(sprintf(
            "The server BAM listing is truncated to the first %d file(s).",
            length(listing$choices)
        ))
    }
    ""
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
#' @param paths BAM paths.
#' @param session_dir Private per-session directory.
#' @param bam_root Optional `server.bam_root`. Every source must be under this
#'   root or in session temporary storage (uploads); the check runs on the
#'   source, before a small file is copied, so a copy cannot bypass D12.
#'
#' @return A data frame with `path` (staged), `source` and `dir`.
#' @noRd
bam_stage_paths <- function(paths, session_dir, bam_root = NULL) {
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
        contained <- tryCatch(
            {
                .check_bam_containment(
                    source, bam_root, session_dir = .session_temp_dir()
                )
                TRUE
            },
            error = function(e) e
        )
        if (inherits(contained, "condition")) {
            unlink(stage_dir, recursive = TRUE, force = TRUE)
            stop(contained)
        }
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
bam_preflight_failure <- function(e, paths, config) {
    list(
        ok = FALSE,
        empty = FALSE,
        table = NULL,
        staged = NULL,
        assembly = NULL,
        error = bam_error_info(e),
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
        stop("At least one BAM is required", call. = FALSE)
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
        stop(
            "Max fit windows must be a positive whole number",
            call. = FALSE
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
            stop(
                sprintf(
                    "Reference sample '%s' is not one of: %s",
                    reference, paste(samples, collapse = ", ")
                ),
                call. = FALSE
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
    label <- .bam_progress_labels[[stage]]
    if (is.null(label) || is.na(label)) {
        label <- stage
    }
    sample <- detail$sample
    if (!is.null(sample) && length(sample) == 1L && !is.na(sample) &&
        nzchar(sample)) {
        label <- paste0(label, " for ", sample)
    }
    label
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

## UI ------------------------------------------------------------------------

#' @noRd
bam_launcher_ui <- function(state, bam_paths = NULL) {
    if (!isTRUE(state$ok)) {
        return(bslib::page_sidebar(
            title = "MaMaMIA",
            bslib::card(
                bslib::card_header("Configuration error"),
                shiny::p(state$error$message),
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
    defaults <- config$defaults
    show_picker <- !is.null(config$server$bam_root) ||
        length(bam_paths) > 0L
    assembly_keys <- names(config$assemblies)
    assembly_labels <- vapply(
        config$assemblies, function(x) as.character(x$label), character(1)
    )
    assembly_choices <- c(
        "Auto-detect" = "",
        stats::setNames(
            assembly_keys,
            sprintf("%s (%s)", assembly_labels, assembly_keys)
        )
    )

    bslib::page_sidebar(
        title = "MaMaMIA",
        sidebar = bslib::sidebar(
            width = 400,
            shiny::h6("Configuration"),
            shiny::p(class = "text-muted", bam_config_summary(config)),
            shiny::fileInput(
                "bams", "BAM files",
                multiple = TRUE, accept = ".bam"
            ),
            if (show_picker) {
                shiny::tagList(
                    shinyWidgets::pickerInput(
                        "server_bams", "Server-side BAMs",
                        choices = NULL, selected = NULL,
                        multiple = TRUE,
                        options = list(
                            `actions-box` = TRUE,
                            `live-search` = TRUE
                        )
                    ),
                    shiny::textOutput("server_note")
                )
            },
            shiny::radioButtons(
                "resolution", "Resolution",
                choices = bam_resolution_choices(defaults),
                selected = bam_default_resolution(defaults)
            ),
            shiny::uiOutput("reference_ui"),
            shiny::selectInput(
                "assembly", "Assembly",
                choices = assembly_choices, selected = ""
            ),
            shiny::tags$details(
                shiny::tags$summary("Advanced"),
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
            shiny::uiOutput("run_ui"),
            shiny::p(
                class = "text-muted",
                paste(
                    "The analysis runs in this session and blocks until it",
                    "finishes; the progress bar reports each stage."
                )
            ),
            shiny::hr(),
            shiny::h6("Explorer filters"),
            explorer_controls_ui("explorer")
        ),
        bslib::card(
            bslib::card_header("Preflight"),
            shiny::tableOutput("preflight_table")
        ),
        bslib::card(
            bslib::card_header("Detected assembly"),
            shiny::uiOutput("assembly_summary")
        ),
        bslib::card(
            bslib::card_header("Problems"),
            shiny::uiOutput("preflight_error"),
            shiny::uiOutput("candidates_ui")
        ),
        shiny::uiOutput("cohort_ui"),
        explorer_results_ui("explorer")
    )
}

## Server --------------------------------------------------------------------

#' @noRd
bam_launcher_server <- function(state, bam_paths = NULL, cache = NULL) {
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

        listing <- bam_server_choices(config, bam_paths)
        if (length(listing$choices) > 0L) {
            shinyWidgets::updatePickerInput(
                session, "server_bams",
                choices = listing$choices,
                selected = listing$selected
            )
        }
        output$server_note <- shiny::renderText(bam_listing_note(listing))

        result <- shiny::reactiveVal(NULL)
        selection_map <- shiny::reactiveVal(list())

        uploaded_paths <- shiny::reactive({
            uploads <- input$bams
            if (is.data.frame(uploads) && "datapath" %in% names(uploads)) {
                paths <- as.character(uploads$datapath)
                paths <- paths[!is.na(paths) & nzchar(paths)]
                return(paths)
            }
            character()
        })

        selected_paths <- shiny::reactive({
            paths <- c(uploaded_paths(), as.character(input$server_bams))
            paths <- paths[!is.na(paths) & nzchar(paths)]
            unique(paths)
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
                value <- bam_stage_paths(
                    paths, session_dir, config$server$bam_root
                )
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

        assembly_override <- shiny::reactive({
            value <- input$assembly
            if (is.null(value) || length(value) != 1L || is.na(value) ||
                !nzchar(value)) {
                NULL
            } else {
                as.character(value)
            }
        })

        preflight_state <- shiny::reactive({
            paths <- selected_paths()
            if (length(paths) == 0L) {
                return(list(
                    ok = FALSE, empty = TRUE, table = NULL, staged = NULL,
                    assembly = NULL, error = NULL, candidates = NULL
                ))
            }
            staged <- tryCatch(staged_bams(), error = function(e) e)
            if (inherits(staged, "condition")) {
                return(bam_preflight_failure(
                    staged, paths, config
                ))
            }
            table <- tryCatch(
                preflightBams(
                    staged$path, config,
                    assembly = assembly_override()
                ),
                error = function(e) e
            )
            if (inherits(table, "condition")) {
                return(bam_preflight_failure(
                    table, staged$path, config
                ))
            }
            table$source <- staged$source
            list(
                ok = TRUE, empty = FALSE, table = table,
                staged = staged$path, assembly = table$assembly[[1L]],
                error = NULL, candidates = NULL
            )
        })

        output$preflight_table <- shiny::renderTable({
            pf <- preflight_state()
            if (isTRUE(pf$empty)) {
                return(data.frame(
                    status = "Select one or more BAM files to see the preflight.",
                    stringsAsFactors = FALSE
                ))
            }
            if (!isTRUE(pf$ok)) {
                return(data.frame(
                    status = "Preflight failed; see the problems card below.",
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
                    "Select one or more BAM files to detect an assembly."
                ))
            }
            if (!isTRUE(pf$ok)) {
                return(shiny::p(
                    class = "text-muted",
                    "No assembly was detected. See the problems card below."
                ))
            }
            assembly <- config$assemblies[[pf$assembly]]
            bam_assembly_summary_ui(
                assembly, config$defaults$allowed_windows
            )
        })

        output$preflight_error <- shiny::renderUI({
            pf <- preflight_state()
            if (isTRUE(pf$empty)) {
                return(NULL)
            }
            if (isTRUE(pf$ok)) {
                return(shiny::p(
                    class = "text-muted", "No preflight problems detected."
                ))
            }
            info <- pf$error
            shiny::tagList(
                shiny::h6(paste0("Stage: ", info$stage)),
                shiny::p(info$message),
                shiny::p(
                    class = "text-muted",
                    paste0("Suggested fix: ", info$suggestion)
                )
            )
        })

        output$candidates_ui <- shiny::renderUI({
            pf <- preflight_state()
            if (is.null(pf$candidates)) {
                return(NULL)
            }
            shiny::tagList(
                shiny::h6("Assembly candidates"),
                shiny::tableOutput("candidate_table")
            )
        })

        output$candidate_table <- shiny::renderTable({
            pf <- preflight_state()
            bam_candidate_display(pf$candidates)
        },
        striped = TRUE, bordered = TRUE, spacing = "xs",
        rownames = FALSE
        )

        output$reference_ui <- shiny::renderUI({
            pf <- preflight_state()
            if (!isTRUE(pf$ok) || nrow(pf$table) < 2L) {
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

        output$run_ui <- shiny::renderUI({
            ready <- isTRUE(preflight_state()$ok)
            shiny::tags$button(
                id = "run", type = "button",
                class = paste(
                    "btn btn-primary action-button",
                    if (!ready) "disabled"
                ),
                disabled = if (ready) NULL else "",
                shiny::span(class = "action-label", "Run analysis")
            )
        })

        ## Changing any input that affects the result invalidates the previous
        ## one, so a stale BEDPE can never be downloaded for a different input.
        run_key <- shiny::reactive({
            list(
                paths = selected_paths(),
                resolution = input$resolution,
                assembly = input$assembly,
                reference = input$reference_sample,
                mapq = input$mapq,
                exclude_flags = input$exclude_flags,
                max_fit_windows = input$max_fit_windows
            )
        })
        shiny::observeEvent(run_key(), {
            result(NULL)
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
                shiny::showNotification(
                    "Fix the preflight problems before running the analysis.",
                    type = "warning", duration = 8
                )
                return()
            }
            samples <- unique(pf$table$sample)
            params <- tryCatch(
                bam_run_parameters(input, config, samples),
                error = function(e) e
            )
            if (inherits(params, "condition")) {
                shiny::showNotification(
                    conditionMessage(params),
                    type = "error", duration = 12
                )
                return()
            }

            previous <- result()
            new_result <- tryCatch(
                shiny::withProgress(
                    message = "MaMaMIA analysis", value = 0,
                    {
                        analyzeBamInput(
                            pf$staged,
                            config,
                            assembly = assembly_override(),
                            windows_per_chromosome = params$resolution,
                            reference_sample = params$reference_sample,
                            normalization = params$normalization,
                            mapq = params$mapq,
                            exclude_flags = params$exclude_flags,
                            cores = params$cores,
                            max_fit_windows = params$max_fit_windows,
                            cache = cache,
                            progress = bam_progress_callback(
                                length(samples), params$normalization
                            )
                        )
                    }
                ),
                error = function(e) {
                    shiny::showNotification(
                        bam_error_text(bam_error_info(e)),
                        type = "error", duration = 20
                    )
                    NULL
                }
            )
            if (is.null(new_result)) {
                result(previous)
                return()
            }
            result(new_result)
            shiny::showNotification(
                sprintf(
                    "%d sample(s), %d segment(s), %d putative introgression(s)",
                    length(new_result$isa),
                    sum(new_result$samples$n_segments),
                    sum(new_result$samples$n_putative_introgressions)
                ),
                type = "message", duration = 6
            )
        })

        output$cohort_ui <- shiny::renderUI({
            res <- result()
            if (is.null(res)) {
                return(NULL)
            }
            samples <- bam_sample_choices(res)
            warnings <- bam_coverage_warnings(res)
            bslib::card(
                bslib::card_header("Cohort"),
                if (length(samples) > 1L) {
                    shiny::selectInput(
                        "result_sample", "Sample",
                        choices = samples,
                        selected = shiny::isolate(result_sample())
                    )
                },
                shiny::tableOutput("cohort_table"),
                if (length(warnings) > 0L) {
                    shiny::div(
                        class = "alert alert-warning",
                        shiny::tags$ul(lapply(warnings, shiny::tags$li))
                    )
                },
                shiny::downloadButton(
                    "download_cohort",
                    "Download selected BEDPE (all samples)"
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
bam_launcher_app <- function(config = NULL, bam_paths = NULL, cache = NULL) {
    require_dashboard_packages()
    if (!is.null(bam_paths) &&
        (!is.character(bam_paths) || anyNA(bam_paths))) {
        stop(
            "`bam_paths` must be NULL or a character vector of BAM paths",
            call. = FALSE
        )
    }
    state <- bam_launcher_state(config)
    ui <- bam_launcher_ui(state, bam_paths)
    server <- bam_launcher_server(state, bam_paths, cache)
    shiny::shinyApp(ui, server)
}
