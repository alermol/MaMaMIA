## Shiny application for BAM input (T6.1-T6.4).
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

.bamhusk_url <- "https://github.com/alermol/bamhusk"
.bamhusk_command <- "husk --drop-seq -@ 8 input.bam output.bam"

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
        "Check that the score table and window grid use the same windows;",
        "cohort recurrence requires at least two samples."
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
    upload = paste(
        "Reduce the BAM locally with bamhusk",
        "(https://github.com/alermol/bamhusk), then upload the smaller file:",
        .bamhusk_command
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
    if (inherits(e, "MaMaMIBamTooLargeError")) {
        return("upload")
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
#' @return A list with `stage`, `class`, `message`, `files`, `suggestion`,
#'   and optional `size`/`limit` fields for upload-limit errors.
#' @noRd
bam_error_info <- function(e, files = NULL, table = NULL) {
    stage <- bam_error_stage(e)
    suggestion <- .bam_error_suggestions[[stage]]
    if (is.null(suggestion) || is.na(suggestion)) {
        suggestion <- .bam_error_suggestions[["analysis"]]
    }
    data <- if (inherits(e, "condition")) e$data else NULL
    list(
        stage = stage,
        class = class(e)[[1L]],
        message = conditionMessage(e),
        files = bam_error_files(e, files = files, table = table),
        suggestion = suggestion,
        size = if (is.list(data)) data$size else NULL,
        limit = if (is.list(data)) data$limit else NULL
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

#' Apply the upload policy for the current deployment
#'
#' Only the ShinyProxy deployment enforces the registry size cap. A direct
#' container run removes it, so uploads are bounded only by local disk space.
#'
#' @param config A validated configuration list.
#' @return The configuration with the effective upload cap.
#' @noRd
mamamia_apply_upload_policy <- function(config) {
    if (!is.list(config)) {
        return(config)
    }
    if (identical(mamamia_deployment_mode(), "direct")) {
        if (is.null(config$server) || !is.list(config$server)) {
            config$server <- list()
        }
        config$server$max_total_bytes <- Inf
    }
    config
}

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
                        "MAMAMIA_CONFIG to the mounted registry."
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
    parsed <- mamamia_apply_upload_policy(parsed)
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
            shiny::uiOutput("resolution_param_label", inline = TRUE),
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
    cores <- mamamia_available_cores()
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

#' Format a byte count for user-facing messages
#' @noRd
bam_format_size <- function(bytes) {
    if (!is.numeric(bytes) || length(bytes) != 1L || is.na(bytes) ||
        !is.finite(bytes) || bytes < 0) {
        return("unknown size")
    }
    units <- c("B", "KiB", "MiB", "GiB", "TiB", "PiB")
    value <- as.numeric(bytes)
    unit <- 1L
    while (value >= 1024 && unit < length(units)) {
        value <- value / 1024
        unit <- unit + 1L
    }
    if (unit == 1L) {
        sprintf(
            "%s B",
            format(round(value), big.mark = ",", scientific = FALSE)
        )
    } else {
        sprintf("%.1f %s", value, units[[unit]])
    }
}

#' Reject an oversized upload before it is staged or parsed
#'
#' Uses the size reported by Shiny's file input, so the server can stop before
#' copying the file into the session directory or reading its header. The
#' generic preflight size check remains as a fallback for non-UI callers.
#'
#' @param uploads The data frame Shiny stores in `input$bams`.
#' @param max_bytes Maximum accepted size in bytes.
#'
#' @return The total size invisibly, or a `MaMaMIBamTooLargeError`.
#' @noRd
bam_upload_limit <- function(uploads, max_bytes) {
    if (!is.data.frame(uploads) || nrow(uploads) == 0L ||
        !all(c("name", "size") %in% names(uploads))) {
        return(invisible(NULL))
    }
    sizes <- suppressWarnings(as.numeric(uploads$size))
    if (length(sizes) == 0L || anyNA(sizes) ||
        any(!is.finite(sizes)) || any(sizes < 0)) {
        return(invisible(NULL))
    }
    total <- sum(sizes)
    if (!is.numeric(max_bytes) || length(max_bytes) != 1L ||
        is.na(max_bytes) || !is.finite(max_bytes) || max_bytes <= 0) {
        return(invisible(total))
    }
    if (total <= max_bytes) {
        return(invisible(total))
    }
    labels <- as.character(uploads$name)
    labels[is.na(labels) | !nzchar(labels)] <- "upload.bam"
    label <- if (length(labels) == 1L) {
        sprintf("'%s'", labels[[1L]])
    } else {
        paste(sprintf("'%s'", labels), collapse = ", ")
    }
    .mamamia_error(
        "MaMaMIBamTooLargeError",
        "BAM %s is %s; this server accepts at most %s.",
        label,
        bam_format_size(total),
        bam_format_size(max_bytes),
        data = list(path = labels, size = total, limit = max_bytes)
    )
}

#' Build error details from the browser's oversized-upload event
#'
#' The client sends `name` and `size` after it clears the file input. The
#' server trusts its own configured limit rather than the client value.
#'
#' @param payload The list Shiny receives in `input$oversized_upload`.
#' @param max_bytes Server limit in bytes.
#'
#' @return A `bam_error_info()` list, or `NULL` when the payload is invalid or
#'   the file is not actually over the limit.
#' @noRd
bam_upload_payload_info <- function(payload, max_bytes) {
    if (!is.list(payload)) {
        return(NULL)
    }
    first <- function(value) {
        value <- unlist(value, use.names = FALSE)
        if (length(value) == 0L) NA else value[[1L]]
    }
    name <- as.character(first(payload$name))
    if (length(name) != 1L || is.na(name) || !nzchar(name)) {
        name <- "upload.bam"
    }
    size <- suppressWarnings(as.numeric(first(payload$size)))
    if (length(size) != 1L || is.na(size) || !is.finite(size) ||
        size <= 0) {
        return(NULL)
    }
    uploads <- data.frame(
        name = name, size = size, stringsAsFactors = FALSE
    )
    error <- tryCatch(
        bam_upload_limit(uploads, max_bytes),
        error = function(e) e
    )
    if (!inherits(error, "condition")) {
        return(NULL)
    }
    bam_error_info(error)
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

#' Score every configured species against one BAM header
#'
#' Used by the dashboard to pick the species that the uploaded BAM actually
#' comes from. A species is an exact match only when every configured contig is
#' present with the expected length; extra header contigs are ignored.
#'
#' @noRd
bam_species_candidate_table <- function(targets, config) {
    species <- config$species
    if (is.null(species) || length(species) == 0L) {
        return(NULL)
    }
    targets <- data.frame(
        name = as.character(targets$name),
        length = as.numeric(targets$length),
        stringsAsFactors = FALSE
    )
    rows <- lapply(names(species), function(key) {
        entry <- species[[key]]
        contigs <- entry$contigs
        if (is.null(contigs) || nrow(contigs) == 0L ||
            !all(c("id", "length") %in% names(contigs))) {
            return(NULL)
        }
        ids <- as.character(contigs$id)
        expected <- as.numeric(contigs$length)
        position <- match(ids, targets$name)
        present <- !is.na(position)
        observed <- rep(NA_real_, length(ids))
        observed[present] <- targets$length[position[present]]
        mismatch <- present & is.finite(expected) & is.finite(observed) &
            observed != expected
        n_targets <- length(ids)
        n_matched <- sum(present & !mismatch)
        label <- entry$label
        if (is.null(label) || length(label) != 1L || is.na(label) ||
            !nzchar(label)) {
            label <- key
        }
        data.frame(
            key = key,
            label = as.character(label),
            n_targets = n_targets,
            n_matched = n_matched,
            n_missing = sum(!present),
            n_length_mismatches = sum(mismatch),
            score = if (n_targets > 0L) n_matched / n_targets else 0,
            exact = n_targets > 0L && !any(!present) && !any(mismatch),
            stringsAsFactors = FALSE
        )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) {
        return(NULL)
    }
    out <- do.call(rbind, rows)
    out[order(-out$score, -out$n_targets, out$key), , drop = FALSE]
}

#' @noRd
bam_species_candidates <- function(path, config) {
    header <- tryCatch(bamHeaderInfo(path), error = function(e) NULL)
    if (is.null(header)) {
        return(NULL)
    }
    bam_species_candidate_table(header$targets, config)
}

#' Detect the most specific exact species match for a BAM header
#' @noRd
bam_detect_species <- function(path, config) {
    candidates <- bam_species_candidates(path, config)
    if (is.null(candidates) || !any(candidates$exact)) {
        return(NULL)
    }
    exact <- candidates[candidates$exact, , drop = FALSE]
    best <- exact[1L, , drop = FALSE]
    if (nrow(exact) > 1L &&
        exact$n_targets[[1L]] == exact$n_targets[[2L]]) {
        return(NULL)
    }
    list(
        species = as.character(best$key[[1L]]),
        label = as.character(best$label[[1L]]),
        candidates = candidates
    )
}

#' Score every ordered species pair against one BAM header
#'
#' The BAM is usually mapped to the combined recipient + donor reference, so a
#' single species rarely matches the whole header. This evaluates the assembly
#' that `resolveSpeciesAssembly()` builds for each ordered pair.
#'
#' @noRd
bam_species_pair_candidate_table <- function(targets, config) {
    species <- config$species
    if (is.null(species) || length(species) < 2L) {
        return(NULL)
    }
    targets <- data.frame(
        name = as.character(targets$name),
        length = as.numeric(targets$length),
        stringsAsFactors = FALSE
    )
    keys <- names(species)
    ## Only pairs whose species contribute at least one header target can be
    ## exact. This keeps the upload-time scan small even for a large registry.
    present <- vapply(keys, function(key) {
        contigs <- species[[key]]$contigs
        is.data.frame(contigs) && "id" %in% names(contigs) &&
            any(as.character(contigs$id) %in% targets$name)
    }, logical(1))
    keys <- keys[present]
    if (length(keys) < 2L) {
        return(NULL)
    }
    rows <- list()
    for (recipient in keys) {
        for (donor in keys) {
            if (identical(recipient, donor)) {
                next
            }
            pair <- tryCatch(
                resolveSpeciesAssembly(
                    config,
                    donor = donor,
                    recipient = recipient,
                    donor_components = bam_species_components(
                        species[[donor]]
                    ),
                    recipient_components = bam_species_components(
                        species[[recipient]]
                    )
                ),
                error = function(e) NULL
            )
            if (is.null(pair) || is.null(pair$assembly)) {
                next
            }
            candidate <- tryCatch(
                .match_candidate(
                    targets, pair$assembly, pair$key, "<header>"
                ),
                error = function(e) NULL
            )
            if (is.null(candidate)) {
                next
            }
            label <- function(key) {
                value <- species[[key]]$label
                if (is.null(value) || length(value) != 1L || is.na(value) ||
                    !nzchar(value)) {
                    key
                } else {
                    as.character(value)
                }
            }
            rows[[length(rows) + 1L]] <- data.frame(
                recipient = recipient,
                donor = donor,
                key = pair$key,
                recipient_label = label(recipient),
                donor_label = label(donor),
                n_targets = candidate$n_targets,
                n_matched = candidate$n_matched,
                n_missing = length(candidate$missing),
                n_length_mismatches = nrow(candidate$length_mismatches),
                score = candidate$score,
                exact = candidate$exact,
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(rows) == 0L) {
        return(NULL)
    }
    out <- do.call(rbind, rows)
    out[order(-out$score, -out$n_targets, out$key), , drop = FALSE]
}

#' @noRd
bam_species_pair_candidates <- function(path, config) {
    header <- tryCatch(bamHeaderInfo(path), error = function(e) NULL)
    if (is.null(header)) {
        return(NULL)
    }
    bam_species_pair_candidate_table(header$targets, config)
}

#' Detect the most specific exact species pair for a BAM header
#'
#' When both orientations of a pair are exact, keep the current donor or the
#' current recipient if that resolves the role ambiguity. This matches the
#' dashboard's default where a newly uploaded BAM keeps the selected donor and
#' changes the recipient.
#'
#' @noRd
bam_choose_species_pair <- function(exact, current_recipient = NULL,
                                    current_donor = NULL) {
    if (!is.data.frame(exact) || nrow(exact) == 0L) {
        return(NULL)
    }
    best <- NULL
    if (!is.null(current_donor) && current_donor %in% exact$donor) {
        best <- exact[exact$donor == current_donor, , drop = FALSE]
    }
    if (is.null(best) && !is.null(current_recipient) &&
        current_recipient %in% exact$recipient) {
        best <- exact[exact$recipient == current_recipient, , drop = FALSE]
    }
    if (is.null(best)) {
        max_targets <- max(exact$n_targets)
        best <- exact[exact$n_targets == max_targets, , drop = FALSE]
    }
    best <- best[order(-best$n_targets, best$key), , drop = FALSE]
    best[1L, , drop = FALSE]
}

bam_detect_species_pair <- function(path, config,
                                    current_recipient = NULL,
                                    current_donor = NULL) {
    candidates <- bam_species_pair_candidates(path, config)
    if (is.null(candidates) || !any(candidates$exact)) {
        return(NULL)
    }
    exact <- candidates[candidates$exact, , drop = FALSE]
    chosen <- bam_choose_species_pair(
        exact,
        current_recipient = current_recipient,
        current_donor = current_donor
    )
    if (is.null(chosen)) {
        return(NULL)
    }
    list(
        recipient = as.character(chosen$recipient[[1L]]),
        donor = as.character(chosen$donor[[1L]]),
        recipient_label = as.character(chosen$recipient_label[[1L]]),
        donor_label = as.character(chosen$donor_label[[1L]]),
        key = as.character(chosen$key[[1L]]),
        candidates = candidates
    )
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

#' Flat default parameters for the dashboard controls
#'
#' A named profile supplies every value. Without a profile the registry
#' defaults are used, so the Custom form still opens with valid controls.
#' @noRd
bam_profile_config_parameters <- function(config) {
    defaults <- config$defaults
    calibration <- resolveCalibrationParameters(
        defaults, mode = "single"
    )
    list(
        windows_per_chromosome = as.integer(
            bam_default_resolution(defaults)
        ),
        metric = as.character(defaults$metric %||% "breadth"),
        score = as.character(defaults$score %||% "paired_log_ratio"),
        donor_aggregation = as.character(
            defaults$donor_aggregation %||% "max"
        ),
        mapq = as.integer(defaults$mapq %||% 20L),
        calibration_fdr = as.numeric(calibration$fdr),
        calibration_floors = as.numeric(calibration$floors),
        calibration_min_mb = calibration$segmentation$min_mb,
        calibration_merge_gap = as.integer(
            calibration$segmentation$merge_gap
        )
    )
}

#' @noRd
bam_profile_default_parameters <- function(config, profile = NULL) {
    spec <- mamamia_analysis_profile(profile)
    if (is.null(spec)) {
        return(bam_profile_config_parameters(config))
    }
    mamamia_profile_flat_parameters(spec$parameters)
}

#' @noRd
bam_profile_resolution_allowed <- function(parameters,
                                           defaults = .config_defaults()) {
    value <- suppressWarnings(as.integer(parameters$windows_per_chromosome))
    if (length(value) != 1L || is.na(value)) {
        return(FALSE)
    }
    range <- bam_resolution_range(defaults)
    value >= range$min && value <= range$max &&
        (value - range$min) %% range$step == 0L
}

#' @noRd
bam_analysis_profile_choices <- function(config) {
    choices <- mamamia_analysis_profile_choices()
    ids <- unname(choices)
    labels <- names(choices)
    keep <- ids == "custom"
    for (i in seq_along(ids)) {
        if (keep[[i]]) {
            next
        }
        profile <- mamamia_analysis_profile(ids[[i]])
        keep[[i]] <- !is.null(profile) &&
            bam_profile_resolution_allowed(profile$parameters, config$defaults)
    }
    if (any(keep & ids != "custom")) {
        choices <- choices[keep]
        names(choices) <- labels[keep]
    } else {
        ## No named profile fits this registry; Advanced is still usable.
        choices <- c("Custom" = "custom")
    }
    choices
}

#' @noRd
bam_default_profile <- function(config) {
    choices <- bam_analysis_profile_choices(config)
    ids <- unname(choices)
    ids <- ids[ids != "custom"]
    if ("balanced" %in% ids) {
        "balanced"
    } else if (length(ids) > 0L) {
        ids[[1L]]
    } else {
        "custom"
    }
}

#' @noRd
bam_profile_has_value <- function(value) {
    if (is.null(value) || length(value) == 0L) {
        return(FALSE)
    }
    !all(is.na(value)) && any(nzchar(as.character(value)))
}

#' @noRd
bam_profile_input_value <- function(input, id, fallback) {
    value <- input[[id]]
    if (!bam_profile_has_value(value) || length(value) != 1L) {
        return(fallback)
    }
    value
}

#' @noRd
bam_profile_whole_number <- function(value, label, minimum, maximum) {
    numeric <- suppressWarnings(as.numeric(value))
    if (length(numeric) != 1L || is.na(numeric) || !is.finite(numeric) ||
        numeric != floor(numeric) ||
        numeric < minimum || numeric > maximum) {
        .mamamia_error(
            "MaMaMIAAnalysisError",
            "%s must be a whole number between %s and %s",
            label,
            format(minimum, scientific = FALSE),
            format(maximum, scientific = FALSE)
        )
    }
    as.integer(numeric)
}

#' @noRd
bam_profile_input_resolution <- function(input, config, resolution = NULL,
                                         fallback = NULL) {
    if (!is.null(resolution) && length(resolution) == 1L &&
        !is.na(resolution)) {
        return(as.integer(resolution))
    }
    for (id in c("resolution_number", "resolution")) {
        value <- input[[id]]
        if (length(value) != 1L || is.na(value)) {
            next
        }
        parsed <- suppressWarnings(as.integer(value))
        if (!is.na(parsed)) {
            return(parsed)
        }
    }
    if (!is.null(fallback) && length(fallback) == 1L &&
        !is.na(fallback)) {
        return(as.integer(fallback))
    }
    as.integer(bam_default_resolution(config$defaults))
}

#' Current Advanced control values as one flat parameter list
#' @noRd
bam_profile_current_parameters <- function(input, config, resolution = NULL,
                                           fallback = NULL) {
    if (is.null(fallback)) {
        fallback <- bam_profile_config_parameters(config)
    }
    n_floors <- bam_calibration_floor_count(config)
    floors <- rep(
        as.numeric(fallback$calibration_floors), length.out = n_floors
    )
    floor_values <- vapply(seq_len(n_floors), function(i) {
        as.numeric(bam_profile_input_value(
            input, paste0("calibration_floor_", i), floors[[i]]
        ))
    }, numeric(1))
    list(
        windows_per_chromosome = bam_profile_input_resolution(
            input, config, resolution,
            fallback = fallback$windows_per_chromosome
        ),
        metric = as.character(bam_profile_input_value(
            input, "metric", fallback$metric
        )),
        score = as.character(bam_profile_input_value(
            input, "score", fallback$score
        )),
        donor_aggregation = as.character(bam_profile_input_value(
            input, "donor_aggregation", fallback$donor_aggregation
        )),
        mapq = bam_profile_whole_number(
            bam_profile_input_value(input, "mapq", fallback$mapq),
            "MAPQ", 0L, 255L
        ),
        calibration_fdr = as.numeric(bam_profile_input_value(
            input, "calibration_fdr", fallback$calibration_fdr
        )),
        calibration_floors = floor_values,
        calibration_min_mb = as.numeric(bam_profile_input_value(
            input, "calibration_min_mb", fallback$calibration_min_mb
        )),
        calibration_merge_gap = bam_profile_whole_number(
            bam_profile_input_value(
                input, "calibration_merge_gap",
                fallback$calibration_merge_gap
            ),
            "Calibration merge gap (windows)", 0L, 10000L
        )
    )
}

#' @noRd
bam_profile_value_changed <- function(current, default) {
    if (is.numeric(current) && is.numeric(default)) {
        return(!isTRUE(all.equal(
            as.numeric(current), as.numeric(default), tolerance = 1e-8,
            check.attributes = FALSE
        )))
    }
    !identical(as.character(current), as.character(default))
}

#' Parameters that differ between two flat profile lists
#' @noRd
bam_profile_changed_parameters <- function(current, defaults) {
    fields <- names(defaults)
    fields[vapply(fields, function(field) {
        bam_profile_value_changed(current[[field]], defaults[[field]])
    }, logical(1))]
}

#' Push every profile default into the Advanced controls
#' @noRd
bam_update_profile_inputs <- function(session, parameters) {
    shiny::updateSliderInput(
        session, "resolution",
        value = as.integer(parameters$windows_per_chromosome)
    )
    shiny::updateNumericInput(
        session, "resolution_number",
        value = as.integer(parameters$windows_per_chromosome)
    )
    shiny::updateSelectInput(
        session, "metric", selected = parameters$metric
    )
    shiny::updateSelectInput(
        session, "score", selected = parameters$score
    )
    shiny::updateSelectInput(
        session, "donor_aggregation",
        selected = parameters$donor_aggregation
    )
    shiny::updateNumericInput(
        session, "mapq", value = as.integer(parameters$mapq)
    )
    shiny::updateNumericInput(
        session, "calibration_fdr",
        value = parameters$calibration_fdr
    )
    floors <- as.numeric(parameters$calibration_floors)
    for (i in seq_along(floors)) {
        shiny::updateNumericInput(
            session, paste0("calibration_floor_", i),
            value = floors[[i]]
        )
    }
    shiny::updateNumericInput(
        session, "calibration_min_mb",
        value = parameters$calibration_min_mb
    )
    shiny::updateNumericInput(
        session, "calibration_merge_gap",
        value = as.integer(parameters$calibration_merge_gap)
    )
    invisible(NULL)
}

#' Short help text for every Advanced parameter
#' @noRd
bam_parameter_help <- function(id, label = NULL) {
    if (grepl("^calibration_floor_", id)) {
        bin <- sub("^Floor ", "", sub(" Mb$", "", label %||% ""))
        return(sprintf(
            paste(
                "Minimum mean score for a call in the %s length bin.",
                "Raise it to reject weak calls in this size range."
            ),
            bin
        ))
    }
    help <- c(
        resolution = paste(
            "Windows per chromosome. More windows give finer breakpoints",
            "but make the run slower."
        ),
        metric = paste(
            "Window metric used by the score. Breadth uses covered bases,",
            "count uses read depth, and combined uses both."
        ),
        score = paste(
            "Paired log-ratio compares donor and recipient coverage.",
            "Donor only ignores the recipient and tolerates copy-number",
            "changes in the recipient."
        ),
        donor_aggregation = paste(
            "How several donor components are combined in one window:",
            "max, mean, or sum."
        ),
        mapq = paste(
            "Minimum mapping quality for a read to count. Higher values",
            "exclude more ambiguous alignments."
        ),
        calibration_fdr = paste(
            "Target false discovery rate for segment thresholds.",
            "Lower values make calls more conservative."
        ),
        calibration_min_mb = paste(
            "Minimum call width in megabases. Shorter segments are",
            "excluded."
        ),
        calibration_merge_gap = paste(
            "Largest gap, in windows, between positive segments that may",
            "be merged into one call."
        ),
        control_samples = paste(
            "Samples used only to calibrate thresholds. Leave empty to",
            "pool every sample in the cohort."
        )
    )
    help[[id]] %||% ""
}

#' @noRd
bam_profile_parameter_label <- function(id, label, changed, help = NULL) {
    shiny::tags$label(
        class = c(
            "control-label",
            if (isTRUE(changed)) {
                "mamamia-param-changed"
            } else {
                "mamamia-param-default"
            }
        ),
        `for` = id,
        shiny::span(class = "mamamia-param-label-text", label),
        if (!is.null(help) && length(help) == 1L && nzchar(help)) {
            bslib::tooltip(
                shiny::span(
                    class = "mamamia-param-help",
                    shiny::icon("circle-question")
                ),
                help,
                placement = "right"
            )
        }
    )
}

#' @noRd
bam_number <- function(value, label, minimum, maximum) {
    if (length(value) != 1L || is.na(value) || !is.numeric(value) ||
        !is.finite(value) || value < minimum || value > maximum) {
        .mamamia_error(
            "MaMaMIAAnalysisError",
            "%s must be a number between %s and %s",
            label,
            format(minimum, scientific = FALSE),
            format(maximum, scientific = FALSE)
        )
    }
    as.numeric(value)
}

#' Validate and resolve the run controls
#'
#' Advanced controls start from the selected profile and may override every
#' profile parameter. Any override makes the effective profile Custom.
#'
#' @param input The Shiny input list (or a list with the same fields).
#' @param config A validated configuration.
#' @param samples Sample names from preflight.
#' @param resolution Optional resolved resolution; overrides `input$resolution`.
#' @param profile_base Named profile the Custom values were based on.
#'
#' @return A list with the resolved `analyzeBamInput()` arguments.
#' @noRd
bam_run_parameters <- function(input, config, samples, resolution = NULL,
                               profile_base = NULL) {
    if (!is.character(samples) || length(samples) == 0L) {
        .mamamia_error("MaMaMIAAnalysisError", "At least one sample is required")
    }
    defaults <- config$defaults
    requested <- input$analysis_profile
    if (is.null(requested) || length(requested) != 1L ||
        is.na(requested) || !nzchar(requested)) {
        requested <- "custom"
    }
    requested <- as.character(requested)
    profile <- mamamia_analysis_profile(requested)
    if (!is.null(profile) &&
        !bam_profile_resolution_allowed(profile$parameters, defaults)) {
        range <- bam_resolution_range(defaults)
        .mamamia_error(
            "MaMaMIAProfileError",
            paste0(
                "Analysis profile '%s' uses %s windows per chromosome; ",
                "this configuration allows %s-%s in steps of %s"
            ),
            profile$label,
            as.integer(profile$parameters$windows_per_chromosome),
            format(range$min, big.mark = ",", scientific = FALSE),
            format(range$max, big.mark = ",", scientific = FALSE),
            range$step
        )
    }
    profile_name <- if (is.null(profile)) "custom" else requested
    template <- if (!is.null(profile)) {
        profile_name
    } else if (!is.null(profile_base) && length(profile_base) == 1L &&
        !is.na(profile_base) && nzchar(profile_base) &&
        !identical(profile_base, "custom")) {
        as.character(profile_base)
    } else {
        NULL
    }
    fallback <- bam_profile_default_parameters(config, template)
    current <- bam_profile_current_parameters(
        input, config, resolution = resolution, fallback = fallback
    )
    resolution <- bam_resolution_value(
        current$windows_per_chromosome, defaults
    )
    mapq <- bam_integer(current$mapq, "MAPQ", 0L, 255L)
    exclude_flags <- bam_integer(
        defaults$exclude_flags, "Exclude flags", 0L, 4095L
    )
    available_cores <- min(256L, bam_max_cores())
    cores <- bam_integer(
        mamamia_default_cores(available_cores),
        "Cores", 1L, available_cores
    )
    donor_aggregation <- as.character(current$donor_aggregation)
    if (!donor_aggregation %in% c("mean", "max", "sum")) {
        .mamamia_error(
            "MaMaMIAAnalysisError",
            "`donor_aggregation` must be one of: mean, max, sum"
        )
    }
    metric <- as.character(current$metric)
    if (!metric %in% c("breadth", "count", "combined")) {
        .mamamia_error(
            "MaMaMIAAnalysisError",
            "`metric` must be one of: breadth, count, combined"
        )
    }
    score <- as.character(current$score)
    if (!score %in% c("paired_log_ratio", "donor_only")) {
        .mamamia_error(
            "MaMaMIAAnalysisError",
            "`score` must be one of: paired_log_ratio, donor_only"
        )
    }
    n_floors <- bam_calibration_floor_count(config)
    calibration_ids <- c(
        "calibration_fdr", paste0("calibration_floor_", seq_len(n_floors)),
        "calibration_min_mb", "calibration_merge_gap"
    )
    has_calibration <- !is.null(profile) ||
        any(vapply(calibration_ids, function(id) {
            bam_profile_has_value(input[[id]])
        }, logical(1)))
    calibration <- if (isTRUE(has_calibration)) {
        floors <- as.numeric(current$calibration_floors)
        if (length(floors) != n_floors || anyNA(floors) || any(floors < 0)) {
            .mamamia_error(
                "MaMaMIAAnalysisError",
                "Each calibration floor must be a non-negative number"
            )
        }
        list(
            fdr = bam_number(current$calibration_fdr, "Calibration FDR", 1e-6, 0.999999),
            floors = floors,
            min_mb = bam_number(
                current$calibration_min_mb,
                "Calibration minimum length (Mb)", 0, Inf
            ),
            merge_gap = bam_integer(
                current$calibration_merge_gap,
                "Calibration merge gap (windows)", 0L, 10000L
            )
        )
    } else {
        list()
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
    if (length(control_samples) > 0L) {
        calibration$controls <- control_samples
    }
    if (length(calibration) == 0L) {
        calibration <- NULL
    }
    changed <- bam_profile_changed_parameters(current, fallback)
    effective <- if (is.null(profile) || length(changed) > 0L) {
        "custom"
    } else {
        profile_name
    }
    base <- if (!is.null(profile)) {
        profile_name
    } else if (!is.null(profile_base) && length(profile_base) == 1L &&
        !is.na(profile_base) && nzchar(profile_base) &&
        !identical(profile_base, "custom")) {
        as.character(profile_base)
    } else {
        "custom"
    }
    list(
        resolution = resolution,
        mapq = mapq,
        exclude_flags = exclude_flags,
        cores = cores,
        metric = metric,
        score = score,
        donor_aggregation = donor_aggregation,
        control_samples = control_samples,
        calibration = calibration,
        profile_name = effective,
        profile_base = base,
        profile_parameters = current,
        profile_changed = changed
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
    mamamia_natural_sort(unique(as.character(result$samples$sample)))
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
#' `set()` evicts the oldest entry once the current limit is reached. The
#' limit stays fixed for the session so old analyses cannot accumulate.
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
                             controls = NULL, profile = NULL,
                             parameters = NULL) {
    if (!isTRUE(preflight$ok) || is.null(preflight$table) ||
        !all(c("path", "sample") %in% names(preflight$table))) {
        return(NULL)
    }
    defaults <- preflight$config$defaults
    if (is.null(defaults)) {
        defaults <- list()
    }
    profile_name <- profile
    if (is.null(profile_name) || length(profile_name) != 1L ||
        is.na(profile_name) || !nzchar(profile_name)) {
        profile_name <- "custom"
    }
    profile_name <- as.character(profile_name)
    profile_spec <- mamamia_analysis_profile(profile_name)
    if (is.null(profile_spec)) {
        profile_name <- "custom"
    }
    samples <- unique(as.character(preflight$table$sample))
    controls <- unique(as.character(controls))
    controls <- controls[nzchar(controls)]
    controls <- sort(intersect(controls, samples))
    if (length(controls) > 0L && setequal(controls, samples)) {
        controls <- character(0)
    }
    if (isTRUE(preflight$example)) {
        profile_name <- "balanced"
        profile_spec <- mamamia_analysis_profile("balanced")
        parameters <- NULL
        controls <- character(0)
    }
    if (!is.null(parameters)) {
        resolution_value <- parameters$windows_per_chromosome
        metric <- parameters$metric
        score <- parameters$score
        donor_aggregation <- parameters$donor_aggregation
        mapq <- parameters$mapq
        calibration <- list(
            fdr = parameters$calibration_fdr,
            floors = parameters$calibration_floors,
            min_mb = parameters$calibration_min_mb,
            merge_gap = parameters$calibration_merge_gap
        )
    } else if (!is.null(profile_spec)) {
        resolution_value <- profile_spec$parameters$windows_per_chromosome
        metric <- profile_spec$parameters$metric
        score <- profile_spec$parameters$score
        donor_aggregation <- profile_spec$parameters$donor_aggregation
        mapq <- profile_spec$parameters$mapq
        calibration <- profile_spec$parameters$calibration
    } else {
        resolution_value <- resolution
        metric <- defaults$metric
        score <- defaults$score
        mapq <- defaults$mapq
        calibration <- list()
    }
    if (is.null(donor_aggregation) || length(donor_aggregation) != 1L ||
        is.na(donor_aggregation) || !nzchar(donor_aggregation)) {
        donor_aggregation <- defaults$donor_aggregation
    }
    if (!is.null(calibration$floors)) {
        calibration$floors <- .calibration_stretch_floors(
            calibration$floors,
            bam_calibration_floor_count(preflight$config)
        )
    }
    if (length(controls) > 0L) {
        calibration$controls <- controls
    }
    if (length(calibration) == 0L) {
        calibration <- NULL
    }
    resolution_key <- if (isTRUE(preflight$example)) {
        as.character(bam_example_default_resolution())
    } else {
        tryCatch(
            as.character(bam_resolution_value(resolution_value, defaults)),
            error = function(e) as.character(resolution_value)
        )
    }
    hash_object(list(
        example = isTRUE(preflight$example),
        paths = unique(as.character(preflight$table$path)),
        samples = samples,
        assembly = as.character(preflight$assembly),
        profile = profile_name,
        resolution = resolution_key,
        mapq = as.integer(mapq),
        exclude_flags = as.integer(defaults$exclude_flags),
        metric = as.character(metric),
        score = as.character(score),
        donor_aggregation = as.character(donor_aggregation),
        calibration = calibration,
        controls = controls
    ))
}

#' @noRd
## UI ------------------------------------------------------------------------

#' @noRd
## Dashboard layout ----------------------------------------------------------


#' Single-file BAM chooser
#'
#' The dashboard runs one BAM per analysis, so the control never offers a
#' batch. `multiple = FALSE` removes multi-select from the native file picker;
#' the launcher's drop guard trims a multi-file drag-and-drop to its first file
#' and the server rejects a multi-file payload. A multi-sample BAM is a single
#' file and remains supported.
#'
#' @param max_bytes Optional maximum upload size in bytes. When supplied, the
#'   upload zone records it as `data-max-bytes` so the client can reject an
#'   oversized file before Shiny starts the upload.
#'
#' @return A `shiny.tag` containing the labelled file input.
#' @noRd
bam_single_file_input <- function(max_bytes = NULL) {
    limit <- if (is.numeric(max_bytes) && length(max_bytes) == 1L &&
        is.finite(max_bytes) && max_bytes > 0) {
        sprintf("%.0f", max_bytes)
    } else {
        NULL
    }
    shiny::div(
        class = "mamamia-upload",
        `data-max-bytes` = limit,
        shiny::div(
            class = "mamamia-upload-head",
            shiny::span(
                class = "mamamia-upload-icon",
                `aria-hidden` = "true",
                shiny::icon("dna")
            ),
            shiny::div(
                class = "mamamia-upload-copy",
                shiny::div(
                    class = "mamamia-upload-title-row",
                    shiny::span(class = "mamamia-upload-title", "BAM file"),
                    if (!is.null(limit)) {
                        shiny::span(
                            class = "mamamia-upload-limit",
                            sprintf("Max %s", bam_format_size(max_bytes)),
                            shiny::actionButton(
                                "bam_upload_help",
                                label = shiny::icon("circle-question"),
                                class = "mamamia-upload-help",
                                title = "About the BAM upload limit",
                                `aria-label` = "About the BAM upload limit"
                            )
                        )
                    }
                ),
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
        ),
        shiny::div(
            class = "mamamia-upload-timer",
            `aria-live` = "off",
            shiny::span(
                class = "mamamia-upload-timer-label", "Time to finish"
            ),
            shiny::span(class = "mamamia-upload-timer-value", "--:--")
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
    too_large <- identical(stage, "upload")
    message <- text_value(info$message)
    suggestion <- text_value(info$suggestion)
    recovery_copy <- if (too_large) {
        list(
            body = paste(
                "This BAM is larger than this server's upload limit,",
                "so it was not analysed."
            ),
            note = NULL
        )
    } else if (assembly) {
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
                shiny::icon(
                    if (too_large) {
                        "triangle-exclamation"
                    } else if (assembly) {
                        "dna"
                    } else {
                        "sliders"
                    }
                )
            ),
            shiny::span(
                class = "mamamia-recovery-title",
                if (too_large) {
                    "BAM is too large"
                } else if (assembly) {
                    "No matching assembly"
                } else {
                    "Setup needs attention"
                }
            )
        ),
        if (nzchar(recovery_copy$body)) {
            shiny::p(class = "mamamia-recovery-message", recovery_copy$body)
        },
        if (too_large && !is.null(info$size) && !is.null(info$limit)) {
            shiny::p(
                class = "mamamia-recovery-message",
                sprintf(
                    "This file is %s; the limit is %s.",
                    bam_format_size(info$size),
                    bam_format_size(info$limit)
                )
            )
        },
        if (!is.null(recovery_copy$note) && nzchar(recovery_copy$note)) {
            shiny::p(class = "mamamia-recovery-note", recovery_copy$note)
        },
        if (too_large) {
            shiny::tagList(
                shiny::p(
                    class = "mamamia-recovery-note",
                    "Slim it locally with ",
                    shiny::tags$a(
                        href = .bamhusk_url,
                        target = "_blank",
                        rel = "noopener noreferrer",
                        "bamhusk"
                    ),
                    ", then upload the smaller BAM."
                ),
                shiny::tags$pre(
                    class = "mamamia-recovery-command",
                    shiny::tags$code(.bamhusk_command)
                )
            )
        },
        if (nzchar(suggestion) && !too_large) {
            shiny::p(
                class = "mamamia-recovery-suggestion",
                shiny::strong("Next step: "),
                suggestion
            )
        },
        if (too_large) {
            shiny::tags$a(
                class = "mamamia-recovery-action",
                href = .bamhusk_url,
                target = "_blank",
                rel = "noopener noreferrer",
                shiny::icon("arrow-up-right-from-square"),
                " Open bamhusk"
            )
        } else if (isTRUE(show_action)) {
            shiny::actionButton(
                "unlock_setup",
                if (assembly) "Choose assembly" else "Review setup",
                icon = shiny::icon("arrow-right"),
                class = "mamamia-recovery-action"
            )
        }
    )
}

#' Modal shown once when an oversized upload is rejected
#' @noRd
bam_upload_too_large_modal <- function(info) {
    if (!is.list(info)) {
        info <- list(message = "The BAM is too large for this server.")
    }
    detail <- if (!is.null(info$size) && !is.null(info$limit)) {
        sprintf(
            "This BAM is %s; this server accepts at most %s.",
            bam_format_size(info$size),
            bam_format_size(info$limit)
        )
    } else {
        as.character(info$message)[[1L]]
    }
    shiny::modalDialog(
        title = "BAM is too large to upload",
        easyClose = TRUE,
        shiny::p(detail),
        if (!is.null(info$files) && length(info$files) > 0L) {
            shiny::p(
                class = "text-muted small",
                basename(as.character(info$files[[1L]]))
            )
        },
        shiny::p(
            "Slim it locally with ",
            shiny::tags$a(
                href = .bamhusk_url,
                target = "_blank",
                rel = "noopener noreferrer",
                "bamhusk"
            ),
            ", then upload the smaller BAM."
        ),
        shiny::tags$pre(
            class = "mamamia-recovery-command",
            shiny::tags$code(.bamhusk_command)
        ),
        shiny::p(
            class = "text-muted small mb-0",
            paste(
                "The upload was rejected before staging or analysis;",
                "no BAM data was processed."
            )
        ),
        footer = shiny::tagList(
            shiny::modalButton("Close"),
            shiny::actionButton(
                "clear_oversized_upload",
                "Choose another file",
                class = "btn-primary"
            )
        )
    )
}

#' Upload-limit help shown by the `?` next to the size hint
#' @noRd
bam_upload_help_ui <- function(max_bytes) {
    limit <- bam_format_size(max_bytes)
    section <- function(title, ...) {
        shiny::div(
            class = "mamamia-help-section",
            shiny::tags$h6(title),
            ...
        )
    }
    shiny::div(
        class = "mamamia-tooltip-help",
        shiny::p(
            class = "mamamia-help-lead",
            "This server accepts BAM files up to ",
            shiny::strong(limit),
            ". The size is checked in your browser before upload."
        ),
        shiny::div(
            class = "mamamia-help-facts",
            shiny::div(
                shiny::span("Limit"),
                shiny::strong(limit)
            ),
            shiny::div(
                shiny::span("Check"),
                shiny::strong("Before upload")
            ),
            shiny::div(
                shiny::span("Typical saving"),
                shiny::strong("50–90%")
            )
        ),
        shiny::div(
            class = "mamamia-help-callout",
            paste(
                "If the BAM is larger, slim it locally with bamhusk.",
                "The command writes a new file and leaves the original",
                "untouched."
            )
        ),
        section(
            "Slim it locally",
            shiny::tags$ol(
                class = "mamamia-help-steps",
                shiny::tags$li(
                    shiny::strong("Get bamhusk"),
                    " from ",
                    shiny::tags$a(
                        href = .bamhusk_url,
                        target = "_blank",
                        rel = "noopener noreferrer",
                        "github.com/alermol/bamhusk"
                    ),
                    "."
                ),
                shiny::tags$li(
                    shiny::strong("Run"),
                    " the command below on the large BAM."
                ),
                shiny::tags$li(
                    shiny::strong("Upload"),
                    " the new, smaller BAM here."
                )
            ),
            shiny::tags$code(
                class = "mamamia-help-command",
                tabindex = "0",
                `aria-label` = "Command to slim the BAM",
                .bamhusk_command
            )
        ),
        section(
            "What changes",
            shiny::tags$dl(
                class = "mamamia-help-list",
                shiny::tags$dt("Kept"),
                shiny::tags$dd(
                    paste(
                        "Everything MaMaMIA reads: FLAG, RNAME, POS,",
                        "MAPQ, CIGAR, @HD/@SQ/@RG and RG tags."
                    )
                ),
                shiny::tags$dt("Removed"),
                shiny::tags$dd(
                    paste(
                        "QUAL, read names, optional tags except RG,",
                        "and @PG/@CO header lines."
                    )
                ),
                shiny::tags$dt("Optional"),
                shiny::tags$dd(
                    paste(
                        "SEQ with --drop-seq. Read length still comes",
                        "from CIGAR, so the analysis is unchanged."
                    )
                )
            )
        )
    )
}

#' Modal for the `?` next to the upload size hint
#' @noRd
bam_upload_help_modal <- function(max_bytes) {
    dialog <- shiny::modalDialog(
        title = "BAM upload limit",
        bam_upload_help_ui(max_bytes),
        size = "l",
        easyClose = TRUE,
        footer = shiny::tagList(
            shiny::tags$a(
                href = .bamhusk_url,
                target = "_blank",
                rel = "noopener noreferrer",
                class = "btn btn-outline-secondary",
                shiny::icon("arrow-up-right-from-square"),
                " Open bamhusk"
            ),
            shiny::modalButton("Close")
        )
    )
    htmltools::tagAppendAttributes(
        dialog,
        class = "mamamia-tooltip-modal modal-dialog-centered"
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
            multiple = TRUE,
            open = "Analysis",
            bslib::accordion_panel(
                "Analysis",
                shiny::uiOutput("analysis_profile_ui")
            ),
            bslib::accordion_panel(
                "Assembly",
                shiny::uiOutput("species_ui")
            ),
            bslib::accordion_panel(
                "Advanced",
                shiny::uiOutput("advanced_parameters_ui"),
                shiny::uiOutput("control_samples_ui")
            )
        ),
        shiny::div(
            class = "mamamia-run-panel",
            shiny::uiOutput("run_ui")
        )
    )
}

#' Wrapper for one Advanced parameter so its label can be restyled alone
#' @noRd
bam_parameter_control_ui <- function(id, input, extra_class = NULL) {
    shiny::div(
        class = c("mamamia-param-control", extra_class),
        shiny::uiOutput(
            paste0("param_label_", id),
            inline = TRUE
        ),
        input
    )
}

#' @noRd
bam_calibration_floor_labels <- function(config) {
    labels <- resolveCalibrationParameters(config$defaults)$labels
    paste0("Floor ", labels, " Mb")
}

#' Number of floor controls declared by the registry bins
#' @noRd
bam_calibration_floor_count <- function(config) {
    length(bam_calibration_floor_labels(config))
}

#' Warning shown when cohort recurrence cannot be computed
#' @noRd
bam_cohort_mask_warning_ui <- function(n_samples) {
    n_samples <- suppressWarnings(as.integer(n_samples))
    if (length(n_samples) != 1L || is.na(n_samples) ||
        n_samples != 1L) {
        return(NULL)
    }
    shiny::div(
        class = "alert alert-warning py-2 px-3 small mb-2",
        role = "status",
        shiny::strong("Cohort recurrence unavailable."),
        paste(
            "This run has one sample. Cohort recurrence masking requires",
            "at least two samples, so no windows will be masked."
        )
    )
}

#' Every profile parameter is editable under Advanced
#' @noRd
bam_advanced_parameters_ui <- function(config, defaults, locked = FALSE) {
    if (isTRUE(locked)) {
        return(shiny::tagList(
            shiny::p(
                class = "mamamia-fixed-note",
                paste(
                    "Analysis parameters are fixed for the bundled example.",
                    "Upload a BAM to edit every profile parameter."
                )
            )
        ))
    }
    n_floors <- bam_calibration_floor_count(config)
    floor_values <- rep(
        as.numeric(defaults$calibration_floors), length.out = n_floors
    )
    floor_controls <- lapply(seq_len(n_floors), function(i) {
        id <- paste0("calibration_floor_", i)
        bam_parameter_control_ui(
            id,
            shiny::numericInput(
                id, NULL,
                value = floor_values[[i]],
                min = 0,
                step = 0.1,
                width = "100%"
            ),
            extra_class = "mamamia-param-floor"
        )
    })
    shiny::tagList(
        shiny::div(
            class = "mamamia-param-actions",
            shiny::actionButton(
                "reset_profile_parameters",
                "Reset all to defaults",
                icon = shiny::icon("rotate-left"),
                class = "btn-sm btn-outline-secondary",
                title = paste(
                    "Restore every profile parameter and control sample",
                    "to the selected profile's defaults"
                )
            )
        ),
        shiny::div(class = "mamamia-param-heading", "Resolution"),
        shiny::uiOutput("resolution_ui"),
        shiny::div(class = "mamamia-param-heading", "Scoring"),
        bam_parameter_control_ui(
            "metric",
            shiny::selectInput(
                "metric", NULL,
                choices = c(
                    "Breadth" = "breadth",
                    "Count" = "count",
                    "Combined" = "combined"
                ),
                selected = defaults$metric
            )
        ),
        bam_parameter_control_ui(
            "score",
            shiny::selectInput(
                "score", NULL,
                choices = c(
                    "Paired log-ratio" = "paired_log_ratio",
                    "Donor only" = "donor_only"
                ),
                selected = defaults$score
            )
        ),
        bam_parameter_control_ui(
            "donor_aggregation",
            shiny::selectInput(
                "donor_aggregation", NULL,
                choices = c("Max" = "max", "Mean" = "mean", "Sum" = "sum"),
                selected = defaults$donor_aggregation
            )
        ),
        bam_parameter_control_ui(
            "mapq",
            shiny::numericInput(
                "mapq", NULL,
                value = defaults$mapq,
                min = 0, max = 255, step = 1,
                width = "100%"
            )
        ),
        shiny::div(class = "mamamia-param-heading", "Masking"),
        shiny::p(
            class = "text-muted small mb-2",
            paste(
                "Cohort recurrence cutoff is derived automatically from",
                "the fitted cohort mixture. Single-sample runs are unmasked."
            )
        ),
        shiny::div(class = "mamamia-param-heading", "Calibration"),
        bam_parameter_control_ui(
            "calibration_fdr",
            shiny::numericInput(
                "calibration_fdr", NULL,
                value = defaults$calibration_fdr,
                min = 0, max = 1, step = 0.001,
                width = "100%"
            )
        ),
        shiny::div(
            class = "mamamia-param-floors",
            floor_controls
        ),
        bam_parameter_control_ui(
            "calibration_min_mb",
            shiny::numericInput(
                "calibration_min_mb", NULL,
                value = defaults$calibration_min_mb,
                min = 0, step = 0.5,
                width = "100%"
            )
        ),
        bam_parameter_control_ui(
            "calibration_merge_gap",
            shiny::numericInput(
                "calibration_merge_gap", NULL,
                value = defaults$calibration_merge_gap,
                min = 0, step = 10,
                width = "100%"
            )
        ),
        shiny::p(
            class = "mamamia-param-legend",
            shiny::strong("Bold"),
            " values are the selected profile's defaults;",
            shiny::span(class = "mamamia-param-changed", "highlighted"),
            " values were changed and make the run Custom."
        )
    )
}

## The result area. The view filters stay visible above all tabs;
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
                shiny::p(
                    class = "mamamia-call-note",
                    paste(
                        "Brief view of passing calls. Hover a column heading",
                        "for its description; use Additional columns to add",
                        "diagnostics. Shift-click headings to sort by several",
                        "columns. Turn off Passing only to include calls that",
                        "failed calibration. The full table, BEDPE files and",
                        "images are in the results bundle."
                    )
                ),
                DT::DTOutput(ns("calls"))
            ),
            bslib::nav_panel(
                "Export",
                shiny::div(
                    class = "mamamia-export-panel",
                    explorer_export_ui("explorer"),
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
                    "var mamamiaAnalysisStartedAt = null;",
                    "var mamamiaAnalysisTimerId = null;",
                    "function mamamiaFormatElapsed(ms) {",
                    "  var total = Math.max(0, Math.floor(ms / 1000));",
                    "  var hours = Math.floor(total / 3600);",
                    "  var minutes = Math.floor((total % 3600) / 60);",
                    "  var seconds = total % 60;",
                    "  function pad(value) {",
                    "    return (value < 10 ? '0' : '') + value;",
                    "  }",
                    "  if (hours > 0) {",
                    "    return hours + ':' + pad(minutes) + ':' + pad(seconds);",
                    "  }",
                    "  return pad(minutes) + ':' + pad(seconds);",
                    "}",
                    "function mamamiaWriteAnalysisTimer(text) {",
                    "  var value = document.querySelectorAll(",
                    "    '.mamamia-analysis-timer-value');",
                    "  Array.prototype.forEach.call(value, function(node) {",
                    "    node.textContent = text;",
                    "  });",
                    "}",
                    "function mamamiaTickAnalysisTimer() {",
                    "  if (mamamiaAnalysisStartedAt === null) return;",
                    "  mamamiaWriteAnalysisTimer(",
                    "    mamamiaFormatElapsed(Date.now() - mamamiaAnalysisStartedAt));",
                    "}",
                    "function mamamiaStartAnalysisTimer() {",
                    "  if (mamamiaAnalysisTimerId !== null) return;",
                    "  mamamiaAnalysisStartedAt = Date.now();",
                    "  mamamiaWriteAnalysisTimer('00:00');",
                    "  mamamiaAnalysisTimerId = window.setInterval(",
                    "    mamamiaTickAnalysisTimer, 250);",
                    "}",
                    "function mamamiaStopAnalysisTimer() {",
                    "  if (mamamiaAnalysisTimerId !== null) {",
                    "    window.clearInterval(mamamiaAnalysisTimerId);",
                    "    mamamiaAnalysisTimerId = null;",
                    "  }",
                    "  mamamiaAnalysisStartedAt = null;",
                    "}",
                    "function mamamiaUploadTimerNode() {",
                    "  var input = document.getElementById('bams');",
                    "  if (!input) return null;",
                    "  var zone = input.closest('.mamamia-upload');",
                    "  return zone ? zone.querySelector('.mamamia-upload-timer') : null;",
                    "}",
                    "function mamamiaWriteUploadTimer(node, text) {",
                    "  if (!node) return;",
                    "  var value = node.querySelector('.mamamia-upload-timer-value');",
                    "  if (value) value.textContent = text;",
                    "}",
                    "function mamamiaUploadProgress() {",
                    "  var bar = document.querySelector('#bams_progress .progress-bar');",
                    "  if (!bar || !bar.style) return null;",
                    "  var width = parseFloat(bar.style.width);",
                    "  if (!isFinite(width)) return null;",
                    "  return Math.max(0, Math.min(1, width / 100));",
                    "}",
                    "function mamamiaUploadRemaining(node) {",
                    "  if (!node || node.mamamiaUploadStartedAt === null ||",
                    "      node.mamamiaUploadStartedAt === undefined) return null;",
                    "  var fraction = mamamiaUploadProgress();",
                    "  if (fraction === null) return null;",
                    "  if (fraction >= 1) return 0;",
                    "  /* Below one percent the estimate is pure noise. */",
                    "  if (fraction < 0.01) return null;",
                    "  var elapsed = Date.now() - node.mamamiaUploadStartedAt;",
                    "  if (elapsed <= 0) return null;",
                    "  return Math.max(0, elapsed * (1 - fraction) / fraction);",
                    "}",
                    "function mamamiaStopUploadTimer(node, failed) {",
                    "  if (node) {",
                    "    if (node.mamamiaUploadTimerId) {",
                    "      window.clearInterval(node.mamamiaUploadTimerId);",
                    "      node.mamamiaUploadTimerId = null;",
                    "    }",
                    "    if (failed) {",
                    "      mamamiaWriteUploadTimer(node, '--:--');",
                    "    }",
                    "    node.mamamiaUploadStartedAt = null;",
                    "    node.classList.remove('is-running');",
                    "    node.classList.toggle('is-error', !!failed);",
                    "  }",
                    "}",
                    "function mamamiaStartUploadTimer(node) {",
                    "  if (!node) return;",
                    "  mamamiaStopUploadTimer(node, false);",
                    "  node.mamamiaUploadStartedAt = Date.now();",
                    "  node.classList.add('is-visible', 'is-running');",
                    "  node.classList.remove('is-error');",
                    "  mamamiaWriteUploadTimer(node, '--:--');",
                    "  node.mamamiaUploadTimerId = window.setInterval(function() {",
                    "    if (node.mamamiaUploadStartedAt === null ||",
                    "        node.mamamiaUploadStartedAt === undefined) return;",
                    "    var bar = document.querySelector(",
                    "      '#bams_progress .progress-bar');",
                    "    if (bar && bar.classList.contains('progress-bar-danger')) {",
                    "      mamamiaStopUploadTimer(node, true);",
                    "      return;",
                    "    }",
                    "    var remaining = mamamiaUploadRemaining(node);",
                    "    if (remaining === null) {",
                    "      mamamiaWriteUploadTimer(node, '--:--');",
                    "      return;",
                    "    }",
                    "    mamamiaWriteUploadTimer(node, mamamiaFormatElapsed(remaining));",
                    "    if (remaining <= 0) {",
                    "      mamamiaStopUploadTimer(node, false);",
                    "    }",
                    "  }, 250);",
                    "}",
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
                    "  if (x.running) {",
                    "    mamamiaStartAnalysisTimer();",
                    "  } else {",
                    "    mamamiaStopAnalysisTimer();",
                    "  }",
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
                    "  function uploadLimit(zone) {",
                    "    if (!zone) return null;",
                    "    var raw = zone.getAttribute('data-max-bytes');",
                    "    if (!raw) return null;",
                    "    var bytes = Number(raw);",
                    "    return isFinite(bytes) && bytes > 0 ? bytes : null;",
                    "  }",
                    "  function clearInput(input) {",
                    "    try {",
                    "      input.files = new DataTransfer().files;",
                    "    } catch (error) {",
                    "      /* Fall back to clearing the value below. */",
                    "    }",
                    "    input.value = '';",
                    "  }",
                    "  function rejectOversized(input, file, zone, event) {",
                    "    var limit = uploadLimit(zone);",
                    "    if (!limit || !file || file.size <= limit) return false;",
                    "    if (event) {",
                    "      event.preventDefault();",
                    "      event.stopImmediatePropagation();",
                    "    }",
                    "    clearInput(input);",
                    "    if (window.Shiny && window.Shiny.setInputValue) {",
                    "      window.Shiny.setInputValue('oversized_upload', {",
                    "        name: file.name || 'upload.bam',",
                    "        size: file.size,",
                    "        limit: limit",
                    "      }, { priority: 'event' });",
                    "    }",
                    "    return true;",
                    "  }",
                    "  document.addEventListener('change', function(event) {",
                    "    var input = event.target;",
                    "    if (!input || input.id !== id || !input.files) return;",
                    "    if (input.files.length > 1) {",
                    "      input.files = onlyFirst(input.files);",
                    "    }",
                    "    if (rejectOversized(",
                    "      input, input.files[0], uploadZone(input), event",
                    "    )) return;",
                    "    if (input.files[0]) {",
                    "      mamamiaStartUploadTimer(mamamiaUploadTimerNode());",
                    "    }",
                    "  }, true);",
                    "  document.addEventListener('shiny:inputchanged', function(event) {",
                    "    var name = event.detail && event.detail.name;",
                    "    var type = event.detail && event.detail.inputType;",
                    "    if (name !== id || type !== 'shiny.fileupload') return;",
                    "    var node = mamamiaUploadTimerNode();",
                    "    if (node) {",
                    "      mamamiaWriteUploadTimer(node, '00:00');",
                    "      mamamiaStopUploadTimer(node, false);",
                    "    }",
                    "  });",
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
                    "      'input[type=\"file\"], .btn-file, label, button, a');",
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
                    "    if (rejectOversized(",
                    "      input, files[0], uploadZone(input), event",
                    "    )) return;",
                    "    input.files = onlyFirst(files);",
                    "    input.dispatchEvent(new Event('change', { bubbles: true }));",
                    "  }, true);",
                    "})();"
                )
            )),
            shiny::tags$script(shiny::HTML(
                paste(
                    "(function() {",
                    "  function filenameFromHeader(header, fallback) {",
                    "    if (!header) return fallback;",
                    "    var match = /filename\\*?=(?:UTF-8'')?[\"']?([^\"';]+)[\"']?/i.exec(header);",
                    "    if (!match) return fallback;",
                    "    try { return decodeURIComponent(match[1]); }",
                    "    catch (error) { return match[1]; }",
                    "  }",
                    "  function statusNode(row) {",
                    "    return row.querySelector('.mamamia-download-status');",
                    "  }",
                    "  function setBusy(row) {",
                    "    var status = statusNode(row);",
                    "    row.classList.add('is-busy');",
                    "    row.setAttribute('aria-busy', 'true');",
                    "    if (!status) return;",
                    "    status.classList.remove('is-error');",
                    "    status.classList.add('is-visible');",
                    "    status.setAttribute('aria-hidden', 'false');",
                    "    var text = status.querySelector('.mamamia-download-status-text');",
                    "    if (text) {",
                    "      text.textContent = row.getAttribute('data-busy-text') ||",
                    "        'Preparing download...';",
                    "    }",
                    "  }",
                    "  function setIdle(row) {",
                    "    var status = statusNode(row);",
                    "    row.classList.remove('is-busy');",
                    "    row.removeAttribute('aria-busy');",
                    "    if (!status) return;",
                    "    status.classList.remove('is-visible', 'is-error');",
                    "    status.setAttribute('aria-hidden', 'true');",
                    "  }",
                    "  function setError(row) {",
                    "    var status = statusNode(row);",
                    "    row.classList.remove('is-busy');",
                    "    row.removeAttribute('aria-busy');",
                    "    if (!status) return;",
                    "    status.classList.add('is-visible', 'is-error');",
                    "    status.setAttribute('aria-hidden', 'false');",
                    "    var text = status.querySelector('.mamamia-download-status-text');",
                    "    if (text) text.textContent = 'Download failed. Try again.';",
                    "  }",
                    "  function saveBlob(blob, filename) {",
                    "    var url = URL.createObjectURL(blob);",
                    "    var link = document.createElement('a');",
                    "    link.href = url;",
                    "    link.download = filename || 'mamamia-download';",
                    "    link.style.display = 'none';",
                    "    document.body.appendChild(link);",
                    "    link.click();",
                    "    setTimeout(function() {",
                    "      link.remove();",
                    "      URL.revokeObjectURL(url);",
                    "    }, 1500);",
                    "  }",
                    "  document.addEventListener('click', function(event) {",
                    "    var link = event.target.closest ?",
                    "      event.target.closest('.mamamia-download-row .shiny-download-link') :",
                    "      null;",
                    "    if (!link) return;",
                    "    if (event.button !== 0 || event.metaKey || event.ctrlKey ||",
                    "        event.shiftKey || event.altKey) return;",
                    "    var row = link.closest('.mamamia-download-row');",
                    "    if (!row || row.getAttribute('data-mamamia-busy') === 'true') {",
                    "      event.preventDefault();",
                    "      return;",
                    "    }",
                    "    var href = link.getAttribute('href');",
                    "    if (!href || href === '#' || /^javascript:/i.test(href)) return;",
                    "    if (typeof window.fetch !== 'function') return;",
                    "    event.preventDefault();",
                    "    event.stopPropagation();",
                    "    row.setAttribute('data-mamamia-busy', 'true');",
                    "    var timer = setTimeout(function() { setBusy(row); }, 120);",
                    "    fetch(href, { credentials: 'same-origin' })",
                    "      .then(function(response) {",
                    "        if (!response.ok) {",
                    "          throw new Error('Download failed with ' + response.status);",
                    "        }",
                    "        var header = response.headers.get('Content-Disposition') || '';",
                    "        var fallback = link.getAttribute('download') ||",
                    "          'mamamia-download';",
                    "        var filename = filenameFromHeader(header, fallback);",
                    "        return response.blob().then(function(blob) {",
                    "          return { blob: blob, filename: filename };",
                    "        });",
                    "      })",
                    "      .then(function(result) {",
                    "        clearTimeout(timer);",
                    "        if (!result.blob || result.blob.size < 1) {",
                    "          throw new Error('Empty download');",
                    "        }",
                    "        saveBlob(result.blob, result.filename);",
                    "        setIdle(row);",
                    "      })",
                    "      .catch(function() {",
                    "        clearTimeout(timer);",
                    "        setError(row);",
                    "        setTimeout(function() { setIdle(row); }, 6000);",
                    "      })",
                    "      .finally(function() {",
                    "        row.setAttribute('data-mamamia-busy', 'false');",
                    "      });",
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
                                    mamamia_app_version()
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
        auto_species <- shiny::reactiveVal(NULL)
        species_manual <- shiny::reactiveVal(FALSE)

        shiny::observeEvent(example_mode(), {
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
            bam_single_file_input(config$server$max_total_bytes)
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
            detected <- auto_species()
            if (!is.null(detected$recipient) &&
                detected$recipient %in% species_keys) {
                return(detected$recipient)
            }
            value <- input$recipient_species
            if (!is.null(value) && length(value) == 1L && !is.na(value) &&
                value %in% species_keys) {
                return(as.character(value))
            }
            default_recipient
        })
        donor_species <- shiny::reactive({
            if (isTRUE(example_mode())) {
                return(bam_example_species_selection()$donor)
            }
            detected <- auto_species()
            if (!is.null(detected$donor) && detected$donor %in% species_keys) {
                return(detected$donor)
            }
            value <- input$donor_species
            if (!is.null(value) && length(value) == 1L && !is.na(value) &&
                value %in% species_keys) {
                return(as.character(value))
            }
            default_donor
        })

        shiny::observeEvent(input$recipient_species, {
            value <- input$recipient_species
            detected <- auto_species()
            if (is.null(detected) || is.null(detected$recipient) ||
                is.null(value) || length(value) != 1L || is.na(value) ||
                identical(as.character(value), detected$recipient)) {
                return()
            }
            species_manual(TRUE)
            detected$recipient <- NULL
            auto_species(detected)
        }, ignoreInit = TRUE)
        shiny::observeEvent(input$donor_species, {
            value <- input$donor_species
            detected <- auto_species()
            if (is.null(detected) || is.null(detected$donor) ||
                is.null(value) || length(value) != 1L || is.na(value) ||
                identical(as.character(value), detected$donor)) {
                return()
            }
            species_manual(TRUE)
            detected$donor <- NULL
            auto_species(detected)
        }, ignoreInit = TRUE)

        shiny::observeEvent(uploaded_paths(), {
            paths <- uploaded_paths()
            if (length(paths) == 0L || isTRUE(example_mode())) {
                auto_species(NULL)
                species_manual(FALSE)
                return()
            }
            species_manual(FALSE)
            current_recipient <- recipient_species()
            current_donor <- donor_species()
            pair <- tryCatch(
                bam_detect_species_pair(
                    paths[[1L]], config,
                    current_recipient = current_recipient,
                    current_donor = current_donor
                ),
                error = function(e) NULL
            )
            if (!is.null(pair)) {
                recipient <- pair$recipient
                donor <- pair$donor
                detected_label <- sprintf(
                    "%s recipient, %s donor",
                    pair$recipient_label, pair$donor_label
                )
            } else {
                detected <- tryCatch(
                    bam_detect_species(paths[[1L]], config),
                    error = function(e) NULL
                )
                if (is.null(detected)) {
                    auto_species(NULL)
                    return()
                }
                recipient <- detected$species
                donor <- current_donor
                if (identical(recipient, donor)) {
                    donor <- current_recipient
                }
                if (is.null(donor) || !donor %in% species_keys ||
                    identical(donor, recipient)) {
                    alternatives <- setdiff(species_keys, recipient)
                    donor <- if (length(alternatives) > 0L) {
                        alternatives[[1L]]
                    } else {
                        recipient
                    }
                }
                detected_label <- detected$label
            }
            auto_species(list(recipient = recipient, donor = donor))
            shiny::updateSelectInput(
                session, "recipient_species", selected = recipient
            )
            shiny::updateSelectInput(
                session, "donor_species", selected = donor
            )
            log_line(sprintf(
                "Detected species from the BAM header: %s.",
                detected_label
            ))
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$swap_species, {
            if (!species_available) {
                return()
            }
            recipient_key <- recipient_species()
            donor_key <- donor_species()
            if (is.null(recipient_key) || is.null(donor_key)) {
                return()
            }
            species_manual(TRUE)
            auto_species(list(recipient = donor_key, donor = recipient_key))
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
            ## Resolve the header match before the first preflight, not only in
            ## the observer that updates the selectors. Otherwise the default
            ## pair can log a transient assembly failure before detection.
            if (!isTRUE(species_manual())) {
                paths <- uploaded_paths()
                if (length(paths) > 0L) {
                    detected <- tryCatch(
                        bam_detect_species_pair(
                            paths[[1L]], config,
                            current_recipient = recipient_key,
                            current_donor = donor_key
                        ),
                        error = function(e) NULL
                    )
                    if (!is.null(detected)) {
                        recipient_key <- detected$recipient
                        donor_key <- detected$donor
                    }
                }
            }
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
            server <- config$server %||% .config_server_defaults()
            limit <- server$max_total_bytes %||%
                .config_server_defaults()$max_total_bytes
            too_large <- tryCatch(
                bam_upload_limit(input$bams, limit),
                error = function(e) e
            )
            if (inherits(too_large, "condition")) {
                return(bam_preflight_failure(too_large, character(), config))
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

        oversized_upload_key <- shiny::reactiveVal(NULL)
        shiny::observeEvent(preflight_state(), {
            pf <- preflight_state()
            if (isTRUE(pf$ok) || isTRUE(pf$empty) ||
                is.null(pf$error) ||
                !identical(pf$error$stage, "upload")) {
                oversized_upload_key(NULL)
                return()
            }
            key <- paste(
                c(pf$error$files, pf$error$size, pf$error$limit),
                collapse = "\n"
            )
            if (identical(oversized_upload_key(), key)) {
                return()
            }
            oversized_upload_key(key)
            shiny::showModal(bam_upload_too_large_modal(pf$error))
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$bam_upload_help, {
            server <- config$server %||% .config_server_defaults()
            limit <- server$max_total_bytes %||%
                .config_server_defaults()$max_total_bytes
            shiny::showModal(bam_upload_help_modal(limit))
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$oversized_upload, {
            server <- config$server %||% .config_server_defaults()
            limit <- server$max_total_bytes %||%
                .config_server_defaults()$max_total_bytes
            info <- bam_upload_payload_info(
                input$oversized_upload, limit
            )
            if (is.null(info)) {
                return()
            }
            upload_epoch(upload_epoch() + 1L)
            shiny::showModal(bam_upload_too_large_modal(info))
        }, ignoreInit = TRUE)

        shiny::observeEvent(input$clear_oversized_upload, {
            shiny::removeModal()
            upload_epoch(upload_epoch() + 1L)
        }, ignoreInit = TRUE)

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

        ## Profile state -------------------------------------------------------
        ## The named profile is the template. Advanced edits never change the
        ## template; they only make the effective profile Custom.
        profile_base <- shiny::reactiveVal(bam_default_profile(config))
        profile_defaults <- shiny::reactive({
            bam_profile_default_parameters(config, profile_base())
        })
        profile_current <- shiny::reactive({
            bam_profile_current_parameters(
                input, config,
                resolution = resolution_value(),
                fallback = profile_defaults()
            )
        })
        profile_changed <- shiny::reactive({
            bam_profile_changed_parameters(
                profile_current(), profile_defaults()
            )
        })
        profile_custom <- shiny::reactive({
            length(profile_changed()) > 0L
        })

        output$analysis_profile_ui <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                return(shiny::tagList(
                    bam_fixed_control_ui(
                        "Analysis profile",
                        mamamia_analysis_profiles()$balanced$label
                    ),
                    shiny::p(
                        class = "mamamia-fixed-note",
                        "Fixed for the bundled example."
                    )
                ))
            }
            choices <- bam_analysis_profile_choices(config)
            base <- profile_base()
            custom <- isTRUE(profile_custom())
            selected <- if (custom) "custom" else base
            profile <- mamamia_analysis_profile(base)
            caption <- if (custom) {
                paste0(
                    "Custom settings based on ",
                    if (is.null(profile)) base else profile$label,
                    ". Highlighted Advanced values differ from the ",
                    "profile defaults."
                )
            } else if (is.null(profile)) {
                "Custom analysis settings."
            } else {
                profile$short_description
            }
            shiny::tagList(
                shiny::selectInput(
                    "analysis_profile", "Analysis profile",
                    choices = choices,
                    selected = selected
                ),
                shiny::p(
                    class = "mamamia-profile-caption",
                    caption
                )
            )
        })
        shiny::outputOptions(
            output, "analysis_profile_ui", suspendWhenHidden = FALSE
        )

        shiny::observeEvent(input$analysis_profile, {
            if (isTRUE(example_mode())) {
                return()
            }
            value <- input$analysis_profile
            if (is.null(value) || length(value) != 1L ||
                is.na(value) || !nzchar(value)) {
                return()
            }
            value <- as.character(value)
            if (identical(value, "custom")) {
                ## Custom is a status the form switches to when a value
                ## differs; it is not a template that resets the controls.
                return()
            }
            profile <- tryCatch(
                mamamia_analysis_profile(value),
                error = function(e) NULL
            )
            if (is.null(profile)) {
                return()
            }
            if (!identical(profile_base(), value)) {
                profile_base(value)
            }
            bam_update_profile_inputs(
                session,
                mamamia_profile_flat_parameters(profile$parameters)
            )
        }, ignoreInit = FALSE)

        shiny::observeEvent(input$reset_profile_parameters, {
            if (isTRUE(example_mode())) {
                return()
            }
            defaults <- profile_defaults()
            bam_update_profile_inputs(session, defaults)
            shiny::updateSelectizeInput(
                session, "control_samples", selected = character(0)
            )
            profile <- mamamia_analysis_profile(profile_base())
            label <- if (is.null(profile)) profile_base() else profile$label
            log_line(sprintf(
                "Reset every Advanced parameter to the %s defaults.", label
            ))
            shiny::showNotification(
                sprintf("Advanced parameters reset to %s defaults.", label),
                type = "message", duration = 4
            )
        }, ignoreInit = TRUE)

        output$advanced_parameters_ui <- shiny::renderUI({
            bam_advanced_parameters_ui(
                config, profile_defaults(),
                locked = isTRUE(example_mode())
            )
        })

        output$mask_recurrence_warning_ui <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                return(NULL)
            }
            pf <- preflight_state()
            if (!isTRUE(pf$ok) || is.null(pf$table) ||
                nrow(pf$table) == 0L) {
                return(NULL)
            }
            bam_cohort_mask_warning_ui(
                length(unique(as.character(pf$table$sample)))
            )
        })
        shiny::outputOptions(
            output, "mask_recurrence_warning_ui", suspendWhenHidden = FALSE
        )

        ## Parameter labels are separate outputs so editing a value never
        ## rebuilds the input and steals keyboard focus.
        output$resolution_param_label <- shiny::renderUI({
            current <- profile_current()$windows_per_chromosome
            default <- profile_defaults()$windows_per_chromosome
            bam_profile_parameter_label(
                "resolution", "Resolution",
                bam_profile_value_changed(current, default),
                help = bam_parameter_help("resolution")
            )
        })
        parameter_labels <- list(
            metric = "Metric",
            score = "Score",
            donor_aggregation = "Donor aggregation",
            mapq = "Minimum MAPQ",
            calibration_fdr = "Calibration FDR",
            calibration_min_mb = "Calibration minimum length (Mb)",
            calibration_merge_gap = "Calibration merge gap (windows)"
        )
        for (parameter_id in names(parameter_labels)) {
            local({
                id <- parameter_id
                label <- parameter_labels[[id]]
                output[[paste0("param_label_", id)]] <- shiny::renderUI({
                    bam_profile_parameter_label(
                        id, label,
                        bam_profile_value_changed(
                            profile_current()[[id]],
                            profile_defaults()[[id]]
                        ),
                        help = bam_parameter_help(id, label)
                    )
                })
            })
        }
        floor_labels <- bam_calibration_floor_labels(config)
        for (i in seq_along(floor_labels)) {
            local({
                index <- i
                id <- paste0("calibration_floor_", index)
                label <- floor_labels[[index]]
                output[[paste0("param_label_", id)]] <- shiny::renderUI({
                    current <- profile_current()$calibration_floors[[index]]
                    default <- profile_defaults()$calibration_floors[[index]]
                    bam_profile_parameter_label(
                        id, label,
                        bam_profile_value_changed(current, default),
                        help = bam_parameter_help(id, label)
                    )
                })
            })
        }

        output$resolution_ui <- shiny::renderUI({
            if (isTRUE(example_mode())) {
                return(bam_example_resolution_slider_ui(config$defaults))
            }
            display <- config$defaults
            display$windows_per_chromosome <-
                profile_defaults()$windows_per_chromosome
            bam_resolution_slider_ui(display)
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
                error = function(e) {
                    fallback <- tryCatch(
                        profile_defaults()$windows_per_chromosome,
                        error = function(e) NULL
                    )
                    if (is.null(fallback) || length(fallback) != 1L ||
                        is.na(fallback)) {
                        fallback <- bam_default_resolution(defaults)
                    }
                    as.integer(fallback)
                }
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

        output$control_samples_ui <- shiny::renderUI({
            pf <- preflight_state()
            if (!isTRUE(pf$ok) || is.null(pf$table) ||
                !"sample" %in% names(pf$table)) {
                return(NULL)
            }
            samples <- unique(as.character(pf$table$sample))
            samples <- samples[!is.na(samples) & nzchar(samples)]
            samples <- mamamia_natural_sort(samples)
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
                    "control_samples",
                    shiny::tagList(
                        "Calibration controls",
                        bslib::tooltip(
                            shiny::span(
                                class = "mamamia-param-help",
                                shiny::icon("circle-question")
                            ),
                            bam_parameter_help("control_samples"),
                            placement = "right"
                        )
                    ),
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
                upload_error <- identical(as.character(stage), "upload")
                label <- if (upload_error) {
                    "BAM too large"
                } else if (assembly_error) {
                    "Needs assembly"
                } else {
                    "Preflight failed"
                }
                badge <- if (upload_error || assembly_error) {
                    "text-bg-warning"
                } else {
                    "text-bg-danger"
                }
                detail <- if (is.null(pf$error)) {
                    "See the calculation log."
                } else if (upload_error) {
                    paste(
                        "Slim it locally with bamhusk, then upload",
                        "the smaller BAM."
                    )
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
            bam_sidebar_species_ui(
                config,
                selection = list(
                    recipient = recipient_species(),
                    donor = donor_species()
                )
            )
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
                },
                shiny::uiOutput("mask_recurrence_warning_ui")
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
                input$control_samples,
                input$analysis_profile,
                parameters = profile_current()
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
            mamamia_natural_sort(unique(as.character(pf$table$sample)))
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
                    resolution = resolution_value(),
                    profile_base = profile_base()
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
                params$profile_name <- "balanced"
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
                    progress = progress,
                    profile = params$profile_name,
                    profile_base = params$profile_base
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
            bam_cohort_export_ui(result())
        })

        ## The summary depends only on the result; rebuilding it whenever the
        ## column picker changes rescans every call and metric table.
        cohort_summary <- shiny::reactive({
            res <- result()
            if (is.null(res)) {
                return(NULL)
            }
            bam_cohort_summary(res)
        })

        output$cohort_table <- DT::renderDT({
            summary <- cohort_summary()
            shiny::validate(
                shiny::need(
                    !is.null(summary),
                    "Run the analysis to see the cohort summary"
                )
            )
            bam_cohort_table(
                summary,
                input$cohort_columns
            )
        })

        output$download_cohort <- shiny::downloadHandler(
            filename = function() {
                "mamamia-cohort-introgressions.bedpe"
            },
            content = function(file) {
                bam_cohort_write_bedpe(result(), file)
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
