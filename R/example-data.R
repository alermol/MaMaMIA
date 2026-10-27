## Bundled example data ------------------------------------------------------
##
## The application offers a small synthetic example without requiring a BAM
## upload or a server-side BAM directory. The data file contains the complete
## `BamAnalysisResult` for every supported resolution, plus a preflight table
## that looks exactly like a real single-BAM preflight. The app replays the
## stored stage details through the same progress callback as
## `analyzeBamInput()`, so the user sees the normal pipeline messages.

#' @noRd
bam_example_dir <- function() {
    dir <- system.file("extdata", "example", package = "MaMaMIA")
    if (!nzchar(dir) || !dir.exists(dir)) {
        candidate <- file.path("inst", "extdata", "example")
        if (dir.exists(candidate)) {
            dir <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
        }
    }
    dir
}

#' @noRd
bam_example_available <- function() {
    dir <- bam_example_dir()
    if (!nzchar(dir) ||
        !file.exists(file.path(dir, "example.rds")) ||
        !file.exists(file.path(dir, "assemblies.yml"))) {
        return(FALSE)
    }
    ## The registry validates the bundled gcrank index; a broken install must
    ## not offer an example that cannot produce a preflight state.
    isTRUE(tryCatch(
        {
            bam_example_config()
            TRUE
        },
        error = function(e) FALSE
    ))
}

.bam_example_cache <- new.env(parent = emptyenv())

#' @noRd
bam_example_data <- function() {
    if (!is.null(.bam_example_cache$data)) {
        return(.bam_example_cache$data)
    }
    dir <- bam_example_dir()
    path <- file.path(dir, "example.rds")
    if (!nzchar(dir) || !file.exists(path)) {
        .mamamia_error(
            "MaMaMIAConfigError",
            paste(
                "The bundled example data is missing;",
                "reinstall MaMaMIA with inst/extdata/example"
            )
        )
    }
    data <- readRDS(path)
    .bam_example_cache$data <- data
    data
}

#' @noRd
bam_example_config <- function() {
    if (!is.null(.bam_example_cache$config)) {
        return(.bam_example_cache$config)
    }
    dir <- bam_example_dir()
    path <- file.path(dir, "assemblies.yml")
    if (!nzchar(dir) || !file.exists(path)) {
        .mamamia_error(
            "MaMaMIAConfigError",
            paste(
                "The bundled example registry is missing;",
                "reinstall MaMaMIA with inst/extdata/example"
            )
        )
    }
    config <- readMaMaMIAConfig(path)
    .bam_example_cache$config <- config
    config
}

#' Supported example resolutions, in ascending order
#' @noRd
bam_example_resolutions <- function() {
    data <- bam_example_data()
    sort(unique(as.integer(names(data$results))))
}

#' @noRd
bam_example_resolution_choices <- function() {
    resolutions <- bam_example_resolutions()
    if (length(resolutions) == 0L) {
        return(character())
    }
    bam_resolution_choices(list(allowed_windows = resolutions))
}

#' @noRd
bam_example_default_resolution <- function() {
    resolutions <- bam_example_resolutions()
    if (length(resolutions) == 0L) {
        return(NULL)
    }
    as.character(if (1000L %in% resolutions) 1000L else resolutions[[1L]])
}

#' @noRd
bam_example_preflight <- function() {
    data <- bam_example_data()
    table <- data$preflight
    if (!is.data.frame(table) || nrow(table) == 0L) {
        .mamamia_error(
            "MaMaMIAConfigError", "The bundled example has no preflight data"
        )
    }
    table
}

#' Resolve the pre-computed result for a resolution
#' @noRd
bam_example_result <- function(resolution) {
    value <- suppressWarnings(as.integer(resolution))
    resolutions <- bam_example_resolutions()
    if (length(value) != 1L || is.na(value) || !value %in% resolutions) {
        .mamamia_error(
            "MaMaMIAWindowError",
            "The bundled example is available at resolution(s): %s",
            paste(resolutions, collapse = ", ")
        )
    }
    data <- bam_example_data()
    value <- data$results[[as.character(value)]]
    if (is.null(value)) {
        .mamamia_error(
            "MaMaMIAWindowError",
            "The bundled example has no result for resolution %s",
            as.character(resolution)
        )
    }
    value
}

#' Build the preflight state for the bundled example
#'
#' Mirrors the successful state returned for uploaded BAMs, but without a
#' staged path: the run loads a pre-computed result instead of reading a BAM.
#' @noRd
bam_example_preflight_state <- function() {
    config <- bam_example_config()
    data <- bam_example_data()
    selection <- data$species_selection
    resolved <- tryCatch(
        resolveSpeciesAssembly(
            config,
            donor = selection$donor,
            recipient = selection$recipient,
            donor_subgenome = selection$donor_subgenome,
            recipient_subgenome = selection$recipient_subgenome
        ),
        error = function(e) e
    )
    if (inherits(resolved, "condition")) {
        .mamamia_error(
            "MaMaMIAConfigError",
            "The bundled example registry cannot resolve its assembly: %s",
            conditionMessage(resolved)
        )
    }
    list(
        ok = TRUE, empty = FALSE, table = bam_example_preflight(),
        staged = NULL, assembly = resolved$key,
        assembly_config = resolved$assembly, config = resolved$config,
        key = resolved$key, unpaired = resolved$unpaired,
        error = NULL, candidates = NULL, example = TRUE
    )
}

#' @noRd
bam_example_species_selection <- function() {
    data <- bam_example_data()
    selection <- data$species_selection
    list(
        recipient = selection$recipient,
        recipient_subgenome = selection$recipient_subgenome,
        donor = selection$donor,
        donor_subgenome = selection$donor_subgenome
    )
}

#' @noRd
bam_example_note <- function() {
    selection <- bam_example_species_selection()
    config <- bam_example_config()
    recipient <- config$species[[selection$recipient]]
    donor <- config$species[[selection$donor]]
    recipient_label <- if (is.null(recipient)) {
        selection$recipient
    } else {
        recipient$label
    }
    donor_label <- if (is.null(donor)) {
        selection$donor
    } else {
        donor$label
    }
    sprintf(
        "Example dataset: %s recipient and %s donor; one sample.",
        recipient_label, donor_label
    )
}

#' Replay the stored pipeline stages through a progress callback
#'
#' The stage order and detail fields are the same as `analyzeBamInput()`, so
#' `bam_progress_label()` renders identical messages for the example and a real
#' run. A small pause keeps the stages visible when the result is loaded
#' instantly from disk.
#'
#' @param result A pre-computed `BamAnalysisResult`.
#' @param progress A `function(stage, detail)` progress callback.
#' @param pause Seconds to wait before each stage.
#' @noRd
bam_example_progress <- function(result, progress, pause = 0) {
    if (!is.function(progress)) {
        return(invisible(NULL))
    }
    report <- function(stage, detail) {
        if (is.numeric(pause) && length(pause) == 1L && is.finite(pause) &&
            pause > 0) {
            Sys.sleep(pause)
        }
        progress(stage, detail)
    }

    samples <- bam_sample_choices(result)
    n_bams <- length(unique(result$bams$path))
    n_samples <- length(samples)
    n_windows <- nrow(result$windows$table)
    assembly_key <- as.character(result$samples$assembly[[1L]])

    report("preflight", list(
        n_bams = n_bams, n_samples = n_samples, assembly = assembly_key
    ))
    report("index", list(
        n_bams = n_bams,
        staged = sum(as.logical(result$provenance$staged), na.rm = TRUE)
    ))
    report("windows", list(
        n_per_chromosome = result$windows$n_per_chromosome,
        n_windows = n_windows
    ))
    report("gc", list(
        n_windows = nrow(result$gc),
        n_indexes = length(result$assembly$gc_indexes)
    ))
    for (sample in samples) {
        total <- sum(result$metrics$count[result$metrics$sample == sample])
        report("metrics", list(sample = sample, total = total))
    }
    report("score", list(
        mode = result$samples$mode[[1L]],
        metric = result$provenance$metric,
        n_windows = nrow(result$scores)
    ))
    report("mask", list(
        masked = sum(result$masks$mask, na.rm = TRUE),
        masked_fraction = mean(result$masks$mask, na.rm = TRUE),
        source = unique(result$masks$source[result$masks$mask])
    ))
    for (sample in samples) {
        report("segmentation", list(
            sample = sample,
            n_segments = sum(result$calls$sample == sample, na.rm = TRUE)
        ))
    }
    report("calibration", list(
        n_calls = sum(result$calls$pass, na.rm = TRUE),
        fdr = result$provenance$fdr
    ))
    invisible(NULL)
}
