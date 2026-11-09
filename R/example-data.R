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
    dir <- mamamia_path("data", "example")
    if (!nzchar(dir) || !dir.exists(dir)) {
        return("")
    }
    normalizePath(dir, winslash = "/", mustWork = TRUE)
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
.bam_example_sanitize_paths <- function(data) {
    bam <- "example/example.bam"
    index <- "example/example.bam.csi"
    rewrite_bam <- function(x) {
        if (is.data.frame(x)) {
            if ("path" %in% names(x)) {
                x$path <- rep(bam, nrow(x))
            }
            if ("index_path" %in% names(x)) {
                x$index_path <- rep(index, nrow(x))
            }
        }
        x
    }
    rewrite_assembly <- function(assembly) {
        if (!is.list(assembly) || !is.list(assembly$contigs)) {
            return(assembly)
        }
        for (role in names(assembly$contigs)) {
            contig <- assembly$contigs[[role]]
            if (is.data.frame(contig) && "index_path" %in% names(contig)) {
                contig$index_path <- rep("example/demo.gcrank", nrow(contig))
                assembly$contigs[[role]] <- contig
            }
        }
        assembly
    }
    if (is.data.frame(data$preflight)) {
        data$preflight <- rewrite_bam(data$preflight)
    }
    if (is.list(data$results)) {
        for (key in names(data$results)) {
            result <- data$results[[key]]
            if (!is.list(result)) {
                next
            }
            result$bams <- rewrite_bam(result$bams)
            result$samples <- rewrite_bam(result$samples)
            result$assembly <- rewrite_assembly(result$assembly)
            provenance <- result$provenance
            if (is.list(provenance)) {
                ## Results written before the app-only layout stored the
                ## application version under an older field name.
                if (is.null(provenance$app)) {
                    provenance$app <- "MaMaMIA"
                }
                if (is.null(provenance$app_version)) {
                    provenance$app_version <- if (
                        !is.null(provenance$package_version)
                    ) {
                        provenance$package_version
                    } else {
                        mamamia_version()
                    }
                }
                if (is.null(provenance$analysis_id) ||
                    length(provenance$analysis_id) != 1L ||
                    is.na(provenance$analysis_id) ||
                    !nzchar(provenance$analysis_id)) {
                    provenance$analysis_id <- paste0("example-", key)
                }
                provenance$bams <- rewrite_bam(provenance$bams)
                if ("config_path" %in% names(provenance)) {
                    provenance$config_path <- "example/assemblies.yml"
                }
                if ("uvr_lock" %in% names(provenance)) {
                    provenance$uvr_lock <- "uvr.lock"
                }
                result$provenance <- provenance
            }
            data$results[[key]] <- result
        }
    }
    data
}

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
                "The bundled example data is missing from data/example"
            )
        )
    }
    data <- .bam_example_sanitize_paths(readRDS(path))
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
                "The bundled example registry is missing from data/example"
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
    recipient <- config$species[[selection$recipient]]
    donor <- config$species[[selection$donor]]
    resolved <- tryCatch(
        resolveSpeciesAssembly(
            config,
            donor = selection$donor,
            recipient = selection$recipient,
            donor_components = unique(as.character(donor$contigs$component)),
            recipient_components = unique(
                as.character(recipient$contigs$component)
            )
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
        recipient_component = selection$recipient_component,
        donor = selection$donor,
        donor_component = selection$donor_component
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
        "Example dataset: %s recipient and %s donor; events on eight chromosomes.",
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
    metric_totals <- tapply(
        as.numeric(result$metrics$count),
        as.character(result$metrics$sample),
        sum, na.rm = TRUE
    )
    for (sample in samples) {
        total <- unname(metric_totals[sample])
        if (length(total) == 0L || is.na(total)) {
            total <- 0
        }
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
    call_counts <- table(as.character(result$calls$sample))
    for (sample in samples) {
        n_segments <- unname(call_counts[sample])
        if (length(n_segments) == 0L || is.na(n_segments)) {
            n_segments <- 0L
        }
        report("segmentation", list(
            sample = sample,
            n_segments = n_segments
        ))
    }
    report("calibration", list(
        n_calls = sum(result$calls$pass, na.rm = TRUE),
        fdr = result$provenance$fdr
    ))
    invisible(NULL)
}
