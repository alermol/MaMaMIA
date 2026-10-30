## Cohort pipeline orchestration.
##
## `analyzeBamInput()` is the only BAM-first entrypoint of the package. It runs
## cheap checks first (preflight), indexes, generates group-aware windows,
## extracts GC, measures counts and breadth in one htslib pass, builds the
## WindowMetrics table, applies a window mask, scores donor/recipient
## homoeologs, segments with CBS, calibrates thresholds by segment length and
## classifies the calls.
##
## The result is the `BamAnalysisResult` contract: `assembly`, `config`,
## `bams`, `samples`, `windows`, `gc`, `metrics`, `scores`, `masks`,
## `calibration`, `calls` and `provenance`.

## Synchronized with DESCRIPTION by tests/testthat/test-pipeline.R.
.MaMaMIAPipelineVersion <- "0.99.2"

## Synchronized with tools/vendor-gcrank.txt by tests/testthat/test-pipeline.R.
.MaMaMIAgcrankCommit <- "aa3a32b4cd748ead4c933225d0c990689499f92f"

.analysis_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAAnalysisError", fmt, ...)
}

.pipeline_integer <- function(x, what, minimum = 1L,
                              maximum = .Machine$integer.max) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
        !is.finite(x) || x != floor(x) || x < minimum || x > maximum) {
        .analysis_error(
            "`%s` must be a single integer between %s and %s",
            what,
            format(minimum, scientific = FALSE),
            format(maximum, scientific = FALSE)
        )
    }
    as.integer(x)
}

.pipeline_config <- function(config) {
    if (is.character(config)) {
        config <- readMaMaMIAConfig(config)
    }
    if (!is.list(config) || is.null(config$assemblies)) {
        .analysis_error("`config` must be a validated MaMaMIA configuration")
    }
    if (length(config$assemblies) == 0L) {
        .analysis_error(
            "`config` has no assembly: pair the species first with resolveSpeciesAssembly()"
        )
    }
    if (is.null(config$defaults) || is.null(config$server)) {
        config <- validateMaMaMIAConfig(config)
    }
    if (!identical(as.integer(config$version), 2L)) {
        .analysis_error(
            "MaMaMIA v1.0 requires a `version: 2` configuration; see MIGRATION.md"
        )
    }
    config
}

.pipeline_progress <- function(progress, stage, detail = list()) {
    if (is.null(progress)) {
        return(invisible(NULL))
    }
    progress(stage, detail)
    invisible(NULL)
}

.validate_progress <- function(progress) {
    if (!is.null(progress) && !is.function(progress)) {
        .analysis_error(
            "`progress` must be NULL or a function(stage, detail)"
        )
    }
    invisible(progress)
}

.pipeline_contigs <- function(assembly) {
    indexes <- assembly$gc_indexes
    if (is.null(indexes) || length(indexes) == 0L) {
        return(NULL)
    }
    pieces <- lapply(
        indexes,
        function(index) readContigsTsv(file.path(index, "contigs.tsv"))
    )
    contigs <- do.call(rbind, lapply(pieces, function(x) {
        x[, c("name", "length"), drop = FALSE]
    }))
    rownames(contigs) <- NULL
    if (anyDuplicated(contigs$name)) {
        duplicate <- unique(contigs$name[duplicated(contigs$name)])[[1L]]
        .analysis_error(
            "Contig '%s' appears in more than one gc_index", duplicate
        )
    }
    contigs
}

.pipeline_index <- function(preflight) {
    paths <- unique(as.character(preflight$path))
    representative <- match(paths, preflight$path)
    n <- length(paths)
    working <- character(n)
    index_path <- character(n)
    staged <- logical(n)

    for (i in seq_len(n)) {
        original <- paths[[i]]
        if (identical(preflight$index_action[[representative[[i]]]], "stage")) {
            working[[i]] <- stageBam(original, .session_temp_dir())
            staged[[i]] <- TRUE
        } else {
            working[[i]] <- original
        }
        index_path[[i]] <- ensureBamIndex(working[[i]])
    }

    list(
        path = working,
        original_path = paths,
        index_path = index_path,
        staged = staged
    )
}

.sample_stage_error <- function(e, sample, stage) {
    class <- if (inherits(e, "MaMaMIAError")) {
        class(e)[[1L]]
    } else {
        "MaMaMIASampleError"
    }
    .mamamia_error(
        class,
        "Sample '%s' (%s): %s",
        sample, stage, conditionMessage(e)
    )
}

.safe_package_version <- function(package) {
    tryCatch(
        as.character(utils::packageVersion(package)),
        error = function(e) NA_character_
    )
}

.htslib_version <- function() {
    version <- htslib_version_cpp()
    if (length(version) != 1L || is.na(version) || !nzchar(version) ||
        identical(version, "unknown")) {
        return(NA_character_)
    }
    version
}

.find_upwards <- function(name, start = getwd(), depth = 6L) {
    directory <- tryCatch(
        normalizePath(start, winslash = "/", mustWork = FALSE),
        error = function(e) start
    )
    for (i in seq_len(depth)) {
        candidate <- file.path(directory, name)
        if (file.exists(candidate)) {
            return(normalizePath(candidate, winslash = "/", mustWork = FALSE))
        }
        parent <- dirname(directory)
        if (identical(parent, directory)) {
            break
        }
        directory <- parent
    }
    NA_character_
}

.uvr_version <- function() {
    executable <- Sys.which("uvr")
    if (!nzchar(executable)) {
        return(NA_character_)
    }
    output <- suppressWarnings(tryCatch(
        system2(executable, "--version", stdout = TRUE, stderr = TRUE),
        error = function(e) character(0)
    ))
    status <- attr(output, "status")
    if (length(output) == 0L || (!is.null(status) && status != 0L)) {
        return(NA_character_)
    }
    trimws(output[[1L]])
}

.uvr_lock_path <- function() {
    candidates <- c(
        Sys.getenv("MAMAMIA_UVRLOCK", unset = ""),
        getOption("MaMaMIA.uvr_lock", ""),
        system.file("uvr.lock", package = "MaMaMIA"),
        .find_upwards("uvr.lock")
    )
    candidates <- candidates[
        nzchar(candidates) & !is.na(candidates) & file.exists(candidates)
    ]
    if (length(candidates) == 0L) {
        return(NA_character_)
    }
    normalizePath(candidates[[1L]], winslash = "/", mustWork = FALSE)
}

.pipeline_config_path <- function(config) {
    path <- attr(config, "config_path", exact = TRUE)
    if (is.null(path) || !is.character(path) || length(path) != 1L ||
        is.na(path) || !nzchar(path)) {
        return(NA_character_)
    }
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.pipeline_provenance <- function(config, preflight, indexed, windows, gc,
                                 n_per, mapq, exclude_flags, metric, score,
                                 donor_aggregation, sample_composition,
                                 mask, calibration, started_at) {
    config_path <- .pipeline_config_path(config)
    config_hash <- if (!is.na(config_path) && file.exists(config_path)) {
        unname(tools::md5sum(config_path))
    } else {
        plain <- config
        attr(plain, "config_path") <- NULL
        hash_object(plain)
    }
    uvr_lock <- .uvr_lock_path()
    bams <- data.frame(
        path = preflight$path,
        sample = preflight$sample,
        sample_source = preflight$sample_source,
        file_size = preflight$file_size,
        mtime = as.character(file.info(preflight$path)$mtime),
        assembly = preflight$assembly,
        stringsAsFactors = FALSE
    )
    mask_sources <- unique(mask$source[mask$mask %in% TRUE])
    list(
        pipeline_version = .MaMaMIAPipelineVersion,
        package = "MaMaMIA",
        package_version = .safe_package_version("MaMaMIA"),
        r_version = as.character(getRversion()),
        dependency_versions = c(
            Rhtslib = .safe_package_version("Rhtslib"),
            GenomicRanges = .safe_package_version("GenomicRanges"),
            DNAcopy = .safe_package_version("DNAcopy")
        ),
        htslib_version = .htslib_version(),
        gcrank_commit = .MaMaMIAgcrankCommit,
        uvr_version = .uvr_version(),
        uvr_lock = uvr_lock,
        uvr_lock_hash = if (is.na(uvr_lock)) {
            NA_character_
        } else {
            unname(tools::md5sum(uvr_lock))
        },
        config_path = config_path,
        config_hash = config_hash,
        assembly = preflight$assembly[[1L]],
        bams = bams,
        staged = indexed$staged,
        n_bams = length(unique(preflight$path)),
        window_rule = "exact_equal_windows",
        windows_per_chromosome = n_per,
        n_samples = nrow(preflight),
        total_windows = nrow(windows$table) * nrow(preflight),
        n_windows = nrow(windows$table),
        n_gc_windows = nrow(gc),
        mapq = mapq,
        exclude_flags = exclude_flags,
        metric = metric,
        score = score,
        donor_aggregation = donor_aggregation,
        sample_composition = as.list(sample_composition),
        breadth_backend = "htslib",
        count_backend = "htslib",
        count_semantics = "aligned_base_overlap",
        gc_source = "gcrank",
        mask_source = if (length(mask_sources) == 0L) {
            "none"
        } else {
            paste(mask_sources, collapse = "+")
        },
        mask_version = {
            versions <- unique(stats::na.omit(mask$version))
            if (length(versions) == 0L) NA_character_ else versions[[1L]]
        },
        calibration_version = {
            versions <- unique(stats::na.omit(calibration$version))
            if (length(versions) == 0L) NA_character_ else versions[[1L]]
        },
        fdr = if (nrow(calibration) == 0L) NA_real_ else calibration$fdr[[1L]],
        started_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
        duration_secs = as.numeric(difftime(
            Sys.time(), started_at, units = "secs"
        ))
    )
}

.pipeline_measure <- function(indexed, preflight, windows, headers,
                              samples, mapq, exclude_flags) {
    paths <- as.character(indexed$original_path)
    selected_paths <- as.character(preflight$path)
    count_one <- function(bam_index) tryCatch({
        rows <- which(selected_paths == paths[[bam_index]])
        if (length(rows) == 0L) {
            return(NULL)
        }
        bam_samples <- samples[rows]
        measured <- tryCatch(
            countWindowMetrics(
                indexed$path[[bam_index]],
                windows,
                samples = bam_samples,
                read_groups = headers[[bam_index]]$read_groups,
                mapq = mapq,
                exclude_flags = exclude_flags
            ),
            error = function(e) e
        )
        if (inherits(measured, "condition")) {
            if (length(rows) == 1L) {
                return(.sample_stage_error(
                    measured, bam_samples[[1L]], "metrics"
                ))
            }
            return(.mamamia_error(
                class(measured)[[1L]],
                "BAM '%s' (%d samples): %s",
                basename(paths[[bam_index]]), length(rows),
                conditionMessage(measured)
            ))
        }
        if (!is.finite(measured$read_length) || measured$read_length <= 0) {
            return(.mamamia_error(
                class = "MaMaMIASampleError",
                "Sample%s %s has no usable reads in the analysis windows",
                if (length(bam_samples) == 1L) "" else "s",
                paste(sprintf("'%s'", bam_samples), collapse = ", ")
            ))
        }
        measured$metrics$read_length <- as.numeric(measured$read_length)
        measured$metrics
    }, error = function(e) e)

    n_bams <- length(paths)
    cores <- min(30L, max(1L, getOption("MaMaMIA.cores", 1L)))
    if (.Platform$OS.type == "unix" && cores > 1L && n_bams > 1L) {
        pieces <- parallel::mclapply(
            seq_len(n_bams), count_one,
            mc.cores = min(cores, n_bams),
            mc.preschedule = FALSE
        )
    } else {
        pieces <- lapply(seq_len(n_bams), count_one)
    }
    failed <- vapply(pieces, function(x) {
        inherits(x, "condition") || inherits(x, "try-error")
    }, logical(1))
    if (any(failed)) {
        piece <- pieces[[which(failed)[[1L]]]]
        if (inherits(piece, "try-error")) {
            .analysis_error("%s", as.character(piece))
        }
        stop(piece)
    }
    do.call(rbind, pieces)
}

#' Analyze one or more BAM files with the MaMaMIA workflow
#'
#' Reads BAM headers, matches them to a configured assembly, ensures CSI
#' indexes, generates group-aware exact-N windows, extracts GC from the
#' configured gcrank indexes, measures breadth and counts in one htslib pass,
#' scores donor/recipient homoeologs, applies a window mask, segments with CBS,
#' calibrates thresholds by segment length and classifies the calls.
#'
#' A single sample is scored and calibrated on its own (single-sample mode).
#' Several samples share a cohort median-polish correction and a cohort
#' calibration. The result supports any number of donor and recipient
#' components.
#'
#' @param bam_paths Character vector of coordinate-sorted BAM paths.
#' @param config A validated configuration list from [readMaMaMIAConfig()] or a
#'   path to a YAML configuration.
#' @param assembly Optional assembly key. Without it, the BAM headers select the
#'   most specific exact match.
#' @param sample_composition Optional named recipient copy-number vector, for
#'   example `c(A = 2, B = 2, D = 2)`.
#' @param metric `"breadth"`, `"count"` or `"combined"`.
#' @param score `"paired_log_ratio"` or `"donor_only"`.
#' @param donor_aggregation `"mean"`, `"max"` or `"sum"`. Defaults to the
#'   registry value; the dashboard exposes it under **Advanced**.
#' @param mask Optional panel `MaskTable`/TSV path, or a list with `panel`,
#'   `blacklist`, `physical_fallback_mb` and `recurrence_threshold`.
#' @param calibration Optional list with `length_bins_mb`, `fdr`, `floors`,
#'   `controls`, `alpha`, `undo_SD`, `min_width`, `merge_gap`, `merge_floor`
#'   and `min_mb`.
#' @param windows_per_chromosome Number of windows per target contig.
#' @param mapq Minimum mapping quality.
#' @param exclude_flags SAM flag mask; a read is excluded when any bit is set.
#' @param cores Number of cores for parallel BAM measurement.
#' @param progress Optional `function(stage, detail)` progress callback.
#'
#' @return A `BamAnalysisResult` list with `assembly`, `config`, `bams`,
#'   `samples`, `windows`, `gc`, `metrics`, `scores`, `masks`, `calibration`,
#'   `calls` and `provenance`.
#' @seealso [preflightBams()], [countWindowMetrics()], [buildWindowMetrics()],
#'   [scoreIntrogression()], [segmentScores()], [calibrateSegments()]
#' @examples
#' demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
#' config <- readMaMaMIAConfig(file.path(demo, "assemblies.yml"))
#' pair <- resolveSpeciesAssembly(
#'     config, donor = "secale_like", recipient = "triticum_like"
#' )
#' result <- analyzeBamInput(
#'     file.path(demo, "demo.bam"), pair$config,
#'     assembly = pair$key,
#'     windows_per_chromosome = 100L,
#'     cores = 1L
#' )
#' result$samples[, c("sample", "n_windows", "n_calls")]
#'
#' @export
analyzeBamInput <- function(bam_paths,
                            config,
                            assembly = NULL,
                            sample_composition = NULL,
                            metric = NULL,
                            score = NULL,
                            donor_aggregation = NULL,
                            mask = NULL,
                            calibration = NULL,
                            windows_per_chromosome = 1000L,
                            mapq = NULL,
                            exclude_flags = NULL,
                            cores = 1L,
                            progress = NULL) {
    started_at <- Sys.time()
    config <- .pipeline_config(config)
    .validate_progress(progress)
    defaults <- config$defaults
    n_per <- .pipeline_integer(
        windows_per_chromosome, "windows_per_chromosome"
    )
    mapq <- .pipeline_integer(
        if (is.null(mapq)) defaults$mapq else mapq, "mapq", 0L, 255L
    )
    exclude_flags <- .pipeline_integer(
        if (is.null(exclude_flags)) defaults$exclude_flags else exclude_flags,
        "exclude_flags", 0L, 4095L
    )
    metric <- if (is.null(metric)) defaults$metric else {
        match.arg(metric, c("breadth", "count", "combined"))
    }
    score <- if (is.null(score)) defaults$score else {
        match.arg(score, c("paired_log_ratio", "donor_only"))
    }
    donor_aggregation <- if (is.null(donor_aggregation)) {
        defaults$donor_aggregation
    } else {
        donor_aggregation
    }
    if (!is.character(donor_aggregation) ||
        length(donor_aggregation) != 1L ||
        is.na(donor_aggregation) || !nzchar(donor_aggregation) ||
        !donor_aggregation %in% c("mean", "max", "sum")) {
        .analysis_error(
            "`donor_aggregation` must be one of: mean, max, sum"
        )
    }
    cores <- .pipeline_integer(cores, "cores", 1L, 256L)
    old_cores <- getOption("MaMaMIA.cores")
    options(MaMaMIA.cores = cores)
    on.exit(options(MaMaMIA.cores = old_cores), add = TRUE)

    if (!is.character(bam_paths) || length(bam_paths) == 0L ||
        anyNA(bam_paths) || any(!nzchar(bam_paths))) {
        .analysis_error(
            "`bam_paths` must be a non-empty character vector of BAM paths"
        )
    }
    calibration_args <- if (is.null(calibration)) list() else calibration
    if (!is.list(calibration_args)) {
        .analysis_error("`calibration` must be NULL or a named list")
    }
    mask_args <- if (is.null(mask)) list() else {
        if (is.list(mask) && !is.data.frame(mask) &&
            any(c("panel", "blacklist", "physical_fallback_mb",
                  "recurrence_threshold") %in% names(mask))) {
            mask
        } else {
            list(panel = mask)
        }
    }

    preflight <- preflightBams(bam_paths, config, assembly)
    n_bams <- length(unique(preflight$path))
    n_samples <- nrow(preflight)
    assembly_key <- preflight$assembly[[1L]]
    assembly_config <- config$assemblies[[assembly_key]]
    if (is.null(assembly_config)) {
        .analysis_error(
            "Assembly '%s' selected by preflight is not in the configuration",
            assembly_key
        )
    }
    .pipeline_progress(progress, "preflight", list(
        n_bams = n_bams,
        n_samples = n_samples,
        assembly = assembly_key
    ))

    indexed <- .pipeline_index(preflight)
    .pipeline_progress(progress, "index", list(
        n_bams = n_bams,
        staged = sum(indexed$staged)
    ))

    contigs <- .pipeline_contigs(assembly_config)
    windows <- makeAnalysisWindows(assembly_config, contigs, n_per)
    .pipeline_progress(progress, "windows", list(
        n_per_chromosome = n_per,
        n_windows = nrow(windows$table)
    ))

    gc <- extractWindowGc(
        windows,
        assembly_config,
        mmap = defaults$gc_mmap,
        validate_crc = defaults$gc_validate_crc
    )
    .pipeline_progress(progress, "gc", list(
        n_windows = nrow(gc),
        n_indexes = length(assembly_config$gc_indexes)
    ))

    samples <- as.character(preflight$sample)
    headers <- lapply(indexed$original_path, bamHeaderInfo)
    measured <- .pipeline_measure(
        indexed, preflight, windows, headers, samples, mapq, exclude_flags
    )
    for (sample in samples) {
        total <- sum(measured$count[measured$sample == sample])
        .pipeline_progress(progress, "metrics", list(
            sample = sample,
            total = total
        ))
        if (!is.finite(total) || total <= 0) {
            .mamamia_error(
                "MaMaMIASampleError",
                "Sample '%s' has zero reads in all analysis windows",
                sample
            )
        }
    }
    metrics <- buildWindowMetrics(
        measured, windows, gc, read_length = NULL
    )
    mode <- if (n_samples == 1L) "single" else "multi"
    composition <- sample_composition
    if (is.null(composition)) {
        composition <- assembly_config$sample_composition
    }
    scores <- scoreIntrogression(
        metrics,
        mode = mode,
        metric = metric,
        donor_aggregation = donor_aggregation,
        recipient_aggregation = defaults$recipient_aggregation,
        sample_composition = composition,
        gc_correction = TRUE,
        score_type = score
    )
    .pipeline_progress(progress, "score", list(
        mode = mode,
        metric = metric,
        score = score,
        n_windows = nrow(scores)
    ))

    mask_spec <- defaults$mask
    mask_spec[names(mask_args)] <- mask_args
    masks <- buildWindowMask(
        metrics,
        scores = scores,
        mode = mode,
        panel_mask = mask_spec$panel,
        blacklist = mask_spec$blacklist,
        physical_fallback_mb = mask_spec$physical_fallback_mb,
        recurrence_threshold = mask_spec$recurrence_threshold
    )
    masked_fraction <- mean(masks$mask, na.rm = TRUE)
    .pipeline_progress(progress, "mask", list(
        masked = sum(masks$mask, na.rm = TRUE),
        masked_fraction = masked_fraction,
        source = unique(masks$source[masks$mask])
    ))

    segment_args <- list(
        alpha = 0.001,
        undo_SD = 1,
        min_width = 5L,
        merge_gap = 100L,
        merge_floor = -0.5,
        min_mb = 5
    )
    segment_keep <- intersect(names(calibration_args), names(segment_args))
    segment_args[segment_keep] <- calibration_args[segment_keep]
    donor_lengths <- tapply(
        windows$table$stop[windows$table$role == "donor"],
        windows$table$group[windows$table$role == "donor"],
        max, na.rm = TRUE
    )
    segment_args$lengths <- donor_lengths
    calls <- do.call(
        segmentScores,
        c(list(scores = scores, mask = masks), segment_args)
    )
    .pipeline_progress(progress, "segmentation", list(
        n_segments = nrow(calls),
        n_samples = length(unique(calls$sample))
    ))

    cal <- defaults$calibration
    cal[names(calibration_args)] <- calibration_args
    cal_floors <- cal$floors
    if (is.null(cal_floors)) {
        cal_floors <- if (identical(mode, "single")) {
            cal$single_sample_floors
        } else {
            rep(0.25, max(1L, length(cal$length_bins_mb) - 1L))
        }
    }
    calibrated <- calibrateSegments(
        calls,
        controls = cal$controls,
        mode = mode,
        length_bins_mb = cal$length_bins_mb,
        fdr = cal$fdr,
        floors = cal_floors
    )
    class_thresholds <- .call_class_thresholds(
        calibrated$calls, cal$controls
    )
    calls <- classifyCalls(
        calibrated$calls,
        controls = cal$controls,
        thresholds = class_thresholds
    )
    calls <- assignCallPositions(
        calls,
        metrics = metrics,
        windows = windows,
        scores = scores,
        assembly = assembly_config,
        masks = masks,
        class_thresholds = class_thresholds,
        controls = cal$controls
    )
    .pipeline_progress(progress, "calibration", list(
        n_calls = sum(calls$pass, na.rm = TRUE),
        fdr = cal$fdr
    ))

    sample_bam <- match(preflight$path, indexed$original_path)
    samples_meta <- data.frame(
        sample = samples,
        path = preflight$path,
        assembly = preflight$assembly,
        index_path = indexed$index_path[sample_bam],
        mode = mode,
        n_windows = rep(nrow(windows$table), length(samples)),
        n_calls = vapply(samples, function(s) {
            sum(calls$pass[calls$sample == s], na.rm = TRUE)
        }, integer(1L)),
        n_high_confidence = vapply(samples, function(s) {
            sum(
                calls$pass[calls$sample == s] &
                    calls$confidence[calls$sample == s] == "high",
                na.rm = TRUE
            )
        }, integer(1L)),
        masked_fraction = rep(masked_fraction, length(samples)),
        stringsAsFactors = FALSE
    )

    bams <- data.frame(
        path = preflight$path,
        sample = samples,
        index_path = indexed$index_path[sample_bam],
        n_targets = preflight$n_targets,
        stringsAsFactors = FALSE
    )
    provenance <- .pipeline_provenance(
        config = config,
        preflight = preflight,
        indexed = indexed,
        windows = windows,
        gc = gc,
        n_per = n_per,
        mapq = mapq,
        exclude_flags = exclude_flags,
        metric = metric,
        score = score,
        donor_aggregation = donor_aggregation,
        sample_composition = composition,
        mask = masks,
        calibration = calibrated$calibration,
        started_at = started_at
    )

    result <- structure(
        list(
            assembly = assembly_config,
            config = config,
            bams = bams,
            samples = samples_meta,
            windows = windows,
            gc = gc,
            metrics = metrics,
            scores = scores,
            masks = masks,
            calibration = calibrated$calibration,
            calls = calls,
            provenance = provenance
        ),
        class = c("BamAnalysisResult", "list")
    )
    validateBamAnalysisResult(result)
    result
}
