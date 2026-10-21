## Cohort pipeline orchestration (T4.3, formerly T3.3).
##
## `analyzeBamInput()` is the BAM-first entrypoint of the package. It runs the
## cheap checks first (preflight), then indexes, generates windows, extracts
## GC, counts reads, builds the per-sample tables, corrects GC bias and
## normalizes the cohort. Only then does it segment the pairwise difference
## profiles, so every sample shares one normalization baseline. Every
## expensive stage is preceded by a validation that names the BAM, sample or
## contig that failed.
##
## The result is the `BamAnalysisResult` contract: a list with `assembly`,
## `config`, `bams`, `samples`, `windows`, `gc`, `counts`, `rca`, `isa`,
## `normalization` and `provenance`. `counts$donor` and `counts$recipient` are
## integer matrices with one column per sample.

## Synchronized with tools/vendor-gcrank.txt by tests/testthat/test-pipeline-bam.R.
.MaMaMIAgcrankCommit <- "8fdbc32557973c7737496959ecb91acd9bed3ee5"

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
    if (!is.list(config) || is.null(config$assemblies) ||
        length(config$assemblies) == 0L) {
        .analysis_error(
            "`config` must be a validated MaMaMIA configuration with at least one assembly; resolve a species pair first with resolveSpeciesAssembly()"
        )
    }
    if (is.null(config$defaults) || is.null(config$server)) {
        config <- validateMaMaMIAConfig(config)
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

.pipeline_normalization <- function(normalization, n_bams) {
    mode <- tryCatch(
        match.arg(normalization, c("median", "reference", "none")),
        error = function(e) .analysis_error("%s", conditionMessage(e))
    )
    ## There is nothing to normalize a single sample against; the requested
    ## mode is still recorded in provenance.
    if (n_bams == 1L) {
        return("none")
    }
    mode
}

.pipeline_cache <- function(cache) {
    if (is.null(cache)) {
        ## A per-process fit cache in session temp: repeated runs in one Shiny
        ## session reuse the fit, while separate processes never share a file.
        return(file.path(
            .session_temp_dir(),
            sprintf("mamamia-fit-cache-%d.rds", Sys.getpid())
        ))
    }
    if (identical(cache, FALSE)) {
        return(NULL)
    }
    if (!is.character(cache) || length(cache) != 1L || is.na(cache) ||
        !nzchar(cache)) {
        .analysis_error(
            "`cache` must be NULL, FALSE or a single .rds file path"
        )
    }
    cache <- path.expand(cache)
    if (dir.exists(cache)) {
        .analysis_error("`cache` points at a directory: %s", cache)
    }
    cache
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
    n <- nrow(preflight)
    working <- character(n)
    index_path <- character(n)
    staged <- logical(n)

    for (i in seq_len(n)) {
        original <- preflight$path[[i]]
        if (identical(preflight$index_action[[i]], "stage")) {
            working[[i]] <- stageBam(original, .session_temp_dir())
            staged[[i]] <- TRUE
        } else {
            working[[i]] <- original
        }
        index_path[[i]] <- ensureBamIndex(working[[i]])
    }

    list(
        path = working,
        original_path = preflight$path,
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
    header <- system.file("include", "htslib", "hts.h", package = "Rhtslib")
    if (!nzchar(header) || !file.exists(header)) {
        return(NA_character_)
    }
    line <- grep(
        "^#define[[:space:]]+HTS_VERSION[[:space:]]+[0-9]+",
        readLines(header, warn = FALSE),
        value = TRUE
    )
    if (length(line) == 0L) {
        return(NA_character_)
    }
    code <- suppressWarnings(as.numeric(
        sub("^.*[[:space:]]([0-9]+)[[:space:]]*$", "\\1", line[[1L]])
    ))
    if (is.na(code)) {
        return(NA_character_)
    }
    ## htslib encodes 1.10 as 101000 and 1.18 as 101800.
    sprintf(
        "%d.%d.%d",
        code %/% 100000,
        (code %/% 1000) %% 100,
        (code %/% 10) %% 100
    )
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
                                 n_per, mapq, exclude_flags, max_fit_windows,
                                 normalization, cache, started_at,
                                 requested_normalization) {
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
    list(
        package = "MaMaMIA",
        package_version = .safe_package_version("MaMaMIA"),
        r_version = as.character(getRversion()),
        dependency_versions = c(
            bamsignals = .safe_package_version("bamsignals"),
            Rhtslib = .safe_package_version("Rhtslib"),
            GenomicRanges = .safe_package_version("GenomicRanges"),
            DNAcopy = .safe_package_version("DNAcopy"),
            glmmTMB = .safe_package_version("glmmTMB")
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
        window_rule = "exact_equal_windows",
        windows_per_chromosome = n_per,
        n_samples = nrow(preflight),
        total_windows = nrow(windows$table) * nrow(preflight),
        n_windows = nrow(windows$table),
        n_gc_windows = nrow(gc),
        mapq = mapq,
        exclude_flags = exclude_flags,
        max_fit_windows = max_fit_windows,
        max_fit_windows_enforced = TRUE,
        normalization_requested = requested_normalization,
        normalization_mode = normalization$mode,
        normalization_by_side = normalization$by_side,
        pseudocount = normalization$pseudocount,
        count_semantics = "bamsignals_start_based",
        gc_source = "gcrank",
        cache = if (is.null(cache)) NA_character_ else cache,
        started_at = format(started_at, "%Y-%m-%dT%H:%M:%S%z"),
        duration_secs = as.numeric(difftime(
            Sys.time(), started_at, units = "secs"
        ))
    )
}

#' Analyze one or more BAM files with the BAM-first workflow
#'
#' Reads BAM headers, matches them to a configured assembly, ensures CSI
#' indexes, generates exact-N analysis windows, extracts GC content from the
#' configured gcrank indexes, counts reads with [countReads()], builds the
#' per-sample tables with [buildSampleTables()], corrects GC bias with
#' [correctReadCounts()] and segments the pairwise coverage differences with
#' [segments()].
#'
#' For a single BAM the normalization mode is always `"none"`; segmentation
#' therefore consumes `cor.gc` (D7). For several BAMs, `normalization = "median"`
#' (the default) or `"reference"` calls [normalizeReadCounts()] after every
#' sample has been corrected and before any segmentation. `"reference"` requires
#' `reference_sample`; the reference sample's dose is exactly 1.
#'
#' The optional `progress` callback is called as `progress(stage, detail)` after
#' each completed stage. Stages are `"preflight"`, `"index"`, `"windows"`,
#' `"gc"`, `"counts"`, `"correction"` and `"segmentation"`. When a cohort is
#' normalized, one `"normalization"` stage is reported between the corrections
#' and the segmentations. `detail` is a list with stage-specific fields such as
#' `sample`, `n_windows`, `n_ideal` or `n_segments`. An error raised by the
#' callback aborts the analysis, which lets a Shiny session cancel a long run.
#'
#' A warning is emitted when the total analysis window count across samples
#' (`length(bam_paths) * nrow(windows$table)`) exceeds `max_fit_windows`; this
#' is a performance guard, not an error. The fit-window cap itself is applied
#' per sample by [correctReadCounts()] and recorded in provenance and in each
#' `RCA$fit_sampling`.
#'
#' BAMs whose directory is not writable are staged in session temp before
#' indexing; the result records the staging decision in `provenance$staged`.
#'
#' @param bam_paths Character vector of coordinate-sorted BAM paths.
#' @param config A validated configuration list from [readMaMaMIAConfig()] or a
#'   path to a YAML configuration.
#' @param assembly Optional assembly key. Without it, the BAM headers select the
#'   most specific exact match.
#' @param windows_per_chromosome Number of windows per target contig.
#' @param reference_sample Sample used as the denominator for
#'   `normalization = "reference"`. Required in that mode and ignored by
#'   `"median"` and `"none"`.
#' @param normalization Normalization mode across samples: `"median"`
#'   (default), `"reference"` or `"none"`. A single BAM is always analyzed with
#'   `"none"`.
#' @param mapq Minimum mapping quality.
#' @param exclude_flags SAM flag mask; a read is excluded when any bit is set.
#' @param cores Number of cores for the GC-bias fit.
#' @param max_fit_windows Maximum number of ideal windows used for each
#'   sample's GC-bias fit. Larger ideal sets are subsampled deterministically by
#'   subgenome and GC decile with at least ten windows per non-empty stratum;
#'   predictions still cover every window. `Inf` disables subsampling. A cohort
#'   whose total window count exceeds this value also emits a performance
#'   warning.
#' @param cache Optional fit-cache `.rds` path. `NULL` uses a per-process file
#'   in session temp; `FALSE` disables the cache.
#' @param progress Optional `function(stage, detail)` progress callback.
#'
#' @return A `BamAnalysisResult` list with `assembly`, `config`, `bams`,
#'   `samples`, `windows`, `gc`, `counts`, `rca`, `isa`, `normalization` and
#'   `provenance`. `samples` is a data frame with one row per sample (`sample`,
#'   `path`, `assembly`, `index_path`, `reference`, `normalization`, `n_windows`,
#'   `n_segments`, `n_putative_introgressions`).
#' @seealso [preflightBams()], [countReads()], [buildSampleTables()]
#' @export
analyzeBamInput <- function(bam_paths,
                            config,
                            assembly = NULL,
                            windows_per_chromosome = 1000L,
                            reference_sample = NULL,
                            normalization = c("median", "reference", "none"),
                            mapq = 20L,
                            exclude_flags = 2308L,
                            cores = 1L,
                            max_fit_windows = 50000L,
                            cache = NULL,
                            progress = NULL) {
    started_at <- Sys.time()
    config <- .pipeline_config(config)
    .validate_progress(progress)
    n_per <- .pipeline_integer(
        windows_per_chromosome, "windows_per_chromosome"
    )
    mapq <- .pipeline_integer(mapq, "mapq", 0L, 255L)
    exclude_flags <- .pipeline_integer(exclude_flags, "exclude_flags", 0L, 4095L)
    if (length(max_fit_windows) != 1L || is.na(max_fit_windows) ||
        !is.numeric(max_fit_windows) || max_fit_windows <= 0 ||
        (!is.infinite(max_fit_windows) && max_fit_windows != floor(max_fit_windows))) {
        .analysis_error(
            "`max_fit_windows` must be a positive integer or Inf"
        )
    }
    if (is.finite(max_fit_windows) &&
        max_fit_windows > .Machine$integer.max) {
        max_fit_windows <- .Machine$integer.max
    }
    cores <- check_cores(cores)

    if (!is.character(bam_paths) || length(bam_paths) == 0L ||
        anyNA(bam_paths) || any(!nzchar(bam_paths))) {
        .analysis_error(
            "`bam_paths` must be a non-empty character vector of BAM paths"
        )
    }
    requested_mode <- tryCatch(
        match.arg(normalization, c("median", "reference", "none")),
        error = function(e) .analysis_error("%s", conditionMessage(e))
    )
    mode <- .pipeline_normalization(requested_mode, length(bam_paths))
    if (!is.null(reference_sample)) {
        if (!is.character(reference_sample) || length(reference_sample) != 1L ||
            is.na(reference_sample) || !nzchar(reference_sample)) {
            .analysis_error(
                "`reference_sample` must be NULL or a single sample name"
            )
        }
    }
    cache <- .pipeline_cache(cache)

    preflight <- preflightBams(bam_paths, config, assembly)
    assembly_key <- preflight$assembly[[1L]]
    assembly_config <- config$assemblies[[assembly_key]]
    if (is.null(assembly_config)) {
        .analysis_error(
            "Assembly '%s' selected by preflight is not in the configuration",
            assembly_key
        )
    }
    if (identical(mode, "reference") && is.null(reference_sample)) {
        .analysis_error(
            "`normalization = 'reference'` requires `reference_sample`"
        )
    }
    if (!is.null(reference_sample) &&
        !reference_sample %in% preflight$sample) {
        .analysis_error(
            "`reference_sample` '%s' is not one of the input samples: %s",
            reference_sample, paste(preflight$sample, collapse = ", ")
        )
    }
    .pipeline_progress(progress, "preflight", list(
        n_bams = nrow(preflight),
        assembly = assembly_key
    ))

    indexed <- .pipeline_index(preflight)
    .pipeline_progress(progress, "index", list(
        n_bams = nrow(preflight),
        staged = sum(indexed$staged)
    ))

    contigs <- .pipeline_contigs(assembly_config)
    windows <- makeAnalysisWindows(assembly_config, contigs, n_per)
    .pipeline_progress(progress, "windows", list(
        n_per_chromosome = n_per,
        n_windows = nrow(windows$table)
    ))

    total_windows <- nrow(preflight) * nrow(windows$table)
    if (is.finite(max_fit_windows) && total_windows > max_fit_windows) {
        warning(
            sprintf(
                paste0(
                    "Total analysis windows across %d sample(s) is %d, above ",
                    "max_fit_windows (%s); the run may be slow and memory-hungry"
                ),
                nrow(preflight), total_windows,
                format(max_fit_windows, scientific = FALSE)
            ),
            call. = FALSE
        )
    }

    gc <- extractWindowGc(
        windows,
        assembly_config,
        mmap = config$defaults$gc_mmap,
        validate_crc = config$defaults$gc_validate_crc
    )
    .pipeline_progress(progress, "gc", list(
        n_windows = nrow(gc),
        n_indexes = length(assembly_config$gc_indexes)
    ))

    samples <- as.character(preflight$sample)
    n_donor <- sum(windows$table$role == "donor")
    n_recipient <- sum(windows$table$role == "recipient")
    donor_counts <- matrix(
        0L,
        nrow = n_donor, ncol = length(samples),
        dimnames = list(NULL, samples)
    )
    recipient_counts <- matrix(
        0L,
        nrow = n_recipient, ncol = length(samples),
        dimnames = list(NULL, samples)
    )

    rcas <- vector("list", length(samples))
    names(rcas) <- samples
    tables_list <- vector("list", length(samples))

    for (i in seq_along(samples)) {
        sample <- samples[[i]]
        original_path <- preflight$path[[i]]
        counts <- tryCatch(
            countReads(
                indexed$path[[i]],
                windows$gr,
                mapq = mapq,
                exclude_flags = exclude_flags
            ),
            error = function(e) .sample_stage_error(e, sample, "counts")
        )
        total <- sum(as.numeric(counts))
        if (total <= 0) {
            .mamamia_error(
                "MaMaMIASampleError",
                "Sample '%s' has zero reads in all analysis windows; cannot estimate GC bias",
                sample
            )
        }
        donor_counts[, i] <- counts[windows$table$role == "donor"]
        recipient_counts[, i] <- counts[windows$table$role == "recipient"]
        .pipeline_progress(progress, "counts", list(
            sample = sample,
            path = original_path,
            total = total
        ))

        tables <- buildSampleTables(
            counts, gc, windows, assembly_config, sample
        )
        tables_list[[i]] <- tables
        rca <- tryCatch(
            do.call(RCA, tables[c("dcov", "dgc", "rcov", "rgc", "meta")]),
            error = function(e) .sample_stage_error(e, sample, "RCA")
        )
        corrected <- tryCatch(
            correctReadCounts(
                rca,
                cores = cores,
                verbose = FALSE,
                cache = cache,
                max_fit_windows = max_fit_windows
            ),
            error = function(e) {
                .sample_stage_error(e, sample, "GC correction")
            }
        )
        .pipeline_progress(progress, "correction", list(
            sample = sample,
            n_ideal = sum(corrected$data$ideal)
        ))

        rcas[[sample]] <- corrected
    }

    ## Normalize only after every sample has been corrected, and before any
    ## CBS (D6). `cor.gc` is preserved; the dose lives in `RCA$data$dose`.
    normalization_out <- list(
        mode = mode,
        reference = NULL,
        by_side = config$defaults$normalization_by_side,
        pseudocount = 0.5
    )
    if (!identical(mode, "none")) {
        reference_used <- if (identical(mode, "reference")) {
            reference_sample
        } else {
            NULL
        }
        rcas <- tryCatch(
            normalizeReadCounts(
                rcas,
                reference = reference_used,
                by_side = config$defaults$normalization_by_side,
                pseudocount = 0.5
            ),
            error = function(e) {
                if (inherits(e, "MaMaMIAError")) {
                    stop(e)
                }
                .analysis_error(
                    "Normalization failed: %s", conditionMessage(e)
                )
            }
        )
        normalization_out <- list(
            mode = mode,
            reference = reference_used,
            by_side = config$defaults$normalization_by_side,
            pseudocount = 0.5
        )
        .pipeline_progress(progress, "normalization", list(
            mode = mode,
            reference = reference_used,
            n_samples = length(samples)
        ))
    }

    ## All samples share one assembly, so their target pairs must agree. Using
    ## a single vector guarantees that every sample is segmented on the same
    ## donor-recipient pairs.
    target_pairs <- tables_list[[1L]]$target_pairs
    if (length(tables_list) > 1L) {
        for (i in seq_along(tables_list)[-1L]) {
            if (!identical(tables_list[[i]]$target_pairs, target_pairs)) {
                .analysis_error(
                    "Samples in one run must share the same donor-recipient target pairs"
                )
            }
        }
    }

    isas <- vector("list", length(samples))
    names(isas) <- samples
    for (i in seq_along(samples)) {
        sample <- samples[[i]]
        isa <- tryCatch(
            segments(
                rcas[[i]],
                target_pairs = target_pairs,
                verbose = FALSE
            ),
            error = function(e) .sample_stage_error(e, sample, "segmentation")
        )
        .pipeline_progress(progress, "segmentation", list(
            sample = sample,
            n_segments = nrow(isa$segments),
            value_col = isa$meta$value_col
        ))
        isas[[sample]] <- isa
    }

    reference_used <- normalization_out$reference
    is_reference <- if (is.null(reference_used)) {
        rep(FALSE, length(samples))
    } else {
        samples == reference_used
    }
    samples_meta <- data.frame(
        sample = samples,
        path = preflight$path,
        assembly = preflight$assembly,
        index_path = indexed$index_path,
        reference = is_reference,
        normalization = rep(normalization_out$mode, length(samples)),
        n_windows = rep(nrow(windows$table), length(samples)),
        n_segments = vapply(isas, function(x) nrow(x$segments), integer(1L)),
        n_putative_introgressions = vapply(
            isas,
            function(x) sum(x$segments$putative_introgression, na.rm = TRUE),
            integer(1L)
        ),
        stringsAsFactors = FALSE
    )

    bams <- data.frame(
        path = preflight$path,
        sample = samples,
        index_path = indexed$index_path,
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
        max_fit_windows = max_fit_windows,
        normalization = normalization_out,
        cache = cache,
        started_at = started_at,
        requested_normalization = requested_mode
    )

    structure(
        list(
            assembly = assembly_config,
            config = config,
            bams = bams,
            samples = samples_meta,
            windows = windows,
            gc = gc,
            counts = list(
                donor = donor_counts,
                recipient = recipient_counts
            ),
            rca = rcas,
            isa = isas,
            normalization = normalization_out,
            provenance = provenance
        ),
        class = c("BamAnalysisResult", "list")
    )
}
