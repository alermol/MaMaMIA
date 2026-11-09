## Multi-sample result bundle: one ZIP with a folder per sample and
## run-wide summaries at the top level.

#' Sample names in a result, in natural order
#' @noRd
bam_bundle_samples <- function(result) {
    samples <- if (is.data.frame(result$samples) &&
        "sample" %in% names(result$samples)) {
        as.character(result$samples$sample)
    } else if (is.data.frame(result$calls) &&
        "sample" %in% names(result$calls)) {
        as.character(result$calls$sample)
    } else {
        character()
    }
    mamamia_natural_sort(unique(samples[!is.na(samples) & nzchar(samples)]))
}

#' Safe, unique folder names for the samples in a cohort bundle
#' @noRd
bam_bundle_sample_folders <- function(samples) {
    samples <- as.character(samples)
    stems <- vapply(samples, bam_bundle_stem, character(1))
    ## Sample folders sit next to the run-wide files, so they must not
    ## shadow a top-level entry such as run/, qc/ or report.html.
    reserved <- c(
        "README.txt", "report.html", "cohort-summary.csv",
        "calls-all.csv", "calls-passing.csv", "calls-confident.csv",
        "calls-summary.csv", "calls-by-status.csv",
        "calls-by-copy-state.csv", "calls-all.bedpe",
        "calls-passing.bedpe", "calls-confident.bedpe",
        "manifest.json", "checksums.md5", "run", "qc", "tracks", "samples"
    )
    used <- tolower(reserved)
    out <- character(length(stems))
    for (i in seq_along(stems)) {
        candidate <- stems[[i]]
        suffix <- 1L
        while (tolower(candidate) %in% used) {
            suffix <- suffix + 1L
            candidate <- paste0(stems[[i]], "-", suffix)
        }
        out[[i]] <- candidate
        used <- c(used, tolower(candidate))
    }
    stats::setNames(out, samples)
}

#' Write the run-wide QC tables for a cohort bundle
#' @noRd
bam_bundle_write_cohort_qc <- function(result, bundle) {
    qc_dir <- file.path(bundle, "qc")
    dir.create(qc_dir, showWarnings = FALSE)
    masks <- bam_bundle_rename_columns(
        result$masks, bam_bundle_mask_columns()
    )
    calibration <- bam_bundle_rename_columns(
        result$calibration, bam_bundle_calibration_columns()
    )
    bam_bundle_validate_qc_tables(
        list(masks = masks, calibration = calibration),
        required_tables = character()
    )
    masks_summary <- bam_bundle_masks_summary(masks)
    bam_bundle_write_csv(
        masks, file.path(qc_dir, "masks-windows.csv.gz"), gzip = TRUE
    )
    bam_bundle_write_csv(
        masks_summary, file.path(qc_dir, "masks-summary.csv")
    )
    bam_bundle_write_csv(
        calibration, file.path(qc_dir, "calibration.csv")
    )
    list(
        masks = masks,
        calibration = calibration,
        masks_summary = masks_summary
    )
}

#' Write cohort-level plots for a multi-sample bundle
#' @noRd
bam_bundle_cohort_plot_files <- function(bundle, cohort, calls_summary,
                                         masks_summary = NULL) {
    plot_dir <- file.path(bundle, "qc", "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    specs <- list(
        list(
            file = "cohort-calls.png",
            title = "Calls per sample in this run",
            caption = paste(
                "All calls and confident calls for every sample.",
                "Raw values: cohort-summary.csv."
            ),
            plot = bam_bundle_plot_cohort(cohort),
            width = bam_bundle_cohort_plot_width(nrow(cohort)),
            height = 5.2
        ),
        list(
            file = "calls-by-class.png",
            title = "Calls by class across the run",
            caption = paste(
                "All calls, passing calls and confident calls by class.",
                "Raw values: calls-summary.csv."
            ),
            plot = bam_bundle_plot_calls(calls_summary),
            width = 9, height = 5.5
        ),
        list(
            file = "masked-fraction-by-group.png",
            title = "Masked windows by chromosome group",
            caption = paste(
                "Masked fraction per group across the run.",
                "Raw values: qc/masks-summary.csv."
            ),
            plot = bam_bundle_plot_masks(masks_summary),
            width = 9, height = 5.0
        )
    )
    rows <- lapply(specs, function(spec) {
        path <- file.path(plot_dir, spec$file)
        if (!bam_bundle_save_plot(
            spec$plot, path, width = spec$width, height = spec$height
        )) {
            return(NULL)
        }
        data.frame(
            path = file.path("qc", "plots", spec$file),
            title = spec$title,
            caption = spec$caption,
            stringsAsFactors = FALSE
        )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) {
        return(data.frame(
            path = character(), title = character(), caption = character(),
            stringsAsFactors = FALSE
        ))
    }
    do.call(rbind, rows)
}

#' HTML table with pre-rendered cells
#' @noRd
bam_bundle_html_table_cells <- function(x, raw_columns = character(),
                                        descriptions = NULL) {
    if (!is.data.frame(x) || nrow(x) == 0L || ncol(x) == 0L) {
        return("<p class=\"muted\">No data for this section.</p>")
    }
    rows <- vapply(seq_len(nrow(x)), function(i) {
        cells <- vapply(names(x), function(column) {
            value <- x[[column]][[i]]
            if (column %in% raw_columns) {
                as.character(value)
            } else {
                bam_bundle_html_value(value)
            }
        }, character(1))
        paste0(
            "<tr><td>", paste(cells, collapse = "</td><td>"),
            "</td></tr>"
        )
    }, character(1))
    bam_bundle_html_table_shell(x, rows, descriptions)
}

#' Self-contained HTML report for a cohort bundle
#' @noRd
bam_bundle_cohort_report_html <- function(result, bundle, samples, folders,
                                          calls_raw, calls_summary, cohort,
                                          plots, generated, version) {
    provenance <- result$provenance
    assembly_label <- if (!is.null(result$assembly$label)) {
        as.character(result$assembly$label)
    } else {
        as.character(result$assembly$key)
    }
    bam_paths <- if (is.data.frame(result$bams) &&
        "path" %in% names(result$bams)) {
        unique(as.character(result$bams$path))
    } else {
        character()
    }
    bam_paths <- bam_paths[!is.na(bam_paths) & nzchar(bam_paths)]
    passing <- sum(isTRUE_vec(calls_raw$pass), na.rm = TRUE)
    confident <- nrow(bam_confident_calls(calls_raw))
    overview <- data.frame(
        Item = c(
            "Samples", "Assembly", "BAM", "Analysis profile",
            "Calls", "Passing calibration", "Confident calls",
            "Duration", "Config hash", "Calibration version",
            "Generated", "MaMaMIA"
        ),
        Value = c(
            paste(samples, collapse = ", "),
            assembly_label,
            paste(bam_paths, collapse = ", "),
            if (!is.null(provenance$analysis_profile) &&
                length(provenance$analysis_profile) == 1L &&
                !is.na(provenance$analysis_profile) &&
                nzchar(provenance$analysis_profile)) {
                provenance$analysis_profile
            } else {
                NA_character_
            },
            format(nrow(calls_raw), big.mark = ",", trim = TRUE),
            format(passing, big.mark = ",", trim = TRUE),
            format(confident, big.mark = ",", trim = TRUE),
            if (!is.null(provenance$duration_secs)) {
                paste0(format(provenance$duration_secs, digits = 4,
                    trim = TRUE), " s")
            } else {
                NA_character_
            },
            bam_bundle_text(provenance$config_hash),
            bam_bundle_text(provenance$calibration_version),
            format(generated, "%Y-%m-%d %H:%M:%S %Z"),
            version
        ),
        stringsAsFactors = FALSE
    )

    summary_columns <- intersect(
        c(
            "sample", "mode", "n_windows", "reads_total", "n_all_calls",
            "n_high_confidence", "pass_rate", "masked_fraction",
            "total_width_mb"
        ),
        names(cohort)
    )
    cohort_table <- cohort[, summary_columns, drop = FALSE]

    sample_table <- data.frame(
        Sample = samples,
        Report = sprintf(
            "<a href=\"%s/report.html\">%s/report.html</a>",
            folders, folders
        ),
        Calls = vapply(samples, function(sample) {
            row <- cohort[as.character(cohort$sample) == sample, ,
                drop = FALSE]
            if (nrow(row) > 0L && "n_all_calls" %in% names(row)) {
                row$n_all_calls[[1L]]
            } else {
                NA_integer_
            }
        }, integer(1)),
        Passing = vapply(samples, function(sample) {
            calls <- calls_raw[as.character(calls_raw$sample) == sample, ,
                drop = FALSE]
            sum(isTRUE_vec(calls$pass), na.rm = TRUE)
        }, integer(1)),
        Confident = vapply(samples, function(sample) {
            calls <- calls_raw[as.character(calls_raw$sample) == sample, ,
                drop = FALSE]
            nrow(bam_confident_calls(calls))
        }, integer(1)),
        stringsAsFactors = FALSE
    )

    c(
        "<!DOCTYPE html>",
        "<html lang=\"en\">",
        "<head>",
        "<meta charset=\"utf-8\">",
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
        paste0(
            "<title>MaMaMIA ",
            if (length(samples) > 1L) "cohort " else "",
            "results</title>"
        ),
        paste0("<style>", bam_bundle_report_css(), "</style>"),
        "</head>",
        "<body><main>",
        "<header>",
        paste0(
            "<p class=\"eyebrow\">MaMaMIA results bundle",
            if (length(samples) > 1L) " for a cohort" else "",
            "</p>"
        ),
        paste0(
            "<h1>",
            if (length(samples) > 1L) "Cohort results" else "Results",
            "</h1>"
        ),
        paste0("<p class=\"lead\">", bam_bundle_html_escape(
            paste(samples, collapse = ", ")
        ), "</p>"),
        "</header>",
        "<section><h2>Run overview</h2>",
        bam_bundle_html_table(overview),
        "</section>",
        "<section><h2>Sample summary</h2>",
        bam_bundle_html_table_cells(
            cohort_table, raw_columns = character()
        ),
        "</section>",
        paste0(
            "<section><h2>",
            if (length(samples) > 1L) "Per-sample bundles" else "Sample bundle",
            "</h2>"
        ),
        bam_bundle_html_table_cells(sample_table, raw_columns = "Report"),
        "<p class=\"muted\">Each linked folder contains the sample-specific ",
        "call tables, BEDPE files, karyotypes, tracks and QC. Run-wide ",
        "diagnostics live at the archive root.</p>",
        "</section>",
        "<section><h2>Calls by class</h2>",
        bam_bundle_plot_html(plots, bundle, "calls-by-class.png"),
        bam_bundle_html_table(calls_summary),
        "</section>",
        "<section><h2>Calls per sample</h2>",
        bam_bundle_plot_html(
            plots, bundle, "cohort-calls.png", wide = TRUE
        ),
        "</section>",
        "<section><h2>Masked windows</h2>",
        bam_bundle_plot_html(
            plots, bundle, "masked-fraction-by-group.png"
        ),
        "</section>",
        paste0(
            "<footer>Generated ", bam_bundle_html_escape(
                format(generated, "%Y-%m-%d %H:%M:%S %Z")
            ),
            " by MaMaMIA ", bam_bundle_html_escape(version), ".</footer>"
        ),
        "</main>",
        bam_bundle_report_script(),
        "</body></html>"
    )
}

#' README for a run bundle
#' @noRd
bam_bundle_cohort_readme <- function(samples, folders, calls,
                                     confident_calls, generated, version) {
    timestamp <- if (inherits(generated, "POSIXt")) {
        format(generated, "%Y-%m-%d %H:%M:%S %Z")
    } else {
        as.character(generated)[[1L]]
    }
    n <- nrow(calls)
    passing <- if ("passes_calibration" %in% names(calls)) {
        sum(isTRUE_vec(calls$passes_calibration), na.rm = TRUE)
    } else {
        n
    }
    confident_n <- if (is.null(confident_calls)) 0L else
        nrow(confident_calls)
    sample_lines <- vapply(seq_along(samples), function(i) {
        paste0(folders[[i]], "/  ", samples[[i]])
    }, character(1))
    title <- if (length(samples) > 1L) {
        "MaMaMIA cohort results bundle"
    } else {
        "MaMaMIA results bundle"
    }
    c(
        title,
        "",
        paste0("Samples: ", paste(samples, collapse = ", ")),
        paste0("Generated: ", timestamp),
        paste0("Version: MaMaMIA ", version),
        paste0(
            "Calls: ", n, " total, ", passing,
            " passing calibration, ", confident_n, " confident"
        ),
        "",
        "Layout",
        "------",
        "cohort-summary.csv      One row per sample: reads, windows, calls,",
        "                        pass rate and mask burden.",
        "calls-all.csv           Every call from every sample.",
        "calls-passing.csv       Calls that pass calibration, all samples.",
        "calls-confident.csv     Passing calls with review_status = confident.",
        "calls-summary.csv       Calls, passing, confident, width and score by",
        "                        call class across the run.",
        "calls-by-status.csv     Calls by review status across the run.",
        "calls-by-copy-state.csv Calls by copy-number state across the run.",
        "calls-all.bedpe         Every call in BEDPE format.",
        "calls-passing.bedpe     Calibrated calls in BEDPE format.",
        "calls-confident.bedpe   Confident calls in BEDPE format.",
        "report.html             Self-contained run summary.",
        "manifest.json           File sizes and MD5 checksums.",
        "checksums.md5           MD5 checksums including the manifest.",
        "",
        "Per-sample folders",
        "------------------",
        sample_lines,
        "",
        "Sample folders sit at the archive root; the files listed above are",
        "run-wide. Each sample folder contains calls-full.csv,",
        "calls-passing.csv, calls-confident.csv, calls-review.csv, the call",
        "summaries, BEDPE",
        "files, karyotype SVGs, tracks, sample-specific qc/ tables and plots,",
        "report.html and README.txt. The sample README repeats the full",
        "column dictionary.",
        "",
        "Run-wide files",
        "--------------",
        "run/                    Provenance, config, parameters, versions,",
        "                        assembly, sample and BAM metadata.",
        "qc/calibration.csv      Calibration thresholds for the run.",
        "qc/masks-summary.csv    Masked windows by chromosome group.",
        "qc/masks-windows.csv.gz Every mask decision.",
        "tracks/masks.bed        Masked windows in group coordinates.",
        "qc/plots/               Run-level call and mask figures.",
        "",
        "Confident set",
        "-------------",
        "The confident set contains calls that pass calibration and have",
        "review_status = confident. The full CSV and BEDPE keep every call.",
        "",
        "Diagnostic conventions",
        "----------------------",
        "Sample qc/metrics-windows.csv.gz has one row per window and role.",
        "Sample qc/scores-windows.csv.gz has one row per group window.",
        "Both files are compressed CSV with a header and descriptive column",
        "names. qc/masks-windows.csv.gz keeps every mask decision, including",
        "windows that were not masked.",
        "",
        "BED tracks are zero-based, half-open intervals. BED scores are zero",
        "because the call class and review status are carried in the name",
        "column. masks.bed uses synthetic contig names group_<group> because",
        "masks are computed in homology-group coordinates. breadth.bedGraph.gz",
        "has one value per window for the sample.",
        "",
        "Column conventions",
        "------------------",
        "Columns ending in _bp are base pairs. Fractions, breadths and",
        "overlaps are between 0 and 1. Residuals and scores are on the",
        "log scale unless the description says otherwise. A blank value",
        "means unavailable or not applicable.",
        "",
        "BEDPE files are tab-separated. They use the same descriptive names,",
        "and calls-passing.bedpe contains only rows where",
        "passes_calibration is TRUE. The combined CSV and BEDPE columns",
        "include sample_name.",
        ""
    )
}

#' Worker budget for writing a result bundle
#'
#' Per-sample export is independent, so it is forked on Unix when more than one
#' worker is allowed. `options(MaMaMIA.bundle_workers = n)` overrides the
#' default; `MaMaMIA.cores` is honoured as a fallback. The default is every
#' available core except one, capped by the container memory limit so the
#' workers do not swap, and the same budget drives 7-Zip compression.
#' @noRd
bam_bundle_workers <- function() {
    value <- getOption("MaMaMIA.bundle_workers", NULL)
    if (is.null(value)) {
        value <- getOption("MaMaMIA.cores", NULL)
    }
    if (is.null(value) || length(value) != 1L || is.na(value) ||
        !is.numeric(value) || value < 1L) {
        value <- min(
            mamamia_default_cores(),
            mamamia_memory_worker_cap()
        )
    }
    value <- as.integer(value)
    if (is.na(value) || value < 1L) 1L else value
}

#' Write a run result bundle
#'
#' Builds one ZIP with top-level run summaries and one folder per sample.
#' Sample-specific calls, BEDPE files, tracks, QC and karyotypes live under
#' `<sample>/`; run provenance and run-wide mask/calibration tables live once
#' at the top level.
#'
#' @param result A `BamAnalysisResult`.
#' @param file Output ZIP path.
#' @param samples Optional sample names; defaults to all samples in `result`.
#' @param current_sample Sample whose filtered payload is supplied.
#' @param current_payload Optional filtered karyotype payload for
#'   `current_sample`.
#' @param full_payload Optional all-call payload for `current_sample`.
#' @param confident_payload Optional confident-call payload for
#'   `current_sample`.
#' @param generated Timestamp shown in the README and reports.
#'
#' @return The output path, invisibly.
#' @noRd
bam_write_cohort_bundle <- function(result, file, samples = NULL,
                                    current_sample = NULL,
                                    current_payload = NULL,
                                    full_payload = NULL,
                                    confident_payload = NULL,
                                    generated = Sys.time()) {
    if (!inherits(result, "BamAnalysisResult")) {
        stop("`result` must be a BamAnalysisResult", call. = FALSE)
    }
    if (!is.character(file) || length(file) != 1L || is.na(file) ||
        !nzchar(file)) {
        stop("`file` must be a single output path", call. = FALSE)
    }
    if (is.null(samples)) {
        samples <- bam_bundle_samples(result)
    }
    samples <- unique(as.character(samples))
    samples <- samples[!is.na(samples) & nzchar(samples)]
    if (length(samples) == 0L) {
        stop("`result` contains no samples to export", call. = FALSE)
    }
    folders <- bam_bundle_sample_folders(samples)
    version <- bam_bundle_version()
    workers <- bam_bundle_workers()
    memory_limit <- mamamia_memory_limit()
    started_at <- proc.time()
    message(sprintf(
        "Result bundle: %d samples, %d workers, memory limit %s",
        length(samples), workers,
        if (is.finite(memory_limit)) {
            sprintf("%.1f GiB", memory_limit / 1024^3)
        } else {
            "unlimited"
        }
    ))
    on.exit(
        message(sprintf(
            "Result bundle finished in %.1f s",
            (proc.time() - started_at)[["elapsed"]]
        )),
        add = TRUE
    )
    staging <- tempfile("mamamia-cohort-bundle-")
    dir.create(staging)
    on.exit(unlink(staging, recursive = TRUE, force = TRUE), add = TRUE)
    bundle_name <- "mamamia-results"
    bundle <- file.path(staging, bundle_name)
    dir.create(bundle)

    ## Compute run-wide tables once. The per-sample loop used to recompute the
    ## cohort summary and rescan every window table for each sample.
    cohort <- bam_cohort_summary(result)
    cohort <- cohort[
        as.character(cohort$sample) %in% samples, , drop = FALSE
    ]
    cohort <- cohort[
        order(match(as.character(cohort$sample), samples)), , drop = FALSE
    ]
    row_index <- function(x) {
        if (!is.data.frame(x) || !"sample" %in% names(x) ||
            nrow(x) == 0L) {
            return(NULL)
        }
        split(
            seq_len(nrow(x)),
            factor(as.character(x$sample), levels = samples)
        )
    }
    metrics_rows <- row_index(result$metrics)
    scores_rows <- row_index(result$scores)
    export_sample <- function(i) {
        sample <- samples[[i]]
        sample_calls <- result$calls[
            as.character(result$calls$sample) == sample, , drop = FALSE
        ]
        sample_current <- if (identical(sample, current_sample) &&
            !is.null(current_payload)) {
            current_payload
        } else {
            karyotype_export_payload(
                result, sample = sample, mode = "calls",
                calls = sample_calls[isTRUE_vec(sample_calls$pass), ,
                    drop = FALSE],
                selected_groups = NULL
            )
        }
        sample_full <- if (identical(sample, current_sample) &&
            !is.null(full_payload)) {
            full_payload
        } else {
            bam_full_karyotype_payload(sample_current)
        }
        sample_confident_calls <- bam_confident_calls(sample_calls)
        sample_confident <- if (identical(sample, current_sample) &&
            !is.null(confident_payload)) {
            confident_payload
        } else {
            karyotype_payload(
                result, sample = sample, mode = "calls",
                calls = sample_confident_calls,
                selected_groups = NULL
            )
        }
        sample_actual <- karyotype_payload(
            result, sample = sample, mode = "calls",
            calls = sample_confident_calls,
            selected_groups = NULL,
            hide_deletions = TRUE
        )
        sample_dir <- file.path(bundle, folders[[i]])
        dir.create(sample_dir)
        parts <- bam_bundle_write_sample_files(
            result, sample_dir, sample,
            current_payload = sample_current,
            full_payload = sample_full,
            confident_payload = sample_confident,
            actual_payload = sample_actual,
            cohort = cohort,
            metrics_rows = if (is.null(metrics_rows)) {
                NULL
            } else {
                metrics_rows[[i]]
            },
            scores_rows = if (is.null(scores_rows)) {
                NULL
            } else {
                scores_rows[[i]]
            },
            include_run = FALSE,
            include_masks = FALSE,
            include_calibration = FALSE,
            include_cohort = FALSE,
            include_cohort_plot = FALSE,
            include_cohort_section = FALSE,
            write_report = TRUE,
            generated = generated
        )
        writeLines(
            bam_bundle_readme(
                sample, parts$calls, generated, version,
                bedpe_columns = parts$bedpe_columns,
                confident_calls = parts$confident_calls,
                result = result,
                bundle = sample_dir,
                include_run = FALSE,
                include_shared = FALSE,
                include_manifest = FALSE
            ),
            file.path(sample_dir, "README.txt"),
            useBytes = TRUE
        )
        NULL
    }
    pieces <- .mamamia_lapply(
        seq_along(samples), export_sample, workers = bam_bundle_workers()
    )
    failed <- vapply(pieces, function(x) inherits(x, "try-error"), logical(1))
    if (any(failed)) {
        stop(pieces[[which(failed)[[1L]]]])
    }

    calls_raw <- result$calls[
        as.character(result$calls$sample) %in% samples, , drop = FALSE
    ]
    calls_raw <- calls_raw[
        order(match(as.character(calls_raw$sample), samples)),
        , drop = FALSE
    ]
    calls <- bam_bundle_descriptive_calls(calls_raw)
    passing_calls <- bam_bundle_descriptive_calls(
        explorer_calls(calls_raw, list(pass_only = TRUE))$calls
    )
    confident_raw <- bam_confident_calls(calls_raw)
    confident_calls <- bam_bundle_descriptive_calls(confident_raw)
    calls_summary <- bam_bundle_call_summary(calls_raw)
    calls_by_status <- bam_bundle_call_status_summary(calls_raw)
    calls_by_copy_state <- bam_bundle_call_status_summary(
        calls_raw, column = "copy_state", label = "copy_number_state"
    )
    bedpe_calls <- if (nrow(calls_raw) > 0L &&
        !"chr_id.don" %in% names(calls_raw)) {
        .annotate_call_intervals(calls_raw, result$windows)
    } else {
        calls_raw
    }
    bedpe_confident <- bam_confident_calls(bedpe_calls)

    utils::write.csv(
        calls, file.path(bundle, "calls-all.csv"),
        row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
    utils::write.csv(
        passing_calls,
        file.path(bundle, "calls-passing.csv"),
        row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
    utils::write.csv(
        confident_calls, file.path(bundle, "calls-confident.csv"),
        row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
    bam_bundle_write_csv(
        calls_summary, file.path(bundle, "calls-summary.csv")
    )
    bam_bundle_write_csv(
        calls_by_status, file.path(bundle, "calls-by-status.csv")
    )
    bam_bundle_write_csv(
        calls_by_copy_state, file.path(bundle, "calls-by-copy-state.csv")
    )
    bam_bundle_write_csv(
        cohort, file.path(bundle, "cohort-summary.csv")
    )
    bam_bundle_write_bedpe(
        bedpe_calls, result$windows,
        file.path(bundle, "calls-all.bedpe"),
        pass_only = FALSE
    )
    bam_bundle_write_bedpe(
        bedpe_calls, result$windows,
        file.path(bundle, "calls-passing.bedpe"),
        pass_only = TRUE
    )
    bam_bundle_write_bedpe(
        bedpe_confident, result$windows,
        file.path(bundle, "calls-confident.bedpe"),
        pass_only = FALSE
    )
    qc <- bam_bundle_write_cohort_qc(result, bundle)
    tracks_dir <- file.path(bundle, "tracks")
    dir.create(tracks_dir, showWarnings = FALSE)
    bam_bundle_write_masks_bed(
        qc$masks, file.path(tracks_dir, "masks.bed")
    )
    plots <- bam_bundle_cohort_plot_files(
        bundle, cohort, calls_summary, qc$masks_summary
    )
    bam_bundle_run_files(result, bundle, samples[[1L]])

    writeLines(
        bam_bundle_cohort_readme(
            samples, folders, calls, confident_calls, generated, version
        ),
        file.path(bundle, "README.txt"),
        useBytes = TRUE
    )
    writeLines(
        bam_bundle_cohort_report_html(
            result, bundle, samples, folders, calls_raw, calls_summary,
            cohort, plots, generated, version
        ),
        file.path(bundle, "report.html"),
        useBytes = TRUE
    )
    manifest <- bam_bundle_manifest(
        bundle,
        if (length(samples) == 1L) samples[[1L]] else "cohort",
        generated, version, samples = samples
    )
    bam_bundle_checksums(bundle, manifest)
    bam_bundle_archive(staging, bundle_name, file)
}
