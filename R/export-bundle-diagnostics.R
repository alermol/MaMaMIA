## Diagnostic tables, genome-browser tracks and reproducibility files for the
## downloadable result bundle. These helpers keep the bundle self-describing:
## every window, score, mask and calibration row that supports the calls can
## be inspected without reopening the dashboard.

## Descriptive names for the per-window metric table.
bam_bundle_metric_columns <- function() {
    c(
        sample = "sample_name",
        chr_id = "contig_id",
        chr_name = "contig_name",
        group = "chromosome_group",
        role = "role",
        component = "component",
        copy_number = "copy_number",
        iid = "window_index",
        start = "start_bp",
        stop = "end_bp",
        length = "window_length_bp",
        count = "read_count",
        covered = "covered_bases",
        mapq_high = "high_mapq_read_count",
        mapq_high_fraction = "high_mapq_read_fraction",
        breadth = "breadth",
        gc = "gc_fraction",
        lambda_breadth = "expected_breadth",
        lambda_count = "expected_read_count",
        read_length = "read_length_bp"
    )
}

## Descriptive names for the per-window score table.
bam_bundle_score_columns <- function() {
    c(
        sample = "sample_name",
        group = "chromosome_group",
        iid = "window_index",
        start = "start_bp",
        stop = "end_bp",
        mode = "analysis_mode",
        metric = "coverage_metric",
        score_type = "score_type",
        score = "introgression_score",
        donor_mean = "mean_donor_coverage",
        recipient_mean = "mean_recipient_coverage",
        donor_component = "donor_component",
        n_donor = "n_donor_contigs",
        n_recipient = "n_recipient_contigs",
        donor_copies = "estimated_donor_copies",
        donor_fraction = "donor_coverage_fraction",
        donor_component_spread = "donor_component_spread",
        recipient_component_spread = "recipient_component_spread",
        n_donor_supported = "supporting_donor_contigs",
        mapq_fraction = "high_mapq_read_fraction",
        chr_name = "group_name",
        length = "group_length_bp",
        n_per_chromosome = "windows_per_chromosome",
        gc_donor = "donor_gc_fraction",
        gc_recipient = "recipient_gc_fraction",
        log_donor = "log_donor_coverage",
        log_recipient = "log_recipient_coverage",
        log_total = "log_total_coverage",
        donor_residual = "donor_residual",
        recipient_residual = "recipient_residual",
        total_residual = "total_depth_residual",
        lambda_copy = "expected_copy_coverage",
        score_raw = "raw_introgression_score"
    )
}

## Descriptive names for the mask table.
bam_bundle_mask_columns <- function() {
    c(
        group = "chromosome_group",
        iid = "window_index",
        start = "start_bp",
        stop = "end_bp",
        f = "mask_signal",
        posterior = "mask_posterior",
        mask = "masked",
        source = "mask_source",
        version = "mask_version"
    )
}

## Descriptive names for the calibration table.
bam_bundle_calibration_columns <- function() {
    c(
        mode = "analysis_mode",
        length_bin = "call_length_bin",
        n_control = "control_calls",
        q99 = "null_q99",
        threshold = "calibration_threshold",
        fdr = "target_fdr",
        version = "calibration_version"
    )
}

bam_bundle_rename_columns <- function(x, mapping) {
    if (!is.data.frame(x) || ncol(x) == 0L) {
        return(x)
    }
    current <- names(x)
    hit <- match(current, names(mapping))
    current[!is.na(hit)] <- unname(mapping[hit[!is.na(hit)]])
    names(x) <- current
    x
}

bam_bundle_group_order <- function(x) {
    values <- unique(as.character(x))
    values <- values[!is.na(values) & nzchar(values)]
    if (length(values) > 0L && all(grepl("^[0-9]+$", values))) {
        return(values[order(as.integer(values))])
    }
    sort(values)
}

bam_bundle_text <- function(x) {
    if (is.null(x) || length(x) == 0L) {
        return(NA_character_)
    }
    as.character(x[[1L]])
}

bam_bundle_run_summary_lines <- function(result, sample) {
    if (is.null(result)) {
        return(character())
    }
    provenance <- result$provenance
    sample_meta <- result$samples[
        as.character(result$samples$sample) == sample, , drop = FALSE
    ]
    n_windows <- if ("n_windows" %in% names(sample_meta) &&
        nrow(sample_meta) > 0L) {
        sample_meta$n_windows[[1L]]
    } else {
        nrow(result$windows$table)
    }
    bam_path <- if ("path" %in% names(sample_meta) &&
        nrow(sample_meta) > 0L) {
        as.character(sample_meta$path[[1L]])
    } else {
        NA_character_
    }
    assembly <- if (!is.null(result$assembly$label)) {
        as.character(result$assembly$label)
    } else {
        as.character(result$assembly$key)
    }
    c(
        paste0("Assembly: ", assembly),
        paste0("BAM: ", bam_path),
        paste0(
            "Windows: ", format(n_windows, big.mark = ",", trim = TRUE)
        ),
        paste0(
            "Duration: ",
            if (!is.null(provenance$duration_secs)) {
                paste0(format(provenance$duration_secs, digits = 4,
                    trim = TRUE), " s")
            } else {
                "unavailable"
            }
        ),
        paste0("Config hash: ", bam_bundle_text(provenance$config_hash))
    )
}

bam_bundle_safe_mean <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else mean(x)
}

bam_bundle_safe_median <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else stats::median(x)
}

bam_bundle_safe_sd <- function(x) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (length(x) < 2L) NA_real_ else stats::sd(x)
}

bam_bundle_write_csv <- function(x, path, gzip = FALSE) {
    if (isTRUE(gzip)) {
        connection <- gzfile(path, "w")
        on.exit(close(connection), add = TRUE)
        utils::write.csv(x, connection, row.names = FALSE, na = "")
    } else {
        utils::write.csv(
            x, path, row.names = FALSE, na = "",
            fileEncoding = "UTF-8"
        )
    }
    invisible(path)
}

bam_bundle_metrics_summary <- function(metrics) {
    empty <- data.frame(
        sample_name = character(), chromosome_group = character(),
        role = character(), component = character(), windows = integer(),
        contigs = integer(), mean_breadth = numeric(),
        median_breadth = numeric(), sd_breadth = numeric(),
        mean_read_count = numeric(), mean_high_mapq_read_fraction = numeric(),
        mean_gc_fraction = numeric(), mean_expected_breadth = numeric(),
        stringsAsFactors = FALSE
    )
    if (!is.data.frame(metrics) || nrow(metrics) == 0L) {
        return(empty)
    }
    key <- paste(
        metrics$sample_name, metrics$chromosome_group,
        metrics$role, metrics$component, sep = "\r"
    )
    parts <- split(seq_len(nrow(metrics)), key)
    out <- do.call(rbind, lapply(parts, function(i) {
        data.frame(
            sample_name = as.character(metrics$sample_name[[i[[1L]]]]),
            chromosome_group = as.character(
                metrics$chromosome_group[[i[[1L]]]]
            ),
            role = as.character(metrics$role[[i[[1L]]]]),
            component = as.character(metrics$component[[i[[1L]]]]),
            windows = length(i),
            contigs = length(unique(
                as.character(metrics$contig_id[i])
            )),
            mean_breadth = bam_bundle_safe_mean(metrics$breadth[i]),
            median_breadth = bam_bundle_safe_median(metrics$breadth[i]),
            sd_breadth = bam_bundle_safe_sd(metrics$breadth[i]),
            mean_read_count = bam_bundle_safe_mean(metrics$read_count[i]),
            mean_high_mapq_read_fraction = bam_bundle_safe_mean(
                metrics$high_mapq_read_fraction[i]
            ),
            mean_gc_fraction = bam_bundle_safe_mean(
                metrics$gc_fraction[i]
            ),
            mean_expected_breadth = bam_bundle_safe_mean(
                metrics$expected_breadth[i]
            ),
            stringsAsFactors = FALSE
        )
    }))
    out <- out[order(
        match(
            out$chromosome_group,
            bam_bundle_group_order(out$chromosome_group)
        ),
        out$role, out$component
    ), , drop = FALSE]
    rownames(out) <- NULL
    out
}

bam_bundle_contig_summary <- function(metrics) {
    empty <- data.frame(
        sample_name = character(), contig_id = character(),
        contig_name = character(), chromosome_group = character(),
        role = character(), component = character(),
        copy_number = numeric(), windows = integer(),
        length_bp = numeric(), mean_breadth = numeric(),
        median_breadth = numeric(), mean_read_count = numeric(),
        mean_high_mapq_read_fraction = numeric(),
        mean_gc_fraction = numeric(), stringsAsFactors = FALSE
    )
    if (!is.data.frame(metrics) || nrow(metrics) == 0L) {
        return(empty)
    }
    parts <- split(seq_len(nrow(metrics)), as.character(metrics$contig_id))
    out <- do.call(rbind, lapply(parts, function(i) {
        first <- i[[1L]]
        data.frame(
            sample_name = as.character(metrics$sample_name[[first]]),
            contig_id = as.character(metrics$contig_id[[first]]),
            contig_name = as.character(metrics$contig_name[[first]]),
            chromosome_group = as.character(
                metrics$chromosome_group[[first]]
            ),
            role = as.character(metrics$role[[first]]),
            component = as.character(metrics$component[[first]]),
            copy_number = as.numeric(metrics$copy_number[[first]]),
            windows = length(i),
            length_bp = sum(
                as.numeric(metrics$window_length_bp[i]), na.rm = TRUE
            ),
            mean_breadth = bam_bundle_safe_mean(metrics$breadth[i]),
            median_breadth = bam_bundle_safe_median(metrics$breadth[i]),
            mean_read_count = bam_bundle_safe_mean(
                metrics$read_count[i]
            ),
            mean_high_mapq_read_fraction = bam_bundle_safe_mean(
                metrics$high_mapq_read_fraction[i]
            ),
            mean_gc_fraction = bam_bundle_safe_mean(
                metrics$gc_fraction[i]
            ),
            stringsAsFactors = FALSE
        )
    }))
    out <- out[order(
        match(
            out$chromosome_group,
            bam_bundle_group_order(out$chromosome_group)
        ),
        out$role, out$contig_id
    ), , drop = FALSE]
    rownames(out) <- NULL
    out
}

bam_bundle_scores_summary <- function(scores) {
    empty <- data.frame(
        sample_name = character(), chromosome_group = character(),
        windows = integer(), mean_introgression_score = numeric(),
        median_introgression_score = numeric(),
        sd_introgression_score = numeric(),
        mean_donor_residual = numeric(),
        mean_recipient_residual = numeric(),
        mean_total_depth_residual = numeric(),
        max_introgression_score = numeric(), stringsAsFactors = FALSE
    )
    if (!is.data.frame(scores) || nrow(scores) == 0L) {
        return(empty)
    }
    parts <- split(
        seq_len(nrow(scores)), as.character(scores$chromosome_group)
    )
    out <- do.call(rbind, lapply(parts, function(i) {
        first <- i[[1L]]
        data.frame(
            sample_name = as.character(scores$sample_name[[first]]),
            chromosome_group = as.character(
                scores$chromosome_group[[first]]
            ),
            windows = length(i),
            mean_introgression_score = bam_bundle_safe_mean(
                scores$introgression_score[i]
            ),
            median_introgression_score = bam_bundle_safe_median(
                scores$introgression_score[i]
            ),
            sd_introgression_score = bam_bundle_safe_sd(
                scores$introgression_score[i]
            ),
            mean_donor_residual = bam_bundle_safe_mean(
                scores$donor_residual[i]
            ),
            mean_recipient_residual = bam_bundle_safe_mean(
                scores$recipient_residual[i]
            ),
            mean_total_depth_residual = bam_bundle_safe_mean(
                scores$total_depth_residual[i]
            ),
            max_introgression_score = if (any(is.finite(
                scores$introgression_score[i]
            ))) {
                max(scores$introgression_score[i], na.rm = TRUE)
            } else {
                NA_real_
            },
            stringsAsFactors = FALSE
        )
    }))
    out <- out[order(
        match(
            out$chromosome_group,
            bam_bundle_group_order(out$chromosome_group)
        )
    ), , drop = FALSE]
    rownames(out) <- NULL
    out
}

bam_bundle_masks_summary <- function(masks) {
    empty <- data.frame(
        chromosome_group = character(), windows = integer(),
        masked_windows = integer(), masked_fraction = numeric(),
        mean_mask_signal = numeric(), mean_mask_posterior = numeric(),
        mask_sources = character(), mask_version = character(),
        stringsAsFactors = FALSE
    )
    if (!is.data.frame(masks) || nrow(masks) == 0L) {
        return(empty)
    }
    parts <- split(seq_len(nrow(masks)), as.character(masks$chromosome_group))
    out <- do.call(rbind, lapply(parts, function(i) {
        masked <- isTRUE_vec(masks$masked[i])
        sources <- unique(as.character(masks$mask_source[i][masked]))
        sources <- sources[
            !is.na(sources) & nzchar(sources) & sources != "none"
        ]
        versions <- unique(as.character(masks$mask_version[i][masked]))
        versions <- versions[!is.na(versions) & nzchar(versions)]
        data.frame(
            chromosome_group = as.character(
                masks$chromosome_group[[i[[1L]]]]
            ),
            windows = length(i),
            masked_windows = sum(masked, na.rm = TRUE),
            masked_fraction = bam_bundle_safe_mean(masked),
            mean_mask_signal = bam_bundle_safe_mean(
                masks$mask_signal[i][masked]
            ),
            mean_mask_posterior = bam_bundle_safe_mean(
                masks$mask_posterior[i][masked]
            ),
            mask_sources = paste(sources, collapse = ";"),
            mask_version = paste(versions, collapse = ";"),
            stringsAsFactors = FALSE
        )
    }))
    out <- out[order(
        match(
            out$chromosome_group,
            bam_bundle_group_order(out$chromosome_group)
        )
    ), , drop = FALSE]
    rownames(out) <- NULL
    out
}

bam_bundle_call_summary <- function(calls) {
    empty <- data.frame(
        call_class = character(), calls = integer(),
        passing_calls = integer(), confident_calls = integer(),
        needs_review_calls = integer(), ambiguous_calls = integer(),
        total_width_bp = numeric(), median_width_bp = numeric(),
        mean_introgression_score = numeric(), stringsAsFactors = FALSE
    )
    if (!is.data.frame(calls) || nrow(calls) == 0L) {
        return(empty)
    }
    classes <- bam_bundle_group_order(calls$class)
    if (length(classes) == 0L) {
        return(empty)
    }
    status <- as.character(calls$class_status)
    passing <- isTRUE_vec(calls$pass)
    confident <- passing & !is.na(status) & status == "confident"
    widths <- as.numeric(calls$width_bp)
    scores <- as.numeric(calls$score_mean)
    out <- do.call(rbind, lapply(classes, function(value) {
        selected <- !is.na(calls$class) & as.character(calls$class) == value
        data.frame(
            call_class = value,
            calls = sum(selected),
            passing_calls = sum(selected & passing),
            confident_calls = sum(selected & confident),
            needs_review_calls = sum(
                selected & !is.na(status) & status == "needs_review"
            ),
            ambiguous_calls = sum(
                selected & !is.na(status) & status == "ambiguous"
            ),
            total_width_bp = if (any(is.finite(widths[selected]))) {
                sum(widths[selected], na.rm = TRUE)
            } else {
                NA_real_
            },
            median_width_bp = bam_bundle_safe_median(widths[selected]),
            mean_introgression_score = bam_bundle_safe_mean(
                scores[selected]
            ),
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    out
}

bam_bundle_call_status_summary <- function(calls, column = "class_status",
                                           label = "review_status") {
    empty <- data.frame(
        value = character(), calls = integer(), passing_calls = integer(),
        total_width_bp = numeric(), median_width_bp = numeric(),
        stringsAsFactors = FALSE
    )
    names(empty)[[1L]] <- label
    if (!is.data.frame(calls) || nrow(calls) == 0L ||
        !column %in% names(calls)) {
        return(empty)
    }
    values <- bam_bundle_group_order(calls[[column]])
    if (length(values) == 0L) {
        return(empty)
    }
    passing <- isTRUE_vec(calls$pass)
    widths <- as.numeric(calls$width_bp)
    out <- do.call(rbind, lapply(values, function(value) {
        selected <- !is.na(calls[[column]]) &
            as.character(calls[[column]]) == value
        frame <- data.frame(
            value = value,
            calls = sum(selected),
            passing_calls = sum(selected & passing),
            total_width_bp = if (any(is.finite(widths[selected]))) {
                sum(widths[selected], na.rm = TRUE)
            } else {
                NA_real_
            },
            median_width_bp = bam_bundle_safe_median(widths[selected]),
            stringsAsFactors = FALSE
        )
        names(frame)[[1L]] <- label
        frame
    }))
    rownames(out) <- NULL
    out
}

bam_bundle_call_group_summary <- function(calls, masks_summary = NULL) {
    empty <- data.frame(
        chromosome_group = character(), calls = integer(),
        passing_calls = integer(), confident_calls = integer(),
        total_width_bp = numeric(), median_width_bp = numeric(),
        masked_fraction = numeric(), stringsAsFactors = FALSE
    )
    if (!is.data.frame(calls) || nrow(calls) == 0L) {
        return(empty)
    }
    groups <- bam_bundle_group_order(calls$group)
    if (length(groups) == 0L) {
        return(empty)
    }
    status <- as.character(calls$class_status)
    passing <- isTRUE_vec(calls$pass)
    confident <- passing & !is.na(status) & status == "confident"
    widths <- as.numeric(calls$width_bp)
    mask_lookup <- NULL
    if (is.data.frame(masks_summary) &&
        all(c("chromosome_group", "masked_fraction") %in%
            names(masks_summary))) {
        mask_lookup <- stats::setNames(
            masks_summary$masked_fraction,
            as.character(masks_summary$chromosome_group)
        )
    }
    out <- do.call(rbind, lapply(groups, function(value) {
        selected <- !is.na(calls$group) & as.character(calls$group) == value
        data.frame(
            chromosome_group = value,
            calls = sum(selected),
            passing_calls = sum(selected & passing),
            confident_calls = sum(selected & confident),
            total_width_bp = if (any(is.finite(widths[selected]))) {
                sum(widths[selected], na.rm = TRUE)
            } else {
                NA_real_
            },
            median_width_bp = bam_bundle_safe_median(widths[selected]),
            masked_fraction = if (!is.null(mask_lookup) &&
                value %in% names(mask_lookup)) {
                unname(mask_lookup[[value]])
            } else {
                NA_real_
            },
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    out
}

bam_bundle_assignment_alternatives <- function(calls) {
    columns <- c(
        "sample_name", "chromosome_group", "call_class", "review_status",
        "passes_calibration", "recipient_chromosome",
        "recipient_assignment_status", "alternative_recipient_chromosome",
        "alternative_recipient_score", "recipient_assignment_margin",
        "recipient_assignment_weight", "group_start_bp", "group_end_bp"
    )
    empty <- data.frame(
        sample_name = character(), chromosome_group = character(),
        call_class = character(), review_status = character(),
        passes_calibration = logical(), recipient_chromosome = character(),
        recipient_assignment_status = character(),
        alternative_recipient_chromosome = character(),
        alternative_recipient_score = numeric(),
        recipient_assignment_margin = numeric(),
        recipient_assignment_weight = numeric(),
        group_start_bp = numeric(), group_end_bp = numeric(),
        stringsAsFactors = FALSE
    )
    if (!is.data.frame(calls) ||
        !"alternative_recipient_chromosome" %in% names(calls)) {
        return(empty)
    }
    alternate <- as.character(calls$alternative_recipient_chromosome)
    keep <- !is.na(alternate) & nzchar(alternate)
    if (!any(keep)) {
        return(empty)
    }
    selected <- intersect(columns, names(calls))
    out <- calls[keep, selected, drop = FALSE]
    rownames(out) <- NULL
    out
}

bam_bundle_write_bed <- function(calls, file, side = c("recipient", "donor")) {
    side <- match.arg(side)
    track <- sprintf(
        paste0(
            "track name=\"MaMaMIA %s calls\" ",
            "description=\"MaMaMIA calls on %s contigs\""
        ),
        side, side
    )
    if (!is.data.frame(calls) || nrow(calls) == 0L) {
        writeLines(track, file)
        return(invisible(file))
    }
    if (identical(side, "recipient")) {
        chrom <- calls$recipient_chr_id
        start <- calls$recipient_start_bp
        stop <- calls$recipient_end_bp
    } else {
        chrom <- calls$donor_chr_id
        start <- calls$donor_start_bp
        stop <- calls$donor_end_bp
    }
    keep <- !is.na(chrom) & nzchar(as.character(chrom)) &
        is.finite(start) & is.finite(stop) & stop > start
    if (!any(keep)) {
        writeLines(track, file)
        return(invisible(file))
    }
    name <- paste(
        as.character(calls$class[keep]),
        as.character(calls$class_status[keep]),
        ifelse(isTRUE_vec(calls$pass[keep]), "pass", "fail"),
        sep = "|"
    )
    out <- data.frame(
        chrom = as.character(chrom[keep]),
        start = as.integer(start[keep]),
        end = as.integer(stop[keep]),
        name = name,
        score = 0L,
        strand = ".",
        stringsAsFactors = FALSE
    )
    writeLines(
        c(track, do.call(paste, c(out, sep = "\t"))),
        file
    )
    invisible(file)
}

bam_bundle_write_masks_bed <- function(masks, file) {
    track <- paste0(
        "track name=\"MaMaMIA masked windows\" ",
        "description=\"Windows masked before scoring, in group coordinates\""
    )
    if (!is.data.frame(masks) || nrow(masks) == 0L) {
        writeLines(track, file)
        return(invisible(file))
    }
    keep <- isTRUE_vec(masks$masked)
    if (!any(keep)) {
        writeLines(track, file)
        return(invisible(file))
    }
    source <- as.character(masks$mask_source[keep])
    source[is.na(source) | !nzchar(source)] <- "masked"
    posterior <- as.numeric(masks$mask_posterior[keep])
    score <- ifelse(
        is.finite(posterior),
        pmax(0, pmin(1000, round(posterior * 1000))),
        0L
    )
    out <- data.frame(
        chrom = paste0("group_", as.character(masks$chromosome_group[keep])),
        start = as.integer(masks$start_bp[keep]),
        end = as.integer(masks$end_bp[keep]),
        name = source,
        score = as.integer(score),
        strand = ".",
        stringsAsFactors = FALSE
    )
    writeLines(c(track, do.call(paste, c(out, sep = "\t"))), file)
    invisible(file)
}

bam_bundle_write_breadth_bedgraph <- function(metrics, file) {
    connection <- gzfile(file, "w")
    on.exit(close(connection), add = TRUE)
    writeLines(
        paste0(
            "track type=bedGraph name=\"MaMaMIA breadth\" ",
            "description=\"Window breadth for this sample\""
        ),
        connection
    )
    if (is.data.frame(metrics) && nrow(metrics) > 0L) {
        keep <- is.finite(metrics$breadth) &
            is.finite(metrics$start_bp) & is.finite(metrics$end_bp)
        if (any(keep)) {
            out <- data.frame(
                contig_id = as.character(metrics$contig_id[keep]),
                start_bp = as.integer(metrics$start_bp[keep]),
                end_bp = as.integer(metrics$end_bp[keep]),
                breadth = round(as.numeric(metrics$breadth[keep]), 6),
                stringsAsFactors = FALSE
            )
            utils::write.table(
                out, connection, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = FALSE
            )
        }
    }
    invisible(file)
}

bam_bundle_json_ready <- function(x) {
    if (is.data.frame(x)) {
        return(x)
    }
    if (inherits(x, "POSIXt") || inherits(x, "Date")) {
        return(as.character(x))
    }
    if (is.list(x)) {
        x[] <- lapply(x, bam_bundle_json_ready)
        return(x)
    }
    if (is.factor(x)) {
        return(as.character(x))
    }
    if (is.atomic(x) && !is.null(names(x))) {
        return(as.list(x))
    }
    x
}

bam_bundle_write_json <- function(x, path) {
    writeLines(
        jsonlite::toJSON(
            bam_bundle_json_ready(x),
            auto_unbox = TRUE, pretty = TRUE, dataframe = "rows",
            na = "null", null = "null"
        ),
        path,
        useBytes = TRUE
    )
    invisible(path)
}

bam_bundle_run_files <- function(result, bundle, sample) {
    run_dir <- file.path(bundle, "run")
    dir.create(run_dir, showWarnings = FALSE)
    provenance <- result$provenance
    bam_bundle_write_json(
        provenance, file.path(run_dir, "provenance.json")
    )
    config_lines <- strsplit(yaml::as.yaml(result$config), "\n",
        fixed = TRUE
    )[[1L]]
    config_path <- provenance$config_path
    header <- c(
        "# MaMaMIA resolved configuration",
        if (!is.null(config_path) && length(config_path) == 1L &&
            !is.na(config_path) && nzchar(config_path)) {
            paste0("# Source: ", config_path)
        } else {
            "# Source: not recorded"
        },
        ""
    )
    writeLines(
        c(header, config_lines),
        file.path(run_dir, "config.yml"),
        useBytes = TRUE
    )
    writeLines(
        utils::capture.output(utils::sessionInfo()),
        file.path(run_dir, "session-info.txt"),
        useBytes = TRUE
    )

    contigs <- result$assembly$contigs
    pick <- function(x, name) {
        if (is.data.frame(x) && name %in% names(x)) {
            x[[name]]
        } else {
            rep(NA, if (is.data.frame(x)) nrow(x) else 0L)
        }
    }
    assembly_key <- if (!is.null(result$assembly$key)) {
        as.character(result$assembly$key)
    } else {
        NA_character_
    }
    assembly_label <- if (!is.null(result$assembly$label)) {
        as.character(result$assembly$label)
    } else {
        NA_character_
    }
    assembly_parts <- lapply(names(contigs), function(role) {
        x <- contigs[[role]]
        if (!is.data.frame(x) || nrow(x) == 0L) {
            return(NULL)
        }
        data.frame(
            assembly_key = rep(assembly_key, nrow(x)),
            assembly_label = rep(assembly_label, nrow(x)),
            role = role,
            contig_id = as.character(pick(x, "id")),
            contig_name = as.character(pick(x, "name")),
            component = as.character(pick(x, "component")),
            chromosome_group = as.character(pick(x, "group")),
            copy_number = as.numeric(pick(x, "copy_number")),
            length_bp = as.numeric(pick(x, "length")),
            index_path = as.character(pick(x, "index_path")),
            stringsAsFactors = FALSE
        )
    })
    assembly_parts <- Filter(Negate(is.null), assembly_parts)
    assembly <- if (length(assembly_parts) > 0L) {
        do.call(rbind, assembly_parts)
    } else {
        data.frame(
            assembly_key = character(), assembly_label = character(),
            role = character(), contig_id = character(),
            contig_name = character(), component = character(),
            chromosome_group = character(), copy_number = numeric(),
            length_bp = numeric(), index_path = character(),
            stringsAsFactors = FALSE
        )
    }
    bam_bundle_write_csv(assembly, file.path(run_dir, "assembly.csv"))
    index_paths <- unique(as.character(assembly$index_path))
    index_paths <- index_paths[!is.na(index_paths) & nzchar(index_paths)]
    if (length(index_paths) > 0L) {
        index_info <- file.info(index_paths)
        index_table <- data.frame(
            index_path = index_paths,
            exists = file.exists(index_paths),
            bytes = as.numeric(index_info$size),
            modified = as.character(index_info$mtime),
            stringsAsFactors = FALSE
        )
    } else {
        index_table <- data.frame(
            index_path = character(), exists = logical(),
            bytes = numeric(), modified = character(),
            stringsAsFactors = FALSE
        )
    }
    bam_bundle_write_csv(
        index_table, file.path(run_dir, "reference-indexes.csv")
    )
    bam_bundle_write_csv(result$samples, file.path(run_dir, "samples.csv"))

    bams <- result$bams
    provenance_bams <- provenance$bams
    if (is.data.frame(bams) && is.data.frame(provenance_bams) &&
        "path" %in% names(bams) && "path" %in% names(provenance_bams)) {
        position <- match(
            as.character(bams$path), as.character(provenance_bams$path)
        )
        for (column in c(
            "file_size", "mtime", "sample_source", "assembly"
        )) {
            if (column %in% names(provenance_bams) &&
                !column %in% names(bams)) {
                bams[[column]] <- provenance_bams[[column]][position]
            }
        }
    }
    bam_bundle_write_csv(bams, file.path(run_dir, "bams.csv"))

    scalar <- vapply(provenance, function(value) {
        is.atomic(value) && length(value) == 1L
    }, logical(1))
    parameters <- data.frame(
        parameter = names(provenance)[scalar],
        value = as.character(unlist(provenance[scalar])),
        stringsAsFactors = FALSE
    )
    bam_bundle_write_csv(
        parameters, file.path(run_dir, "parameters.csv")
    )

    dependency <- provenance$dependency_versions
    get_dependency <- function(name) {
        if (!is.null(names(dependency)) && name %in% names(dependency)) {
            as.character(dependency[[name]])
        } else {
            NA_character_
        }
    }
    versions <- data.frame(
        component = c(
            "MaMaMIA", "R", "Rhtslib", "GenomicRanges", "DNAcopy",
            "htslib", "uvr", "uvr_lock_hash", "gcrank_commit",
            "config_hash", "calibration_version"
        ),
        version = c(
            bam_bundle_text(provenance$package_version),
            bam_bundle_text(provenance$r_version),
            get_dependency("Rhtslib"),
            get_dependency("GenomicRanges"),
            get_dependency("DNAcopy"),
            bam_bundle_text(provenance$htslib_version),
            bam_bundle_text(provenance$uvr_version),
            bam_bundle_text(provenance$uvr_lock_hash),
            bam_bundle_text(provenance$gcrank_commit),
            bam_bundle_text(provenance$config_hash),
            bam_bundle_text(provenance$calibration_version)
        ),
        stringsAsFactors = FALSE
    )
    bam_bundle_write_csv(versions, file.path(run_dir, "versions.csv"))
    versions
}

bam_bundle_qc_files <- function(result, bundle, sample) {
    qc_dir <- file.path(bundle, "qc")
    dir.create(qc_dir, showWarnings = FALSE)
    metrics <- result$metrics
    if (is.data.frame(metrics) && "sample" %in% names(metrics)) {
        metrics <- metrics[
            as.character(metrics$sample) == sample, , drop = FALSE
        ]
    }
    scores <- result$scores
    if (is.data.frame(scores) && "sample" %in% names(scores)) {
        scores <- scores[
            as.character(scores$sample) == sample, , drop = FALSE
        ]
    }
    metrics <- bam_bundle_rename_columns(
        metrics, bam_bundle_metric_columns()
    )
    scores <- bam_bundle_rename_columns(
        scores, bam_bundle_score_columns()
    )
    masks <- bam_bundle_rename_columns(
        result$masks, bam_bundle_mask_columns()
    )
    calibration <- bam_bundle_rename_columns(
        result$calibration, bam_bundle_calibration_columns()
    )

    bam_bundle_write_csv(
        metrics, file.path(qc_dir, "metrics-windows.csv.gz"), gzip = TRUE
    )
    bam_bundle_write_csv(
        scores, file.path(qc_dir, "scores-windows.csv.gz"), gzip = TRUE
    )
    bam_bundle_write_csv(
        masks, file.path(qc_dir, "masks-windows.csv.gz"), gzip = TRUE
    )
    bam_bundle_write_csv(
        calibration, file.path(qc_dir, "calibration.csv")
    )

    metrics_summary <- bam_bundle_metrics_summary(metrics)
    contig_summary <- bam_bundle_contig_summary(metrics)
    scores_summary <- bam_bundle_scores_summary(scores)
    masks_summary <- bam_bundle_masks_summary(masks)
    bam_bundle_write_csv(
        metrics_summary, file.path(qc_dir, "metrics-summary.csv")
    )
    bam_bundle_write_csv(
        contig_summary, file.path(qc_dir, "contig-summary.csv")
    )
    bam_bundle_write_csv(
        scores_summary, file.path(qc_dir, "scores-summary.csv")
    )
    bam_bundle_write_csv(
        masks_summary, file.path(qc_dir, "masks-summary.csv")
    )
    cohort <- bam_cohort_summary(result)
    bam_bundle_write_csv(
        cohort, file.path(qc_dir, "cohort-summary.csv")
    )
    list(
        metrics = metrics,
        scores = scores,
        masks = masks,
        calibration = calibration,
        metrics_summary = metrics_summary,
        contig_summary = contig_summary,
        scores_summary = scores_summary,
        masks_summary = masks_summary,
        cohort = cohort
    )
}

bam_bundle_write_diagnostics <- function(result, bundle, sample,
                                         calls_raw, calls_desc) {
    run_dir <- file.path(bundle, "run")
    qc_dir <- file.path(bundle, "qc")
    tracks_dir <- file.path(bundle, "tracks")
    dir.create(run_dir, showWarnings = FALSE)
    dir.create(qc_dir, showWarnings = FALSE)
    dir.create(tracks_dir, showWarnings = FALSE)

    status <- as.character(calls_raw$class_status)
    passing <- isTRUE_vec(calls_raw$pass)
    confident <- passing & !is.na(status) & status == "confident"
    bam_bundle_write_csv(
        calls_desc[passing, , drop = FALSE],
        file.path(bundle, "calls-passing.csv")
    )
    bam_bundle_write_csv(
        calls_desc[!confident, , drop = FALSE],
        file.path(bundle, "calls-review.csv")
    )

    calls_summary <- bam_bundle_call_summary(calls_raw)
    calls_by_status <- bam_bundle_call_status_summary(calls_raw)
    calls_by_copy_state <- bam_bundle_call_status_summary(
        calls_raw, column = "copy_state", label = "copy_number_state"
    )
    masks_for_groups <- bam_bundle_rename_columns(
        result$masks, bam_bundle_mask_columns()
    )
    masks_summary <- bam_bundle_masks_summary(masks_for_groups)
    calls_by_group <- bam_bundle_call_group_summary(
        calls_raw, masks_summary
    )
    alternatives <- bam_bundle_assignment_alternatives(calls_desc)

    bam_bundle_write_csv(
        calls_summary, file.path(bundle, "calls-summary.csv")
    )
    bam_bundle_write_csv(
        calls_by_group, file.path(bundle, "calls-by-group.csv")
    )
    bam_bundle_write_csv(
        calls_by_status, file.path(bundle, "calls-by-status.csv")
    )
    bam_bundle_write_csv(
        calls_by_copy_state, file.path(bundle, "calls-by-copy-state.csv")
    )
    bam_bundle_write_csv(
        alternatives, file.path(bundle, "assignment-alternatives.csv")
    )

    versions <- bam_bundle_run_files(result, bundle, sample)
    qc <- bam_bundle_qc_files(result, bundle, sample)

    bam_bundle_write_bed(
        calls_raw, file.path(tracks_dir, "calls-recipient.bed"),
        side = "recipient"
    )
    bam_bundle_write_bed(
        calls_raw, file.path(tracks_dir, "calls-donor.bed"), side = "donor"
    )
    bam_bundle_write_masks_bed(
        qc$masks, file.path(tracks_dir, "masks.bed")
    )
    bam_bundle_write_breadth_bedgraph(
        qc$metrics, file.path(tracks_dir, "breadth.bedGraph.gz")
    )

    list(
        calls_summary = calls_summary,
        calls_by_group = calls_by_group,
        calls_by_status = calls_by_status,
        calls_by_copy_state = calls_by_copy_state,
        alternatives = alternatives,
        metrics_summary = qc$metrics_summary,
        contig_summary = qc$contig_summary,
        scores_summary = qc$scores_summary,
        masks_summary = qc$masks_summary,
        cohort = qc$cohort,
        calibration = qc$calibration,
        versions = versions
    )
}

bam_bundle_manifest <- function(bundle, sample, generated, version) {
    files <- list.files(bundle, recursive = TRUE, full.names = TRUE)
    files <- files[file.info(files)$isdir %in% FALSE]
    relative <- substring(files, nchar(bundle) + 2L)
    relative <- gsub("\\\\", "/", relative)
    info <- file.info(files)
    manifest <- list(
        bundle = basename(bundle),
        sample = sample,
        generated = format(generated, "%Y-%m-%dT%H:%M:%S%z"),
        package = "MaMaMIA",
        package_version = version,
        files = data.frame(
            path = relative,
            bytes = as.numeric(info$size),
            md5 = unname(tools::md5sum(files)),
            stringsAsFactors = FALSE
        )
    )
    bam_bundle_write_json(
        manifest, file.path(bundle, "manifest.json")
    )
    invisible(manifest)
}

bam_bundle_checksums <- function(bundle) {
    files <- list.files(bundle, recursive = TRUE, full.names = TRUE)
    files <- files[basename(files) != "checksums.md5"]
    files <- files[file.info(files)$isdir %in% FALSE]
    relative <- substring(files, nchar(bundle) + 2L)
    relative <- gsub("\\\\", "/", relative)
    md5 <- unname(tools::md5sum(files))
    order_index <- order(relative)
    writeLines(
        paste(md5[order_index], relative[order_index], sep = "  "),
        file.path(bundle, "checksums.md5"),
        useBytes = TRUE
    )
    invisible(file.path(bundle, "checksums.md5"))
}

bam_bundle_html_escape <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub("\"", "&quot;", x, fixed = TRUE)
}

bam_bundle_html_value <- function(value) {
    if (length(value) == 0L || is.null(value)) {
        return("&mdash;")
    }
    if (is.logical(value)) {
        if (is.na(value)) "&mdash;" else if (value) "yes" else "no"
    } else if (is.numeric(value)) {
        if (!is.finite(value)) {
            "&mdash;"
        } else {
            format(value, digits = 4, trim = TRUE, scientific = FALSE)
        }
    } else {
        text <- as.character(value)
        if (is.na(text) || !nzchar(text)) "&mdash;" else
            bam_bundle_html_escape(text)
    }
}

bam_bundle_html_table <- function(x) {
    if (!is.data.frame(x) || nrow(x) == 0L || ncol(x) == 0L) {
        return("<p class=\"muted\">No data for this section.</p>")
    }
    header <- paste0(
        "<tr>",
        paste0(
            "<th>", bam_bundle_html_escape(names(x)), "</th>",
            collapse = ""
        ),
        "</tr>"
    )
    rows <- vapply(seq_len(nrow(x)), function(i) {
        cells <- vapply(x, function(column) {
            bam_bundle_html_value(column[[i]])
        }, character(1))
        paste0("<tr><td>", paste(cells, collapse = "</td><td>"),
            "</td></tr>")
    }, character(1))
    paste0(
        "<div class=\"scroll\"><table>", header, paste(rows, collapse = ""),
        "</table></div>"
    )
}

bam_bundle_report_html <- function(result, bundle, sample, calls_raw,
                                   diagnostics, generated, version,
                                   has_full_svg = FALSE) {
    provenance <- result$provenance
    sample_meta <- result$samples[
        as.character(result$samples$sample) == sample, , drop = FALSE
    ]
    n_windows <- if ("n_windows" %in% names(sample_meta) &&
        nrow(sample_meta) > 0L) {
        sample_meta$n_windows[[1L]]
    } else {
        nrow(result$windows$table)
    }
    mode <- if ("mode" %in% names(sample_meta) && nrow(sample_meta) > 0L) {
        as.character(sample_meta$mode[[1L]])
    } else {
        NA_character_
    }
    masked_fraction <- if ("masked_fraction" %in% names(sample_meta) &&
        nrow(sample_meta) > 0L) {
        as.numeric(sample_meta$masked_fraction[[1L]])
    } else {
        NA_real_
    }
    passing <- sum(isTRUE_vec(calls_raw$pass), na.rm = TRUE)
    confident <- nrow(bam_confident_calls(calls_raw))
    bam_path <- if ("path" %in% names(sample_meta) &&
        nrow(sample_meta) > 0L) {
        as.character(sample_meta$path[[1L]])
    } else {
        NA_character_
    }
    assembly_label <- if (!is.null(result$assembly$label)) {
        as.character(result$assembly$label)
    } else {
        as.character(result$assembly$key)
    }
    overview <- data.frame(
        Item = c(
            "Sample", "Assembly", "BAM", "Analysis mode", "Windows",
            "Calls", "Passing calibration", "Confident calls",
            "Masked windows", "Duration", "Config hash",
            "Calibration version", "Generated", "MaMaMIA"
        ),
        Value = c(
            sample,
            assembly_label,
            bam_path,
            mode,
            format(n_windows, big.mark = ",", trim = TRUE),
            format(nrow(calls_raw), big.mark = ",", trim = TRUE),
            format(passing, big.mark = ",", trim = TRUE),
            format(confident, big.mark = ",", trim = TRUE),
            if (is.finite(masked_fraction)) {
                paste0(format(masked_fraction * 100, digits = 3,
                    trim = TRUE), "%")
            } else {
                NA_character_
            },
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

    files <- list.files(bundle, recursive = TRUE, full.names = FALSE)
    file_info <- file.info(file.path(bundle, files))
    file_table <- data.frame(
        path = gsub("\\\\", "/", files),
        bytes = as.numeric(file_info$size),
        stringsAsFactors = FALSE
    )
    file_table <- file_table[order(file_table$path), , drop = FALSE]

    cohort <- diagnostics$cohort
    if (is.data.frame(cohort) && nrow(cohort) > 1L) {
        keep <- intersect(
            c(
                "sample", "mode", "n_windows", "n_all_calls",
                "n_class_confident", "pass_rate", "masked_fraction",
                "reads_total", "total_width_mb"
            ),
            names(cohort)
        )
        cohort <- cohort[, keep, drop = FALSE]
    }

    css <- paste(
        "body{margin:0;background:#f4f6f3;color:#1f2a24;",
        "font:14px/1.5 -apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif}",
        "main{max-width:1100px;margin:0 auto;padding:32px 24px 48px}",
        "header{border-bottom:2px solid #1f6f5c;padding-bottom:16px;",
        "margin-bottom:8px}",
        "h1{font-size:28px;line-height:1.2;margin:4px 0 6px}",
        "h2{font-size:17px;margin:28px 0 10px;color:#1f6f5c}",
        "h3{font-size:14px;margin:20px 0 8px;color:#33443a}",
        ".eyebrow{font-size:11px;font-weight:700;letter-spacing:.14em;",
        "text-transform:uppercase;color:#65736a;margin:0}",
        ".lead{color:#4f5f55;margin:0}",
        "table{border-collapse:collapse;width:100%;background:#fff;",
        "font-size:12px}",
        "th,td{border-bottom:1px solid #dfe5e0;padding:6px 8px;",
        "text-align:left;vertical-align:top}",
        "th{background:#eef3ef;font-size:10px;letter-spacing:.06em;",
        "text-transform:uppercase;color:#5b6b60}",
        "tr:nth-child(even) td{background:#fafbfa}",
        ".scroll{overflow-x:auto;border:1px solid #dfe5e0;border-radius:6px}",
        "img,.karyotype svg{max-width:100%;height:auto;",
        "border:1px solid #dfe5e0;background:#fff}",
        ".muted{color:#65736a}",
        "code{background:#eef3ef;border-radius:3px;padding:1px 4px;",
        "font-size:12px}",
        "footer{margin-top:36px;padding-top:14px;",
        "border-top:1px solid #dfe5e0;color:#65736a;font-size:11px}",
        sep = ""
    )

    karyotype <- if (isTRUE(has_full_svg)) {
        svg_path <- file.path(bundle, "karyotype-full.svg")
        svg <- if (file.exists(svg_path)) {
            paste(readLines(svg_path, warn = FALSE), collapse = "\n")
        } else {
            ""
        }
        svg <- sub("^<\\?xml[^>]*>\\s*", "", svg)
        if (nzchar(svg)) {
            c(
                "<section><h2>Karyotype</h2>",
                "<div class=\"karyotype\">",
                svg,
                "</div>",
                "<p class=\"muted\">Full karyotype for this sample. ",
                "The filtered and confident-only views are ",
                "<code>karyotype-current.svg</code> and ",
                "<code>karyotype-confident.svg</code>.</p></section>"
            )
        } else {
            c(
                "<section><h2>Karyotype</h2>",
                "<img src=\"karyotype-full.svg\" ",
                "alt=\"Karyotype with all calls\">",
                "<p class=\"muted\">Full karyotype for this sample. ",
                "The filtered and confident-only views are ",
                "<code>karyotype-current.svg</code> and ",
                "<code>karyotype-confident.svg</code>.</p></section>"
            )
        }
    } else {
        character()
    }

    c(
        "<!DOCTYPE html>",
        "<html lang=\"en\">",
        "<head>",
        "<meta charset=\"utf-8\">",
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">",
        paste0("<title>MaMaMIA results - ", bam_bundle_html_escape(sample),
            "</title>"),
        paste0("<style>", css, "</style>"),
        "</head>",
        "<body><main>",
        "<header>",
        "<p class=\"eyebrow\">MaMaMIA results bundle</p>",
        paste0("<h1>", bam_bundle_html_escape(sample), "</h1>"),
        paste0("<p class=\"lead\">", bam_bundle_html_escape(assembly_label),
            "</p>"),
        "</header>",
        "<section><h2>Run overview</h2>",
        bam_bundle_html_table(overview),
        "</section>",
        karyotype,
        "<section><h2>Calls by class</h2>",
        bam_bundle_html_table(diagnostics$calls_summary),
        "</section>",
        "<section><h2>Calls by chromosome group</h2>",
        bam_bundle_html_table(diagnostics$calls_by_group),
        "</section>",
        "<section><h2>Review status</h2>",
        bam_bundle_html_table(diagnostics$calls_by_status),
        "</section>",
        "<section><h2>Copy-number state</h2>",
        bam_bundle_html_table(diagnostics$calls_by_copy_state),
        "</section>",
        "<section><h2>Coverage and GC by group</h2>",
        bam_bundle_html_table(diagnostics$metrics_summary),
        "</section>",
        "<section><h2>Score summary by group</h2>",
        bam_bundle_html_table(diagnostics$scores_summary),
        "</section>",
        "<section><h2>Masked windows</h2>",
        bam_bundle_html_table(diagnostics$masks_summary),
        "</section>",
        "<section><h2>Calibration</h2>",
        bam_bundle_html_table(diagnostics$calibration),
        "</section>",
        if (is.data.frame(cohort) && nrow(cohort) > 1L) {
            c(
                "<section><h2>All samples in this run</h2>",
                bam_bundle_html_table(cohort),
                "</section>"
            )
        } else {
            character()
        },
        "<section><h2>Files in this bundle</h2>",
        "<p class=\"muted\">Every content file is listed here. ",
        "<code>manifest.json</code> also records sizes and MD5 checksums, ",
        "and <code>checksums.md5</code> covers the manifest as well.</p>",
        bam_bundle_html_table(file_table),
        "</section>",
        paste0(
            "<footer>Generated ", bam_bundle_html_escape(
                format(generated, "%Y-%m-%d %H:%M:%S %Z")
            ),
            " by MaMaMIA ", bam_bundle_html_escape(version), ".</footer>"
        ),
        "</main></body></html>"
    )
}
