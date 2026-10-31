## Downloadable result bundle: full call table, BEDPE files and static
## karyotype images in one ZIP archive. The table shown in the Calls tab is a
## short inspection view; this bundle keeps every column and every call.

#' Expand a filtered karyotype payload to every call
#'
#' `karyotype_export_payload()` already computes the unfiltered chromosome
#' layout and calls. Reuse them instead of rebuilding the full payload for the
#' ZIP.
#' @noRd
bam_full_karyotype_payload <- function(payload) {
    full <- payload
    if (is.list(payload$full_chromosomes) &&
        length(payload$full_chromosomes) > 0L) {
        full$chromosomes <- lapply(
            payload$full_chromosomes,
            function(chromosome) {
                chromosome$selected <- TRUE
                chromosome
            }
        )
    } else {
        calls <- payload$full_calls
        full$chromosomes <- lapply(
            seq_along(payload$chromosomes),
            function(i) {
                chromosome <- payload$chromosomes[[i]]
                if (is.list(calls) && length(calls) >= i &&
                    !is.null(calls[[i]])) {
                    chromosome$calls <- calls[[i]]
                }
                chromosome$selected <- TRUE
                chromosome
            }
        )
    }
    full$full_calls <- NULL
    full$full_chromosomes <- NULL
    full
}

#' Descriptive column dictionary for the result bundle
#'
#' The CSV and BEDPE files in the bundle use readable column names. The
#' dictionary is also rendered into README.txt so every field can be looked up
#' without the package source.
#' @noRd
bam_bundle_column_dictionary <- function() {
    data.frame(
        old = c(
            "sample", "group", "start_iid", "end_iid", "start_bp", "end_bp",
            "width_bp", "score_mean", "donor_mean", "recipient_mean",
            "total_mean", "donor_component_spread",
            "recipient_component_spread", "n_donor_supported",
            "mapq_fraction", "n_windows", "method", "mode", "metric",
            "length_bin", "threshold", "pass", "class", "class_evidence",
            "copy_state", "total_evidence", "class_confidence",
            "class_status", "component_status", "mappability_status",
            "confidence", "class_group", "recipient_component",
            "recipient_chr_id", "recipient_chr_name", "recipient_group",
            "recipient_start_bp", "recipient_end_bp",
            "recipient_run_start_bp", "recipient_run_end_bp",
            "recipient_deletion_width_bp", "recipient_expected_breadth",
            "recipient_observed_breadth", "recipient_mapq_fraction",
            "recipient_copy_number", "recipient_mean_depletion",
            "recipient_depletion_se", "recipient_depletion_lower",
            "recipient_depletion_z", "recipient_width_overlap",
            "recipient_fraction_ratio", "donor_copies", "copy_state_group",
            "affected_fraction", "affected_copies", "dosage_confidence",
            "recipient_weight", "recipient_margin",
            "recipient_assignment", "recipient_cross_group",
            "recipient_reliable", "recipient_strong",
            "recipient_alternative", "recipient_alternative_score",
            "donor_component", "donor_chr_id", "donor_chr_name",
            "donor_start_bp", "donor_end_bp"
        ),
        name = c(
            "sample_name", "chromosome_group", "start_window_index",
            "end_window_index", "group_start_bp", "group_end_bp",
            "call_width_bp", "mean_introgression_score",
            "mean_donor_residual", "mean_recipient_residual",
            "mean_total_depth_residual", "donor_component_spread",
            "recipient_component_spread", "supporting_donor_components",
            "high_mapq_read_fraction", "windows_in_call",
            "scoring_method", "analysis_mode", "coverage_metric",
            "call_length_bin", "calibration_threshold",
            "passes_calibration", "call_class", "call_class_evidence",
            "copy_number_state", "total_depth_evidence", "call_class_confidence",
            "review_status", "component_consistency", "mappability_quality",
            "donor_confidence", "group_level_call_class",
            "recipient_copy_set", "recipient_contig_id",
            "recipient_chromosome", "recipient_chromosome_group",
            "recipient_start_bp", "recipient_end_bp",
            "recipient_deletion_run_start_bp",
            "recipient_deletion_run_end_bp",
            "recipient_deletion_run_width_bp",
            "recipient_expected_breadth", "recipient_observed_breadth",
            "recipient_high_mapq_read_fraction",
            "recipient_component_copy_number",
            "recipient_mean_depletion",
            "recipient_depletion_standard_error",
            "recipient_depletion_lower_bound",
            "recipient_depletion_z_score",
            "recipient_width_reciprocal_overlap",
            "donor_to_recipient_width_ratio", "estimated_donor_copies",
            "group_level_copy_state", "affected_copy_fraction",
            "affected_copy_number", "dosage_estimate_confidence",
            "recipient_assignment_weight", "recipient_assignment_margin",
            "recipient_assignment_status", "recipient_is_cross_group",
            "recipient_candidate_reliable", "recipient_candidate_strong",
            "alternative_recipient_chromosome",
            "alternative_recipient_score", "donor_copy_set",
            "donor_contig_id", "donor_chromosome",
            "donor_start_bp", "donor_end_bp"
        ),
        description = c(
            "Sample name read from the BAM @RG SM tag.",
            "Chromosome homology group used for scoring and segmentation, for example 7.",
            "First analysis window index in the call, zero-based within the group.",
            "Last analysis window index in the call, zero-based and inclusive.",
            "Call start in group-level coordinates, in base pairs, scaled to the donor group length.",
            "Call end in group-level coordinates, in base pairs.",
            "Call width in base pairs, calculated as end_bp minus start_bp.",
            "Mean introgression score across the call; higher values mean stronger donor evidence.",
            "Mean donor residual across the call; positive values mean donor enrichment.",
            "Mean recipient residual across the call; negative values mean recipient loss.",
            "Mean total-depth residual across the call; positive values mean a total gain and negative values a total loss.",
            "Spread of donor residuals across donor copy sets; large values flag component heterogeneity.",
            "Spread of recipient residuals across recipient copy sets; large values flag component heterogeneity.",
            "Number of donor copy sets that support the call.",
            "Mean fraction of high-MAPQ reads across the call; values below 0.5 flag low mappability.",
            "Number of analysis windows covered by the call.",
            "Scoring method: paired_lr for a paired donor/recipient log ratio, or donor_only.",
            "Analysis mode: single for one sample, or multi for a cohort.",
            "Coverage metric used for scoring: breadth, count or combined.",
            "Segment-length bin used to select the calibration threshold.",
            "Calibrated score threshold for this call's mode and length bin.",
            "TRUE when score_mean exceeds the calibrated threshold.",
            "Final call class after recipient assignment: substitution, addition, deletion, complex or other.",
            "Evidence that supported the call: donor, recipient, donor+recipient, total imbalance or none.",
            "Copy-number interpretation: balanced, gain, loss, compound or neutral.",
            "Total-depth evidence: gain, loss, neutral, or blank when unavailable.",
            "Class-level confidence: high, moderate or low.",
            "Review status: confident, ambiguous or needs_review.",
            "Component consistency: consistent or heterogeneous.",
            "Mappability flag: ok or low; low means the high-MAPQ read fraction is below 0.5.",
            "Legacy donor-effect confidence: high, moderate or low.",
            "Group-level call class before recipient assignment.",
            "Copy-set label of the assigned recipient contig.",
            "Assembly contig identifier of the assigned recipient.",
            "Display name of the assigned recipient contig.",
            "Homology group of the assigned recipient contig; differs from chromosome_group for cross-group placements.",
            "Start of the corresponding interval on the assigned recipient contig, in base pairs.",
            "End of that interval on the assigned recipient contig, in base pairs.",
            "Start of the contiguous recipient depletion run used for assignment, in base pairs.",
            "End of that depletion run, in base pairs.",
            "Width of the recipient depletion run, in base pairs; blank when no run was found.",
            "Baseline breadth expected for the recipient contig.",
            "Mean breadth observed inside the call interval on the recipient contig.",
            "Mean high-MAPQ read fraction in the recipient interval.",
            "Configured copy number of the recipient component, default 2.",
            "Mean corrected log-depth depletion in the recipient interval; positive values mean loss.",
            "Standard error of the mean recipient depletion.",
            "Lower confidence bound of the recipient depletion estimate.",
            "Z score of the recipient depletion estimate.",
            "Reciprocal overlap between the donor segment and the recipient depletion run, from 0 to 1.",
            "Donor width fraction divided by recipient run width fraction; values near 1 mean matched widths.",
            "Estimated number of donor copies in the call.",
            "Copy state before recipient assignment.",
            "Estimated fraction of recipient copies affected by the event.",
            "Estimated number of recipient copies affected; blank when the estimate is too uncertain.",
            "Dosage confidence: high, moderate, low or unavailable.",
            "Weight of the chosen recipient among plausible candidates, from 0 to 1.",
            "Score difference between the chosen recipient and the best alternative; larger values are more decisive.",
            "Recipient placement outcome: assigned, ambiguous, cross_group or no_recipient_loss.",
            "TRUE when the assigned recipient contig is outside the call's homology group.",
            "TRUE when the recipient candidate passed the calibrated reliability screen.",
            "TRUE when the recipient candidate passed the stronger depletion screen.",
            "Display name of the next-best recipient contig, when available.",
            "Assignment score of the next-best recipient contig.",
            "Copy-set label of the donor contig that contributed the call.",
            "Assembly contig identifier of the donor contig.",
            "Display name of the donor contig.",
            "Donor segment start in base pairs.",
            "Donor segment end in base pairs."
        ),
        stringsAsFactors = FALSE
    )
}

#' BEDPE-only aliases to descriptive column names
#' @noRd
bam_bundle_bedpe_aliases <- function() {
    c(
        "chr_id.don" = "donor_contig_id",
        "start.don" = "donor_start_bp",
        "stop.don" = "donor_end_bp",
        "chr_id.rec" = "recipient_contig_id",
        "start.rec" = "recipient_start_bp",
        "stop.rec" = "recipient_end_bp"
    )
}

#' Rename call-table columns with the descriptive bundle names
#' @noRd
bam_bundle_descriptive_calls <- function(calls) {
    dictionary <- bam_bundle_column_dictionary()
    aliases <- bam_bundle_bedpe_aliases()
    old <- names(calls)
    aliased <- aliases[old]
    old[!is.na(aliased)] <- aliased[!is.na(aliased)]
    replacement <- dictionary$name[match(old, dictionary$old)]
    replacement[is.na(replacement)] <- old[is.na(replacement)]
    names(calls) <- replacement
    calls
}

#' Confident calls used by the default plot
#'
#' A confident call passes calibration and has class status `confident`.
#' @noRd
bam_confident_calls <- function(calls) {
    if (is.null(calls) || nrow(calls) == 0L) {
        return(calls)
    }
    keep <- rep(TRUE, nrow(calls))
    if ("pass" %in% names(calls)) {
        keep <- keep & isTRUE_vec(calls$pass)
    }
    if ("class_status" %in% names(calls)) {
        keep <- keep &
            as.character(calls$class_status) %in% "confident"
    } else if ("class_confidence" %in% names(calls)) {
        keep <- keep &
            as.character(calls$class_confidence) %in% "high"
    }
    calls[keep, , drop = FALSE]
}

#' Brief call-table columns for the Calls tab
#'
#' The brief view is meant to answer "what did we find and where" without
#' horizontal scrolling. The full table is available in the result bundle.
#' @noRd
bam_brief_call_columns <- function(x = NULL) {
    columns <- c(
        "group", "class", "copy_state", "class_status",
        "recipient_chr_name", "recipient_assignment",
        "start_bp", "end_bp", "width_bp", "pass"
    )
    include_sample <- is.data.frame(x) &&
        "sample" %in% names(x) &&
        length(unique(as.character(x$sample))) > 1L
    if (include_sample) {
        columns <- c("sample", columns)
    }
    columns
}

#' Human-readable labels for the brief call table
#' @noRd
bam_brief_call_labels <- function(columns) {
    labels <- c(
        sample = "Sample",
        group = "Group",
        class = "Class",
        copy_state = "Copy state",
        class_status = "Status",
        recipient_chr_name = "Recipient",
        recipient_assignment = "Placement",
        start_bp = "Group start (bp)",
        end_bp = "Group end (bp)",
        width_bp = "Width (bp)",
        pass = "Pass"
    )
    unname(labels[columns])
}

#' Safe stem for a bundle filename
#' @noRd
bam_bundle_stem <- function(sample) {
    sample <- as.character(sample)
    if (length(sample) != 1L || is.na(sample) || !nzchar(sample)) {
        return("sample")
    }
    safe <- gsub("[^A-Za-z0-9._-]+", "-", sample)
    safe <- gsub("^[.-]+|[.-]+$", "", safe)
    if (!nzchar(safe)) "sample" else safe
}

#' ZIP filename for the result bundle
#' @noRd
bam_bundle_filename <- function(sample = NULL) {
    paste0("mamamia-", bam_bundle_stem(sample), "-results.zip")
}

#' Merge payload palette with the package defaults
#' @noRd
bam_bundle_palette <- function(payload) {
    palette <- karyotype_palette()
    supplied <- payload$palette
    if (is.list(supplied)) {
        for (name in names(supplied)) {
            value <- supplied[[name]]
            if (is.character(value) && length(value) == 1L &&
                !is.na(value) && nzchar(value)) {
                palette[[name]] <- value
            }
        }
    }
    palette
}

#' Escape text for an SVG document
#' @noRd
bam_svg_escape <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub("\"", "&quot;", x, fixed = TRUE)
}

#' Format a finite SVG number without scientific notation
#' @noRd
bam_svg_number <- function(x, digits = 2) {
    x <- suppressWarnings(as.numeric(x))
    if (length(x) != 1L || !is.finite(x)) {
        return("0")
    }
    formatC(x, format = "f", digits = digits)
}

#' SVG text element
#' @noRd
bam_svg_text <- function(text, x, y, size = 11, weight = 400,
                         fill = "#3a3d35", anchor = "start") {
    sprintf(
        paste0(
            '<text x="%s" y="%s" font-family="Inter, system-ui, ',
            'sans-serif" font-size="%s" font-weight="%s" fill="%s" ',
            'text-anchor="%s">%s</text>'
        ),
        bam_svg_number(x), bam_svg_number(y), bam_svg_number(size, 1),
        bam_svg_number(weight, 0), bam_svg_escape(fill), anchor,
        bam_svg_escape(text)
    )
}

#' SVG rectangle element
#' @noRd
bam_svg_rect <- function(x, y, width, height, fill, stroke = "none",
                         stroke_width = 0.7, opacity = NULL,
                         dash = NULL, rx = 2) {
    dash_attr <- if (is.null(dash)) {
        ""
    } else {
        sprintf(' stroke-dasharray="%s"', bam_svg_escape(dash))
    }
    opacity_attr <- if (is.null(opacity)) {
        ""
    } else {
        sprintf(' opacity="%s"', bam_svg_number(opacity, 2))
    }
    sprintf(
        paste0(
            '<rect x="%s" y="%s" width="%s" height="%s" rx="%s" ',
            'fill="%s" stroke="%s" stroke-width="%s"%s%s/>'
        ),
        bam_svg_number(x), bam_svg_number(y), bam_svg_number(width),
        bam_svg_number(height), bam_svg_number(rx, 1),
        bam_svg_escape(fill), bam_svg_escape(stroke),
        bam_svg_number(stroke_width, 2), dash_attr, opacity_attr
    )
}

#' Pick a palette colour by name
#' @noRd
bam_bundle_colour <- function(palette, name, fallback = "other") {
    name <- as.character(name)
    if (length(name) != 1L || is.na(name) || !nzchar(name)) {
        name <- fallback
    }
    value <- palette[[name]]
    if (is.null(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
        value <- palette[[fallback]]
    }
    if (is.null(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
        value <- "#7570b3"
    }
    value
}

#' Static call legend for the bundle image
#' @noRd
bam_svg_legend <- function(palette, x, y, width = NULL) {
    entries <- list(
        c("substitution", "Substitution"),
        c("addition", "Addition"),
        c("deletion", "Deletion"),
        c("complex", "Complex"),
        c("other", "Other"),
        c("deleted", "Deleted recipient"),
        c("added", "Added donor")
    )
    item_widths <- vapply(entries, function(entry) {
        18 + nchar(entry[[2L]], type = "width") * 5.6
    }, numeric(1))
    total_width <- sum(item_widths) + 16 * (length(entries) - 1L)
    if (!is.null(width)) {
        x <- x + max(0, (as.numeric(width)[[1L]] - total_width) / 2)
    }
    parts <- character()
    cursor <- x
    for (i in seq_along(entries)) {
        entry <- entries[[i]]
        key <- entry[[1L]]
        label <- entry[[2L]]
        dash <- if (identical(key, "deleted")) "3 2" else NULL
        stroke <- if (identical(key, "deleted")) {
            palette$negative
        } else {
            palette$text
        }
        parts <- c(
            parts,
            bam_svg_rect(
                cursor, y - 9, 12, 12,
                fill = bam_bundle_colour(palette, key),
                stroke = stroke, stroke_width = 0.6,
                dash = dash, rx = 1.5
            ),
            bam_svg_text(
                label, cursor + 18, y,
                size = 10, weight = 500, fill = palette$text
            )
        )
        cursor <- cursor + item_widths[[i]] + 16
    }
    paste(parts, collapse = "")
}

#' Opacity for a call segment
#' @noRd
bam_svg_call_opacity <- function(call) {
    confidence <- call$class_confidence
    if (is.null(confidence) || length(confidence) != 1L ||
        is.na(confidence) || !nzchar(as.character(confidence))) {
        confidence <- call$confidence
    }
    confidence <- as.character(confidence)
    opacity <- if (identical(confidence, "high")) {
        0.95
    } else if (identical(confidence, "moderate")) {
        0.72
    } else {
        0.45
    }
    if (identical(call$pass, FALSE)) {
        opacity <- opacity * 0.45
    }
    assignment <- as.character(call$recipient_assignment)
    if (identical(assignment, "ambiguous")) {
        opacity <- opacity * 0.55
    }
    if (identical(assignment, "no_recipient_loss")) {
        opacity <- opacity * 0.35
    }
    opacity
}

#' Truncate text to fit a fixed character budget
#' @noRd
bam_svg_truncate <- function(text, max_chars) {
    text <- as.character(text)
    if (length(text) != 1L || is.na(text)) {
        return("")
    }
    max_chars <- max(3L, as.integer(max_chars)[[1L]])
    if (nchar(text, type = "width") <= max_chars) {
        return(text)
    }
    paste0(substr(text, 1L, max_chars - 3L), "...")
}

#' Static SVG karyotype drawn from a karyotype payload
#'
#' The interactive widget exports SVG in the browser. The bundle is assembled
#' on the server, so it uses this compact static rendering of the same payload.
#'
#' @param payload A karyotype payload.
#' @param width,height Image dimensions in pixels.
#' @param title,subtitle Optional headings.
#' @param footer Optional footer text.
#'
#' @return An SVG document as a single string.
#' @noRd
bam_karyotype_static_svg <- function(payload, width = 1600, height = 900,
                                     title = "MaMaMIA karyotype",
                                     subtitle = NULL, footer = NULL) {
    width <- max(640, as.numeric(width)[[1L]])
    height <- max(420, as.numeric(height)[[1L]])
    palette <- bam_bundle_palette(payload)
    chromosomes <- payload$chromosomes
    sample <- if (is.null(payload$sample)) "" else as.character(payload$sample)
    if (is.null(subtitle) || !nzchar(subtitle)) {
        subtitle <- paste0(
            if (nzchar(sample)) paste0("Sample: ", sample) else "Sample",
            " \u00b7 Calls"
        )
    }
    if (is.null(footer)) {
        footer <- "Generated by MaMaMIA"
    }
    max_title_chars <- floor((width - 120) / (18 * 0.62))
    max_subtitle_chars <- floor((width - 120) / (11 * 0.62))
    max_footer_chars <- floor((width - 120) / (9 * 0.62))
    title <- bam_svg_truncate(title, max_title_chars)
    subtitle <- bam_svg_truncate(subtitle, max_subtitle_chars)

    header <- c(
        sprintf(
            paste0(
                '<svg xmlns="http://www.w3.org/2000/svg" width="%s" ',
                'height="%s" viewBox="0 0 %s %s" ',
                'preserveAspectRatio="xMidYMid meet" role="img" ',
                'aria-label="MaMaMIA karyotype">'
            ),
            bam_svg_number(width, 0), bam_svg_number(height, 0),
            bam_svg_number(width, 0), bam_svg_number(height, 0)
        ),
        bam_svg_rect(
            0, 0, width, height, fill = "#ffffff", stroke = "none",
            stroke_width = 0, rx = 0
        ),
        bam_svg_rect(
            0.5, 0.5, width - 1, height - 1, fill = "none",
            stroke = palette$grid, stroke_width = 1, rx = 0
        ),
        bam_svg_text(
            title, width / 2, 36, size = 18, weight = 700,
            fill = palette$text, anchor = "middle"
        ),
        bam_svg_text(
            subtitle, width / 2, 58, size = 11, weight = 500,
            fill = palette$text, anchor = "middle"
        ),
        bam_svg_legend(palette, 0, 86, width = width)
    )

    if (is.null(chromosomes) || length(chromosomes) == 0L) {
        body <- bam_svg_text(
            "No chromosomes in this result.", width / 2, height / 2,
            size = 14, weight = 500, fill = palette$text, anchor = "middle"
        )
        return(paste0(
            '<?xml version="1.0" encoding="UTF-8"?>\n',
            paste(c(header, body, "</svg>"), collapse = ""),
            "\n"
        ))
    }

    margin_left <- 52
    margin_right <- 52
    top <- 118
    bottom <- 92
    plot_width <- max(1, width - margin_left - margin_right)
    plot_height <- max(1, height - top - bottom)
    n <- length(chromosomes)
    lengths <- vapply(chromosomes, function(x) {
        value <- suppressWarnings(as.numeric(x$length))
        if (length(value) != 1L || !is.finite(value) || value <= 0) {
            0
        } else {
            value
        }
    }, numeric(1))
    max_length <- max(c(1, lengths))
    pair_gap <- 3
    ids <- vapply(chromosomes, function(x) {
        value <- x$id
        if (is.null(value) || length(value) != 1L || is.na(value)) {
            ""
        } else {
            as.character(value)
        }
    }, character(1))
    added_flags <- vapply(
        chromosomes, function(x) isTRUE(x$added), logical(1)
    )
    copies_vec <- vapply(chromosomes, function(x) {
        value <- suppressWarnings(as.integer(x$copies))
        if (length(value) != 1L || is.na(value) || value < 1L) {
            1L
        } else {
            value
        }
    }, integer(1))
    same_pair <- c(
        FALSE,
        ids[-1L] == ids[-n] & added_flags[-1L] == added_flags[-n]
    )
    group_index <- cumsum(!same_pair)
    n_groups <- max(group_index)
    slot <- plot_width / max(1, n_groups)
    max_copies <- max(c(1L, copies_vec))
    bar_width <- max(
        7, min(34, (slot - pair_gap * (max_copies - 1L) - 10) / max_copies)
    )
    label_size <- if (bar_width < 14) 8 else 10
    label_max <- max(3L, floor((slot - 4) / (label_size * 0.62)))

    body <- character()
    for (i in seq_along(chromosomes)) {
        chromosome <- chromosomes[[i]]
        added <- isTRUE(chromosome$added)
        deleted <- !added && isTRUE(chromosome$deleted)
        selected <- !identical(chromosome$selected, FALSE)
        calls <- chromosome$calls
        if (is.null(calls)) {
            calls <- list()
        }
        copies <- copies_vec[[i]]
        copy <- suppressWarnings(as.integer(chromosome$copy))
        if (length(copy) != 1L || is.na(copy) || copy < 1L) {
            copy <- 1L
        }
        pair_width <- copies * bar_width + (copies - 1L) * pair_gap
        group_left <- margin_left + (group_index[[i]] - 1L) * slot +
            (slot - pair_width) / 2
        x <- group_left + (copy - 1L) * (bar_width + pair_gap)
        pair_center <- margin_left + (group_index[[i]] - 0.5) * slot
        bar_height <- max(3, plot_height * lengths[[i]] / max_length)
        added_class <- if (added && length(calls) == 1L) {
            calls[[1L]]$class
        } else {
            NULL
        }
        base_fill <- if (deleted) {
            palette$deleted
        } else if (added) {
            bam_bundle_colour(palette, added_class, "added")
        } else {
            "#f4f4f4"
        }
        base_stroke <- if (deleted) {
            palette$negative
        } else if (added) {
            bam_bundle_colour(palette, added_class, "added")
        } else {
            "#666666"
        }
        group <- sprintf(
            '<g%s>',
            if (selected) "" else ' opacity="0.32"'
        )
        group <- paste0(
            group,
            bam_svg_rect(
                x, top, bar_width, bar_height,
                fill = base_fill, stroke = base_stroke,
                stroke_width = if (deleted) 1.4 else if (added) 1.8 else 0.7,
                dash = if (deleted) "4 3" else NULL
            )
        )
        for (call in calls) {
            start <- suppressWarnings(as.numeric(call$start))
            end <- suppressWarnings(as.numeric(call$end))
            if (length(start) != 1L || !is.finite(start)) start <- 0
            if (length(end) != 1L || !is.finite(end)) end <- 1
            start <- max(0, min(1, start))
            end <- max(0, min(1, end))
            if (end <= start) {
                next
            }
            call_y <- top + start * bar_height
            call_height <- max(1.5, (end - start) * bar_height)
            cross <- isTRUE(call$cross_group)
            group <- paste0(
                group,
                bam_svg_rect(
                    x, call_y, bar_width, call_height,
                    fill = bam_bundle_colour(palette, call$class, "other"),
                    stroke = if (cross) "#333333" else "#ffffff",
                    stroke_width = if (cross) 1.3 else 0.5,
                    opacity = bam_svg_call_opacity(call),
                    dash = if (cross) "3 2" else NULL,
                    rx = 0
                )
            )
        }
        if (deleted && copy == 1L) {
            group <- paste0(
                group,
                bam_svg_text(
                    "Deleted", pair_center, top - 8,
                    size = 9, weight = 700, fill = palette$negative,
                    anchor = "middle"
                )
            )
        }
        if (added && copy == 1L) {
            group <- paste0(
                group,
                bam_svg_text(
                    "Added", pair_center, top - 8,
                    size = 9, weight = 700, fill = palette$addition,
                    anchor = "middle"
                )
            )
        }
        if (copy == 1L || copies == 1L) {
            label <- chromosome$label
            if (is.null(label) || !nzchar(as.character(label))) {
                label <- chromosome$name
            }
            if (is.null(label) || !nzchar(as.character(label))) {
                label <- chromosome$component
            }
            label <- bam_svg_truncate(label, label_max)
            group <- paste0(
                group,
                bam_svg_text(
                    label, pair_center, top + plot_height + 22,
                    size = label_size, weight = 600, fill = palette$text,
                    anchor = "middle"
                )
            )
        }
        group <- paste0(group, "</g>")
        body <- c(body, group)
    }

    footer_text <- if (nzchar(sample)) {
        paste0(footer, " \u00b7 ", sample)
    } else {
        footer
    }
    footer_text <- bam_svg_truncate(footer_text, max_footer_chars)
    footer_parts <- bam_svg_text(
        footer_text, width / 2, height - 26,
        size = 9, weight = 400, fill = palette$text,
        anchor = "middle"
    )
    paste0(
        '<?xml version="1.0" encoding="UTF-8"?>\n',
        paste(c(header, body, footer_parts, "</svg>"), collapse = ""),
        "\n"
    )
}

#' Write a BEDPE file with descriptive column names
#' @noRd
bam_bundle_write_bedpe <- function(calls, windows, file, pass_only = FALSE) {
    raw <- tempfile(fileext = ".bedpe")
    on.exit(unlink(raw), add = TRUE)
    writeCallsToBedpe(
        calls, raw, pass_only = pass_only, windows = windows
    )
    table <- utils::read.delim(
        raw, check.names = FALSE, stringsAsFactors = FALSE,
        quote = "", comment.char = ""
    )
    table <- bam_bundle_descriptive_calls(table)
    utils::write.table(
        table, file, sep = "\t", quote = FALSE,
        row.names = FALSE, na = ""
    )
    names(table)
}

#' README for a result bundle
#' @noRd
bam_bundle_readme <- function(sample, calls, generated, version,
                              bedpe_columns = character(),
                              confident_calls = NULL) {
    n <- nrow(calls)
    passing <- if ("passes_calibration" %in% names(calls)) {
        sum(isTRUE_vec(calls$passes_calibration), na.rm = TRUE)
    } else if ("pass" %in% names(calls)) {
        sum(isTRUE_vec(calls$pass), na.rm = TRUE)
    } else {
        n
    }
    confident_n <- if (is.null(confident_calls)) {
        0L
    } else {
        nrow(confident_calls)
    }
    timestamp <- if (inherits(generated, "POSIXt")) {
        format(generated, "%Y-%m-%d %H:%M:%S %Z")
    } else {
        as.character(generated)[[1L]]
    }
    dictionary <- bam_bundle_column_dictionary()
    dictionary$in_csv <- dictionary$name %in% names(calls)
    dictionary$in_bedpe <- dictionary$name %in% bedpe_columns
    rows <- vapply(seq_len(nrow(dictionary)), function(i) {
        sprintf(
            "| `%s` | %s | %s | %s |",
            dictionary$name[[i]],
            if (dictionary$in_csv[[i]]) "yes" else "-",
            if (dictionary$in_bedpe[[i]]) "yes" else "-",
            dictionary$description[[i]]
        )
    }, character(1))
    c(
        "MaMaMIA results bundle",
        "",
        paste0("Sample: ", sample),
        paste0("Generated: ", timestamp),
        paste0("Package: MaMaMIA ", version),
        paste0(
            "Calls: ", n, " total, ", passing,
            " passing calibration, ", confident_n, " confident"
        ),
        "",
        "Contents",
        "--------",
        "calls-full.csv       Every call for this sample with every column.",
        "calls-confident.csv  Passing calls with review_status = confident.",
        "calls-all.bedpe      Every call for this sample in BEDPE format.",
        "calls-passing.bedpe  Calibrated calls (pass = TRUE) in BEDPE format.",
        "calls-confident.bedpe Confident calls in BEDPE format.",
        "karyotype-full.svg   Every call for the sample, without display filters.",
        "karyotype-current.svg Calls for the sample after the active filters.",
        "karyotype-confident.svg Confident calls only.",
        "",
        "Confident set",
        "-------------",
        "The confident set contains calls that pass calibration and have",
        "review_status = confident. The Karyotype and Calls tabs use this",
        "set by default; clear the Status filter to see ambiguous and",
        "needs-review calls. The full CSV and BEDPE keep every call.",
        "",
        "Column conventions",
        "------------------",
        "Columns ending in _bp are base pairs. Fractions, breadths and",
        "overlaps are between 0 and 1. Residuals and scores are on the",
        "log scale unless the description says otherwise. A blank value",
        "means unavailable or not applicable.",
        "",
        "group_start_bp and group_end_bp are in the homology-group",
        "coordinate space used for scoring. Use recipient_start_bp and",
        "recipient_end_bp for the assigned recipient contig, and",
        "donor_start_bp and donor_end_bp for the donor segment.",
        "",
        "Column descriptions",
        "-------------------",
        "BEDPE files are tab-separated. They use the same descriptive names,",
        "and calls-passing.bedpe contains only rows where",
        "passes_calibration is TRUE. The CSV and BEDPE columns below show",
        "whether a field is present in each file.",
        "",
        "| Column | CSV | BEDPE | Description |",
        "|---|---|---|---|",
        rows,
        ""
    )
}

#' Write the full result bundle
#'
#' @param result A `BamAnalysisResult`.
#' @param file Output ZIP path.
#' @param sample Sample to export.
#' @param current_payload Optional karyotype payload for the filtered view.
#' @param full_payload Optional karyotype payload with every call.
#' @param confident_payload Optional karyotype payload with confident calls.
#' @param generated Timestamp shown in the README.
#'
#' @return The output path, invisibly.
#' @noRd
bam_write_result_bundle <- function(result, file, sample,
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
    sample <- as.character(sample)
    if (length(sample) != 1L || is.na(sample) || !nzchar(sample)) {
        sample <- karyotype_sample(result, NULL)
    }
    calls_raw <- result$calls
    calls_raw <- calls_raw[
        as.character(calls_raw$sample) == sample, , drop = FALSE
    ]
    confident_raw <- bam_confident_calls(calls_raw)
    calls <- bam_bundle_descriptive_calls(calls_raw)
    confident_calls <- bam_bundle_descriptive_calls(confident_raw)
    staging <- tempfile("mamamia-bundle-")
    dir.create(staging)
    old <- getwd()
    on.exit({
        setwd(old)
        unlink(staging, recursive = TRUE, force = TRUE)
    }, add = TRUE)
    bundle_name <- paste0("mamamia-", bam_bundle_stem(sample), "-results")
    bundle <- file.path(staging, bundle_name)
    dir.create(bundle)

    utils::write.csv(
        calls, file.path(bundle, "calls-full.csv"),
        row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
    utils::write.csv(
        confident_calls, file.path(bundle, "calls-confident.csv"),
        row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
    bedpe_columns <- bam_bundle_write_bedpe(
        calls_raw, result$windows, file.path(bundle, "calls-all.bedpe"),
        pass_only = FALSE
    )
    bam_bundle_write_bedpe(
        calls_raw, result$windows, file.path(bundle, "calls-passing.bedpe"),
        pass_only = TRUE
    )
    bam_bundle_write_bedpe(
        confident_raw, result$windows,
        file.path(bundle, "calls-confident.bedpe"),
        pass_only = FALSE
    )
    version <- tryCatch(
        as.character(utils::packageVersion("MaMaMIA")),
        error = function(e) ""
    )
    if (!is.null(current_payload)) {
        writeLines(
            bam_karyotype_static_svg(
                current_payload,
                title = "MaMaMIA karyotype",
                subtitle = paste0(
                    "Sample: ", sample, " \u00b7 Filtered calls"
                ),
                footer = paste0("MaMaMIA ", version)
            ),
            file.path(bundle, "karyotype-current.svg"),
            useBytes = TRUE
        )
    }
    if (!is.null(full_payload)) {
        writeLines(
            bam_karyotype_static_svg(
                full_payload,
                title = "MaMaMIA karyotype",
                subtitle = paste0("Sample: ", sample, " \u00b7 All calls"),
                footer = paste0("MaMaMIA ", version)
            ),
            file.path(bundle, "karyotype-full.svg"),
            useBytes = TRUE
        )
    }
    if (!is.null(confident_payload)) {
        writeLines(
            bam_karyotype_static_svg(
                confident_payload,
                title = "MaMaMIA karyotype",
                subtitle = paste0(
                    "Sample: ", sample, " \u00b7 Confident calls"
                ),
                footer = paste0("MaMaMIA ", version)
            ),
            file.path(bundle, "karyotype-confident.svg"),
            useBytes = TRUE
        )
    }
    writeLines(
        bam_bundle_readme(
            sample, calls, generated, version,
            bedpe_columns = bedpe_columns,
            confident_calls = confident_calls
        ),
        file.path(bundle, "README.txt"),
        useBytes = TRUE
    )

    if (file.exists(file)) {
        unlink(file)
    }
    setwd(staging)
    status <- tryCatch(
        utils::zip(
            zipfile = file, files = bundle_name, flags = "-q -r"
        ),
        error = function(e) e
    )
    setwd(old)
    if (inherits(status, "error") || !identical(as.integer(status), 0L)) {
        stop(
            "Could not create the results ZIP archive. ",
            "Check that the `zip` command is available.",
            call. = FALSE
        )
    }
    invisible(file)
}
