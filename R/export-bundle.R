## Downloadable result bundle: full call table, BEDPE files and static
## karyotype images in one ZIP archive. The table shown in the Calls tab is a
## short inspection view; this bundle keeps every column and every call.

#' Return the self-contained full karyotype view
#'
#' The current contract requires `payload$views$full`, a complete payload
#' whose chromosomes own their calls and scores. A payload without that view
#' has no full export; the old positional `full_chromosomes`/`full_calls` side
#' channels are no longer read or reconstructed.
#'
#' @return The full payload, or `NULL` when the payload carries no full view.
#' @noRd
bam_full_karyotype_payload <- function(payload) {
    view <- payload$views$full
    if (!is.list(view) || !is.list(view$chromosomes) ||
        length(view$chromosomes) == 0L) {
        return(NULL)
    }
    full <- view
    full$full_calls <- NULL
    full$full_chromosomes <- NULL
    full$views <- NULL
    ## The full export is evidence, not the hidden live view.
    full$hide_deletions <- FALSE
    full$hidden_deletions <- list(segments = 0L, bases = 0)
    full$scene <- karyotype_scene(full)
    full
}

#' Descriptive column dictionary for the result bundle
#'
#' The CSV and BEDPE files in the bundle use readable column names. The
#' dictionary is also rendered into README.txt so every field can be looked up
#' without the application source.
#' @noRd
bam_bundle_column_dictionary <- function() {
    schema <- call_schema()
    data.frame(
        old = schema$name,
        name = schema$export,
        description = schema$description,
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

#' Confident calls used by the confident export
#'
#' A confident call passes calibration and has class status `confident`.
#' @noRd
bam_confident_calls <- function(calls) {
    if (is.null(calls) || nrow(calls) == 0L) {
        return(calls)
    }
    calls <- explorer_calls(calls, list(pass_only = TRUE))$calls
    if (nrow(calls) == 0L) {
        return(calls)
    }
    keep <- rep(TRUE, nrow(calls))
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
        "recipient_chr_name", "donor_chr_name", "recipient_assignment",
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

#' Human-readable labels for call-table columns
#' @noRd
bam_call_column_labels <- function(columns = NULL) {
    labels <- call_schema_labels()
    if (is.null(columns)) {
        return(labels)
    }
    out <- unname(labels[columns])
    out[is.na(out)] <- columns[is.na(out)]
    out
}

#' Column metadata for the Calls tab
#'
#' Combines short display labels with the bundle dictionary descriptions.
#' @param columns Call-table columns; every dictionary column by default.
#' @noRd
bam_call_column_metadata <- function(columns = NULL) {
    dictionary <- bam_bundle_column_dictionary()
    if (is.null(columns)) {
        columns <- dictionary$old
    }
    index <- match(columns, dictionary$old)
    description <- dictionary$description[index]
    description[is.na(description)] <- "No description available."
    data.frame(
        column = columns,
        label = bam_call_column_labels(columns),
        description = description,
        stringsAsFactors = FALSE
    )
}

#' Selectable extra columns for the Calls tab
#'
#' The brief inspection columns stay visible. Every other diagnostic column
#' can be added from the picker. `sample` is omitted because it appears
#' automatically when the result contains more than one sample.
#' @noRd
bam_call_extra_column_choices <- function() {
    metadata <- bam_call_column_metadata()
    defaults <- bam_brief_call_columns()
    extra <- metadata[
        !metadata$column %in% c(defaults, "sample"),
        ,
        drop = FALSE
    ]
    stats::setNames(extra$column, extra$label)
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
    "mamamia-results.zip"
}

#' ZIP filename for a multi-sample cohort bundle
#' @noRd
bam_cohort_bundle_filename <- function() {
    bam_bundle_filename()
}

#' MaMaMIA version for bundle metadata
#' @noRd
bam_bundle_version <- function() {
    tryCatch(mamamia_version(), error = function(e) "")
}

#' Write one karyotype SVG into a bundle
#' @noRd
bam_bundle_write_karyotype <- function(payload, path, sample, label,
                                       version) {
    if (is.null(payload)) {
        return(invisible(NULL))
    }
    writeLines(
        bam_karyotype_static_svg(
            payload,
            title = "MaMaMIA karyotype",
            subtitle = paste0("Sample: ", sample, " \u00b7 ", label),
            footer = paste0("MaMaMIA ", version)
        ),
        path,
        useBytes = TRUE
    )
    invisible(path)
}

#' Merge payload palette with the application defaults
#' @noRd
bam_bundle_palette <- function(payload) {
    karyotype_scene_palette(payload)
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
    formatC(x, format = "f", digits = digits, decimal.mark = ".")
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
                         dash = NULL, rx = 2, title = NULL) {
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
    close <- if (is.null(title)) {
        "/>"
    } else {
        paste0(
            "><title>", bam_svg_escape(title), "</title></rect>"
        )
    }
    sprintf(
        paste0(
            '<rect x="%s" y="%s" width="%s" height="%s" rx="%s" ',
            'fill="%s" stroke="%s" stroke-width="%s"%s%s%s'
        ),
        bam_svg_number(x), bam_svg_number(y), bam_svg_number(width),
        bam_svg_number(height), bam_svg_number(rx, 1),
        bam_svg_escape(fill), bam_svg_escape(stroke),
        bam_svg_number(stroke_width, 2), dash_attr, opacity_attr, close
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

#' Static legend for the bundle image
#'
#' Consumes the scene legend so the image cannot drift from the browser.
#' @noRd
bam_svg_legend <- function(scene, x, y, width = NULL,
                           show_deletions = TRUE) {
    ## Accept the old `(palette, x, y, width, show_deletions)` call shape by
    ## resolving it through the same scene policy.
    if (is.null(scene$legend)) {
        palette <- bam_bundle_palette(list(palette = scene))
        scene <- list(
            palette = palette,
            legend = list(
                kind = "calls",
                entries = karyotype_scene_legend_entries(
                    palette, !isTRUE(show_deletions)
                )
            )
        )
    }
    palette <- scene$palette
    legend <- scene$legend
    if (identical(legend$kind, "heatmap")) {
        heatmap <- legend$heatmap
        total_width <- 190
        if (!is.null(width)) {
            x <- x + max(0, (as.numeric(width)[[1L]] - total_width) / 2)
        }
        gradient <- sprintf(
            paste0(
                '<defs><linearGradient id="mamamia-karyotype-score" ',
                'x1="0%%" x2="100%%" y1="0%%" y2="0%%">',
                '<stop offset="0%%" stop-color="%s"/>',
                '<stop offset="50%%" stop-color="%s"/>',
                '<stop offset="100%%" stop-color="%s"/>',
                "</linearGradient></defs>"
            ),
            bam_svg_escape(heatmap$negative),
            bam_svg_escape(heatmap$zero),
            bam_svg_escape(heatmap$positive)
        )
        return(paste0(
            gradient,
            bam_svg_text(
                "Score", x, y, size = 10, weight = 600,
                fill = palette$text
            ),
            bam_svg_rect(
                x + 36, y - 9, 126, 10, fill = "url(#mamamia-karyotype-score)",
                stroke = "none", stroke_width = 0, rx = 1.5
            ),
            bam_svg_text(
                heatmap$labels$min, x + 36, y + 13, size = 9,
                weight = 500, fill = palette$text
            ),
            bam_svg_text(
                heatmap$labels$zero, x + 99, y + 13, size = 9,
                weight = 500, fill = palette$text, anchor = "middle"
            ),
            bam_svg_text(
                heatmap$labels$max, x + 162, y + 13, size = 9,
                weight = 500, fill = palette$text, anchor = "end"
            )
        ))
    }
    entries <- legend$entries
    if (length(entries) == 0L) {
        return("")
    }
    item_widths <- vapply(entries, function(entry) {
        18 + nchar(entry$label, type = "width") * 5.6
    }, numeric(1))
    total_width <- sum(item_widths) + 16 * (length(entries) - 1L)
    if (!is.null(width)) {
        x <- x + max(0, (as.numeric(width)[[1L]] - total_width) / 2)
    }
    parts <- character()
    cursor <- x
    for (i in seq_along(entries)) {
        entry <- entries[[i]]
        parts <- c(
            parts,
            bam_svg_rect(
                cursor, y - 9, 12, 12,
                fill = entry$fill,
                stroke = entry$stroke,
                stroke_width = entry$stroke_width,
                dash = entry$dash, rx = 1.5
            ),
            bam_svg_text(
                entry$label, cursor + 18, y,
                size = 10, weight = 500, fill = palette$text
            )
        )
        cursor <- cursor + item_widths[[i]] + 16
    }
    paste(parts, collapse = "")
}

#' Opacity for a call segment
#'
#' Kept as a thin scene accessor for callers that still ask for one call's
#' opacity directly.
#' @noRd
bam_svg_call_opacity <- function(call) {
    karyotype_scene_call(call, 1L, karyotype_palette())$opacity
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
    margin_left <- 52
    margin_right <- 52
    top <- 118
    bottom <- 92
    plot_width <- max(1, width - margin_left - margin_right)
    plot_height <- max(1, height - top - bottom)
    ## New payloads carry the resolved scene, including heatmap bins. Use it
    ## directly so the static image cannot bin the score profile differently
    ## from the browser. Only a payload outside the current contract is
    ## resolved here with the same fixed bin policy.
    scene <- if (is.list(payload$scene)) {
        payload$scene
    } else {
        karyotype_scene(payload)
    }
    palette <- scene$palette
    chromosomes <- payload$chromosomes
    scene_chromosomes <- scene$chromosomes
    sample <- scene$sample
    if (is.null(subtitle) || !nzchar(subtitle)) {
        subtitle <- paste0(
            if (nzchar(sample)) paste0("Sample: ", sample) else "Sample",
            " \u00b7 ",
            if (identical(scene$mode, "heatmap")) "Heatmap" else "Calls"
        )
    }
    if (is.null(footer)) {
        footer <- "Generated by MaMaMIA"
    }
    if (!is.null(scene$omitted_note)) {
        subtitle <- paste0(
            subtitle, " \u00b7 ", scene$omitted_note$text
        )
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
        bam_svg_legend(scene, 0, 86, width = width)
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

    n <- length(scene_chromosomes)
    lengths <- vapply(scene_chromosomes, function(x) {
        value <- suppressWarnings(as.numeric(x$display_length))
        if (length(value) != 1L || !is.finite(value) || value <= 0) 0 else value
    }, numeric(1))
    max_length <- max(c(1, lengths))
    pair_gap <- 3
    ids <- vapply(scene_chromosomes, function(x) {
        as.character(x$id %||% "")
    }, character(1))
    added_flags <- vapply(
        scene_chromosomes, function(x) isTRUE(x$added), logical(1)
    )
    copies_vec <- vapply(scene_chromosomes, function(x) {
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
    for (i in seq_along(scene_chromosomes)) {
        scene_chromosome <- scene_chromosomes[[i]]
        added <- isTRUE(scene_chromosome$added)
        deleted <- isTRUE(scene_chromosome$deleted)
        selected <- isTRUE(scene_chromosome$selected)
        copies <- copies_vec[[i]]
        copy <- scene_chromosome$copy
        pair_width <- copies * bar_width + (copies - 1L) * pair_gap
        group_left <- margin_left + (group_index[[i]] - 1L) * slot +
            (slot - pair_width) / 2
        x <- group_left + (copy - 1L) * (bar_width + pair_gap)
        pair_center <- margin_left + (group_index[[i]] - 0.5) * slot
        bar_height <- max(3, plot_height * lengths[[i]] / max_length)
        bar_label <- scene_chromosome$label
        bar_state <- if (deleted) {
            " (deleted)"
        } else if (added) {
            " (added)"
        } else {
            ""
        }
        call_count <- length(scene_chromosome$calls)
        bar_title <- paste0(
            as.character(bar_label), bar_state,
            if (call_count > 0L) {
                paste0(" \u00b7 ", call_count, " calls")
            } else {
                ""
            }
        )
        base <- scene_chromosome$base
        group <- sprintf(
            '<g%s>',
            if (selected) "" else ' opacity="0.32"'
        )
        group <- paste0(
            group,
            bam_svg_rect(
                x, top, bar_width, bar_height,
                fill = base$fill, stroke = base$stroke,
                stroke_width = base$stroke_width,
                dash = base$dash,
                title = bar_title
            )
        )
        draw_calls <- scene_chromosome$calls
        if (identical(scene$mode, "heatmap") && !deleted) {
            for (bin in scene_chromosome$heatmap_bins) {
                bin_y <- top + bin$start * bar_height
                bin_height <- max(0.5, (bin$end - bin$start) * bar_height + 0.6)
                group <- paste0(
                    group,
                    bam_svg_rect(
                        x, bin_y, bar_width, bin_height,
                        fill = bin$fill, stroke = "none",
                        stroke_width = 0, rx = 0
                    )
                )
            }
            draw_calls <- Filter(
                function(call) isTRUE(call$cross_group),
                scene_chromosome$calls
            )
        }
        for (call in draw_calls) {
            start <- call$start
            end <- call$end
            if (end <= start) {
                next
            }
            call_y <- top + start * bar_height
            call_height <- max(1.5, (end - start) * bar_height)
            call_bp <- {
                from_bp <- call$start_bp
                to_bp <- call$end_bp
                if (length(from_bp) == 1L && length(to_bp) == 1L &&
                    is.finite(from_bp) && is.finite(to_bp)) {
                    paste0(
                        format(from_bp, big.mark = ",", scientific = FALSE,
                            trim = TRUE),
                        "\u2013",
                        format(to_bp, big.mark = ",", scientific = FALSE,
                            trim = TRUE),
                        " bp"
                    )
                } else {
                    paste0(
                        sprintf("%.1f%%", 100 * start), "\u2013",
                        sprintf("%.1f%%", 100 * end)
                    )
                }
            }
            call_state <- as.character(call$class_status)
            if (length(call_state) != 1L || is.na(call_state)) {
                call_state <- ""
            }
            call_score <- call$score
            call_title <- paste0(
                "Call: ", as.character(call$class),
                if (length(call_score) == 1L && is.finite(call_score)) {
                    paste0(" \u00b7 score ",
                        format(call_score, digits = 3, trim = TRUE))
                } else {
                    ""
                },
                " \u00b7 ", call_bp,
                if (nzchar(call_state)) {
                    paste0(" \u00b7 ", call_state)
                } else {
                    ""
                }
            )
            group <- paste0(
                group,
                bam_svg_rect(
                    x, call_y, bar_width, call_height,
                    fill = call$fill,
                    stroke = call$stroke,
                    stroke_width = call$stroke_width,
                    opacity = call$opacity,
                    dash = call$dash,
                    rx = 0,
                    title = call_title
                )
            )
        }
        if (!is.null(scene_chromosome$caption)) {
            group <- paste0(
                group,
                bam_svg_text(
                    scene_chromosome$caption$text, pair_center, top - 8,
                    size = 9, weight = 700,
                    fill = scene_chromosome$caption$fill,
                    anchor = "middle"
                )
            )
        }
        if (isTRUE(scene_chromosome$show_label)) {
            label <- bam_svg_truncate(bar_label, label_max)
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
    table <- .bedpe_table(
        calls, pass_only = pass_only, windows = windows
    )
    table <- bam_bundle_descriptive_calls(table)
    utils::write.table(
        table, file, sep = "\t", quote = FALSE,
        row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
    names(table)
}

#' README for a result bundle
#' @noRd
bam_bundle_readme <- function(sample, calls, generated, version,
                              bedpe_columns = character(),
                              confident_calls = NULL, result = NULL,
                              bundle = NULL, include_run = TRUE,
                              include_shared = TRUE,
                              include_manifest = TRUE) {
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
    run_lines <- bam_bundle_run_summary_lines(result, sample)
    has_file <- function(relative) {
        is.null(bundle) || file.exists(file.path(bundle, relative))
    }
    has_plots <- is.null(bundle) ||
        dir.exists(file.path(bundle, "qc", "plots"))
    karyotype_lines <- c(
        if (has_file("karyotype-full.svg")) {
            c(
                "karyotype-full.svg    Every call for the sample, without display",
                "                      filters."
            )
        },
        if (has_file("karyotype-current.svg")) {
            "karyotype-current.svg Calls for the sample after the active filters."
        },
        if (has_file("karyotype-confident.svg")) {
            "karyotype-confident.svg Confident calls only."
        },
        if (has_file("karyotype-actual.svg")) {
            c(
                "karyotype-actual.svg Confident calls with deleted chromosomes",
                "                      and deletion segments removed."
            )
        }
    )
    plot_lines <- if (has_plots) {
        c(
            "qc/plots/              PNG diagnostics for coverage, scores,",
            "                       residuals, masks and calls.",
            if (isTRUE(include_shared)) {
                c(
                    "                       Cohort plots, when present, are in",
                    "                       qc/plots/."
                )
            } else {
                character()
            }
        )
    } else {
        character()
    }
    c(
        "MaMaMIA results bundle",
        "",
        paste0("Sample: ", sample),
        paste0("Generated: ", timestamp),
        paste0("Version: MaMaMIA ", version),
        paste0(
            "Calls: ", n, " total, ", passing,
            " passing calibration, ", confident_n, " confident"
        ),
        "",
        if (length(run_lines) > 0L) {
            c("Run", "---", run_lines, "")
        } else {
            character()
        },
        "Contents",
        "--------",
        "Primary results",
        "calls-full.csv        Every call for this sample with every column.",
        "calls-passing.csv     Calls that pass calibration.",
        "calls-confident.csv   Passing calls with review_status = confident.",
        "calls-review.csv      Every call outside the confident set.",
        "calls-summary.csv     Calls, passing, confident, width and score by",
        "                      call class.",
        "calls-by-group.csv    Calls and masked fraction by chromosome group.",
        "calls-by-status.csv   Calls by review status.",
        "calls-by-copy-state.csv Calls by copy-number state.",
        "assignment-alternatives.csv Next-best recipient candidates, when",
        "                      present.",
        "calls-all.bedpe       Every call for this sample in BEDPE format.",
        "calls-passing.bedpe   Calibrated calls (pass = TRUE) in BEDPE format.",
        "calls-confident.bedpe Confident calls in BEDPE format.",
        karyotype_lines,
        "report.html           Self-contained summary with tables and the",
        "                      full karyotype.",
        "",
        "Diagnostics",
        "-----------",
        if (isTRUE(include_shared)) {
            "qc/cohort-summary.csv  One row per sample in the run."
        } else {
            character()
        },
        "qc/metrics-summary.csv Coverage, GC and MAPQ by chromosome group and",
        "                       role.",
        "qc/contig-summary.csv  Per-contig coverage and GC.",
        "qc/scores-summary.csv  Score and residual summaries by group.",
        "qc/metrics-windows.csv.gz Every window metric for this sample.",
        "qc/scores-windows.csv.gz Every window score for this sample.",
        if (isTRUE(include_shared)) {
            c(
                "qc/masks-summary.csv   Masked windows by chromosome group.",
                "qc/masks-windows.csv.gz Every mask decision.",
                "qc/calibration.csv     Calibration thresholds by mode and",
                "                       length bin."
            )
        } else {
            character()
        },
        "qc/component-evidence.csv Per-call evidence for every donor component.",
        plot_lines,
        "",
        "Genome-browser tracks",
        "---------------------",
        "tracks/calls-recipient.bed Calls on assigned recipient contigs.",
        "tracks/calls-donor.bed     Calls on donor contigs.",
        if (isTRUE(include_shared)) {
            "tracks/masks.bed           Masked windows in group coordinates."
        } else {
            character()
        },
        "tracks/breadth.bedGraph.gz Window breadth for this sample.",
        "",
        if (isTRUE(include_run)) {
            c(
                "Run and reproducibility",
                "-----------------------",
                "run/provenance.json    Full provenance record.",
                "run/config.yml         Resolved configuration.",
                "run/parameters.csv     Key analysis parameters.",
                "run/profile-parameters.csv Every effective profile parameter with",
                "                       its default and whether it was changed.",
                "run/versions.csv       Software and reference versions.",
                "run/session-info.txt   R session information.",
                "run/assembly.csv       Contigs, groups and roles.",
                "run/reference-indexes.csv Index paths, sizes and timestamps.",
                "run/samples.csv        Sample metadata.",
                "run/bams.csv           BAM paths, sizes, timestamps and target",
                "                       counts.",
                ""
            )
        } else {
            character()
        },
        if (isTRUE(include_manifest)) {
            c(
                "manifest.json          File sizes and MD5 checksums.",
                "checksums.md5          MD5 checksums including the manifest.",
                ""
            )
        } else {
            character()
        },
        "Confident set",
        "-------------",
        "The confident set contains calls that pass calibration and have",
        "review_status = confident. The Karyotype and Calls tabs show",
        "passing calls by default; turn off Passing only to see every call,",
        "or select a Status to narrow by review class. The full CSV and",
        "BEDPE keep every call.",
        "",
        "Diagnostic conventions",
        "----------------------",
        "qc/metrics-windows.csv.gz has one row per sample, window and role.",
        "qc/scores-windows.csv.gz has one row per sample and group window.",
        "Both files are compressed CSV with a header and descriptive column",
        "names. qc/masks-windows.csv.gz keeps every mask decision, including",
        "windows that were not masked.",
        "",
        "BED tracks are zero-based, half-open intervals. BED scores are zero",
        "because the call class and review status are carried in the name",
        "column. masks.bed uses synthetic contig names group_<group> because",
        "masks are computed in homology-group coordinates. breadth.bedGraph.gz",
        "has one value per window for this sample.",
        "",
        "manifest.json lists every payload file with its size and MD5",
        "checksum. checksums.md5 repeats the MD5 sums and also covers the",
        "manifest. report.html is a self-contained human-readable summary.",
        "qc/plots/, when present, holds the same diagnostic figures as PNG",
        "files; each plot is generated from the raw table named in its",
        "caption.",
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

#' Write the files for one sample into a bundle directory
#'
#' The caller adds the README, manifest, checksums and ZIP archive. Keeping
#' the sample payload separate lets the cohort bundle write the same files
#' under `<sample>/` without duplicating the export logic.
#'
#' @param result A `BamAnalysisResult`.
#' @param bundle Directory that receives the sample files.
#' @param sample Sample name.
#' @param current_payload,full_payload,confident_payload Optional karyotype
#'   payloads. `NULL` skips the corresponding SVG.
#' @param actual_payload Optional confident-call payload with deletions removed.
#' @param include_run Write the shared `run/` metadata.
#' @param include_masks,include_calibration Write run-level mask and
#'   calibration tables.
#' @param include_cohort Write the run-wide cohort summary.
#' @param include_cohort_plot Write `cohort-calls.png`.
#' @param include_cohort_section Include the cohort section in the report.
#' @param write_report Write `report.html`.
#' @param generated Timestamp used in the report.
#'
#' @return A list with the descriptive call tables, BEDPE columns,
#'   diagnostics and application version.
#' @noRd
bam_bundle_write_sample_files <- function(result, bundle, sample,
                                          current_payload = NULL,
                                          full_payload = NULL,
                                          confident_payload = NULL,
                                          actual_payload = NULL,
                                          cohort = NULL,
                                          metrics_rows = NULL,
                                          scores_rows = NULL,
                                          include_run = TRUE,
                                          include_masks = TRUE,
                                          include_calibration = TRUE,
                                          include_cohort = TRUE,
                                          include_cohort_plot = TRUE,
                                          include_cohort_section = TRUE,
                                          write_report = TRUE,
                                          generated = Sys.time()) {
    calls_raw <- result$calls
    calls_raw <- calls_raw[
        as.character(calls_raw$sample) == sample, , drop = FALSE
    ]
    confident_raw <- bam_confident_calls(calls_raw)
    calls <- bam_bundle_descriptive_calls(calls_raw)
    confident_calls <- bam_bundle_descriptive_calls(confident_raw)
    ## Annotate once; each BEDPE writer otherwise rescans the windows table
    ## for the same calls.
    bedpe_calls <- if (nrow(calls_raw) > 0L &&
        !"chr_id.don" %in% names(calls_raw)) {
        .annotate_call_intervals(calls_raw, result$windows)
    } else {
        calls_raw
    }
    bedpe_confident <- bam_confident_calls(bedpe_calls)

    utils::write.csv(
        calls, file.path(bundle, "calls-full.csv"),
        row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
    utils::write.csv(
        confident_calls, file.path(bundle, "calls-confident.csv"),
        row.names = FALSE, na = "", fileEncoding = "UTF-8"
    )
    bedpe_columns <- bam_bundle_write_bedpe(
        bedpe_calls, result$windows, file.path(bundle, "calls-all.bedpe"),
        pass_only = FALSE
    )
    bam_bundle_write_bedpe(
        bedpe_calls, result$windows, file.path(bundle, "calls-passing.bedpe"),
        pass_only = TRUE
    )
    bam_bundle_write_bedpe(
        bedpe_confident, result$windows,
        file.path(bundle, "calls-confident.bedpe"),
        pass_only = FALSE
    )
    diagnostics <- bam_bundle_write_diagnostics(
        result, bundle, sample, calls_raw, calls,
        cohort = cohort,
        metrics_rows = metrics_rows,
        scores_rows = scores_rows,
        include_run = include_run,
        include_masks = include_masks,
        include_calibration = include_calibration,
        include_cohort = include_cohort,
        include_cohort_plot = include_cohort_plot,
        include_masks_track = include_masks
    )
    version <- bam_bundle_version()
    bam_bundle_write_karyotype(
        current_payload, file.path(bundle, "karyotype-current.svg"),
        sample, "Filtered calls", version
    )
    bam_bundle_write_karyotype(
        full_payload, file.path(bundle, "karyotype-full.svg"),
        sample, "All calls", version
    )
    bam_bundle_write_karyotype(
        confident_payload, file.path(bundle, "karyotype-confident.svg"),
        sample, "Confident calls", version
    )
    bam_bundle_write_karyotype(
        actual_payload, file.path(bundle, "karyotype-actual.svg"),
        sample, "Confident calls, deletions removed", version
    )
    if (isTRUE(write_report)) {
        writeLines(
            bam_bundle_report_html(
                result, bundle, sample, calls_raw, diagnostics,
                generated, version,
                has_full_svg = !is.null(full_payload),
                show_cohort = include_cohort_section
            ),
            file.path(bundle, "report.html"),
            useBytes = TRUE
        )
    }
    list(
        calls = calls,
        confident_calls = confident_calls,
        calls_raw = calls_raw,
        bedpe_columns = bedpe_columns,
        diagnostics = diagnostics,
        version = version
    )
}

#' Locate the 7-Zip executable
#'
#' Bundles are archived with 7-Zip so compression itself uses the worker
#' budget. Both the `7z` name (Debian's `7zip` package) and the standalone
#' `7zz`/`7za` names are accepted.
#' @noRd
bam_bundle_zip_tool <- function() {
    paths <- Sys.which(c("7z", "7zz", "7za"))
    paths <- paths[!is.na(paths) & nzchar(paths)]
    if (length(paths) == 0L) {
        stop(
            "Could not create the results ZIP archive. ",
            "The `7z` command is required.",
            call. = FALSE
        )
    }
    unname(paths[[1L]])
}

#' 7-Zip arguments for a result bundle
#' @noRd
bam_bundle_archive_args <- function(zipfile, bundle_name,
                                    workers = bam_bundle_workers()) {
    workers <- suppressWarnings(as.integer(workers))
    if (length(workers) != 1L || is.na(workers) || workers < 1L) {
        workers <- 1L
    }
    c(
        "a",
        "-tzip",
        "-mx=9",
        paste0("-mmt=", workers),
        "-bso0",
        "-bsp0",
        shQuote(zipfile),
        shQuote(bundle_name)
    )
}

#' Build a bundle ZIP from a staging directory
#' @noRd
bam_bundle_archive <- function(staging, bundle_name, file) {
    old <- getwd()
    ## 7-Zip runs inside the staging directory, so a relative output path
    ## would land in staging and be deleted with it.
    file <- if (grepl("^(/|[A-Za-z]:[/\\\\])", file)) {
        file
    } else {
        file.path(old, file)
    }
    tool <- bam_bundle_zip_tool()
    ## Build the archive next to the destination first; only replace the
    ## previous bundle once 7-Zip has succeeded.
    tmp_zip <- tempfile(tmpdir = dirname(file), fileext = ".zip")
    on.exit({
        setwd(old)
        unlink(tmp_zip)
    }, add = TRUE)
    setwd(staging)
    status <- tryCatch(
        system2(
            tool,
            bam_bundle_archive_args(tmp_zip, bundle_name),
            stdout = FALSE,
            stderr = FALSE
        ),
        error = function(e) e
    )
    setwd(old)
    if (inherits(status, "error") ||
        !identical(as.integer(status), 0L)) {
        stop(
            "Could not create the results ZIP archive with 7-Zip.",
            call. = FALSE
        )
    }
    if (file.exists(file)) {
        unlink(file)
    }
    if (!file.rename(tmp_zip, file) &&
        !file.copy(tmp_zip, file, overwrite = TRUE)) {
        stop(
            "Could not move the results ZIP archive to '", file, "'",
            call. = FALSE
        )
    }
    invisible(file)
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
    bam_write_cohort_bundle(
        result, file,
        samples = sample,
        current_sample = sample,
        current_payload = current_payload,
        full_payload = full_payload,
        confident_payload = confident_payload,
        generated = generated
    )
}
