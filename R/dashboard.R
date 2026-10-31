## The introgression explorer: karyotype, calls, filters and BEDPE export.

dashboard_packages <- c(
    "bslib", "DT", "htmltools", "htmlwidgets", "shiny", "shinyWidgets"
)

require_dashboard_packages <- function() {
    missing <- dashboard_packages[
        !vapply(dashboard_packages, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(missing) == 0L) {
        return(invisible(TRUE))
    }
    stop(
        "Install dashboard dependencies first: ",
        sprintf('install.packages(c(%s))',
            paste(sprintf('"%s"', missing), collapse = ", ")),
        call. = FALSE
    )
}

#' Build a score profile plot
#'
#' One ggplot per result: a score line faceted by chromosome group with a
#' shaded rectangle over every visible call. Split out of the Shiny module so
#' it can be rendered and tested without a session.
#'
#' @param scores A `ScoreTable` with `group`, `start`, `score` and `mode`.
#' @param calls Optional `Calls` table with `start_bp` and `end_bp`.
#'
#' @return A `ggplot` object.
#' @noRd
explorer_score_plot <- function(scores, calls = NULL) {
    modes <- unique(as.character(scores$mode))
    mode_colours <- c(single = "#0f5c53", multi = "#b7791f")
    mode_colours <- mode_colours[intersect(names(mode_colours), modes)]
    p <- ggplot2::ggplot(
        scores,
        ggplot2::aes(
            x = .data$start / 1e6, y = .data$score,
            colour = .data$mode
        )
    ) +
        ggplot2::geom_hline(
            yintercept = 0, colour = "#d8d4c8", linewidth = 0.3
        ) +
        ggplot2::geom_line(linewidth = 0.35) +
        ggplot2::facet_wrap(
            ~group,
            scales = "free_x", ncol = 1,
            labeller = ggplot2::as_labeller(
                function(x) paste("Group", x)
            )
        ) +
        ggplot2::labs(
            x = "Position (Mb)", y = "Introgression score",
            colour = NULL
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(
                fill = "#ffffff", colour = NA
            ),
            panel.background = ggplot2::element_rect(
                fill = "#ffffff", colour = NA
            ),
            panel.grid.major = ggplot2::element_line(
                colour = "#ece9df", linewidth = 0.3
            ),
            panel.grid.minor = ggplot2::element_blank(),
            panel.spacing = ggplot2::unit(8, "pt"),
            axis.text = ggplot2::element_text(
                colour = "#6d7168", size = 8.5
            ),
            axis.title = ggplot2::element_text(
                colour = "#4a4d45", size = 9.5
            ),
            strip.text = ggplot2::element_text(
                colour = "#1b1e1a", face = "bold", size = 8.5,
                hjust = 0, margin = ggplot2::margin(b = 3)
            ),
            strip.background = ggplot2::element_blank(),
            legend.position = "top",
            legend.justification = "right",
            legend.key.width = ggplot2::unit(16, "pt"),
            legend.text = ggplot2::element_text(
                colour = "#6d7168", size = 8
            ),
            plot.margin = ggplot2::margin(8, 12, 6, 10)
        )
    if (length(mode_colours) > 0L) {
        p <- p + ggplot2::scale_colour_manual(values = mode_colours)
    }
    if (!is.null(calls) && nrow(calls) > 0L) {
        rect <- calls[
            as.character(calls$group) %in% unique(as.character(scores$group)),
            ,
            drop = FALSE
        ]
        if (nrow(rect) == 0L) {
            return(p)
        }
        rect$group <- factor(rect$group, levels = unique(scores$group))
        rect$ymin <- 0
        rect$ymax <- max(scores$score, na.rm = TRUE)
        p <- p + ggplot2::geom_rect(
            data = rect,
            ggplot2::aes(
                xmin = .data$start_bp / 1e6,
                xmax = .data$end_bp / 1e6,
                ymin = .data$ymin, ymax = .data$ymax
            ),
            inherit.aes = FALSE, alpha = 0.1, colour = NA,
            fill = "#b23a3a"
        )
    }
    p
}

.dashboard_calls <- function(result, sample = NULL) {
    if (!inherits(result, "BamAnalysisResult")) {
        stop("`result` must be a BamAnalysisResult", call. = FALSE)
    }
    calls <- result$calls
    if (is.null(calls) || nrow(calls) == 0L) {
        return(calls)
    }
    if (!is.null(sample)) {
        calls <- calls[calls$sample %in% sample, , drop = FALSE]
    }
    calls
}

.dashboard_scores <- function(result, sample = NULL) {
    scores <- result$scores
    if (is.null(scores) || nrow(scores) == 0L) {
        return(scores)
    }
    if (!is.null(sample)) {
        scores <- scores[scores$sample %in% sample, , drop = FALSE]
    }
    scores
}

explorer_view_filters_ui <- function(id) {
    ns <- shiny::NS(id)
    picker <- function(inputId, label, title, choices,
                       selected = character(0)) {
        shinyWidgets::pickerInput(
            ns(inputId),
            label = shiny::span(class = "visually-hidden", label),
            choices = choices,
            selected = selected,
            multiple = TRUE,
            width = "fit",
            options = shinyWidgets::pickerOptions(
                actionsBox = FALSE,
                title = title,
                noneSelectedText = title,
                selectedTextFormat = "count > 1",
                countSelectedText = "{0} selected",
                liveSearch = FALSE
            )
        )
    }
    shiny::div(
        class = "mamamia-filter-pills",
        picker(
            "class", "Call class", "Class",
            c(
                "Substitution" = "substitution",
                "Addition" = "addition",
                "Deletion" = "deletion",
                "Complex" = "complex",
                "Other" = "other"
            )
        ),
        picker(
            "class_status", "Class status", "Status",
            c(
                "Confident" = "confident",
                "Ambiguous" = "ambiguous",
                "Needs review" = "needs_review"
            ),
            selected = "confident"
        ),
        picker(
            "copy_state", "Copy state", "Copy",
            c(
                "Balanced" = "balanced",
                "Gain" = "gain",
                "Loss" = "loss",
                "Compound" = "compound",
                "Neutral" = "neutral"
            )
        ),
        picker(
            "confidence", "Confidence", "Confidence",
            c("High" = "high", "Moderate" = "moderate", "Low" = "low")
        ),
        shiny::div(
            class = "mamamia-filter-switch",
            bslib::input_switch(
                ns("pass_only"), "Passing only", value = TRUE
            )
        )
    )
}

explorer_call_filters_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::numericInput(
        ns("min_width"), "Minimum call width (Mb)", value = 5, min = 0, step = 1
    )
}

explorer_selection_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::div(
        class = "d-flex gap-2",
        shiny::actionButton(
            ns("select_all"), "Select visible", class = "btn-sm flex-fill"
        ),
        shiny::actionButton(
            ns("clear_selection"), "Clear", class = "btn-sm flex-fill"
        )
    )
}

explorer_export_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shiny::div(
            class = "mamamia-export-section",
            shiny::h5("Results bundle"),
            shiny::p(
                class = "mamamia-export-lead",
                paste(
                    "Download full and confident call tables, all-call,",
                    "passing and confident BEDPE files, and matching",
                    "karyotype SVGs in one ZIP archive."
                )
            ),
            shiny::downloadButton(
                ns("download_bundle"),
                "Download results bundle (.zip)",
                class = "btn-sm btn-primary"
            )
        ),
        shiny::div(
            class = "mamamia-export-section",
            shiny::h5("Selected calls"),
            shiny::textOutput(ns("selection_status")),
            shiny::downloadButton(
                ns("download_bedpe"),
                "Download selected BEDPE",
                class = "btn-sm"
            )
        )
    )
}

explorer_tooltip_help_ui <- function() {
    row <- function(term, description) {
        if (inherits(description, "shiny.tag.list")) {
            parts <- vapply(
                description,
                function(part) as.character(part),
                character(1)
            )
            rendered <- parts[[1L]]
            if (length(parts) > 1L) {
                for (i in seq.int(2L, length(parts))) {
                    previous <- gsub("<[^>]*>", "", parts[[i - 1L]])
                    following <- gsub("<[^>]*>", "", parts[[i]])
                    previous_end <- substr(
                        previous[[1L]],
                        nchar(previous[[1L]]),
                        nchar(previous[[1L]])
                    )
                    following_start <- substr(following[[1L]], 1L, 1L)
                    add_space <- nzchar(previous_end) &&
                        nzchar(following_start) &&
                        !grepl("[[:space:]]", previous_end) &&
                        !grepl("[[:space:]]", following_start) &&
                        !following_start %in% c(
                            ".", ",", ";", ":", "!", "?", "%", ")", "]", "}"
                        )
                    rendered <- paste0(
                        rendered,
                        if (add_space) " " else "",
                        parts[[i]]
                    )
                }
            }
            description <- shiny::HTML(rendered)
        }
        shiny::tagList(
            shiny::tags$dt(term),
            shiny::tags$dd(description)
        )
    }
    section <- function(title, intro = NULL, ...) {
        shiny::div(
            class = "mamamia-help-section",
            shiny::tags$h6(title),
            if (!is.null(intro)) {
                shiny::p(class = "mamamia-help-intro", intro)
            },
            shiny::tags$dl(class = "mamamia-help-list", ...)
        )
    }
    shiny::div(
        class = "mamamia-tooltip-help",
        shiny::p(
            class = "mamamia-help-lead",
            paste(
                "Hover a chromosome bar to inspect the nearest introgression",
                "call. The card is a call-level summary, not a per-window",
                "trace. It follows the current class, class-status, copy-state,",
                "confidence and width filters, so a filtered-out call cannot",
                "appear."
            )
        ),
        shiny::div(
            class = "mamamia-help-callout",
            paste(
                "Read top to bottom: the header identifies the call and its",
                "class confidence; Placement says where it was matched;",
                "Evidence says why; Dosage and Geometry describe what changed;",
                "Position (heatmap only) adds local context."
            )
        ),
        section(
            "Quick interpretation",
            paste(
                "Use these patterns as a starting point, then check placement",
                "and geometry before treating a call as a confirmed event."
            ),
            row(
                "Balanced substitution",
                shiny::tagList(
                    "Donor excess, matching recipient loss, and a ",
                    shiny::tags$strong("near-neutral total"),
                    ". Usually class ",
                    shiny::tags$strong("substitution"),
                    " with copy state ",
                    shiny::tags$strong("balanced"), "."
                )
            ),
            row(
                "Addition (gain)",
                shiny::tagList(
                    "Donor excess without reciprocal recipient loss, so total",
                    "depth rises. Usually copy state ",
                    shiny::tags$strong("gain"), "."
                )
            ),
            row(
                "Deletion (loss)",
                shiny::tagList(
                    "Recipient loss without donor support, so total depth",
                    "falls. Usually copy state ",
                    shiny::tags$strong("loss"), "."
                )
            ),
            row(
                "Complex",
                shiny::tagList(
                    "Donor and recipient evidence disagree, or the donor and",
                    "recipient spans are geometrically incompatible. Check",
                    "Evidence and Geometry together."
                )
            ),
            row(
                "Other",
                shiny::tagList(
                    "A scored segment with enough signal to report but not",
                    "enough support for a more specific class. Treat it as a",
                    "review lead, not a confirmed event."
                )
            ),
            row(
                "No host",
                shiny::tagList(
                    "The donor segment is assigned by chromosome group but",
                    "there is no recipient deletion to anchor it. The class",
                    "can still be confident; the host is not."
                )
            )
        ),
        section(
            "Header",
            paste(
                "The header always answers three questions: what is the",
                "call, how confident is the class, and which group and donor",
                "span does it describe?"
            ),
            row(
                "Class",
                shiny::tagList(
                    shiny::tags$strong("substitution"),
                    ", ",
                    shiny::tags$strong("addition"),
                    ", ",
                    shiny::tags$strong("deletion"),
                    ", ",
                    shiny::tags$strong("complex"),
                    " or ",
                    shiny::tags$strong("other"),
                    ". The left border and class title use the same colour as",
                    "the legend."
                )
            ),
            row(
                "Status",
                shiny::tagList(
                    "Class confidence: ",
                    shiny::tags$strong("confident"),
                    ", ",
                    shiny::tags$strong("ambiguous"),
                    " or ",
                    shiny::tags$strong("needs_review"),
                    ". This describes the class, not the recipient placement."
                )
            ),
            row(
                "Context",
                paste(
                    "Sample name, chromosome group and donor span. Compare",
                    "the donor span with the recipient span under Placement."
                )
            ),
            row(
                "Accent colour",
                paste(
                    "The card border and class title inherit the call colour",
                    "from the legend; the status pill uses the class-status",
                    "palette."
                )
            )
        ),
        section(
            "Placement",
            paste(
                "Placement answers where the donor segment was matched.",
                "Placement confidence and class confidence are separate: a",
                "confident substitution can still have an ambiguous host."
            ),
            row(
                "Recipient",
                paste(
                    "Assigned recipient chromosome and interval. The interval",
                    "is the deletion run used for placement, not necessarily",
                    "the full donor span."
                )
            ),
            row(
                "Donor",
                "Donor chromosome, group and interval."
            ),
            row(
                "Relation",
                shiny::tagList(
                    shiny::tags$strong("same group"),
                    ", ",
                    shiny::tags$strong("cross-group candidate"),
                    " or ",
                    shiny::tags$strong("not placed"),
                    ". A cross-group call is a translocation candidate, not a",
                    "confirmed host."
                )
            ),
            row(
                "Assignment",
                shiny::tagList(
                    shiny::tags$strong("assigned"),
                    ", ",
                    shiny::tags$strong("cross-group candidate"),
                    ", ",
                    shiny::tags$strong("ambiguous"),
                    ", ",
                    shiny::tags$strong("no host"),
                    " or ",
                    shiny::tags$strong("group only"),
                    ". This is the placement state used in the Calls table."
                )
            ),
            row(
                "Margin",
                paste(
                    "Lead of the chosen host over the next candidate. A small",
                    "margin is a reason to inspect the alternative placement."
                )
            ),
            row(
                "Group class",
                paste(
                    "Shown only when recipient placement changes the call",
                    "class, for example addition to substitution."
                )
            ),
            row(
                "Alternative",
                paste(
                    "Another possible host shown when the assignment is",
                    "ambiguous. Compare its group, interval and margin before",
                    "choosing a host."
                )
            ),
            row(
                "Reading placement",
                paste(
                    "Strongest: same group and assigned. Cross-group is a",
                    "candidate only. No host means there is no recipient loss",
                    "to anchor the call. Group only is a fallback when no",
                    "contig-level host is available."
                )
            )
        ),
        section(
            "Evidence",
            paste(
                "Evidence compares donor and recipient breadth in the same",
                "windows. This is the core reason for the call; the numbers",
                "are signed residuals, not raw read counts."
            ),
            row(
                "Score",
                paste(
                    "Paired donor/recipient statistic for the segment. Larger",
                    "positive values mean stronger donor excess. Thresholds are",
                    "calibrated for the window count and segment width."
                )
            ),
            row(
                "Donor",
                "Mean donor-side residual. Positive means donor gain; negative means donor loss."
            ),
            row(
                "Recipient",
                "Mean recipient-side residual. Positive means recipient gain; negative means recipient loss."
            ),
            row(
                "Total",
                paste(
                    "Relative total-depth residual. Near zero means balanced;",
                    "positive means total gain; negative means total loss. This",
                    "is what separates a replacement from an addition or",
                    "deletion."
                )
            ),
            row(
                "Colour",
                paste(
                    "Metrics use two decimals: green is positive, red is",
                    "negative and muted is missing. The sign matters more than",
                    "the absolute size for interpretation."
                )
            ),
            row(
                "Copy state",
                shiny::tagList(
                    shiny::tags$strong("balanced"),
                    ", ",
                    shiny::tags$strong("gain"),
                    ", ",
                    shiny::tags$strong("loss"),
                    ", ",
                    shiny::tags$strong("compound"),
                    " or ",
                    shiny::tags$strong("neutral"),
                    ". This is a hypothesis about the copy configuration, not",
                    "an absolute copy number."
                )
            ),
            row(
                "Evidence source",
                shiny::tagList(
                    shiny::tags$strong("donor"),
                    ", ",
                    shiny::tags$strong("recipient"),
                    ", ",
                    shiny::tags$strong("donor+recipient"),
                    " or ",
                    shiny::tags$strong("none"),
                    ". It says which side drove the class."
                )
            ),
            row(
                "Typical pattern",
                paste(
                    "Substitution: donor positive, recipient negative, total",
                    "near zero. Addition: donor positive, recipient near zero,",
                    "total positive. Deletion: donor near zero, recipient",
                    "negative, total negative. Complex: conflicting signs or a",
                    "width mismatch."
                )
            )
        ),
        section(
            "Status chips",
            paste(
                "Chips appear only when their condition holds. They are",
                "shortcuts into the call table and the full result columns."
            ),
            row(
                "class high / moderate / low",
                "Class confidence reported by the classifier."
            ),
            row(
                "donor / recipient / donor+recipient / none",
                "Which evidence sides supported the class."
            ),
            row(
                "balanced / gain / loss / compound / neutral",
                "The copy-state hypothesis for the segment."
            ),
            row(
                "total gain / total loss",
                "Total-depth evidence is not neutral."
            ),
            row(
                "components heterogeneous",
                paste(
                    "Component spread or donor-component support is uneven;",
                    "the group-level score may hide a component-specific event."
                )
            ),
            row(
                "low mappability",
                paste(
                    "BAM-derived high-MAPQ fraction is below 0.5; placement",
                    "and dosage are less reliable."
                )
            ),
            row(
                "not passing",
                "The call is shown for review but did not pass the configured filters."
            )
        ),
        section(
            "Dosage",
            paste(
                "Dosage is relative to the registry copy number. It is not an",
                "absolute copy-number measurement and should be read as an",
                "estimate alongside Evidence and Geometry."
            ),
            row(
                "Copies",
                paste(
                    "Estimated affected copies out of the configured copy",
                    "number, for example 2 of 2 copies. Shows unavailable when",
                    "no copy count is estimated."
                )
            ),
            row(
                "Estimated fraction",
                "Affected fraction of recipient copies, shown as a percentage."
            ),
            row(
                "Estimated donor copies",
                "Relative donor-copy estimate when the signal supports one."
            ),
            row(
                "Confidence",
                shiny::tagList(
                    "Dosage confidence: ",
                    shiny::tags$strong("high"),
                    ", ",
                    shiny::tags$strong("moderate"),
                    ", ",
                    shiny::tags$strong("low"),
                    " or ",
                    shiny::tags$strong("unavailable"),
                    "."
                )
            ),
            row(
                "Caveat",
                paste(
                    "Only Copies is rounded to an integer. The fraction and",
                    "donor-copy values are continuous estimates. An addition",
                    "derives Copies from donor copies; a deletion without donor",
                    "support derives it from the recipient loss fraction at",
                    "lower confidence."
                )
            )
        ),
        section(
            "Geometry",
            paste(
                "Geometry checks whether the donor and recipient spans are",
                "compatible. It can turn a substitution into a complex call",
                "even when the depth evidence looks clean."
            ),
            row("Donor width", "Width of the donor segment."),
            row(
                "Recipient run",
                "Width of the recipient deletion run used for placement."
            ),
            row(
                "Overlap",
                paste(
                    "Reciprocal overlap of donor and recipient intervals in",
                    "contig-fraction space. Low overlap weakens the pairing."
                )
            ),
            row(
                "Width ratio",
                paste(
                    "Donor fraction divided by recipient-run fraction. A value",
                    "near 1 supports a balanced replacement; values far from 1",
                    "suggest an addition, deletion or complex event."
                )
            ),
            row(
                "Interpretation",
                paste(
                    "Overlap below 0.5 or a width ratio far from 1 is a reason",
                    "to prefer complex or review the recipient assignment."
                )
            )
        ),
        section(
            "Position and footer",
            paste(
                "Position is available in Heatmap mode; Calls mode ends with",
                "the footer instead."
            ),
            row(
                "Position",
                "Hovered position along the chromosome, in bp, kb or Mb."
            ),
            row("Score", "Heatmap score at that relative position."),
            row(
                "Mode",
                "Heatmap shows the score profile; Calls shows coloured call segments."
            ),
            row(
                "Click",
                paste(
                    "Click a chromosome to filter its group; click the",
                    "background to clear the filter."
                )
            )
        ),
        section(
            "No-call card",
            paste(
                "When the pointer is not over a call, the tooltip describes",
                "the chromosome instead of a call."
            ),
            row(
                "Header",
                "Chromosome state: recipient chromosome, deleted recipient or added donor."
            ),
            row(
                "Chromosome block",
                paste(
                    "Component, state, homolog number (when more than one",
                    "copy is configured) and configured length."
                )
            ),
            row(
                "Call block",
                paste(
                    "Calls mode only: says that no assigned donor segment is",
                    "present at the hovered position."
                )
            ),
            row(
                "Position block",
                "Heatmap mode only: the hovered position and score."
            )
        ),
        section(
            "Caveats",
            NULL,
            row(
                "Relative, not absolute",
                paste(
                    "Copy state and dosage describe relative read-depth",
                    "changes. They are not absolute copy numbers."
                )
            ),
            row(
                "Class vs placement",
                paste(
                    "Class status describes confidence in the label; the",
                    "Assignment row describes confidence in the host."
                )
            ),
            row(
                "Cross-group",
                paste(
                    "A cross-group result is a translocation candidate, not a",
                    "confirmed donor-recipient pair."
                )
            ),
            row(
                "Filters",
                paste(
                    "The card follows the visible calls. Clear a filter to",
                    "inspect calls that are currently hidden."
                )
            ),
            row(
                "Heatmap score",
                paste(
                    "A heatmap value is local context, not call-level evidence.",
                    "Use the Score, Donor, Recipient and Total metrics for",
                    "the call itself."
                )
            )
        )
    )
}

#' @noRd
mamamia_analysis_loader_ui <- function() {
    shiny::div(
        class = "mamamia-analysis-loading",
        style = "display: none;",
        role = "status",
        `aria-live` = "polite",
        `aria-hidden` = "true",
        shiny::div(
            class = "mamamia-analysis-loader",
            `aria-hidden` = "true",
            shiny::HTML(paste0(
                '<svg viewBox="0 0 80 80" ',
                'xmlns="http://www.w3.org/2000/svg">',
                '<circle class="mamamia-analysis-orbit" ',
                'cx="40" cy="40" r="31"></circle>',
                '<path class="mamamia-analysis-strand" ',
                'd="M28 16 C 28 28, 34 34, 40 40 ',
                'C 46 46, 52 52, 52 64"></path>',
                '<path class="mamamia-analysis-strand ',
                'mamamia-analysis-strand-warm" ',
                'd="M52 16 C 52 28, 46 34, 40 40 ',
                'C 34 46, 28 52, 28 64"></path>',
                '<circle class="mamamia-analysis-chiasma" ',
                'cx="40" cy="40" r="3.2"></circle>',
                "</svg>"
            ))
        ),
        shiny::div(class = "mamamia-analysis-eyebrow", "Analysis"),
        shiny::div(
            class = "mamamia-analysis-stage", "Preparing analysis"
        ),
        shiny::div(
            class = "mamamia-analysis-track",
            shiny::tags$span(class = "mamamia-analysis-fill")
        )
    )
}

explorer_karyotype_ui <- function(id, height = "calc(100vh - 360px)") {
    ns <- shiny::NS(id)
    shiny::div(
        class = "mamamia-karyotype",
        shiny::div(
            style = paste(
                "display:flex; flex-wrap:wrap; align-items:center;",
                "gap:.75rem; padding:.2rem .3rem .1rem;"
            ),
            shiny::radioButtons(
                ns("karyotype_mode"), "Plot",
                choices = c("Heatmap" = "heatmap", "Calls" = "calls"),
                selected = "heatmap", inline = TRUE
            ),
            shiny::div(
                class = paste(
                    "ms-auto d-flex flex-wrap align-items-center gap-2"
                ),
                shiny::tags$button(
                    id = ns("download_karyotype_current"),
                    type = "button",
                    class = "btn btn-sm btn-outline-secondary",
                    disabled = "disabled",
                    title = paste(
                        "Download the karyotype exactly as filtered",
                        "and displayed (SVG)."
                    ),
                    shiny::icon("download"),
                    "Current view"
                ),
                shiny::tags$button(
                    id = ns("download_karyotype_full"),
                    type = "button",
                    class = "btn btn-sm btn-outline-secondary",
                    disabled = "disabled",
                    title = paste(
                        "Download every call for this sample,",
                        "ignoring the active filters (SVG)."
                    ),
                    shiny::icon("download"),
                    "Full view"
                ),
                shiny::actionButton(
                    ns("tooltip_help"), "Tooltip help",
                    icon = shiny::icon("circle-question"),
                    class = "btn-sm btn-outline-secondary"
                )
            )
        ),
        shiny::div(
            class = "mamamia-karyotype-plot",
            mamamia_karyotype_output(ns("karyotype"), height = height),
            mamamia_analysis_loader_ui()
        ),
        shiny::tags$script(shiny::HTML(
            "$(document).off('click.mamamiaKaryotypeCurrent').on(",
            "'click.mamamiaKaryotypeCurrent', ",
            "'[id$=\"-download_karyotype_current\"]', function() {",
            "  var id = this.id.replace(",
            "    /-download_karyotype_current$/, '-karyotype');",
            "  var el = document.getElementById(id);",
            "  if (el && el.mamamiaDownloadKaryotype) {",
            "    el.mamamiaDownloadKaryotype('current');",
            "  }",
            "});",
            "$(document).off('click.mamamiaKaryotypeFull').on(",
            "'click.mamamiaKaryotypeFull', ",
            "'[id$=\"-download_karyotype_full\"]', function() {",
            "  var id = this.id.replace(",
            "    /-download_karyotype_full$/, '-karyotype');",
            "  var el = document.getElementById(id);",
            "  if (el && el.mamamiaDownloadKaryotype) {",
            "    el.mamamiaDownloadKaryotype('full');",
            "  }",
            "});"
        ))
    )
}

#' Build the karyotype payload with an unfiltered export copy
#'
#' The dashboard renders the current sample, mode and filters. The full
#' export ignores the call filters and group dimming, so the widget receives
#' the per-chromosome calls for every call in the same sample.
#'
#' @param result A `BamAnalysisResult`.
#' @param sample Sample to display; defaults to the first.
#' @param mode `"heatmap"` or `"calls"`.
#' @param calls Calls to overlay; defaults to every call in the result.
#' @param selected_groups Groups currently selected in the explorer.
#'
#' @return A karyotype payload with `full_chromosomes` and `full_calls`
#'   entries for the unfiltered export.
#' @noRd
karyotype_export_payload <- function(result, sample = NULL,
                                     mode = c("heatmap", "calls"),
                                     calls = NULL, selected_groups = NULL) {
    mode <- match.arg(mode)
    payload <- karyotype_payload(
        result,
        sample = sample,
        mode = mode,
        calls = calls,
        selected_groups = selected_groups
    )
    full <- karyotype_payload(
        result,
        sample = sample,
        mode = mode,
        calls = NULL,
        selected_groups = NULL
    )
    payload$full_chromosomes <- lapply(
        full$chromosomes,
        function(chromosome) {
            chromosome$scores <- numeric()
            chromosome
        }
    )
    payload$full_calls <- lapply(
        full$chromosomes,
        function(chromosome) chromosome$calls
    )
    payload
}

explorer_controls_ui <- function(id) {
    shiny::tagList(
        explorer_view_filters_ui(id),
        explorer_call_filters_ui(id),
        explorer_selection_ui(id),
        explorer_export_ui(id)
    )
}

explorer_results_ui <- function(id) {
    ns <- shiny::NS(id)
    bslib::navset_card_tab(
        id = ns("tabs"),
        bslib::nav_panel(
            "Karyotype",
            explorer_karyotype_ui(id, height = "calc(100vh - 360px)")
        ),
        bslib::nav_panel(
            "Calls",
            DT::DTOutput(ns("calls"))
        )
    )
}

#' Server for the explorer
#'
#' @param id Module id.
#' @param result Reactive returning a `BamAnalysisResult`, or the result itself.
#' @param sample Optional sample name or reactive.
#' @noRd
explorer_server <- function(id, result, sample = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        selected_group <- shiny::reactiveVal(NULL)
        result_value <- function() {
            if (shiny::is.reactive(result)) result() else result
        }
        sample_value <- function() {
            if (is.null(sample)) {
                return(NULL)
            }
            value <- if (shiny::is.reactive(sample)) sample() else sample
            if (is.null(value) || length(value) != 1L || is.na(value) ||
                !nzchar(value)) {
                return(NULL)
            }
            as.character(value)
        }

        scores <- shiny::reactive({
            res <- result_value()
            if (is.null(res)) {
                return(NULL)
            }
            .dashboard_scores(res, sample_value())
        })

        calls <- shiny::reactive({
            res <- result_value()
            if (is.null(res)) {
                return(NULL)
            }
            .dashboard_calls(res, sample_value())
        })

        group_filter <- selected_group

        shiny::observeEvent(scores(), {
            selected_group(NULL)
        }, ignoreNULL = TRUE)

        shiny::observeEvent(input$karyotype_click, {
            clicked <- input$karyotype_click
            group <- clicked$group
            x <- scores()
            if (is.null(x) || nrow(x) == 0L) {
                return()
            }
            groups <- sort(unique(as.character(x$group)))
            if (is.null(group) || length(group) != 1L || is.na(group) ||
                !nzchar(group)) {
                selected_group(NULL)
            } else if (group %in% groups) {
                selected_group(group)
            }
        }, ignoreInit = TRUE)

        call_display_group <- function(x) {
            groups <- as.character(x$group)
            if (all(c("recipient_group", "recipient_cross_group") %in%
                names(x))) {
                cross <- as.logical(x$recipient_cross_group) %in% TRUE
                recipient <- as.character(x$recipient_group)
                move <- cross & !is.na(recipient) & nzchar(recipient)
                groups[move] <- recipient[move]
            }
            groups
        }

        visible_calls <- shiny::reactive({
            x <- calls()
            if (is.null(x) || nrow(x) == 0L) {
                return(x)
            }
            keep <- rep(TRUE, nrow(x))
            if (isTRUE(input$pass_only) && "pass" %in% names(x)) {
                keep <- keep & isTRUE_vec(x$pass)
            }
            groups <- group_filter()
            if (length(groups) > 0L) {
                keep <- keep & (
                    as.character(x$group) %in% groups |
                        call_display_group(x) %in% groups
                )
            }
            if (length(input$class) > 0L) {
                keep <- keep & x$class %in% input$class
            }
            if (length(input$class_status) > 0L) {
                keep <- keep & x$class_status %in% input$class_status
            }
            if (length(input$copy_state) > 0L &&
                "copy_state" %in% names(x)) {
                keep <- keep & x$copy_state %in% input$copy_state
            }
            if (length(input$confidence) > 0L) {
                keep <- keep & x$confidence %in% input$confidence
            }
            width <- suppressWarnings(as.numeric(input$min_width))
            if (length(width) == 1L && is.finite(width) && width > 0) {
                keep <- keep & x$width_bp >= width * 1e6
            }
            x[keep, , drop = FALSE]
        })

        selected <- shiny::reactive({
            x <- visible_calls()
            rows <- input$calls_rows_selected
            if (is.null(x) || length(rows) == 0L) {
                return(integer())
            }
            as.integer(rows)
        })

        output$selection_status <- shiny::renderText({
            if (is.null(result_value())) {
                return("No analysis loaded yet")
            }
            paste(length(selected()), "call(s) selected for BEDPE export")
        })

        proxy <- DT::dataTableProxy("calls", session = session)
        shiny::observeEvent(input$select_all, {
            x <- visible_calls()
            if (!is.null(x) && nrow(x) > 0L) {
                DT::selectRows(proxy, seq_len(nrow(x)))
            }
        })
        shiny::observeEvent(input$clear_selection, {
            DT::selectRows(proxy, NULL)
        })

        shiny::observeEvent(input$tooltip_help, {
            dialog <- shiny::modalDialog(
                title = "Tooltip guide",
                explorer_tooltip_help_ui(),
                size = "l",
                easyClose = TRUE,
                footer = shiny::modalButton("Close")
            )
            dialog <- htmltools::tagAppendAttributes(
                dialog,
                class = "mamamia-tooltip-modal modal-dialog-centered"
            )
            shiny::showModal(dialog)
        })

        output$karyotype <- render_mamamia_karyotype({
            mode <- input$karyotype_mode
            if (length(mode) != 1L || is.na(mode) ||
                !mode %in% c("heatmap", "calls")) {
                mode <- "heatmap"
            }
            res <- result_value()
            if (is.null(res)) {
                return(mamamia_karyotype(list(
                    mode = mode,
                    sample = "",
                    score_limit = 0.5,
                    palette = karyotype_palette(),
                    groups = list(),
                    message = "Run the analysis to plot",
                    hint = paste(
                        "Load a BAM or the bundled example data, then run",
                        "the analysis to map introgressions."
                    )
                )))
            }
            mamamia_karyotype(karyotype_export_payload(
                res,
                sample = sample_value(),
                mode = mode,
                calls = visible_calls(),
                selected_groups = group_filter()
            ))
        })

        output$calls <- DT::renderDT({
            shiny::validate(
                shiny::need(!is.null(result_value()), "Run the analysis to list calls")
            )
            x <- visible_calls()
            columns <- intersect(bam_brief_call_columns(x), names(x))
            display <- x[, columns, drop = FALSE]
            humanise <- function(values, labels) {
                values <- as.character(values)
                out <- unname(labels[values])
                out[is.na(out)] <- values[is.na(out)]
                out
            }
            if ("class" %in% names(display)) {
                display$class <- humanise(display$class, c(
                    substitution = "Substitution",
                    addition = "Addition",
                    deletion = "Deletion",
                    complex = "Complex",
                    other = "Other"
                ))
            }
            if ("copy_state" %in% names(display)) {
                display$copy_state <- humanise(display$copy_state, c(
                    balanced = "Balanced",
                    gain = "Gain",
                    loss = "Loss",
                    compound = "Compound",
                    neutral = "Neutral"
                ))
            }
            if ("class_status" %in% names(display)) {
                display$class_status <- humanise(display$class_status, c(
                    confident = "Confident",
                    ambiguous = "Ambiguous",
                    needs_review = "Needs review"
                ))
            }
            if ("recipient_assignment" %in% names(display)) {
                display$recipient_assignment <- humanise(
                    display$recipient_assignment,
                    c(
                        assigned = "Assigned",
                        ambiguous = "Ambiguous",
                        cross_group = "Cross-group",
                        no_recipient_loss = "No host"
                    )
                )
            }
            if ("pass" %in% names(display)) {
                display$pass <- ifelse(
                    is.na(display$pass), NA_character_,
                    ifelse(isTRUE_vec(display$pass), "Yes", "No")
                )
            }
            DT::datatable(
                display,
                colnames = bam_brief_call_labels(columns),
                rownames = FALSE,
                selection = "multiple",
                options = list(pageLength = 10, scrollX = TRUE)
            )
        })

        output$download_bundle <- shiny::downloadHandler(
            filename = function() {
                bam_bundle_filename(sample_value())
            },
            content = function(file) {
                res <- result_value()
                if (is.null(res)) {
                    stop(
                        "Run the analysis before downloading",
                        call. = FALSE
                    )
                }
                sample <- sample_value()
                if (is.null(sample)) {
                    sample <- karyotype_sample(res, NULL)
                }
                mode <- input$karyotype_mode
                if (length(mode) != 1L || is.na(mode) ||
                    !mode %in% c("heatmap", "calls")) {
                    mode <- "heatmap"
                }
                shown <- visible_calls()
                if (!is.null(shown) && nrow(shown) > 0L) {
                    shown <- shown[
                        as.character(shown$sample) == sample, ,
                        drop = FALSE
                    ]
                }
                current_payload <- karyotype_export_payload(
                    res,
                    sample = sample,
                    mode = mode,
                    calls = shown,
                    selected_groups = group_filter()
                )
                full_payload <- bam_full_karyotype_payload(current_payload)
                confident_calls <- bam_confident_calls(
                    res$calls[
                        as.character(res$calls$sample) == sample, ,
                        drop = FALSE
                    ]
                )
                confident_payload <- karyotype_payload(
                    res,
                    sample = sample,
                    mode = "calls",
                    calls = confident_calls,
                    selected_groups = NULL
                )
                bam_write_result_bundle(
                    res, file,
                    sample = sample,
                    current_payload = current_payload,
                    full_payload = full_payload,
                    confident_payload = confident_payload
                )
            }
        )

        output$download_bedpe <- shiny::downloadHandler(
            filename = function() {
                sample_name <- sample_value()
                if (is.null(sample_name)) {
                    "mamamia-selected-calls.bedpe"
                } else {
                    paste0(
                        "mamamia-",
                        gsub("[^A-Za-z0-9._-]+", "-", sample_name),
                        "-selected-calls.bedpe"
                    )
                }
            },
            content = function(file) {
                res <- result_value()
                x <- visible_calls()
                rows <- selected()
                if (is.null(res) || is.null(x) || length(rows) == 0L) {
                    stop(
                        "Select at least one call before downloading",
                        call. = FALSE
                    )
                }
                selected_calls <- x[rows, , drop = FALSE]
                writeCallsToBedpe(
                    selected_calls, file,
                    samples = sample_value(),
                    pass_only = FALSE,
                    windows = res$windows
                )
            }
        )
    })
}

#' Explore introgression calls in an interactive dashboard
#'
#' Opens a Shiny application for a `BamAnalysisResult` from
#' [analyzeBamInput()]. The dashboard offers class, class-status, copy-state,
#' confidence and width filters, a vertical karyotype with a score heatmap or
#' a call map, a brief sortable call table, a selected-call BEDPE export and a
#' results bundle with the full call table, BEDPE files and karyotype images.
#'
#' @param result A `BamAnalysisResult`.
#' @param host Host passed to [shiny::runApp()].
#' @param port Port passed to [shiny::runApp()]. Use `NULL` to let Shiny choose.
#' @param launch.browser Logical passed to [shiny::runApp()].
#'
#' @return The return value of [shiny::runApp()].
#' @export
exploreIntrogressions <- function(result,
                                  host = "127.0.0.1", port = NULL,
                                  launch.browser = interactive()) {
    if (!inherits(result, "BamAnalysisResult")) {
        stop(
            "`result` must be a BamAnalysisResult from analyzeBamInput().",
            call. = FALSE
        )
    }
    require_dashboard_packages()
    ui <- htmltools::attachDependencies(
        bslib::page_sidebar(
            title = "MaMaMIA introgression explorer",
            fillable = FALSE,
            theme = bam_dashboard_theme(),
            sidebar = bslib::sidebar(explorer_controls_ui("explorer")),
            explorer_results_ui("explorer")
        ),
        mamamia_font_dependency()
    )
    server <- function(input, output, session) {
        explorer_server("explorer", shiny::reactive(result))
    }
    app <- shiny::shinyApp(ui, server)
    args <- list(appDir = app, host = host, launch.browser = launch.browser)
    if (!is.null(port)) {
        args$port <- port
    }
    do.call(shiny::runApp, args)
}
