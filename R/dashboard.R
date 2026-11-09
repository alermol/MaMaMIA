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
            )
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
    shiny::tagList(
        shiny::numericInput(
            ns("min_width"),
            "Minimum call width (Mb)",
            value = 0,
            min = 0,
            step = 1
        ),
        shiny::div(
            class = "mamamia-call-columns",
            shinyWidgets::pickerInput(
                ns("call_columns"),
                "Additional columns",
                choices = bam_call_extra_column_choices(),
                selected = character(0),
                multiple = TRUE,
                width = "fit",
                options = shinyWidgets::pickerOptions(
                    actionsBox = TRUE,
                    liveSearch = TRUE,
                    title = "Add columns",
                    noneSelectedText = "Add columns",
                    selectedTextFormat = "count > 0",
                    countSelectedText = "{0} added"
                )
            )
        )
    )
}

bam_table_header <- function(metadata) {
    shiny::tags$table(
        class = "display",
        shiny::tags$thead(
            shiny::tags$tr(
                lapply(seq_len(nrow(metadata)), function(i) {
                    description <- metadata$description[[i]]
                    shiny::tags$th(
                        title = description,
                        `aria-label` = paste0(
                            metadata$label[[i]], ". ", description
                        ),
                        shiny::tags$span(
                            class = "mamamia-call-label",
                            metadata$label[[i]]
                        )
                    )
                })
            )
        )
    )
}

bam_call_table_header <- function(columns) {
    bam_table_header(bam_call_column_metadata(columns))
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

explorer_download_button <- function(id, label, class = "btn-sm",
                                     busy = "Preparing download...") {
    shiny::div(
        class = "mamamia-download-row",
        `data-busy-text` = busy,
        shiny::downloadButton(id, label, class = class),
        shiny::tags$span(
            class = "mamamia-download-status",
            role = "status",
            `aria-live` = "polite",
            `aria-hidden` = "true",
            shiny::tags$span(
                class = "mamamia-download-spinner",
                `aria-hidden` = "true"
            ),
            shiny::tags$span(
                class = "mamamia-download-status-text", busy
            )
        )
    )
}

explorer_export_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::div(
        class = "mamamia-export-actions",
        shiny::div(
            class = "mamamia-export-section mamamia-export-card",
            shiny::h5("Results bundle"),
            shiny::p(
                class = "mamamia-export-lead",
                paste(
                    "Download full and confident call tables, all-call,",
                    "passing and confident BEDPE files, and matching",
                    "karyotype SVGs in one ZIP archive. The archive has one",
                    "folder per sample plus run-wide summaries at the top",
                    "level."
                )
            ),
            explorer_download_button(
                ns("download_bundle"),
                "Download results bundle (.zip)",
                class = "btn-sm btn-primary",
                busy = "Preparing bundle..."
            )
        ),
        shiny::div(
            class = "mamamia-export-section mamamia-export-card",
            shiny::h5("All calls"),
            shiny::p(
                class = "mamamia-export-lead",
                paste(
                    "Download every call in the current run as one BEDPE",
                    "file, without the rest of the bundle."
                )
            ),
            explorer_download_button(
                ns("download_all_bedpe"),
                "Download all BEDPE",
                class = "btn-sm",
                busy = "Preparing BEDPE..."
            )
        ),
        shiny::div(
            class = "mamamia-export-section mamamia-export-card",
            shiny::h5("Selected calls"),
            shiny::textOutput(ns("selection_status")),
            explorer_download_button(
                ns("download_bedpe"),
                "Download selected BEDPE",
                class = "btn-sm",
                busy = "Preparing BEDPE..."
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
                "Evidence says why; Geometry describes what changed;",
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
                    ". This is a relative read-depth direction, not a",
                    "quantitative measure."
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
                    "is less reliable."
                )
            ),
            row(
                "not passing",
                "The call is shown for review but did not pass the configured filters."
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
                "Component, state and configured length."
            ),
            row(
                "Partial addition",
                paste(
                    "A donor chromosome bar appears when a partial addition",
                    "covers at least 20% of that donor. Only the added part",
                    "is blue; the remainder is an empty dashed donor deletion."
                )
            ),
            row(
                "Added donor fill",
                paste(
                    "An added donor chromosome uses the light-blue",
                    "chromosome-state colour. A call drawn on it keeps its",
                    "class colour. In heatmap mode the score bins cover the",
                    "call interval and the rest of the bar stays light blue."
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
            "Hide deletions",
            paste(
                "The switch beside the Plot buttons removes deletions from",
                "the karyotype instead of drawing them as dashed or red",
                "segments. The counter beside it reports how many segments",
                "were omitted in the current view."
            ),
            row(
                "Partial deletion",
                paste(
                    "The deleted interval is spliced out and the remaining",
                    "recipient chromosome closes the gap and becomes shorter."
                )
            ),
            row(
                "Complete deletion",
                "The recipient chromosome is not drawn at all."
            ),
            row(
                "Substitution",
                paste(
                    "The added donor chromosome takes the place of the",
                    "deleted recipient chromosome."
                )
            )
        ),
        section(
            "Caveats",
            NULL,
            row(
                "Relative, not absolute",
                paste(
                    "Copy state describes a relative read-depth direction.",
                    "It is not a quantitative measure."
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
            class = "mamamia-analysis-timer",
            `aria-live` = "off",
            shiny::span(
                class = "mamamia-analysis-timer-label", "Processing time"
            ),
            shiny::span(
                class = "mamamia-analysis-timer-value", "00:00"
            )
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
            shiny::div(
                class = "mamamia-karyotype-plot-controls",
                shiny::radioButtons(
                    ns("karyotype_mode"), "Plot",
                    choices = c("Heatmap" = "heatmap", "Calls" = "calls"),
                    selected = "heatmap", inline = TRUE
                ),
                shiny::div(
                    class = "mamamia-karyotype-deletions",
                    bslib::input_switch(
                        ns("hide_deletions"), "Hide deletions",
                        value = FALSE
                    ),
                    shiny::textOutput(ns("deletion_info"), inline = TRUE)
                )
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
#' @param hide_deletions Remove deleted chromosomes and deletion segments
#'   from the current view only.
#' @param full Optional unfiltered payload for the full export.
#'
#' @return A self-contained karyotype payload. The top level is the current
#'   view; `views$full` is a complete unfiltered payload with each chromosome
#'   carrying its own calls and scores.
#' @noRd
karyotype_export_payload <- function(result, sample = NULL,
                                     mode = c("heatmap", "calls"),
                                     calls = NULL, selected_groups = NULL,
                                     hide_deletions = FALSE,
                                     full = NULL) {
    mode <- match.arg(mode)
    payload <- karyotype_payload(
        result,
        sample = sample,
        mode = mode,
        calls = calls,
        selected_groups = selected_groups,
        hide_deletions = hide_deletions
    )
    if (is.null(full)) {
        full <- karyotype_payload(
            result,
            sample = sample,
            mode = mode,
            calls = NULL,
            selected_groups = NULL,
            hide_deletions = FALSE
        )
    }
    ## Every view is a complete payload. No positional `full_calls` or
    ## `full_chromosomes` pairing is needed by the exporter.
    payload$views <- list(full = full)
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

        visible_calls <- shiny::reactive({
            explorer_calls(calls(), list(
                pass_only = input$pass_only,
                groups = group_filter(),
                class = input$class,
                class_status = input$class_status,
                copy_state = input$copy_state,
                confidence = input$confidence,
                min_width_bp = suppressWarnings(
                    as.numeric(input$min_width)
                ) * 1e6
            ))$calls
        })

        selected <- shiny::reactive({
            x <- visible_calls()
            rows <- input$calls_rows_selected
            if (is.null(x) || length(rows) == 0L) {
                return(integer())
            }
            rows <- suppressWarnings(as.integer(rows))
            rows <- rows[!is.na(rows) & rows >= 1L & rows <= nrow(x)]
            unique(rows)
        })

        output$selection_status <- shiny::renderText({
            if (is.null(result_value())) {
                return("No analysis loaded yet")
            }
            paste(length(selected()), "call(s) selected for BEDPE export")
        })

        proxy <- DT::dataTableProxy("calls", session = session)
        ## A new sample or filter produces a different table; row indices from
        ## the previous table would silently select different calls.
        shiny::observeEvent(visible_calls(), {
            DT::selectRows(proxy, NULL)
        }, ignoreNULL = FALSE)
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

        ## The unfiltered export depends only on result, sample and mode, not
        ## on the call filters, so build it once per view instead of on every
        ## filter change.
        full_karyotype <- shiny::reactive({
            res <- result_value()
            mode <- input$karyotype_mode
            if (is.null(res) || length(mode) != 1L || is.na(mode) ||
                !mode %in% c("heatmap", "calls")) {
                return(NULL)
            }
            karyotype_payload(
                res,
                sample = sample_value(),
                mode = mode,
                calls = NULL,
                selected_groups = NULL,
                hide_deletions = FALSE
            )
        })

        ## The live view and its counter share one payload so the "Current
        ## view" export cannot disagree with the number shown by the switch.
        current_karyotype <- shiny::reactive({
            res <- result_value()
            if (is.null(res)) {
                return(NULL)
            }
            mode <- input$karyotype_mode
            if (length(mode) != 1L || is.na(mode) ||
                !mode %in% c("heatmap", "calls")) {
                mode <- "heatmap"
            }
            karyotype_export_payload(
                res,
                sample = sample_value(),
                mode = mode,
                calls = visible_calls(),
                selected_groups = group_filter(),
                hide_deletions = isTRUE(input$hide_deletions),
                full = full_karyotype()
            )
        })

        output$deletion_info <- shiny::renderText({
            payload <- current_karyotype()
            if (is.null(payload)) {
                return("")
            }
            if (!isTRUE(payload$hide_deletions)) {
                return("Deletions shown")
            }
            segments <- payload$hidden_deletions$segments
            if (is.null(segments) || length(segments) != 1L ||
                is.na(segments) || segments <= 0L) {
                return("No deletions omitted")
            }
            sprintf(
                "%d segment%s omitted",
                as.integer(segments),
                if (as.integer(segments) == 1L) "" else "s"
            )
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
            payload <- current_karyotype()
            if (is.null(payload)) {
                return(NULL)
            }
            mamamia_karyotype(payload)
        })

        output$calls <- DT::renderDT({
            shiny::validate(
                shiny::need(!is.null(result_value()), "Run the analysis to list calls")
            )
            x <- visible_calls()
            defaults <- intersect(bam_brief_call_columns(x), names(x))
            selected_columns <- input$call_columns
            if (is.null(selected_columns)) {
                selected_columns <- character(0)
            }
            extras <- names(x)[
                names(x) %in% as.character(selected_columns)
            ]
            columns <- unique(c(defaults, extras))
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
            logical_columns <- extras[vapply(
                x[extras], is.logical, logical(1)
            )]
            for (column in logical_columns) {
                values <- x[[column]]
                display[[column]] <- ifelse(
                    is.na(values), NA_character_,
                    ifelse(values, "Yes", "No")
                )
            }
            number_columns <- extras[vapply(
                x[extras], is.numeric, logical(1)
            )]
            table <- DT::datatable(
                display,
                container = bam_call_table_header(columns),
                rownames = FALSE,
                selection = "multiple",
                options = list(
                    pageLength = 10,
                    scrollX = TRUE,
                    orderMulti = TRUE,
                    columnDefs = list(list(
                        orderSequence = c("asc", "desc", ""),
                        targets = "_all"
                    )),
                    initComplete = DT::JS(
                        "function(settings) {
                            var api = this.api();
                            api.on('order.dt', function() {
                                var order = api.order();
                                var unsorted = order.some(function(item) {
                                    return item[1] === '';
                                });
                                if (unsorted) {
                                    api.order([]).draw();
                                }
                            });
                        }"
                    )
                )
            )
            if (length(number_columns) > 0L) {
                table <- DT::formatSignif(
                    table, number_columns, digits = 3
                )
            }
            table
        }, server = TRUE)

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
                bam_write_cohort_bundle(
                    res, file,
                    current_sample = sample,
                    current_payload = current_payload,
                    full_payload = full_payload,
                    confident_payload = confident_payload
                )
            }
        )

        output$download_all_bedpe <- shiny::downloadHandler(
            filename = function() {
                sample <- sample_value()
                if (is.null(sample)) {
                    "mamamia-all-calls.bedpe"
                } else {
                    paste0(
                        "mamamia-", bam_bundle_stem(sample),
                        "-all-calls.bedpe"
                    )
                }
            },
            content = function(file) {
                res <- result_value()
                if (is.null(res)) {
                    stop(
                        "Run the analysis before downloading",
                        call. = FALSE
                    )
                }
                calls <- res$calls
                sample <- sample_value()
                if (!is.null(sample)) {
                    calls <- calls[
                        as.character(calls$sample) == sample, ,
                        drop = FALSE
                    ]
                }
                bam_bundle_write_bedpe(
                    calls, res$windows, file, pass_only = FALSE
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
#' The results bundle contains one folder per sample plus run-wide summaries
#' at the top level.
#'
#' @param result A `BamAnalysisResult`.
#' @param host Host passed to [shiny::runApp()].
#' @param port Port passed to [shiny::runApp()]. Use `NULL` to let Shiny choose.
#' @param launch.browser Logical passed to [shiny::runApp()].
#'
#' @return The return value of [shiny::runApp()].
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
