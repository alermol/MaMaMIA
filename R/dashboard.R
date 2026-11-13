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

explorer_view_filters_ui <- function(id, lang = mamamia_language_default(),
                                     selected = list()) {
    ns <- shiny::NS(id)
    selected_value <- function(name) {
        value <- selected[[name]]
        if (is.null(value)) character(0) else as.character(value)
    }
    picker <- function(inputId, label, title, choices, selected) {
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
                countSelectedText = mamamia_tr("results.picker_selected", lang),
                liveSearch = FALSE
            )
        )
    }
    class_choices <- c(
        stats::setNames(
            "substitution", mamamia_tr("results.class_substitution", lang)
        ),
        stats::setNames(
            "addition", mamamia_tr("results.class_addition", lang)
        ),
        stats::setNames(
            "deletion", mamamia_tr("results.class_deletion", lang)
        ),
        stats::setNames("complex", mamamia_tr("results.class_complex", lang)),
        stats::setNames("other", mamamia_tr("results.class_other", lang))
    )
    status_choices <- c(
        stats::setNames(
            "confident", mamamia_tr("results.status_confident", lang)
        ),
        stats::setNames(
            "ambiguous", mamamia_tr("results.status_ambiguous", lang)
        ),
        stats::setNames(
            "needs_review", mamamia_tr("results.status_needs_review", lang)
        )
    )
    copy_choices <- c(
        stats::setNames("balanced", mamamia_tr("results.copy_balanced", lang)),
        stats::setNames("gain", mamamia_tr("results.copy_gain", lang)),
        stats::setNames("loss", mamamia_tr("results.copy_loss", lang)),
        stats::setNames("compound", mamamia_tr("results.copy_compound", lang)),
        stats::setNames("neutral", mamamia_tr("results.copy_neutral", lang))
    )
    confidence_choices <- c(
        stats::setNames("high", mamamia_tr("results.confidence_high", lang)),
        stats::setNames(
            "moderate", mamamia_tr("results.confidence_moderate", lang)
        ),
        stats::setNames("low", mamamia_tr("results.confidence_low", lang))
    )
    pass_only <- selected[["pass_only"]]
    if (is.null(pass_only) || length(pass_only) != 1L || is.na(pass_only)) {
        pass_only <- TRUE
    }
    shiny::div(
        class = "mamamia-filter-pills",
        picker(
            "class", mamamia_tr("results.filter_class", lang),
            mamamia_tr("results.filter_class_title", lang),
            class_choices, selected_value("class")
        ),
        picker(
            "class_status", mamamia_tr("results.filter_status", lang),
            mamamia_tr("results.filter_status_title", lang),
            status_choices, selected_value("class_status")
        ),
        picker(
            "copy_state", mamamia_tr("results.filter_copy", lang),
            mamamia_tr("results.filter_copy_title", lang),
            copy_choices, selected_value("copy_state")
        ),
        picker(
            "confidence", mamamia_tr("results.filter_confidence", lang),
            mamamia_tr("results.filter_confidence", lang),
            confidence_choices, selected_value("confidence")
        ),
        shiny::div(
            class = "mamamia-filter-switch",
            bslib::input_switch(
                ns("pass_only"), mamamia_tr("results.filter_passing", lang),
                value = isTRUE(pass_only)
            )
        )
    )
}

explorer_call_filters_ui <- function(id, lang = mamamia_language_default(),
                                     selected = list()) {
    ns <- shiny::NS(id)
    minimum <- selected[["min_width"]]
    if (is.null(minimum) || length(minimum) != 1L || is.na(minimum)) {
        minimum <- 0
    }
    columns <- selected[["call_columns"]]
    if (is.null(columns)) {
        columns <- character(0)
    }
    shiny::tagList(
        shiny::numericInput(
            ns("min_width"),
            mamamia_tr("results.min_width", lang),
            value = minimum,
            min = 0,
            step = 1
        ),
        shiny::div(
            class = "mamamia-call-columns",
            shinyWidgets::pickerInput(
                ns("call_columns"),
                mamamia_tr("results.additional_columns", lang),
                choices = bam_call_extra_column_choices(),
                selected = as.character(columns),
                multiple = TRUE,
                width = "fit",
                options = shinyWidgets::pickerOptions(
                    actionsBox = TRUE,
                    liveSearch = TRUE,
                    title = mamamia_tr("results.add_columns", lang),
                    noneSelectedText = mamamia_tr(
                        "results.add_columns", lang
                    ),
                    selectedTextFormat = "count > 0",
                    countSelectedText = mamamia_tr(
                        "results.picker_added", lang
                    )
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

bam_call_table_header <- function(columns,
                                  lang = mamamia_language_default()) {
    metadata <- bam_call_column_metadata(columns)
    keys <- mamamia_i18n_keys(lang)
    for (i in seq_len(nrow(metadata))) {
        column <- metadata$column[[i]]
        label_key <- paste0("call_table.label.", column)
        desc_key <- paste0("call_table.desc.", column)
        if (label_key %in% keys) {
            metadata$label[[i]] <- mamamia_tr(label_key, lang)
        }
        if (desc_key %in% keys) {
            metadata$description[[i]] <- mamamia_tr(desc_key, lang)
        }
    }
    bam_table_header(metadata)
}

explorer_selection_ui <- function(id, lang = mamamia_language_default()) {
    ns <- shiny::NS(id)
    shiny::div(
        class = "d-flex gap-2",
        shiny::actionButton(
            ns("select_all"), mamamia_tr("results.select_visible", lang),
            class = "btn-sm flex-fill"
        ),
        shiny::actionButton(
            ns("clear_selection"), mamamia_tr("results.clear", lang),
            class = "btn-sm flex-fill"
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

explorer_export_ui <- function(id, lang = mamamia_language_default()) {
    ns <- shiny::NS(id)
    shiny::div(
        class = "mamamia-export-actions",
        shiny::div(
            class = "mamamia-export-section mamamia-export-card",
            shiny::h5(mamamia_tr("export.results_bundle", lang)),
            shiny::p(
                class = "mamamia-export-lead",
                mamamia_tr("export.results_lead", lang)
            ),
            explorer_download_button(
                ns("download_bundle"),
                mamamia_tr("export.download_bundle", lang),
                class = "btn-sm btn-primary",
                busy = mamamia_tr("export.preparing_bundle", lang)
            )
        ),
        shiny::div(
            class = "mamamia-export-section mamamia-export-card",
            shiny::h5(mamamia_tr("export.all_calls", lang)),
            shiny::p(
                class = "mamamia-export-lead",
                mamamia_tr("export.all_calls_lead", lang)
            ),
            explorer_download_button(
                ns("download_all_bedpe"),
                mamamia_tr("export.download_all_bedpe", lang),
                class = "btn-sm",
                busy = mamamia_tr("export.preparing_bedpe", lang)
            )
        ),
        shiny::div(
            class = "mamamia-export-section mamamia-export-card",
            shiny::h5(mamamia_tr("export.selected_calls", lang)),
            shiny::textOutput(ns("selection_status")),
            explorer_download_button(
                ns("download_bedpe"),
                mamamia_tr("export.download_selected_bedpe", lang),
                class = "btn-sm",
                busy = mamamia_tr("export.preparing_bedpe", lang)
            )
        )
    )
}

explorer_tooltip_help_ui <- function(lang = mamamia_language_default()) {
    mark_strong <- function(text, tokens) {
        if (length(tokens) == 0L) {
            return(text)
        }
        marks <- paste0("\u0001", seq_along(tokens), "\u0002")
        for (i in seq_along(tokens)) {
            text <- gsub(tokens[[i]], marks[[i]], text, fixed = TRUE)
        }
        for (i in seq_along(tokens)) {
            text <- gsub(
                marks[[i]],
                paste0("<strong>", tokens[[i]], "</strong>"),
                text,
                fixed = TRUE
            )
        }
        text
    }
    row <- function(term, description, strong = character()) {
        shiny::tagList(
            shiny::tags$dt(term),
            shiny::tags$dd(shiny::HTML(
                mark_strong(as.character(description), strong)
            ))
        )
    }
    row_key <- function(term_key, desc_key = NULL, strong = character()) {
        desc <- if (is.null(desc_key)) "" else mamamia_tr(desc_key, lang)
        row(mamamia_tr(term_key, lang), desc, strong)
    }
    section <- function(title, intro = NULL, ...) {
        shiny::div(
            class = "mamamia-help-section",
            shiny::tags$h6(title),
            if (!is.null(intro) && nzchar(intro)) {
                shiny::p(class = "mamamia-help-intro", intro)
            },
            shiny::tags$dl(class = "mamamia-help-list", ...)
        )
    }
    tr <- function(key) mamamia_tr(key, lang)
    ru <- identical(lang, "ru")
    shiny::div(
        class = "mamamia-tooltip-help",
        shiny::p(class = "mamamia-help-lead", tr("tooltip_guide.lead")),
        shiny::div(
            class = "mamamia-help-callout",
            tr("tooltip_guide.callout")
        ),
        section(
            tr("tooltip_guide.quick_interpretation.title"),
            tr("tooltip_guide.quick_interpretation.intro"),
            row_key(
                "tooltip_guide.quick_interpretation.balanced_substitution",
                "tooltip_guide.quick_interpretation.balanced_substitution.desc",
                if (ru) {
                    c("почти нейтральная сумма", "замещение", "сбалансированное")
                } else {
                    c("near-neutral total", "substitution", "balanced")
                }
            ),
            row_key(
                "tooltip_guide.quick_interpretation.addition",
                "tooltip_guide.quick_interpretation.addition.desc",
                if (ru) "прирост" else "gain"
            ),
            row_key(
                "tooltip_guide.quick_interpretation.deletion",
                "tooltip_guide.quick_interpretation.deletion.desc",
                if (ru) "потеря" else "loss"
            ),
            row_key(
                "tooltip_guide.quick_interpretation.complex",
                "tooltip_guide.quick_interpretation.complex.desc"
            ),
            row_key(
                "tooltip_guide.quick_interpretation.other",
                "tooltip_guide.quick_interpretation.other.desc"
            ),
            row_key(
                "tooltip_guide.quick_interpretation.no_host",
                "tooltip_guide.quick_interpretation.no_host.desc"
            )
        ),
        section(
            tr("tooltip_guide.header.title"),
            tr("tooltip_guide.header.intro"),
            row_key(
                "tooltip_guide.header.class",
                "tooltip_guide.header.class.desc",
                if (ru) {
                    c("замещение", "добавление", "делеция", "сложный", "прочее")
                } else {
                    c("substitution", "addition", "deletion", "complex", "other")
                }
            ),
            row_key(
                "tooltip_guide.header.status",
                "tooltip_guide.header.status.desc",
                if (ru) {
                    c("уверенный", "неоднозначный", "требует проверки")
                } else {
                    c("confident", "ambiguous", "needs_review")
                }
            ),
            row_key(
                "tooltip_guide.header.context",
                "tooltip_guide.header.context.desc"
            ),
            row_key(
                "tooltip_guide.header.accent_colour",
                "tooltip_guide.header.accent_colour.desc"
            )
        ),
        section(
            tr("tooltip_guide.placement.title"),
            tr("tooltip_guide.placement.intro"),
            row_key(
                "tooltip_guide.placement.recipient",
                "tooltip_guide.placement.recipient.desc"
            ),
            row_key(
                "tooltip_guide.placement.donor",
                "tooltip_guide.placement.donor.desc"
            ),
            row_key(
                "tooltip_guide.placement.relation",
                "tooltip_guide.placement.relation.desc",
                if (ru) {
                    c("та же группа", "межгрупповой кандидат", "не размещено")
                } else {
                    c("same group", "cross-group candidate", "not placed")
                }
            ),
            row_key(
                "tooltip_guide.placement.assignment",
                "tooltip_guide.placement.assignment.desc",
                if (ru) {
                    c(
                        "размещён", "межгрупповой кандидат", "неоднозначный",
                        "нет реципиента", "только группа"
                    )
                } else {
                    c(
                        "assigned", "cross-group candidate", "ambiguous",
                        "no host", "group only"
                    )
                }
            ),
            row_key(
                "tooltip_guide.placement.margin",
                "tooltip_guide.placement.margin.desc"
            ),
            row_key(
                "tooltip_guide.placement.group_class",
                "tooltip_guide.placement.group_class.desc"
            ),
            row_key(
                "tooltip_guide.placement.alternative",
                "tooltip_guide.placement.alternative.desc"
            ),
            row_key(
                "tooltip_guide.placement.reading",
                "tooltip_guide.placement.reading.desc"
            )
        ),
        section(
            tr("tooltip_guide.evidence.title"),
            tr("tooltip_guide.evidence.intro"),
            row_key(
                "tooltip_guide.evidence.score",
                "tooltip_guide.evidence.score.desc"
            ),
            row_key(
                "tooltip_guide.evidence.donor",
                "tooltip_guide.evidence.donor.desc"
            ),
            row_key(
                "tooltip_guide.evidence.recipient",
                "tooltip_guide.evidence.recipient.desc"
            ),
            row_key(
                "tooltip_guide.evidence.total",
                "tooltip_guide.evidence.total.desc"
            ),
            row_key(
                "tooltip_guide.evidence.colour",
                "tooltip_guide.evidence.colour.desc"
            ),
            row_key(
                "tooltip_guide.evidence.copy_state",
                "tooltip_guide.evidence.copy_state.desc",
                if (ru) {
                    c(
                        "сбалансированное", "прирост", "потеря", "составное",
                        "нейтральное"
                    )
                } else {
                    c("balanced", "gain", "loss", "compound", "neutral")
                }
            ),
            row_key(
                "tooltip_guide.evidence.evidence_source",
                "tooltip_guide.evidence.evidence_source.desc",
                if (ru) {
                    c("донор+реципиент", "донор", "реципиент", "нет")
                } else {
                    c("donor+recipient", "donor", "recipient", "none")
                }
            ),
            row_key(
                "tooltip_guide.evidence.typical_pattern",
                "tooltip_guide.evidence.typical_pattern.desc"
            )
        ),
        section(
            tr("tooltip_guide.status.title"),
            tr("tooltip_guide.status.intro"),
            row_key(
                "tooltip_guide.status.class_confidence",
                "tooltip_guide.status.class_confidence.desc"
            ),
            row_key(
                "tooltip_guide.status.evidence_sides",
                "tooltip_guide.status.evidence_sides.desc"
            ),
            row_key(
                "tooltip_guide.status.copy_states",
                "tooltip_guide.status.copy_states.desc"
            ),
            row_key(
                "tooltip_guide.status.total",
                "tooltip_guide.status.total.desc"
            ),
            row_key(
                "tooltip_guide.status.components",
                "tooltip_guide.status.components.desc"
            ),
            row_key(
                "tooltip_guide.status.low_mappability",
                "tooltip_guide.status.low_mappability.desc"
            ),
            row_key(
                "tooltip_guide.status.not_passing",
                "tooltip_guide.status.not_passing.desc"
            )
        ),
        section(
            tr("tooltip_guide.geometry.title"),
            tr("tooltip_guide.geometry.intro"),
            row_key(
                "tooltip_guide.geometry.donor_width",
                "tooltip_guide.geometry.donor_width.desc"
            ),
            row_key(
                "tooltip_guide.geometry.recipient_run",
                "tooltip_guide.geometry.recipient_run.desc"
            ),
            row_key(
                "tooltip_guide.geometry.overlap",
                "tooltip_guide.geometry.overlap.desc"
            ),
            row_key(
                "tooltip_guide.geometry.width_ratio",
                "tooltip_guide.geometry.width_ratio.desc"
            ),
            row_key(
                "tooltip_guide.geometry.interpretation",
                "tooltip_guide.geometry.interpretation.desc"
            )
        ),
        section(
            tr("tooltip_guide.position.title"),
            tr("tooltip_guide.position.intro"),
            row_key(
                "tooltip_guide.position.position",
                "tooltip_guide.position.position.desc"
            ),
            row_key(
                "tooltip_guide.position.score",
                "tooltip_guide.position.score.desc"
            ),
            row_key(
                "tooltip_guide.position.mode",
                "tooltip_guide.position.mode.desc"
            ),
            row_key(
                "tooltip_guide.position.click",
                "tooltip_guide.position.click.desc"
            )
        ),
        section(
            tr("tooltip_guide.no_call.title"),
            tr("tooltip_guide.no_call.intro"),
            row_key(
                "tooltip_guide.no_call.header",
                "tooltip_guide.no_call.header.desc"
            ),
            row_key(
                "tooltip_guide.no_call.chromosome_block",
                "tooltip_guide.no_call.chromosome_block.desc"
            ),
            row_key(
                "tooltip_guide.no_call.partial_addition",
                "tooltip_guide.no_call.partial_addition.desc"
            ),
            row_key(
                "tooltip_guide.no_call.added_donor_fill",
                "tooltip_guide.no_call.added_donor_fill.desc"
            ),
            row_key(
                "tooltip_guide.no_call.call_block",
                "tooltip_guide.no_call.call_block.desc"
            ),
            row_key(
                "tooltip_guide.no_call.position_block",
                "tooltip_guide.no_call.position_block.desc"
            )
        ),
        section(
            tr("tooltip_guide.hide_deletions.title"),
            tr("tooltip_guide.hide_deletions.intro"),
            row_key(
                "tooltip_guide.hide_deletions.partial_deletion",
                "tooltip_guide.hide_deletions.partial_deletion.desc"
            ),
            row_key(
                "tooltip_guide.hide_deletions.complete_deletion",
                "tooltip_guide.hide_deletions.complete_deletion.desc"
            ),
            row_key(
                "tooltip_guide.hide_deletions.substitution",
                "tooltip_guide.hide_deletions.substitution.desc"
            )
        ),
        section(
            tr("tooltip_guide.unplaced.title"),
            tr("tooltip_guide.unplaced.intro")
        ),
        section(
            tr("tooltip_guide.caveats.title"),
            NULL,
            row_key(
                "tooltip_guide.caveats.relative",
                "tooltip_guide.caveats.relative.desc"
            ),
            row_key(
                "tooltip_guide.caveats.class_vs_placement",
                "tooltip_guide.caveats.class_vs_placement.desc"
            ),
            row_key(
                "tooltip_guide.caveats.cross_group",
                "tooltip_guide.caveats.cross_group.desc"
            ),
            row_key(
                "tooltip_guide.caveats.filters",
                "tooltip_guide.caveats.filters.desc"
            ),
            row_key(
                "tooltip_guide.caveats.heatmap_score",
                "tooltip_guide.caveats.heatmap_score.desc"
            )
        )
    )
}

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
        shiny::div(
            class = "mamamia-analysis-eyebrow",
            `data-mamamia-i18n` = "run.analysis_eyebrow",
            "Analysis"
        ),
        shiny::div(
            class = "mamamia-analysis-stage", "Preparing analysis"
        ),
        shiny::div(
            class = "mamamia-analysis-timer",
            `aria-live` = "off",
            shiny::span(
                class = "mamamia-analysis-timer-label",
                `data-mamamia-i18n` = "run.processing_time",
                "Processing time"
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

explorer_karyotype_controls_ui <- function(id,
                                             lang = mamamia_language_default(),
                                             selected = list()) {
    ns <- shiny::NS(id)
    mode <- selected[["mode"]]
    if (is.null(mode) || length(mode) != 1L || is.na(mode) ||
        !mode %in% c("heatmap", "calls")) {
        mode <- "heatmap"
    }
    hide_deletions <- selected[["hide_deletions"]]
    if (is.null(hide_deletions) || length(hide_deletions) != 1L ||
        is.na(hide_deletions)) {
        hide_deletions <- FALSE
    }
    hide_unplaced <- selected[["hide_unplaced"]]
    if (is.null(hide_unplaced) || length(hide_unplaced) != 1L ||
        is.na(hide_unplaced)) {
        hide_unplaced <- TRUE
    }
    shiny::div(
        style = paste(
            "display:flex; flex-wrap:wrap; align-items:center;",
            "gap:.75rem; padding:.2rem .3rem .1rem;"
        ),
        shiny::div(
            class = "mamamia-karyotype-plot-controls",
            shiny::radioButtons(
                ns("karyotype_mode"), mamamia_tr("karyotype.plot", lang),
                choices = stats::setNames(
                    c("heatmap", "calls"),
                    c(
                        mamamia_tr("karyotype.heatmap", lang),
                        mamamia_tr("karyotype.calls", lang)
                    )
                ),
                selected = mode, inline = TRUE
            ),
            shiny::div(
                class = "mamamia-karyotype-deletions",
                bslib::input_switch(
                    ns("hide_deletions"),
                    mamamia_tr("karyotype.hide_deletions", lang),
                    value = isTRUE(hide_deletions)
                ),
                shiny::textOutput(ns("deletion_info"), inline = TRUE)
            ),
            shiny::div(
                class = "mamamia-karyotype-deletions",
                bslib::input_switch(
                    ns("hide_unplaced"),
                    mamamia_tr("karyotype.hide_unplaced", lang),
                    value = isTRUE(hide_unplaced)
                ),
                shiny::textOutput(ns("unplaced_info"), inline = TRUE)
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
                title = mamamia_tr("karyotype.download_current_title", lang),
                shiny::icon("download"),
                mamamia_tr("karyotype.current_view", lang)
            ),
            shiny::tags$button(
                id = ns("download_karyotype_full"),
                type = "button",
                class = "btn btn-sm btn-outline-secondary",
                disabled = "disabled",
                title = mamamia_tr("karyotype.download_full_title", lang),
                shiny::icon("download"),
                mamamia_tr("karyotype.full_view", lang)
            ),
            shiny::actionButton(
                ns("tooltip_help"), mamamia_tr("karyotype.tooltip_help", lang),
                icon = shiny::icon("circle-question"),
                class = "btn-sm btn-outline-secondary"
            )
        )
    )
}

explorer_karyotype_ui <- function(id, height = "calc(100vh - 360px)",
                                  lang = mamamia_language_default(),
                                  dynamic = FALSE) {
    ns <- shiny::NS(id)
    controls <- if (isTRUE(dynamic)) {
        shiny::uiOutput("karyotype_controls_ui")
    } else {
        explorer_karyotype_controls_ui(id, lang = lang)
    }
    shiny::div(
        class = "mamamia-karyotype",
        controls,
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
#' @param hide_unplaced Remove unplaced calls from the current view only.
#' @param show_unplaced Draw neutral unplaced donor bars when unplaced calls
#'   are not hidden.
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
                                     hide_unplaced = FALSE,
                                     show_unplaced = FALSE,
                                     full = NULL,
                                     lang = mamamia_language_default()) {
    mode <- match.arg(mode)
    payload <- karyotype_payload(
        result,
        sample = sample,
        mode = mode,
        calls = calls,
        selected_groups = selected_groups,
        hide_deletions = hide_deletions,
        hide_unplaced = hide_unplaced,
        show_unplaced = show_unplaced,
        lang = lang
    )
    if (is.null(full)) {
        full <- karyotype_payload(
            result,
            sample = sample,
            mode = mode,
            calls = NULL,
            selected_groups = NULL,
            hide_deletions = FALSE,
            hide_unplaced = FALSE,
            show_unplaced = FALSE,
            lang = lang
        )
    }
    ## Every view is a complete payload. No positional `full_calls` or
    ## `full_chromosomes` pairing is needed by the exporter.
    payload$views <- list(full = full)
    payload
}

explorer_controls_ui <- function(id, lang = mamamia_language_default()) {
    shiny::tagList(
        explorer_view_filters_ui(id, lang = lang),
        explorer_call_filters_ui(id, lang = lang),
        explorer_selection_ui(id, lang = lang),
        explorer_export_ui(id, lang = lang)
    )
}

explorer_results_ui <- function(id, lang = mamamia_language_default()) {
    ns <- shiny::NS(id)
    bslib::navset_card_tab(
        id = ns("tabs"),
        bslib::nav_panel(
            mamamia_tr("tabs.karyotype", lang),
            value = "karyotype",
            explorer_karyotype_ui(
                id, height = "calc(100vh - 360px)", lang = lang
            )
        ),
        bslib::nav_panel(
            mamamia_tr("tabs.calls", lang),
            value = "calls",
            shiny::div(
                class = "small text-muted mb-2",
                mamamia_tr("results.unplaced_note", lang)
            ),
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
explorer_server <- function(id, result, sample = NULL,
                            lang = shiny::reactive(
                                mamamia_language_default()
                            )) {
    shiny::moduleServer(id, function(input, output, session) {
        lang_value <- function() {
            if (shiny::is.reactive(lang)) lang() else lang
        }
        selected_group <- shiny::reactiveVal(NULL)
        bundle_cache <- bundle_cache_state()
        session$onSessionEnded(function() bundle_cache_clear(bundle_cache))
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

        visible_call_state <- shiny::reactive({
            raw <- calls()
            if (!is.null(raw)) {
                raw$.call_id <- seq_len(nrow(raw))
            }
            explorer_calls(raw, list(
                pass_only = input$pass_only,
                groups = group_filter(),
                class = input$class,
                class_status = input$class_status,
                copy_state = input$copy_state,
                confidence = input$confidence,
                min_width_bp = suppressWarnings(
                    as.numeric(input$min_width)
                ) * 1e6
            ))
        })

        visible_calls <- shiny::reactive({
            visible_call_state()$calls
        })

        ## The table follows exactly what the current karyotype payload placed
        ## on a chromosome. `current_karyotype()` is defined below; reactive
        ## bodies resolve it when they run.
        rendered_calls <- shiny::reactive({
            payload <- current_karyotype()
            explorer_rendered_calls(
                visible_calls(),
                if (is.null(payload)) NULL else payload$rendered_call_ids
            )
        })

        output$unplaced_info <- shiny::renderText({
            if (is.null(result_value())) {
                return("")
            }
            payload <- current_karyotype()
            if (is.null(payload)) {
                return("")
            }
            ids <- payload$unplaced_call_ids
            count <- if (is.null(ids)) 0L else length(ids)
            if (count <= 0L) {
                return(mamamia_tr("results.unplaced_none", lang_value()))
            }
            if (isTRUE(input$hide_unplaced)) {
                sprintf(
                    mamamia_tr("results.unplaced_hidden", lang_value()),
                    as.integer(count)
                )
            } else {
                sprintf(
                    mamamia_tr("results.unplaced_shown", lang_value()),
                    as.integer(count)
                )
            }
        })

        selected <- shiny::reactive({
            x <- rendered_calls()
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
                return(mamamia_tr("results.no_analysis", lang_value()))
            }
            n <- length(selected())
            if (identical(lang_value(), "ru")) {
                sprintf(
                    mamamia_tr("results.selection_status", "ru"),
                    as.integer(n)
                )
            } else {
                paste(n, mamamia_tr("results.selection_status", "en"))
            }
        })

        proxy <- DT::dataTableProxy("calls", session = session)
        ## A new sample or filter produces a different table; row indices from
        ## the previous table would silently select different calls.
        shiny::observeEvent(rendered_calls(), {
            DT::selectRows(proxy, NULL)
        }, ignoreNULL = FALSE)
        shiny::observeEvent(input$select_all, {
            x <- rendered_calls()
            if (!is.null(x) && nrow(x) > 0L) {
                DT::selectRows(proxy, seq_len(nrow(x)))
            }
        })
        shiny::observeEvent(input$clear_selection, {
            DT::selectRows(proxy, NULL)
        })

        shiny::observeEvent(input$tooltip_help, {
            dialog <- shiny::modalDialog(
                title = mamamia_tr("export.tooltip_guide", lang_value()),
                explorer_tooltip_help_ui(lang_value()),
                size = "l",
                easyClose = TRUE,
                footer = shiny::modalButton(
                    mamamia_tr("errors.modal_close", lang_value())
                )
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
                hide_deletions = FALSE,
                hide_unplaced = FALSE,
                show_unplaced = FALSE,
                lang = lang_value()
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
                hide_unplaced = isTRUE(input$hide_unplaced),
                show_unplaced = !isTRUE(input$hide_unplaced),
                full = full_karyotype(),
                lang = lang_value()
            )
        })

        bundle_request_value <- shiny::reactive({
            res <- result_value()
            payload <- current_karyotype()
            if (is.null(res) || is.null(payload)) {
                return(NULL)
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
            bundle_request(
                res,
                sample = sample,
                mode = mode,
                calls = rendered_calls(),
                groups = group_filter(),
                hide_deletions = isTRUE(input$hide_deletions),
                hide_unplaced = isTRUE(input$hide_unplaced),
                current_payload = payload
            )
        })

        output$deletion_info <- shiny::renderText({
            payload <- current_karyotype()
            if (is.null(payload)) {
                return("")
            }
            if (!isTRUE(payload$hide_deletions)) {
                return(mamamia_tr("results.deletions_shown", lang_value()))
            }
            segments <- payload$hidden_deletions$segments
            if (is.null(segments) || length(segments) != 1L ||
                is.na(segments) || segments <= 0L) {
                return(mamamia_tr("results.deletions_none", lang_value()))
            }
            if (identical(lang_value(), "ru")) {
                sprintf(
                    mamamia_tr("results.deletions_omitted", "ru"),
                    as.integer(segments)
                )
            } else {
                sprintf(
                    mamamia_tr("results.deletions_omitted", "en"),
                    as.integer(segments),
                    if (as.integer(segments) == 1L) "" else "s"
                )
            }
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
                    message = mamamia_tr("results.run_to_plot", lang_value()),
                    hint = mamamia_tr("results.run_hint", lang_value()),
                    labels = mamamia_karyotype_labels(lang_value())
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
                shiny::need(
                    !is.null(result_value()),
                    mamamia_tr("results.run_to_list", lang_value())
                )
            )
            x <- rendered_calls()
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
                    substitution = mamamia_tr(
                        "results.class_substitution", lang_value()
                    ),
                    addition = mamamia_tr(
                        "results.class_addition", lang_value()
                    ),
                    deletion = mamamia_tr(
                        "results.class_deletion", lang_value()
                    ),
                    complex = mamamia_tr(
                        "results.class_complex", lang_value()
                    ),
                    other = mamamia_tr(
                        "results.class_other", lang_value()
                    )
                ))
            }
            if ("copy_state" %in% names(display)) {
                display$copy_state <- humanise(display$copy_state, c(
                    balanced = mamamia_tr(
                        "results.copy_balanced", lang_value()
                    ),
                    gain = mamamia_tr("results.copy_gain", lang_value()),
                    loss = mamamia_tr("results.copy_loss", lang_value()),
                    compound = mamamia_tr(
                        "results.copy_compound", lang_value()
                    ),
                    neutral = mamamia_tr(
                        "results.copy_neutral", lang_value()
                    )
                ))
            }
            if ("class_status" %in% names(display)) {
                display$class_status <- humanise(display$class_status, c(
                    confident = mamamia_tr(
                        "results.status_confident", lang_value()
                    ),
                    ambiguous = mamamia_tr(
                        "results.status_ambiguous", lang_value()
                    ),
                    needs_review = mamamia_tr(
                        "results.status_needs_review", lang_value()
                    )
                ))
            }
            if ("recipient_assignment" %in% names(display)) {
                display$recipient_assignment <- humanise(
                    display$recipient_assignment,
                    c(
                        assigned = mamamia_tr(
                            "results.assignment_assigned", lang_value()
                        ),
                        ambiguous = mamamia_tr(
                            "results.status_ambiguous", lang_value()
                        ),
                        cross_group = mamamia_tr(
                            "results.assignment_cross_group", lang_value()
                        ),
                        no_recipient_loss = mamamia_tr(
                            "results.assignment_no_recipient_loss",
                            lang_value()
                        )
                    )
                )
            }
            if ("pass" %in% names(display)) {
                display$pass <- ifelse(
                    is.na(display$pass), NA_character_,
                    ifelse(
                        isTRUE_vec(display$pass),
                        mamamia_tr("common.yes", lang_value()),
                        mamamia_tr("common.no", lang_value())
                    )
                )
            }
            logical_columns <- extras[vapply(
                x[extras], is.logical, logical(1)
            )]
            for (column in logical_columns) {
                values <- x[[column]]
                display[[column]] <- ifelse(
                    is.na(values), NA_character_,
                    ifelse(
                        values,
                        mamamia_tr("common.yes", lang_value()),
                        mamamia_tr("common.no", lang_value())
                    )
                )
            }
            number_columns <- extras[vapply(
                x[extras], is.numeric, logical(1)
            )]
            table <- DT::datatable(
                display,
                container = bam_call_table_header(
                    columns, lang = lang_value()
                ),
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
                request <- bundle_request_value()
                if (is.null(request)) {
                    stop(
                        mamamia_tr(
                            "results.run_before_download", lang_value()
                        ),
                        call. = FALSE
                    )
                }
                key <- bundle_request_key(request)
                bundle_cache_serve(bundle_cache, key, file, function(path) {
                    res <- request$result
                    sample <- request$sample
                    current_payload <- request$current_payload
                    full_payload <- bam_full_karyotype_payload(
                        current_payload
                    )
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
                        res, path,
                        current_sample = sample,
                        current_payload = current_payload,
                        full_payload = full_payload,
                        confident_payload = confident_payload
                    )
                })
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
                        mamamia_tr(
                            "results.run_before_download", lang_value()
                        ),
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
                x <- rendered_calls()
                rows <- selected()
                if (is.null(res) || is.null(x) || length(rows) == 0L) {
                    stop(
                        mamamia_tr(
                            "results.select_before_download", lang_value()
                        ),
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
