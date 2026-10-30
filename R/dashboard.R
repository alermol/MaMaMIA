## The introgression explorer: score profiles, calls, filters and BEDPE
## export.

dashboard_packages <- c(
    "bslib", "DT", "htmltools", "htmlwidgets", "shiny"
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

#' Build the explorer score profile
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
    p <- ggplot2::ggplot(
        scores,
        ggplot2::aes(
            x = .data$start / 1e6, y = .data$score,
            colour = .data$mode
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
        ggplot2::geom_line(linewidth = 0.3) +
        ggplot2::facet_wrap(~group, scales = "free_x", ncol = 1) +
        ggplot2::labs(
            x = "Position (Mb)", y = "Introgression score",
            colour = "Mode"
        ) +
        ggplot2::theme_minimal(base_size = 12)
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
            inherit.aes = FALSE, alpha = 0.15,
            fill = "firebrick"
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
    shiny::tagList(
        shiny::checkboxGroupInput(
            ns("class"), "Call class",
            choices = c(
                "substitution", "addition", "deletion", "complex", "other"
            ),
            selected = c(
                "substitution", "addition", "deletion", "complex", "other"
            ), inline = TRUE
        ),
        shiny::checkboxGroupInput(
            ns("class_status"), "Class status",
            choices = c("confident", "ambiguous", "needs_review"),
            selected = c("confident", "ambiguous", "needs_review"),
            inline = TRUE
        ),
        shiny::checkboxGroupInput(
            ns("copy_state"), "Copy state",
            choices = c(
                "balanced", "gain", "loss", "compound", "neutral"
            ),
            selected = c(
                "balanced", "gain", "loss", "compound", "neutral"
            ), inline = TRUE
        ),
        shiny::checkboxGroupInput(
            ns("confidence"), "Confidence",
            choices = c("high", "moderate", "low"),
            selected = c("high", "moderate", "low"), inline = TRUE
        ),
        shiny::checkboxInput(ns("pass_only"), "Show calibrated calls only", TRUE)
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
        shiny::textOutput(ns("selection_status")),
        shiny::downloadButton(ns("download_bedpe"), "Download selected BEDPE")
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
            )
        ),
        mamamia_karyotype_output(ns("karyotype"), height = height)
    )
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
            "Score profile",
            shiny::plotOutput(ns("score_plot"), height = "calc(100vh - 280px)")
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

        visible_scores <- shiny::reactive({
            x <- scores()
            if (is.null(x) || nrow(x) == 0L) {
                return(x)
            }
            groups <- group_filter()
            if (length(groups) > 0L) {
                x <- x[as.character(x$group) %in% groups, , drop = FALSE]
            }
            x
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

        output$score_plot <- shiny::renderPlot({
            x <- visible_scores()
            shiny::validate(
                shiny::need(!is.null(result_value()), "Run the analysis to plot"),
                shiny::need(
                    !is.null(x) && nrow(x) > 0L,
                    "Select at least one chromosome group"
                )
            )
            explorer_score_plot(x, visible_calls())
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
                    message = "Run the analysis to plot"
                )))
            }
            mamamia_karyotype(karyotype_payload(
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
            columns <- c(
                "sample", "group", "class", "copy_state", "class_group",
                "class_evidence", "class_confidence", "class_status",
                "total_evidence", "component_status", "mappability_status",
                "confidence", "recipient_group", "recipient_cross_group",
                "recipient_component", "recipient_chr_name",
                "recipient_start_bp", "recipient_end_bp",
                "recipient_deletion_width_bp", "recipient_assignment",
                "recipient_weight", "recipient_margin",
                "recipient_expected_breadth", "recipient_observed_breadth",
                "recipient_copy_number", "affected_fraction",
                "affected_copies", "dosage_confidence",
                "recipient_mean_depletion", "recipient_depletion_se",
                "recipient_depletion_lower", "recipient_width_overlap",
                "recipient_fraction_ratio", "recipient_mapq_fraction",
                "recipient_alternative", "recipient_alternative_score",
                "donor_component", "donor_chr_name",
                "start_bp", "end_bp", "width_bp",
                "score_mean", "donor_mean", "recipient_mean",
                "method", "mode", "pass"
            )
            DT::datatable(
                x[, intersect(columns, names(x)), drop = FALSE],
                rownames = FALSE,
                selection = "multiple",
                filter = "top",
                options = list(pageLength = 10, scrollX = TRUE)
            )
        })

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
#' [analyzeBamInput()]. The dashboard offers chromosome-group, class,
#' confidence and width filters, a vertical karyotype with a score heatmap or
#' a call map, a score profile with called intervals, a sortable call table and
#' a BEDPE export of the selected calls.
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
    ui <- bslib::page_sidebar(
        title = "MaMaMIA introgression explorer",
        fillable = FALSE,
        sidebar = bslib::sidebar(explorer_controls_ui("explorer")),
        explorer_results_ui("explorer")
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
