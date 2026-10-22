## The introgression explorer: the Shiny module with the chromosome and pair
## filters, the segment table and the BEDPE export, plus the
## `exploreIntrogressions()` launcher around it.

dashboard_packages <- c("bslib", "DT", "gridlayout", "shiny", "shinyWidgets")

## The mirror plot draws two points per window, so a high-resolution run has
## hundreds of thousands of them. The dashboard thins the drawn points to this
## cap; segments and the BEDPE export stay complete.
dashboard_max_plot_points <- 20000L

require_dashboard_packages <- function() {
    missing <- dashboard_packages[!vapply(dashboard_packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) == 0L) {
        return(invisible(TRUE))
    }
    cran <- setdiff(missing, "gridlayout")
    instructions <- character()
    if (length(cran) > 0L) {
        instructions <- c(
            instructions,
            sprintf(
                "install.packages(c(%s))",
                paste(sprintf('"%s"', cran), collapse = ", ")
            )
        )
    }
    if ("gridlayout" %in% missing) {
        instructions <- c(
            instructions,
            "remotes::install_github(\"rstudio/gridlayout\")"
        )
    }
    stop(
        "Install dashboard dependencies first: ",
        paste(instructions, collapse = "; "),
        call. = FALSE
    )
}

isa_dashboard_segments <- function(ISA) {
    segments <- ISA$segments
    segments$.row_id <- seq_len(nrow(segments))
    if (nrow(segments) == 0L) {
        return(segments)
    }
    annotate_segment_coords(segments, ISA$out)
}

## Recipients first, then donors in reverse order; this is the order the
## explorer uses for its chromosome pickers.
isa_dashboard_genome <- function(ISA, chr_ids = NULL) {
    lengths <- stats::setNames(ISA$data$stop, ISA$data$chr_id)
    lengths <- tapply(lengths, names(lengths), max, na.rm = TRUE)
    lengths <- lengths[is.finite(lengths)]
    if (!is.null(chr_ids)) {
        lengths <- lengths[intersect(chr_ids, names(lengths))]
    }
    is_donor <- names(lengths) %in% ISA$meta$don_chr_ids
    as.list(lengths[c(names(lengths)[!is_donor], rev(names(lengths)[is_donor]))])
}


isa_dashboard_pair_labels <- function(segments) {
    if (nrow(segments) == 0L) {
        return(character())
    }
    paste(segments$chr_id.don, segments$chr_id.rec, sep = " -> ")
}

filter_dashboard_segments <- function(segments, chr_ids, pair_ids, putative_only) {
    keep <- rep(TRUE, nrow(segments))
    if (!is.null(chr_ids)) {
        keep <- keep & segments$chr_id.don %in% chr_ids & segments$chr_id.rec %in% chr_ids
    }
    if (!is.null(pair_ids)) {
        keep <- keep & isa_dashboard_pair_labels(segments) %in% pair_ids
    }
    if (isTRUE(putative_only)) {
        keep <- keep & segments$putative_introgression
    }
    segments[keep, , drop = FALSE]
}

selected_isa_for_bedpe <- function(ISA, selected_rows) {
    out <- ISA
    out$segments$putative_introgression <- seq_len(nrow(out$segments)) %in% selected_rows
    out
}

## Subset an ISA to the explorer's chromosome, pair and width filters so that
## the plot follows the table. Only the pieces `plot.ISA()` reads are subset.
## Returns `NULL` when no window survives the filters.
filter_dashboard_isa <- function(ISA, chr_ids = NULL, pair_ids = NULL,
                                 min_width = NULL) {
    out <- ISA$out
    keep <- rep(TRUE, nrow(out))
    if (!is.null(chr_ids)) {
        keep <- keep & out$chr_id.don %in% chr_ids & out$chr_id.rec %in% chr_ids
    }
    pair_labels <- paste(out$chr_id.don, out$chr_id.rec, sep = " -> ")
    if (!is.null(pair_ids)) {
        keep <- keep & pair_labels %in% pair_ids
    }
    out <- out[keep, , drop = FALSE]
    if (nrow(out) == 0L) {
        return(NULL)
    }

    segments <- ISA$segments
    segment_keep <- isa_dashboard_pair_labels(segments) %in%
        paste(out$chr_id.don, out$chr_id.rec, sep = " -> ")
    if (!is.null(min_width)) {
        segment_keep <- segment_keep &
            (segments$loc.end - segments$loc.start + 1L) >= min_width
    }

    ISA$out <- out
    ISA$segments <- segments[segment_keep, , drop = FALSE]
    ISA
}


isa_dashboard_table <- function(segments) {
    cols <- c(
        ".row_id", "chr_id.don", "start.don", "stop.don",
        "chr_id.rec", "start.rec", "stop.rec", "num.mark",
        "seg.mean", "seg.median", "putative_introgression",
        "chr_name.don", "chr_name.rec", "subgenome.don", "subgenome.rec"
    )
    out <- segments[intersect(cols, names(segments))]
    names(out)[names(out) == ".row_id"] <- "row"
    out
}


## Reusable explorer UI: the filter controls can live in a page sidebar
## (`exploreIntrogressions()`) or in a card (`runMaMaMIA()`).
explorer_controls_ui <- function(id) {
    ns <- shiny::NS(id)
    shiny::tagList(
        shinyWidgets::pickerInput(
            ns("chr_ids"), "Chromosomes",
            choices = NULL, selected = NULL,
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `live-search` = TRUE)
        ),
        shinyWidgets::pickerInput(
            ns("pair_ids"), "Chromosome pairs",
            choices = NULL, selected = NULL,
            multiple = TRUE,
            options = list(`actions-box` = TRUE, `live-search` = TRUE)
        ),
        shiny::checkboxInput(ns("putative_only"), "Show putative introgressions only", TRUE),
        shiny::radioButtons(
            ns("threshold_type"), "Call threshold metric",
            choices = c("Median" = "median", "Mean" = "mean"),
            inline = TRUE
        ),
        shiny::numericInput(ns("threshold_value"), "Minimum segment metric", value = 0),
        shiny::numericInput(ns("min_width"), "Minimum segment width (windows)", value = 2, min = 2, step = 1),
        shiny::div(
            class = "d-flex gap-2",
            shiny::actionButton(ns("select_all"), "Select visible", class = "btn-sm flex-fill"),
            shiny::actionButton(ns("clear_selection"), "Clear", class = "btn-sm flex-fill")
        ),
        shiny::textOutput(ns("selection_status")),
        shiny::downloadButton(ns("download_bedpe"), "Download selected BEDPE")
    )
}

explorer_results_ui <- function(id) {
    ns <- shiny::NS(id)
    bslib::navset_card_tab(
        id = ns("tabs"),
        bslib::nav_panel(
            "Introgression circos",
            isa_dashboard_circos_output(
                ns("circos"), height = "calc(100vh - 265px)"
            )
        ),
        bslib::nav_panel(
            "Coverage / dosage",
            shiny::radioButtons(
                ns("plot_type"), "Plot",
                choices = c("Mirror" = "mirror", "Difference" = "diff"),
                inline = TRUE
            ),
            shiny::plotOutput(ns("isa_plot"), height = "calc(100vh - 335px)")
        ),
        bslib::nav_panel(
            "Segments",
            DT::DTOutput(ns("segments"))
        )
    )
}

## Server for the explorer. `ISA` is a reactive returning an `ISA` object, or
## `NULL` while no analysis has run yet; the explorer stays idle in that case.
## `chr_ids` optionally limits the chromosomes selected on arrival. `sample` is
## an optional sample name (or reactive returning one); when supplied, the BEDPE
## download gains the `sample` column, so a per-sample file can be merged into
## a cohort export without losing its identity.
##
## `selection` optionally supplies the row ids (the `.row_id` values) to
## preselect when the table is rendered, so the parent app can restore the rows
## chosen for another sample when the user switches back. `on_selection` is
## called with `(rows, sample)` whenever the current selection changes; the
## parent uses it to keep a per-sample selection map for the combined export.
explorer_server <- function(id, ISA, chr_ids = NULL, sample = NULL,
                            selection = NULL, on_selection = NULL) {
    shiny::moduleServer(id, function(input, output, session) {
        ## Changing the threshold re-flags the segments without re-running CBS
        ## (see `subset.ISA`); the segment rows are kept in their original order,
        ## so `.row_id` remains the stable key for the table and the export.
        isa_thresholded <- shiny::reactive({
            isa <- ISA()
            if (is.null(isa)) {
                return(NULL)
            }
            threshold_type <- input$threshold_type
            if (!identical(threshold_type, "mean") && !identical(threshold_type, "median")) {
                threshold_type <- if (identical(isa$meta$metrics, "seg.mean")) "mean" else "median"
            }
            threshold_value <- input$threshold_value
            if (!is.numeric(threshold_value) || length(threshold_value) != 1L || is.na(threshold_value)) {
                threshold_value <- isa$meta$metrics_threshold
            }
            tryCatch(
                if (identical(threshold_type, "mean")) {
                    subset(isa, threshold_type = "mean", min_segment_mean = threshold_value)
                } else {
                    subset(isa, threshold_type = "median", min_segment_median = threshold_value)
                },
                error = function(e) {
                    shiny::showNotification(
                        paste("Could not apply the call threshold:", conditionMessage(e)),
                        type = "warning", duration = 8
                    )
                    isa
                }
            )
        })

        current_segments <- shiny::reactive({
            isa <- isa_thresholded()
            if (is.null(isa)) {
                return(NULL)
            }
            isa_dashboard_segments(isa)
        })

        min_width <- shiny::reactive({
            value <- input$min_width
            if (!is.numeric(value) || length(value) != 1L || is.na(value)) {
                return(2L)
            }
            max(2L, as.integer(value))
        })

        visible_segments <- shiny::reactive({
            segments <- current_segments()
            if (is.null(segments)) {
                return(NULL)
            }
            ## Same width criterion as `subset.ISA(min_width = )`, applied to the
            ## annotated rows so that `.row_id` stays stable for the export.
            segments <- segments[(segments$loc.end - segments$loc.start + 1L) >= min_width(), , drop = FALSE]
            filter_dashboard_segments(
                segments,
                chr_ids = input$chr_ids,
                pair_ids = input$pair_ids,
                putative_only = input$putative_only
            )
        })

        ## The plot follows the table filters: the thresholded ISA is narrowed to
        ## the selected chromosomes and pairs before `plot.ISA()` sees it.
        plot_isa <- shiny::reactive({
            isa <- isa_thresholded()
            if (is.null(isa)) {
                return(NULL)
            }
            filter_dashboard_isa(
                isa,
                chr_ids = input$chr_ids,
                pair_ids = input$pair_ids,
                min_width = min_width()
            )
        })

        ## Sample name for the BEDPE download; the parent app can pass a reactive
        ## that follows its sample selector.
        current_sample <- function() {
            if (is.null(sample)) {
                return(NULL)
            }
            value <- if (shiny::is.reactive(sample)) sample() else sample
            if (is.null(value) || length(value) != 1L ||
                is.na(value) || !nzchar(value)) {
                return(NULL)
            }
            as.character(value)
        }

        ## Rows the parent asked to preselect. Read inside `renderDT()` with
        ## `isolate()` so a selection callback cannot re-render the table in a
        ## loop; a sample switch invalidates `visible_segments()` instead.
        preselected_rows <- function() {
            if (is.null(selection)) {
                return(integer())
            }
            value <- if (shiny::is.reactive(selection)) selection() else selection
            if (is.null(value)) {
                return(integer())
            }
            suppressWarnings(as.integer(value))
        }

        ## A new ISA object resets the pickers to the full set and syncs the
        ## threshold controls with the metrics the object was segmented with.
        shiny::observeEvent(ISA(), {
            isa <- ISA()
            if (is.null(isa)) {
                return()
            }
            segments <- isa_dashboard_segments(isa)
            chr_ids_all <- names(isa_dashboard_genome(isa))
            initial_chr_ids <- if (is.null(chr_ids)) chr_ids_all else intersect(chr_ids, chr_ids_all)
            pairs <- sort(unique(isa_dashboard_pair_labels(segments)))
            shinyWidgets::updatePickerInput(session, "chr_ids", choices = chr_ids_all, selected = initial_chr_ids)
            shinyWidgets::updatePickerInput(session, "pair_ids", choices = pairs, selected = pairs)
            shiny::updateRadioButtons(
                session, "threshold_type",
                selected = if (identical(isa$meta$metrics, "seg.mean")) "mean" else "median"
            )
            shiny::updateNumericInput(session, "threshold_value", value = isa$meta$metrics_threshold)
            shiny::updateNumericInput(session, "min_width", value = 2L)
        })

        selected_rows <- shiny::reactive({
            visible <- visible_segments()
            selected <- input$segments_rows_selected
            if (is.null(visible) || length(selected) == 0L || nrow(visible) == 0L) {
                return(integer())
            }
            visible$.row_id[selected]
        })

        ## Report the selection to the parent so it can keep a per-sample map for
        ## the combined cohort export.
        shiny::observeEvent(selected_rows(), {
            if (is.null(on_selection)) {
                return()
            }
            on_selection(selected_rows(), current_sample())
        }, ignoreNULL = FALSE)

        output$selection_status <- shiny::renderText({
            if (is.null(ISA())) {
                return("No analysis loaded yet")
            }
            paste(length(selected_rows()), "segment(s) selected for BEDPE export")
        })

        segments_proxy <- DT::dataTableProxy("segments", session = session)

        shiny::observeEvent(input$select_all, {
            visible <- visible_segments()
            if (is.null(visible) || nrow(visible) == 0L) {
                return()
            }
            DT::selectRows(segments_proxy, seq_len(nrow(visible)))
        })

        shiny::observeEvent(input$clear_selection, {
            DT::selectRows(segments_proxy, NULL)
        })

        output$circos <- isa_dashboard_render_circos({
            isa <- plot_isa()
            if (is.null(isa)) {
                return(NULL)
            }
            genome <- isa_dashboard_genome(isa, input$chr_ids)
            if (length(genome) == 0L) {
                return(NULL)
            }
            isa_dashboard_circos(isa, visible_segments(), genome)
        })

        output$isa_plot <- shiny::renderPlot({
            isa <- plot_isa()
            shiny::validate(
                shiny::need(!is.null(ISA()), "Run the analysis to plot"),
                shiny::need(!is.null(isa), "Select at least one chromosome and pair to plot")
            )
            plot_type <- input$plot_type
            if (!identical(plot_type, "mirror") && !identical(plot_type, "diff")) {
                plot_type <- "mirror"
            }
            tryCatch(
                plot(
                    isa,
                    plot.type = plot_type,
                    all_segments = !isTRUE(input$putative_only),
                    max_points = dashboard_max_plot_points
                ),
                error = function(e) {
                    shiny::showNotification(
                        paste("Could not draw the plot:", conditionMessage(e)),
                        type = "warning", duration = 8
                    )
                    shiny::validate(shiny::need(FALSE, "Could not draw the plot"))
                }
            )
        })

        output$segments <- DT::renderDT({
            shiny::validate(
                shiny::need(!is.null(ISA()), "Run the analysis to list segments")
            )
            visible <- visible_segments()
            selected <- shiny::isolate(preselected_rows())
            positions <- if (is.null(visible)) {
                integer()
            } else {
                match(selected, visible$.row_id)
            }
            positions <- positions[!is.na(positions)]
            DT::datatable(
                isa_dashboard_table(visible),
                rownames = FALSE,
                selection = list(
                    mode = "multiple",
                    selected = as.integer(positions)
                ),
                filter = "top",
                options = list(pageLength = 10, scrollX = TRUE)
            )
        })

        output$download_bedpe <- shiny::downloadHandler(
            filename = function() {
                sample_name <- current_sample()
                if (is.null(sample_name)) {
                    "mamamia-selected-introgressions.bedpe"
                } else {
                    safe <- gsub("[^A-Za-z0-9._-]+", "-", sample_name)
                    paste0("mamamia-", safe, "-selected-introgressions.bedpe")
                }
            },
            content = function(file) {
                isa <- ISA()
                rows <- selected_rows()
                if (is.null(isa) || length(rows) == 0L) {
                    stop("Select at least one segment before downloading", call. = FALSE)
                }
                writeToBedpe(
                    selected_isa_for_bedpe(isa, rows),
                    file,
                    sample = current_sample()
                )
            }
        )
    })
}

#' Explore introgression calls in an interactive dashboard
#'
#' Opens a small Shiny application for interactively inspecting an `ISA` object
#' returned by [segments()]. The app offers searchable chromosome and
#' chromosome-pair filters, a circos overview built with `interacCircos`, a
#' mirror or difference plot of the donor-recipient difference, a sortable
#' `DT` segment table, and exports the currently selected table rows with
#' [writeToBedpe()].
#'
#' The circos view stacks the chromosome ideograms, the per-window
#' donor-recipient difference, the called segment spans and links between the
#' paired intervals. It follows the same chromosome, pair, width and putative
#' filters as the plot and the table.
#'
#' Call and threshold filters re-flag the segments on the fly through
#' [subset.ISA], so the call threshold and the minimum segment width can be
#' revisited without re-running segmentation.
#'
#' To run the whole pipeline (input files to dashboard) instead, use
#' [runMaMaMIA()].
#'
#' Dashboard packages are optional MaMaMIA dependencies. Install them with
#' `install.packages(c("shiny", "bslib", "shinyWidgets", "DT"))` and
#' `remotes::install_github("rstudio/gridlayout")` before calling this function.
#'
#' @param ISA An object of class `ISA`, as returned by [segments()].
#' @param chr_ids Optional chromosome IDs to select initially. Defaults to
#'   all chromosomes present in `ISA`.
#' @param host Host passed to [shiny::runApp()].
#' @param port Port passed to [shiny::runApp()]. Use `NULL` to let Shiny choose.
#' @param launch.browser Logical passed to [shiny::runApp()].
#'
#' @return The return value of [shiny::runApp()]. Called for its side effect of
#'   launching the dashboard.
#' @seealso [runMaMaMIA()]
#' @examples
#' \dontrun{
#' data(triticum)
#' rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
#' rca <- correctReadCounts(rca, cores = 4L, verbose = FALSE)
#' isa <- segments(rca, target_pairs = c("OY997261.1" = "CP155614.1"), verbose = FALSE)
#' exploreIntrogressions(isa)
#' }
#' @export
exploreIntrogressions <- function(ISA,
                                  chr_ids = NULL,
                                  host = "127.0.0.1",
                                  port = NULL,
                                  launch.browser = interactive()) {
    stopifnot("Input is not ISA object" = inherits(ISA, "ISA"))
    require_dashboard_packages()

    all_chr_ids <- names(isa_dashboard_genome(ISA))
    if (is.null(chr_ids)) {
        chr_ids <- all_chr_ids
    }
    unknown <- setdiff(chr_ids, all_chr_ids)
    if (length(unknown) > 0L) {
        stop("Unknown chromosome IDs: ", paste(unknown, collapse = ", "), call. = FALSE)
    }

    ui <- bslib::page_sidebar(
        title = "MaMaMIA introgression explorer",
        ## See bam_launcher_ui(): a fillable page turns the result cards into
        ## flex items and collapses the plot output to zero height.
        fillable = FALSE,
        sidebar = bslib::sidebar(explorer_controls_ui("explorer")),
        explorer_results_ui("explorer")
    )

    server <- function(input, output, session) {
        explorer_server("explorer", shiny::reactive(ISA), chr_ids = chr_ids)
    }

    app <- shiny::shinyApp(ui, server)
    args <- list(appDir = app, host = host, launch.browser = launch.browser)
    if (!is.null(port)) {
        args$port <- port
    }
    do.call(shiny::runApp, args)
}
