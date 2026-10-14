## The introgression explorer: the Shiny module with the chromosome and pair
## filters, the segment table and the BEDPE export, plus the
## `exploreIntrogressions()` launcher around it. The disabled BioCircos circos
## view lives in circos.R.

dashboard_packages <- c(
    ## "BioCircos" is temporarily removed together with the circos view; add it
    ## back here and in DESCRIPTION when the circos code is reinstated.
    "bslib", "DT", "shiny", "shinyWidgets"
)

require_dashboard_packages <- function() {
    missing <- dashboard_packages[!vapply(dashboard_packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing) > 0L) {
        stop(
            "Install dashboard dependencies first: install.packages(c(",
            paste(sprintf('"%s"', missing), collapse = ", "),
            "))",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

isa_dashboard_segments <- function(ISA) {
    segments <- ISA$segments
    segments$.row_id <- seq_len(nrow(segments))
    if (nrow(segments) == 0L) {
        return(segments)
    }
    annotate_segment_coords(segments, ISA$out)
}

## Recipients first, then donors in reverse order: the layout the circos view
## expects, with every donor across the 12 o'clock line from its recipient.
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

## Add a clicked call to the selection, or remove it when it was already picked.
## Shared by the explorer table and the parked circos view.
isa_dashboard_toggle_selection <- function(picked, target) {
    if (length(target) != 1L || is.na(target)) {
        return(picked)
    }
    if (target %in% picked) setdiff(picked, target) else c(picked, target)
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
    shiny::tagList(
        ## The circos card is temporarily disabled with the BioCircos dependency:
        # bslib::card(
        #     bslib::card_header("Introgression circos"),
        #     BioCircos::BioCircosOutput(ns("circos"), height = "650px")
        # ),
        bslib::card(
            bslib::card_header("Segments"),
            shiny::p(
                class = "text-muted",
                "The circos view is temporarily unavailable while the BioCircos",
                "dependency is removed. Use the table below to filter calls and",
                "export the selected rows as BEDPE."
            ),
            DT::DTOutput(ns("segments"))
        )
    )
}

## Server for the explorer. `ISA` is a reactive returning an `ISA` object, or
## `NULL` while no analysis has run yet; the explorer stays idle in that case.
## `chr_ids` optionally limits the chromosomes selected on arrival.
explorer_server <- function(id, ISA, chr_ids = NULL) {
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

        output$selection_status <- shiny::renderText({
            if (is.null(ISA())) {
                return("No analysis loaded yet")
            }
            paste(length(selected_rows()), "segment(s) selected for BEDPE export")
        })

        segments_proxy <- DT::dataTableProxy("segments", session = session)

        ## Clicking a circos band toggles the same selection the table drives,
        ## so the plot and the export agree on what is selected. Temporarily
        ## disabled with the BioCircos dependency:
        # shiny::observeEvent(input$circos_clicked_row, {
        #     visible <- visible_segments()
        #     target <- input$circos_clicked_row
        #     if (is.null(visible) || is.null(target) || nrow(visible) == 0L) {
        #         return()
        #     }
        #     picked <- visible$.row_id[input$segments_rows_selected]
        #     keep <- isa_dashboard_toggle_selection(picked, target)
        #     rows <- which(visible$.row_id %in% keep)
        #     DT::selectRows(segments_proxy, if (length(rows) > 0L) rows else NULL)
        # }, ignoreNULL = TRUE)

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

        ## Highlight changes were pushed to the circos widget instead of
        ## re-rendering it, because a rebuild costs thousands of arcs. Disabled
        ## with the BioCircos dependency:
        # shiny::observeEvent(selected_rows(), {
        #     session$sendCustomMessage("mamamia-selection", as.list(selected_rows()))
        # }, ignoreNULL = FALSE)

        ## Rendered with BioCircos::renderBioCircos() once the circos view is
        ## reinstated; temporarily disabled with the BioCircos dependency:
        # output$circos <- BioCircos::renderBioCircos({
        #     isa <- isa_thresholded()
        #     shiny::validate(
        #         shiny::need(!is.null(isa), "Run the analysis to draw the circos")
        #     )
        #     genome <- isa_dashboard_genome(isa, input$chr_ids)
        #     shiny::validate(
        #         shiny::need(length(genome) > 0L, "Select at least one chromosome")
        #     )
        #     isa_dashboard_circos(isa, visible_segments(), genome)
        # })

        output$segments <- DT::renderDT({
            shiny::validate(
                shiny::need(!is.null(ISA()), "Run the analysis to list segments")
            )
            DT::datatable(
                isa_dashboard_table(visible_segments()),
                rownames = FALSE,
                selection = "multiple",
                filter = "top",
                options = list(pageLength = 10, scrollX = TRUE)
            )
        })

        output$download_bedpe <- shiny::downloadHandler(
            filename = function() "mamamia-selected-introgressions.bedpe",
            content = function(file) {
                isa <- ISA()
                rows <- selected_rows()
                if (is.null(isa) || length(rows) == 0L) {
                    stop("Select at least one segment before downloading", call. = FALSE)
                }
                writeToBedpe(selected_isa_for_bedpe(isa, rows), file)
            }
        )
    })
}

#' Explore introgression calls in an interactive dashboard
#'
#' Opens a small Shiny application for interactively inspecting an `ISA` object
#' returned by [segments()]. The app offers searchable chromosome and
#' chromosome-pair filters, a sortable `DT` segment table, and exports the
#' currently selected table rows with [writeToBedpe()].
#'
#' The BioCircos genomic circos plot is temporarily unavailable: the BioCircos
#' dependency has been removed and the circos code is commented out. The
#' explorer shows the segment table and the BEDPE export only.
#'
#' Call and threshold filters re-flag the segments on the fly through
#' [subset.ISA], so the call threshold and the minimum segment width can be
#' revisited without re-running segmentation.
#'
#' To run the whole pipeline (input files to dashboard) instead, use
#' [runMaMaMIA()].
#'
#' Dashboard packages are optional MaMaMIA dependencies. Install them with
#' `install.packages(c("shiny", "bslib", "shinyWidgets", "DT"))`
#' before calling this function.
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
