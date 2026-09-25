dashboard_packages <- c("BioCircos", "bslib", "DT", "htmltools", "htmlwidgets", "shiny", "shinyWidgets")

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

## Recipients first, then donors in reverse order, so every donor sits across the
## 12 o'clock line from its recipient instead of opposite it.
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

## For the mirrored donor half a coordinate c becomes length - c, so a feature
## keeps its distance from the shared 12 o'clock edge. Chromosomes absent from
## `don_len` are left alone, which is how the recipient side stays untouched.
isa_dashboard_mirror <- function(chr, start, stop, don_len) {
    len <- unname(don_len[chr])
    list(
        start = ifelse(is.na(len), start, len - stop),
        stop = ifelse(is.na(len), stop, len - start)
    )
}

isa_dashboard_mirror_donor <- function(data, don_len) {
    if (nrow(data) == 0L) {
        return(data)
    }
    mirrored <- isa_dashboard_mirror(data$chr, data$start, data$stop, don_len)
    data$start <- mirrored$start
    data$stop <- mirrored$stop
    data
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

## Shared by the table and the plot: a click adds a call to the selection, or
## removes it when it was already picked.
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

## The circos mirrors tools/illustrations/render-illustrations.R::fig_circos(),
## which stacks three concentric layers: the chromosome ideogram, the per-window
## donor - recipient coverage difference, then the ribbons for the calls.
##
## BioCircos radial units: every track radius is scaled by 0.7 * maxRadius and
## the built-in ideogram band spans 0.70-0.80 of maxRadius, that is 1.000-1.143
## in track units, so every layer added here must stay below 1.0. CNV arcs place
## a constant-thickness mark at minRadius + (value - range[1]) / (range[0] -
## range[1]) * (maxRadius - minRadius), while CNV `width` is in pixels.
## Links are single quadratic curves and gene1Starts/gene1Ends only supply the
## midpoint, so segment size is shown by splitting calls into size classes.
dashboard_palette <- c(
    donor = "#737373",
    recipient = "#c9c9c9",
    diff_high = "#b2182b",
    diff_low = "#2166ac",
    ribbon = "#cb181d"
)
ribbon_selected_color <- "#D95F02"
label_on_donor <- "#ffffff"
label_on_recipient <- "#2f3437"
diff_zero_radius <- 0.82
diff_half_range <- 0.11
diff_mark_width_px <- 4
## The ribbons start at the difference ring, not out at the ideogram, so the
## chunk marks sit at the same radius as the link radius and read as their roots.
span_mark_radius <- 0.70
span_mark_width_px <- 8
link_size_widths <- c(1.5, 3, 5, 7.5)
link_max_radius <- 0.70
diff_max_windows_per_chr <- 400L

isa_dashboard_size_class <- function(size, n_class) {
    usable <- size
    usable[!is.finite(usable) | usable <= 0] <- NA_real_
    breaks <- unique(stats::quantile(usable, probs = seq(0, 1, length.out = n_class + 1L), na.rm = TRUE))
    if (length(breaks) < 2L) {
        return(rep(1L, length(size)))
    }
    class <- cut(usable, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    class[is.na(class)] <- 1L
    pmin(as.integer(class), n_class)
}

## Keep the diff ring light: BioCircos serialises one arc per window, so very
## long chromosomes are thinned to the peak |diff| window of each bin.
isa_dashboard_thin_windows <- function(windows, max_per_chr) {
    if (nrow(windows) <= max_per_chr) {
        return(windows)
    }
    parts <- lapply(split(windows, windows$chr), function(rows) {
        if (nrow(rows) <= max_per_chr) {
            return(rows)
        }
        bin <- cut(seq_len(nrow(rows)), breaks = max_per_chr, labels = FALSE)
        keep <- tapply(seq_len(nrow(rows)), bin, function(i) i[which.max(abs(rows$diff[i]))])
        rows[sort(as.integer(keep)), , drop = FALSE]
    })
    do.call(rbind, unname(parts))
}

isa_dashboard_diff_windows <- function(ISA, chr_ids) {
    out <- ISA$out
    out <- out[out$chr_id.don %in% chr_ids & out$chr_id.rec %in% chr_ids, , drop = FALSE]
    if (nrow(out) == 0L) {
        return(data.frame(chr = character(), start = numeric(), stop = numeric(), diff = numeric()))
    }
    ## fig_circos() draws the difference on both the donor and the recipient side.
    windows <- rbind(
        data.frame(chr = out$chr_id.don, start = out$start.don, stop = out$stop.don, diff = out$diff),
        data.frame(chr = out$chr_id.rec, start = out$start.rec, stop = out$stop.rec, diff = out$diff)
    )
    isa_dashboard_thin_windows(windows, diff_max_windows_per_chr)
}

isa_dashboard_diff_limit <- function(windows) {
    if (nrow(windows) == 0L) {
        return(1)
    }
    lim <- as.numeric(stats::quantile(abs(windows$diff), 0.99, na.rm = TRUE))
    if (!is.finite(lim) || lim <= 0) 1 else lim
}

isa_dashboard_diff_track <- function(name, rows, sign, lim, color) {
    if (sign > 0) {
        min_radius <- diff_zero_radius
        max_radius <- diff_zero_radius + diff_half_range
        range <- c(lim, 0)
    } else {
        min_radius <- diff_zero_radius - diff_half_range
        max_radius <- diff_zero_radius
        range <- c(0, lim)
    }
    BioCircos::BioCircosCNVTrack(
        name,
        rows$chr, rows$start, rows$stop,
        values = abs(rows$diff),
        minRadius = min_radius,
        maxRadius = max_radius,
        width = diff_mark_width_px,
        color = color,
        range = range
    )
}

isa_dashboard_diff_tracks <- function(windows, lim) {
    tracks <- BioCircos::BioCircosTracklist()
    if (nrow(windows) == 0L) {
        return(tracks)
    }
    diff <- windows$diff
    diff[!is.finite(diff)] <- 0
    windows$diff <- pmax(pmin(diff, lim), -lim)
    sides <- list(
        list(name = "diff_positive", sign = 1, keep = windows$diff > 0, color = dashboard_palette[["diff_high"]]),
        list(name = "diff_negative", sign = -1, keep = windows$diff < 0, color = dashboard_palette[["diff_low"]])
    )
    for (side in sides) {
        rows <- windows[side$keep, , drop = FALSE]
        if (nrow(rows) == 0L) {
            next
        }
        tracks <- tracks + isa_dashboard_diff_track(side$name, rows, side$sign, lim, side$color)
    }
    tracks
}

## Constant-thickness marks spanning each called segment on both chromosomes, so
## the chunk extent is visible even though a BioCircos link collapses to a curve.
isa_dashboard_span_track <- function(segments) {
    if (nrow(segments) == 0L) {
        return(BioCircos::BioCircosTracklist())
    }
    spans <- rbind(
        data.frame(chr = segments$chr_id.don, start = segments$start.don, stop = segments$stop.don),
        data.frame(chr = segments$chr_id.rec, start = segments$start.rec, stop = segments$stop.rec)
    )
    BioCircos::BioCircosCNVTrack(
        "introgression_spans",
        spans$chr, spans$start, spans$stop,
        values = rep(1, nrow(spans)),
        minRadius = span_mark_radius,
        maxRadius = span_mark_radius,
        width = span_mark_width_px,
        color = dashboard_palette[["ribbon"]],
        range = c(0, 1)
    )
}

## Band hover text comes from the link `labels`; the hit area is a separate
## near-transparent wide track so a band is easy to point at and click.
link_hit_width_px <- 22
link_hit_color <- "rgba(255,255,255,0.01)"

isa_dashboard_mb <- function(bp) {
    mb <- bp / 1e6
    ifelse(abs(mb) < 10, formatC(mb, format = "f", digits = 2), formatC(mb, format = "f", digits = 1))
}

isa_dashboard_band_label <- function(segments) {
    paste0(
        isa_dashboard_mb(segments$size), " Mb",
        " | donor ", segments$chr_id.don,
        " -> recipient ", segments$chr_id.rec,
        " | segment ", segments$.row_id
    )
}

## One hover string per chromosome, listing the calls that touch it.
isa_dashboard_chromosome_summary <- function(segments, chr_ids) {
    summary <- stats::setNames(vector("list", length(chr_ids)), chr_ids)
    for (chr in chr_ids) {
        rows <- segments[segments$chr_id.don == chr | segments$chr_id.rec == chr, , drop = FALSE]
        if (nrow(rows) == 0L) {
            summary[[chr]] <- paste0(chr, ": no introgressions")
            next
        }
        entries <- vapply(seq_len(nrow(rows)), function(i) {
            row <- rows[i, ]
            as_donor <- row$chr_id.don == chr
            from <- if (as_donor) row$start.don else row$start.rec
            to <- if (as_donor) row$stop.don else row$stop.rec
            sprintf(
                "#%d %s %s-%s Mb (%s Mb)",
                row$.row_id,
                if (as_donor) "donated to" else "received from",
                isa_dashboard_mb(from), isa_dashboard_mb(to), isa_dashboard_mb(to - from)
            )
        }, character(1))
        summary[[chr]] <- paste0(
            chr, ": ", nrow(rows), " introgression(s)",
            "<br>", paste(entries, collapse = "<br>")
        )
    }
    summary
}

## The interactive layer lives in inst/mamamia/ and ships as a widget dependency,
## which is the only channel honoured by a Shiny render call: prependContent()
## warns "Ignoring prepended content" there, and a bare htmlDependency() passed
## to it is dropped even in a static render. The stylesheet matters on its own,
## since it is what keeps the plot from swallowing clicks meant for the table.
isa_dashboard_circos_assets <- function() {
    htmltools::htmlDependency(
        "mamamia-circos",
        "0.1.0",
        src = "mamamia",
        script = "mamamia-circos.js",
        stylesheet = "mamamia-circos.css",
        package = "MaMaMIA"
    )
}

isa_dashboard_circos_hook <- function() {
    "function(el, x) { window.mamamiaCircosEnhance(el, x); }"
}

isa_dashboard_link_track <- function(name, rows, width, color) {
    BioCircos::BioCircosLinkTrack(
        name,
        rows$chr_id.don, rows$start.don, rows$stop.don,
        rows$chr_id.rec, rows$start.rec, rows$stop.rec,
        labels = rows$.label,
        color = color,
        maxRadius = link_max_radius,
        width = paste0(width, "px"),
        displayAxis = FALSE,
        displayLabel = FALSE
    )
}

## Selection is deliberately not part of the widget: a rebuild costs thousands of
## arcs, so the plot is drawn once and the picked bands are restyled in place.
isa_dashboard_circos <- function(ISA, segments, genome) {
    chr_ids <- names(genome)
    lengths <- unlist(genome)
    ## The donor half mirrors the recipient one, so its features are handed over
    ## as length - coordinate and its sectors run in reverse order.
    don_len <- lengths[names(lengths) %in% ISA$meta$don_chr_ids]
    tracks <- BioCircos::BioCircosTracklist()

    windows <- isa_dashboard_mirror_donor(isa_dashboard_diff_windows(ISA, chr_ids), don_len)
    tracks <- tracks + isa_dashboard_diff_tracks(windows, isa_dashboard_diff_limit(windows))

    ribbon_rows <- integer()
    span_rows <- integer()
    segments <- segments[
        segments$chr_id.don %in% chr_ids & segments$chr_id.rec %in% chr_ids,
        ,
        drop = FALSE
    ]
    if (nrow(segments) > 0L) {
        segments$size <- segments$stop.don - segments$start.don
        segments$.label <- isa_dashboard_band_label(segments)
        size_class <- isa_dashboard_size_class(segments$size, length(link_size_widths))

        mirrored <- isa_dashboard_mirror(segments$chr_id.don, segments$start.don, segments$stop.don, don_len)
        segments$start.don <- mirrored$start
        segments$stop.don <- mirrored$stop

        tracks <- tracks + isa_dashboard_span_track(segments)

        ## Row ids in the order the drawn elements appear, so the client can map a
        ## clicked or selected row onto its ribbon and its two chromosome marks.
        ribbon_rows <- unlist(lapply(seq_along(link_size_widths), function(k) {
            segments$.row_id[size_class == k]
        }), use.names = FALSE)
        span_rows <- c(segments$.row_id, segments$.row_id)

        for (k in seq_along(link_size_widths)) {
            rows <- segments[size_class == k, , drop = FALSE]
            if (nrow(rows) == 0L) {
                next
            }
            tracks <- tracks + isa_dashboard_link_track(
                paste0("introgression_size", k),
                rows,
                link_size_widths[k],
                paste0(dashboard_palette[["ribbon"]], "55")
            )
        }

        ## Wide, near-transparent copy of every band, on top, so the whole ribbon
        ## responds to hover and click instead of just the drawn stroke.
        tracks <- tracks + isa_dashboard_link_track(
            "introgression_hit_area",
            segments,
            link_hit_width_px,
            link_hit_color
        )
    }

    widget <- BioCircos::BioCircos(
        tracks,
        genome = genome,
        genomeFillColor = ifelse(
            chr_ids %in% ISA$meta$don_chr_ids,
            dashboard_palette[["donor"]],
            dashboard_palette[["recipient"]]
        ),
        chrPad = 0.03,
        displayGenomeBorder = TRUE,
        genomeBorderColor = "#58636F",
        genomeBorderSize = 0.4,
        genomeTicksDisplay = FALSE,
        genomeLabelTextSize = "5pt",
        ## Dark ink by default, so the names stay readable on the light band even
        ## if the enhancement script never runs; that script whitens the donor
        ## labels and pulls every name onto its own segment.
        genomeLabelTextColor = "#2f3437",
        genomeLabelDy = 0,
        genomeLabelOrientation = 0,
        LINKMouseOverStrokeColor = ribbon_selected_color,
        LINKMouseOverStrokeWidth = 5,
        LINKMouseOverOpacity = 0.95,
        LINKMouseOverTooltipsHtml01 = "",
        LINKMouseOverTooltipsHtml02 = ""
    )
    widget$x$chromosome_summary <- isa_dashboard_chromosome_summary(segments, chr_ids)
    widget$x$chromosome_ids <- chr_ids
    widget$x$donor_chr_ids <- ISA$meta$don_chr_ids
    widget$x$hit_area_rows <- if (nrow(segments) > 0L) as.integer(segments$.row_id) else integer()
    widget$x$ribbon_rows <- as.integer(ribbon_rows)
    widget$x$span_rows <- as.integer(span_rows)
    widget$dependencies <- c(widget$dependencies, list(isa_dashboard_circos_assets()))
    htmlwidgets::onRender(widget, isa_dashboard_circos_hook())
}

#' Explore introgression calls in an interactive dashboard
#'
#' Opens a small Shiny application for interactively inspecting an `ISA` object
#' returned by [segments()]. The app uses a BioCircos genomic circos plot,
#' searchable chromosome and chromosome-pair filters, a sortable `DT` segment
#' table, and exports the currently selected table rows with [writeToBedpe()].
#'
#' Dashboard packages are optional MaMaMIA dependencies. Install them with
#' `install.packages(c("shiny", "bslib", "shinyWidgets", "DT", "BioCircos"))`
#' before calling this function.
#'
#' @param ISA An object of class `ISA`, as returned by [segments()].
#' @param chr_ids Optional chromosome IDs to show initially. Defaults to all
#'   chromosomes present in `ISA`.
#' @param host Host passed to [shiny::runApp()].
#' @param port Port passed to [shiny::runApp()]. Use `NULL` to let Shiny choose.
#' @param launch.browser Logical passed to [shiny::runApp()].
#'
#' @return The return value of [shiny::runApp()]. Called for its side effect of
#'   launching the dashboard.
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

    all_segments <- isa_dashboard_segments(ISA)
    pair_choices <- sort(unique(isa_dashboard_pair_labels(all_segments)))

    ui <- bslib::page_sidebar(
        title = "MaMaMIA introgression explorer",
        sidebar = bslib::sidebar(
            shinyWidgets::pickerInput(
                "chr_ids", "Chromosomes",
                choices = all_chr_ids, selected = chr_ids,
                multiple = TRUE,
                options = list(`actions-box` = TRUE, `live-search` = TRUE)
            ),
            shinyWidgets::pickerInput(
                "pair_ids", "Chromosome pairs",
                choices = pair_choices, selected = pair_choices,
                multiple = TRUE,
                options = list(`actions-box` = TRUE, `live-search` = TRUE)
            ),
            shiny::checkboxInput("putative_only", "Show putative introgressions only", TRUE),
            shiny::textOutput("selection_status"),
            shiny::downloadButton("download_bedpe", "Download selected BEDPE")
        ),
        bslib::card(
            bslib::card_header("Introgression circos"),
            BioCircos::BioCircosOutput("circos", height = "650px")
        ),
        bslib::card(
            bslib::card_header("Segments"),
            DT::DTOutput("segments")
        )
    )

    server <- function(input, output, session) {
        visible_segments <- shiny::reactive({
            filter_dashboard_segments(
                all_segments,
                chr_ids = input$chr_ids,
                pair_ids = input$pair_ids,
                putative_only = input$putative_only
            )
        })

        selected_rows <- shiny::reactive({
            visible <- visible_segments()
            selected <- input$segments_rows_selected
            if (length(selected) == 0L || nrow(visible) == 0L) {
                return(integer())
            }
            visible$.row_id[selected]
        })

        output$selection_status <- shiny::renderText({
            paste(length(selected_rows()), "segment(s) selected for BEDPE export")
        })

        ## Clicking a band toggles the same selection the table drives, so the
        ## plot and the export agree on what is selected.
        segments_proxy <- DT::dataTableProxy("segments")
        shiny::observeEvent(input$circos_clicked_row, {
            visible <- visible_segments()
            target <- input$circos_clicked_row
            if (is.null(target) || nrow(visible) == 0L) {
                return()
            }
            picked <- visible$.row_id[input$segments_rows_selected]
            keep <- isa_dashboard_toggle_selection(picked, target)
            rows <- which(visible$.row_id %in% keep)
            DT::selectRows(segments_proxy, if (length(rows) > 0L) rows else NULL)
        }, ignoreNULL = TRUE)

        ## Highlight changes are pushed to the plot instead of re-rendering it:
        ## rebuilding the widget costs thousands of arcs on every click.
        shiny::observeEvent(selected_rows(), {
            session$sendCustomMessage("mamamia-selection", as.list(selected_rows()))
        }, ignoreNULL = FALSE)

        output$circos <- BioCircos::renderBioCircos({
            genome <- isa_dashboard_genome(ISA, input$chr_ids)
            shiny::validate(
                shiny::need(length(genome) > 0L, "Select at least one chromosome")
            )
            isa_dashboard_circos(ISA, visible_segments(), genome)
        })

        output$segments <- DT::renderDT({
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
                writeToBedpe(selected_isa_for_bedpe(ISA, selected_rows()), file)
            }
        )
    }

    app <- shiny::shinyApp(ui, server)
    args <- list(appDir = app, host = host, launch.browser = launch.browser)
    if (!is.null(port)) {
        args$port <- port
    }
    do.call(shiny::runApp, args)
}
