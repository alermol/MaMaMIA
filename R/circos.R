## The circos view of the introgression explorer, built with interacCircos.
##
## Three layers are stacked inside the genome ideogram:
##   * a CNV track with the per-window donor - recipient difference,
##   * a thin arc track marking every called segment on both chromosomes,
##   * links joining the donor and recipient interval of every call.
##
## The data-preparation helpers are independent of the widget so they can be
## tested without rendering JavaScript.

## High-resolution runs have thousands of windows per chromosome. The CNV
## track is thinned to this many arcs per chromosome; segments and links are
## never thinned.
circos_max_windows_per_chr <- 400L
circos_diff_positive <- "#b2182b"
circos_diff_negative <- "#2166ac"
circos_diff_neutral <- "#d9d9d9"
circos_segment_putative <- "#cb181d"
circos_segment_other <- "#8c8c8c"
circos_link_color <- "rgba(203,24,29,0.35)"
circos_donor_color <- "#737373"
circos_recipient_color <- "#c9c9c9"

## Map chromosome IDs to the names from the registry. Some registries have
## repeated chromosome names across subgenomes, so duplicate labels are
## disambiguated before they reach the widget.
#' @noRd
isa_dashboard_chr_labels <- function(ISA, chr_ids = NULL) {
    out <- ISA$out
    if (is.null(out) || nrow(out) == 0L) {
        return(stats::setNames(character(), character()))
    }
    info <- rbind(
        data.frame(
            chr_id = as.character(out$chr_id.don),
            chr_name = as.character(out$chr_name.don),
            subgenome = as.character(out$subgenome.don),
            stringsAsFactors = FALSE
        ),
        data.frame(
            chr_id = as.character(out$chr_id.rec),
            chr_name = as.character(out$chr_name.rec),
            subgenome = as.character(out$subgenome.rec),
            stringsAsFactors = FALSE
        )
    )
    info <- info[!duplicated(info$chr_id), , drop = FALSE]
    if (!is.null(chr_ids)) {
        position <- match(as.character(chr_ids), info$chr_id)
        info <- info[position, , drop = FALSE]
        info$chr_id <- as.character(chr_ids)
    }
    labels <- info$chr_name
    missing <- is.na(labels) | !nzchar(labels)
    labels[missing] <- info$chr_id[missing]
    duplicated_label <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    if (any(duplicated_label)) {
        suffix <- ifelse(
            !is.na(info$subgenome) & nzchar(info$subgenome),
            paste0(" (", info$subgenome, ")"),
            paste0(" [", info$chr_id, "]")
        )
        labels[duplicated_label] <- paste0(
            labels[duplicated_label], suffix[duplicated_label]
        )
    }
    duplicated_label <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    if (any(duplicated_label)) {
        labels[duplicated_label] <- paste0(
            labels[duplicated_label], " [", info$chr_id[duplicated_label], "]"
        )
    }
    stats::setNames(labels, info$chr_id)
}

#' @noRd
circos_chr_labels <- function(labels, chr_ids) {
    chr_ids <- as.character(chr_ids)
    if (is.null(labels)) {
        return(chr_ids)
    }
    out <- unname(labels[chr_ids])
    missing <- is.na(out) | !nzchar(out)
    out[missing] <- chr_ids[missing]
    out
}

#' @noRd
isa_dashboard_empty_windows <- function() {
    data.frame(
        chr = character(), start = numeric(), stop = numeric(),
        diff = numeric(), stringsAsFactors = FALSE
    )
}

#' Keep the CNV track light: one arc per bin, at the peak absolute difference.
#'
#' @noRd
isa_dashboard_thin_windows <- function(windows,
                                       max_per_chr = circos_max_windows_per_chr) {
    if (nrow(windows) == 0L) {
        return(windows)
    }
    max_per_chr <- as.integer(max_per_chr)[[1L]]
    if (is.na(max_per_chr) || max_per_chr < 1L) {
        max_per_chr <- 1L
    }
    parts <- lapply(split(windows, windows$chr), function(rows) {
        if (nrow(rows) <= max_per_chr) {
            return(rows)
        }
        bin <- floor((seq_len(nrow(rows)) - 1L) * max_per_chr / nrow(rows)) + 1L
        keep <- tapply(
            seq_len(nrow(rows)), bin,
            function(i) i[which.max(abs(rows$diff[i]))]
        )
        rows[sort(as.integer(keep)), , drop = FALSE]
    })
    out <- do.call(rbind, unname(parts))
    rownames(out) <- NULL
    out
}

#' Per-window difference on the donor and the recipient side.
#'
#' @noRd
isa_dashboard_diff_windows <- function(ISA, chr_ids,
                                       max_per_chr = circos_max_windows_per_chr,
                                       labels = NULL) {
    empty <- isa_dashboard_empty_windows()
    out <- ISA$out
    if (is.null(out) || nrow(out) == 0L || length(chr_ids) == 0L) {
        return(empty)
    }
    keep <- out$chr_id.don %in% chr_ids & out$chr_id.rec %in% chr_ids
    out <- out[keep, , drop = FALSE]
    if (nrow(out) == 0L) {
        return(empty)
    }
    windows <- rbind(
        data.frame(
            chr = out$chr_id.don,
            start = out$start.don,
            stop = out$stop.don,
            diff = out$diff,
            stringsAsFactors = FALSE
        ),
        data.frame(
            chr = out$chr_id.rec,
            start = out$start.rec,
            stop = out$stop.rec,
            diff = out$diff,
            stringsAsFactors = FALSE
        )
    )
    windows$chr <- circos_chr_labels(labels, windows$chr)
    isa_dashboard_thin_windows(windows, max_per_chr)
}

#' Symmetric colour limit for the difference track.
#'
#' @noRd
isa_dashboard_diff_limit <- function(windows) {
    if (nrow(windows) == 0L) {
        return(1)
    }
    lim <- as.numeric(stats::quantile(
        abs(windows$diff), 0.99,
        na.rm = TRUE
    ))
    if (!is.finite(lim) || lim <= 0) 1 else lim
}

#' Convert the difference windows to the interacCircos CNV data contract.
#'
#' @noRd
isa_dashboard_circos_cnv <- function(windows, limit) {
    if (nrow(windows) == 0L) {
        return(data.frame(
            chr = character(), start = numeric(), end = numeric(),
            value = numeric(), link = character(), color = character(),
            html = character(), stringsAsFactors = FALSE
        ))
    }
    value <- windows$diff
    value[!is.finite(value)] <- 0
    value <- pmax(pmin(value, limit), -limit)
    color <- ifelse(
        value > 0, circos_diff_positive,
        ifelse(value < 0, circos_diff_negative, circos_diff_neutral)
    )
    data.frame(
        chr = windows$chr,
        start = windows$start,
        end = windows$stop,
        value = value,
        link = NA_character_,
        color = color,
        html = NA_character_,
        stringsAsFactors = FALSE
    )
}

#' @noRd
isa_dashboard_circos_arcs <- function(segments, labels = NULL) {
    if (is.null(segments) || nrow(segments) == 0L) {
        return(data.frame(
            chr = character(), start = numeric(), end = numeric(),
            color = character(), des = character(), link = character(),
            html = character(), stringsAsFactors = FALSE
        ))
    }
    chr_don <- circos_chr_labels(labels, segments$chr_id.don)
    chr_rec <- circos_chr_labels(labels, segments$chr_id.rec)
    size <- segments$stop.don - segments$start.don
    label <- sprintf(
        "#%d %s -> %s | %.2f Mb | %s",
        segments$.row_id,
        chr_don,
        chr_rec,
        size / 1e6,
        ifelse(
            segments$putative_introgression,
            "putative introgression",
            "not flagged"
        )
    )
    color <- ifelse(
        segments$putative_introgression,
        circos_segment_putative,
        circos_segment_other
    )
    data.frame(
        chr = c(chr_don, chr_rec),
        start = c(segments$start.don, segments$start.rec),
        end = c(segments$stop.don, segments$stop.rec),
        color = rep(color, 2L),
        des = rep(label, 2L),
        link = NA_character_,
        html = NA_character_,
        stringsAsFactors = FALSE
    )
}

#' @noRd
isa_dashboard_circos_links <- function(segments, labels = NULL) {
    if (is.null(segments) || nrow(segments) == 0L) {
        return(data.frame(
            g1chr = character(), g1start = numeric(), g1end = numeric(),
            g1name = character(), g2chr = character(), g2start = numeric(),
            g2end = numeric(), g2name = character(), fusion = character(),
            link = character(), html = character(), stringsAsFactors = FALSE
        ))
    }
    chr_don <- circos_chr_labels(labels, segments$chr_id.don)
    chr_rec <- circos_chr_labels(labels, segments$chr_id.rec)
    size <- segments$stop.don - segments$start.don
    label <- sprintf(
        "#%d %s -> %s | %.2f Mb | %s",
        segments$.row_id,
        chr_don,
        chr_rec,
        size / 1e6,
        ifelse(
            segments$putative_introgression,
            "putative introgression",
            "not flagged"
        )
    )
    data.frame(
        g1chr = chr_don,
        g1start = segments$start.don,
        g1end = segments$stop.don,
        g1name = sprintf(
            "%s:%.2f-%.2f Mb",
            chr_don,
            segments$start.don / 1e6,
            segments$stop.don / 1e6
        ),
        g2chr = chr_rec,
        g2start = segments$start.rec,
        g2end = segments$stop.rec,
        g2name = sprintf(
            "%s:%.2f-%.2f Mb",
            chr_rec,
            segments$start.rec / 1e6,
            segments$stop.rec / 1e6
        ),
        fusion = label,
        link = NA_character_,
        html = NA_character_,
        stringsAsFactors = FALSE
    )
}

## Remove tooltips left behind by a previous widget render and hide the current
## tooltips when the pointer leaves the plot or the page scrolls.
#' @noRd
circos_hover_cleanup <- function() {
    htmlwidgets::JS(
        "function(el, x) {",
        "  var selectors = [",
        "    '.NGCircosCNVTooltip',",
        "    '.NGCircosARCTooltip',",
        "    '.NGCircosLINKTooltip'",
        "  ];",
        "  selectors.forEach(function(selector) {",
        "    var nodes = document.querySelectorAll(selector);",
        "    for (var i = 0; i < nodes.length - 1; i++) {",
        "      if (nodes[i].parentNode) nodes[i].parentNode.removeChild(nodes[i]);",
        "    }",
        "  });",
        "  var hide = function() {",
        "    selectors.forEach(function(selector) {",
        "      document.querySelectorAll(selector).forEach(function(node) {",
        "        node.style.opacity = '0';",
        "      });",
        "    });",
        "  };",
        "  if (!el.__mamamiaCircosHover) {",
        "    el.__mamamiaCircosHover = hide;",
        "    el.addEventListener('mouseleave', hide);",
        "  }",
        "  if (!document.__mamamiaCircosHover) {",
        "    document.__mamamiaCircosHover = hide;",
        "    document.addEventListener('scroll', hide, true);",
        "  }",
        "}"
    )
}

#' Build the interacCircos widget for one filtered ISA
#'
#' @param ISA A filtered `ISA` object (the explorer passes `plot_isa()`).
#' @param segments The annotated and filtered segment table (the explorer
#'   passes `visible_segments()`).
#' @param genome Named list of chromosome lengths, typically from
#'   [`isa_dashboard_genome()`].
#'
#' @return An `interacCircos` htmlwidget.
#' @noRd
isa_dashboard_circos <- function(ISA, segments, genome) {
    genome <- unlist(genome, use.names = TRUE)
    genome <- genome[!is.na(genome) & genome > 0]
    if (length(genome) == 0L) {
        stop("No chromosome lengths available for the circos view", call. = FALSE)
    }
    chr_ids <- names(genome)
    labels <- isa_dashboard_chr_labels(ISA, chr_ids)
    display_ids <- circos_chr_labels(labels, chr_ids)
    genome <- stats::setNames(unname(genome), display_ids)
    windows <- isa_dashboard_diff_windows(ISA, chr_ids, labels = labels)
    limit <- isa_dashboard_diff_limit(windows)
    cnv <- isa_dashboard_circos_cnv(windows, limit)
    arcs <- isa_dashboard_circos_arcs(segments, labels = labels)
    links <- isa_dashboard_circos_links(segments, labels = labels)

    fill <- ifelse(
        chr_ids %in% ISA$meta$don_chr_ids,
        circos_donor_color, circos_recipient_color
    )

    modules <- interacCircos::CircosModuleList()
    if (nrow(cnv) > 0L) {
        modules <- modules +
            interacCircos::CircosBackground(
                "cnv_background",
                minRadius = 150, maxRadius = 190,
                fillColors = "#F4F4F4", borderColors = "#D9D9D9",
                axisShow = TRUE, axisColor = "#B8B8B8",
                axisNum = 4, axisWidth = 0.3, borderSize = 0.3
            ) +
            interacCircos::CircosCnv(
                "difference",
                data = cnv,
                minRadius = 150, maxRadius = 190, width = 6,
                color = circos_diff_neutral,
                ValueAxisManualScale = TRUE,
                ValueAxisMaxScale = limit,
                ValueAxisMinScale = -limit,
                strokeColor = "none", strokeWidth = 0,
                opacity = 0.9
            )
    }
    if (nrow(arcs) > 0L) {
        modules <- modules +
            interacCircos::CircosArc(
                "segment_spans",
                data = arcs,
                innerRadius = 198, outerRadius = 206,
                opacity = 0.9
            )
    }
    if (nrow(links) > 0L) {
        modules <- modules +
            interacCircos::CircosLink(
                "introgressions",
                data = links,
                radius = 125, fillColor = circos_link_color,
                width = 3, type = "Q",
                displayLinkAxis = FALSE, displayLinkLabel = FALSE
            )
    }

    widget <- interacCircos::Circos(
        moduleList = modules,
        genome = as.list(genome),
        genomeFillColor = fill,
        width = "100%",
        height = 650,
        innerRadius = 216,
        outerRadius = 240,
        chrPad = 0.03,
        displayGenomeBorder = TRUE,
        genomeBorderColor = "#58636F",
        genomeBorderSize = 0.4,
        genomeTicksDisplay = FALSE,
        genomeLabelDisplay = TRUE,
        genomeLabelTextSize = "8pt",
        genomeLabelTextColor = "#2f3437",
        CNVMouseOverDisplay = TRUE,
        CNVMouseOutDisplay = TRUE,
        CNVMouseOutColor = "none",
        CNVMouseOutArcOpacity = "none",
        CNVMouseOutArcStrokeColor = "none",
        CNVMouseOutArcStrokeWidth = "none",
        ARCMouseOverDisplay = TRUE,
        ARCMouseOutDisplay = TRUE,
        ARCMouseOutColor = "none",
        ARCMouseOutArcOpacity = "none",
        ARCMouseOutArcStrokeColor = "none",
        ARCMouseOutArcStrokeWidth = "none",
        LINKMouseOverDisplay = TRUE,
        LINKMouseOutDisplay = TRUE,
        LINKMouseOutOpacity = "none",
        LINKMouseOutStrokeColor = "none",
        LINKMouseOutStrokeWidth = "none",
        zoom = TRUE
    )
    htmlwidgets::onRender(widget, circos_hover_cleanup())
}

#' Shiny output for the interacCircos widget
#'
#' interacCircos does not ship Shiny bindings, so pair the standard htmlwidgets
#' output and render functions with its widget name.
#'
#' @noRd
isa_dashboard_circos_output <- function(outputId, width = "100%",
                                        height = "650px") {
    htmlwidgets::shinyWidgetOutput(
        outputId, "interacCircos", width, height,
        package = "interacCircos"
    )
}

#' @noRd
isa_dashboard_render_circos <- function(expr, env = parent.frame()) {
    expr <- substitute(expr)
    htmlwidgets::shinyRenderWidget(
        expr, isa_dashboard_circos_output, env, quoted = TRUE
    )
}
