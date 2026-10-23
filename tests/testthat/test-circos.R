## The interacCircos overview: data preparation and widget construction.

test_that("diff windows cover both chromosomes and are thinned per chromosome", {
    isa <- seg_step()
    windows <- isa_dashboard_diff_windows(isa, c("d1", "r1"))

    expect_setequal(names(windows), c("chr", "start", "stop", "diff"))
    expect_setequal(unique(windows$chr), c("d1", "r1"))
    expect_true(all(is.finite(windows$diff)))
    expect_true(all(windows$start < windows$stop))

    thinned <- isa_dashboard_thin_windows(windows, 10L)
    expect_equal(nrow(thinned), 20L)
    expect_true(mean(abs(thinned$diff)) >= mean(abs(windows$diff)))

    empty <- isa_dashboard_diff_windows(isa, "nope")
    expect_equal(nrow(empty), 0L)
    expect_equal(isa_dashboard_diff_limit(empty), 1)
})

test_that("CNV data are clipped and coloured by the sign of the difference", {
    isa <- seg_step()
    windows <- isa_dashboard_diff_windows(isa, c("d1", "r1"))
    limit <- isa_dashboard_diff_limit(windows)
    cnv <- isa_dashboard_circos_cnv(windows, limit)

    expect_equal(nrow(cnv), nrow(windows))
    expect_setequal(names(cnv), c("chr", "start", "end", "value", "link", "color", "html"))
    expect_true(all(cnv$value >= -limit & cnv$value <= limit))
    expect_true(all(cnv$color %in% c(
        circos_diff_positive, circos_diff_negative, circos_diff_neutral
    )))
    expect_identical(unique(cnv$color[cnv$value > 0]), circos_diff_positive)
    expect_identical(unique(cnv$color[cnv$value < 0]), circos_diff_negative)

    zero <- isa_dashboard_circos_cnv(windows[0, , drop = FALSE], limit)
    expect_equal(nrow(zero), 0L)
})

test_that("segment spans and links follow the annotated segment table", {
    isa <- seg_step()
    segments <- isa_dashboard_segments(isa)
    arcs <- isa_dashboard_circos_arcs(segments)
    links <- isa_dashboard_circos_links(segments)

    expect_equal(nrow(arcs), 2L * nrow(segments))
    expect_setequal(names(arcs), c("chr", "start", "end", "color", "des", "link", "html"))
    expect_setequal(unique(arcs$chr), c("d1", "r1"))

    expect_equal(nrow(links), nrow(segments))
    expect_setequal(names(links), c(
        "g1chr", "g1start", "g1end", "g1name", "g2chr", "g2start", "g2end",
        "g2name", "fusion", "link", "html"
    ))
    expect_true(all(links$g1chr == "d1"))
    expect_true(all(links$g2chr == "r1"))
    expect_true(all(nzchar(links$g1name) & nzchar(links$g2name) & nzchar(links$fusion)))

    empty <- segments[0, , drop = FALSE]
    expect_equal(nrow(isa_dashboard_circos_arcs(empty)), 0L)
    expect_equal(nrow(isa_dashboard_circos_links(empty)), 0L)
})

test_that("the interacCircos widget stacks difference, spans and links", {
    testthat::skip_if_not_installed("interacCircos")

    isa <- seg_step()
    segments <- isa_dashboard_segments(isa)
    genome <- isa_dashboard_genome(isa)
    labels <- isa_dashboard_chr_labels(isa, names(genome))
    widget <- isa_dashboard_circos(isa, segments, genome)

    expect_s3_class(widget, "interacCircos")
    expect_s3_class(widget, "htmlwidget")

    modules <- widget$x$moduleList_NGCircos
    names <- vapply(modules, function(module) module[[1L]], character(1))
    expect_true("CNV_difference" %in% names)
    expect_true("ARC_segment_spans" %in% names)
    expect_true("LINK_introgressions" %in% names)
    expect_identical(unname(labels[names(genome)]), c("1A", "1At"))
    expect_equal(widget$x$genome, as.list(stats::setNames(genome, labels[names(genome)])))
    expect_equal(
        widget$x$genomeFillColor[names(genome) %in% isa$meta$don_chr_ids],
        rep(circos_donor_color, sum(names(genome) %in% isa$meta$don_chr_ids))
    )

    ## The plot is zoomable and the hover tooltips are hidden on mouse-out.
    expect_true(widget$x$zoom)
    expect_true(widget$x$CNVMouseOutDisplay)
    expect_true(widget$x$ARCMouseOutDisplay)
    expect_true(widget$x$LINKMouseOutDisplay)
    expect_length(widget$jsHooks$render, 1L)
    hover_js <- paste(widget$jsHooks$render[[1L]]$code, collapse = "\n")
    expect_match(hover_js, "mouseleave", fixed = TRUE)
    expect_match(hover_js, "NGCircosCNVTooltip", fixed = TRUE)

    cnv_data <- modules[[which(names == "CNV_difference")]][[3L]]
    arc_data <- modules[[which(names == "ARC_segment_spans")]][[3L]]
    link_data <- modules[[which(names == "LINK_introgressions")]][[3L]]
    expect_setequal(
        unique(vapply(cnv_data, function(row) row$chr, character(1))),
        c("1A", "1At")
    )
    expect_setequal(
        unique(vapply(arc_data, function(row) row$chr, character(1))),
        c("1A", "1At")
    )
    expect_true(all(vapply(
        link_data, function(row) row$g1chr == "1At", logical(1)
    )))
    expect_true(all(vapply(
        link_data, function(row) row$g2chr == "1A", logical(1)
    )))

    ## Without calls the CNV track still draws the genome and the difference.
    bare <- isa_dashboard_circos(isa, segments[0, , drop = FALSE], genome)
    bare_modules <- vapply(
        bare$x$moduleList_NGCircos,
        function(module) module[[1L]], character(1)
    )
    expect_true("CNV_difference" %in% bare_modules)
    expect_false("LINK_introgressions" %in% bare_modules)
    expect_false("ARC_segment_spans" %in% bare_modules)
})

test_that("circos data use chromosome names instead of registry IDs", {
    isa <- seg_step()
    genome <- isa_dashboard_genome(isa)
    labels <- isa_dashboard_chr_labels(isa, names(genome))

    expect_identical(unname(labels[c("d1", "r1")]), c("1At", "1A"))

    windows <- isa_dashboard_diff_windows(
        isa, names(genome), labels = labels
    )
    expect_setequal(unique(windows$chr), c("1A", "1At"))

    segments <- isa_dashboard_segments(isa)
    arcs <- isa_dashboard_circos_arcs(segments, labels = labels)
    links <- isa_dashboard_circos_links(segments, labels = labels)

    expect_identical(
        unique(isa_dashboard_pair_labels(segments, labels = labels)),
        "1At -> 1A"
    )
    expect_setequal(unique(arcs$chr), c("1A", "1At"))
    expect_false(any(grepl("\\bd1\\b|\\br1\\b", arcs$des)))
    expect_true(all(links$g1chr == "1At"))
    expect_true(all(links$g2chr == "1A"))
    expect_false(any(grepl("\\bd1\\b|\\br1\\b", links$fusion)))
    expect_true(all(grepl("^1At:", links$g1name)))
    expect_true(all(grepl("^1A:", links$g2name)))
})

test_that("duplicate chromosome names are disambiguated", {
    isa <- seg_step()
    isa$out$chr_name.don <- "1"
    isa$out$chr_name.rec <- "1"
    labels <- isa_dashboard_chr_labels(isa, c("d1", "r1"))

    expect_identical(length(unique(labels)), 2L)
    expect_true(all(grepl("^1", labels)))
})

test_that("the explorer renders an interacCircos output", {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("bslib")
    testthat::skip_if_not_installed("DT")

    html <- as.character(explorer_results_ui("explorer"))
    expect_true(any(grepl("explorer-circos", html)))
})
