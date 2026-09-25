test_that("dashboard helpers annotate segment coordinates", {
    isa <- seg_step()
    dash <- isa_dashboard_segments(isa)

    expect_equal(dash$.row_id, seq_len(nrow(isa$segments)))
    expect_true(all(c("start.don", "stop.don", "start.rec", "stop.rec") %in% names(dash)))
    expect_true(all(dash$start.don <= dash$stop.don))
    expect_true(all(dash$start.rec <= dash$stop.rec))
})

test_that("dashboard filtering respects chromosomes, pairs, and putative flag", {
    isa <- seg_step()
    dash <- isa_dashboard_segments(isa)

    expect_equal(nrow(filter_dashboard_segments(dash, "d1", "d1 -> r1", FALSE)), 0L)
    filtered <- filter_dashboard_segments(dash, c("d1", "r1"), "d1 -> r1", TRUE)
    expect_true(all(filtered$putative_introgression))
    expect_equal(nrow(filtered), sum(dash$putative_introgression))
})

test_that("selected_isa_for_bedpe() marks selected source rows only", {
    isa <- seg_step()
    selected <- selected_isa_for_bedpe(isa, c(1L, 3L))

    expect_s3_class(selected, "ISA")
    expect_equal(which(selected$segments$putative_introgression), intersect(c(1L, 3L), seq_len(nrow(isa$segments))))
    expect_false(identical(selected$segments$putative_introgression, isa$segments$putative_introgression))
})

test_that("dashboard genome uses chromosome lengths from ISA windows", {
    isa <- seg_step()
    genome <- isa_dashboard_genome(isa)

    expect_named(genome, c("d1", "r1"), ignore.order = TRUE)
    expect_true(all(unlist(genome) > 0))
    expect_named(isa_dashboard_genome(isa, "d1"), "d1")
})

test_that("diff windows cover both chromosomes with a clipped limit", {
    isa <- seg_step()
    windows <- isa_dashboard_diff_windows(isa, c("d1", "r1"))

    expect_setequal(names(windows), c("chr", "start", "stop", "diff"))
    expect_setequal(unique(windows$chr), c("d1", "r1"))
    expect_true(all(is.finite(windows$diff)))

    lim <- isa_dashboard_diff_limit(windows)
    expect_true(lim > 0)
    expect_true(stats::quantile(abs(windows$diff), 0.99, na.rm = TRUE) <= lim)

    empty <- isa_dashboard_diff_windows(isa, "nope")
    expect_equal(nrow(empty), 0L)
    expect_equal(isa_dashboard_diff_limit(empty), 1)
})

test_that("diff windows are thinned per chromosome around the peaks", {
    windows <- data.frame(
        chr = "c1",
        start = seq_len(1000),
        stop = seq_len(1000) + 1,
        diff = rep(c(1, 9), 500)
    )
    thinned <- isa_dashboard_thin_windows(windows, 100L)

    expect_equal(nrow(thinned), 100L)
    expect_true(mean(abs(thinned$diff)) > mean(abs(windows$diff)))
    expect_equal(nrow(isa_dashboard_thin_windows(windows, 5000L)), 1000L)
})

test_that("size classes are monotone in segment size", {
    size <- c(1, 1, 5, 10, 100, 1000)
    class <- isa_dashboard_size_class(size, 4L)

    expect_equal(length(class), length(size))
    expect_true(all(class >= 1L & class <= 4L))
    expect_true(all(diff(class[order(size)]) >= 0L))
    expect_true(all(isa_dashboard_size_class(c(5, 5, 5), 4L) == 1L))
})

test_that("selection toggles from the table and from a band click", {
    expect_equal(isa_dashboard_toggle_selection(integer(), 3L), 3L)
    expect_equal(isa_dashboard_toggle_selection(c(1L, 3L), 3L), 1L)
    expect_equal(isa_dashboard_toggle_selection(c(1L, 3L), 2L), c(1L, 3L, 2L))
    expect_equal(isa_dashboard_toggle_selection(c(1L, 3L), NA_integer_), c(1L, 3L))
    expect_equal(isa_dashboard_toggle_selection(integer(), integer()), integer())
})

test_that("the circos stacks difference, spans and size-classed links", {
    testthat::skip_if_not_installed("BioCircos")

    isa <- seg_step()
    segments <- isa_dashboard_segments(isa)
    genome <- isa_dashboard_genome(isa)
    circos <- isa_dashboard_circos(isa, segments, genome)

    expect_s3_class(circos, "htmlwidget")
    track_names <- vapply(circos$x$tracklist, function(track) track[[1]], character(1))
    expect_true(any(grepl("CNV_diff_positive", track_names)))
    expect_true(any(grepl("CNV_diff_negative", track_names)))
    expect_true(any(grepl("CNV_introgression_spans", track_names)))
    expect_true(any(grepl("introgression_hit_area", track_names)))
    expect_true(all(grepl("LINK_", track_names[grepl("LINK_", track_names)])))
    link_tracks <- track_names[grepl("LINK_", track_names)]
    expect_true(length(link_tracks) > 0L)
    widths <- vapply(
        circos$x$tracklist[grepl("LINK_", track_names)],
        function(track) as.numeric(sub("px$", "", track[[2]]$LinkWidth)),
        numeric(1)
    )
    expect_true(all(diff(widths) >= 0))

    fill <- circos$x$genomeFillColor
    expect_equal(length(fill), length(genome))
    expect_setequal(unique(fill), c(dashboard_palette[["donor"]], dashboard_palette[["recipient"]]))
    expect_equal(fill[names(genome) %in% isa$meta$don_chr_ids], dashboard_palette[["donor"]])
    expect_equal(circos$x$genomeLabelTextColor, "#2f3437")
    ## The stylesheet and the enhancement script must ride on the widget as a
    ## dependency: prependContent() is ignored in a Shiny render, which is how the
    ## widget once ended up unreadable, selectable, and swallowing table clicks.
    deps <- Filter(function(d) inherits(d, "html_dependency"), circos$dependencies)
    dep_names <- vapply(deps, function(d) d$name, character(1))
    expect_true("mamamia-circos" %in% dep_names)
    dep <- deps[[match("mamamia-circos", dep_names)]]
    expect_equal(dep$script, "mamamia-circos.js")
    expect_equal(dep$stylesheet, "mamamia-circos.css")

    summary <- circos$x$chromosome_summary
    expect_setequal(names(summary), names(genome))
    expect_true(all(nzchar(unlist(summary))))
    expect_setequal(circos$x$donor_chr_ids, isa$meta$don_chr_ids)
    expect_setequal(circos$x$chromosome_ids, names(genome))
    expect_equal(circos$x$hit_area_rows, as.integer(segments$.row_id))
    expect_equal(length(circos$x$ribbon_rows), nrow(segments))
    expect_equal(length(circos$x$span_rows), 2L * nrow(segments))
})
