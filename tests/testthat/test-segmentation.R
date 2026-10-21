test_that("segments() returns and annotates requested calls", {
    isa <- seg_step()
    expect_s3_class(isa, "ISA"); expect_true(nrow(isa$segments) >= 1L)
    expect_true(all(c("chr_id.don", "chr_id.rec", "loc.start", "loc.end", "num.mark", "seg.mean", "seg.sd", "seg.median", "seg.mad", "putative_introgression", "chr_name.don", "subgenome.don", "chr_name.rec", "subgenome.rec") %in% names(isa$segments)))
    expect_type(isa$segments$putative_introgression, "logical")
    expect_true(all(isa$segments$chr_id.don == "d1" & isa$segments$chr_id.rec == "r1"))
    expect_true(all(isa$segments$loc.start >= 1 & isa$segments$loc.end <= 60))
    expect_true(any(isa$segments$seg.median > 0))
    expect_equal(isa$meta$metrics, "seg.median")
    expect_setequal(isa$meta$don_chr_ids, "d1"); expect_setequal(isa$meta$rec_chr_ids, "r1")
})

test_that("segments() is unchanged when dose is absent", {
    plain <- seg_step()
    with_dose <- step_rca()
    with_dose$data$dose <- with_dose$data$cor.gc
    normalized <- segments(with_dose, target_pairs = tp, verbose = FALSE)

    expect_identical(plain$meta$value_col, "cor.gc")
    expect_identical(normalized$meta$value_col, "dose")
    expect_identical(plain$segments, normalized$segments)
    expect_equal(
        normalized$out$diff,
        normalized$out$dose.don - normalized$out$dose.rec
    )
})

test_that("segments(value_col =) selects the difference column", {
    x <- step_rca()
    x$data$dose <- x$data$cor.gc * 10
    dose_isa <- segments(x, target_pairs = tp, verbose = FALSE)
    cor_isa <- segments(
        x,
        target_pairs = tp,
        value_col = "cor.gc",
        verbose = FALSE
    )

    expect_identical(dose_isa$meta$value_col, "dose")
    expect_identical(cor_isa$meta$value_col, "cor.gc")
    expect_equal(dose_isa$out$diff, 10 * cor_isa$out$diff)
    expect_error(
        segments(x, target_pairs = tp, value_col = "nope", verbose = FALSE),
        "value_col"
    )
    expect_error(
        segments(x, target_pairs = tp, value_col = 1, verbose = FALSE),
        "value_col"
    )
})

test_that("plot.ISA() uses the normalized value column", {
    x <- step_rca()
    x$data$dose <- x$data$cor.gc
    isa <- segments(x, target_pairs = tp, verbose = FALSE)

    mirror <- plot(isa, plot.type = "mirror")
    difference <- plot(isa, plot.type = "diff")
    expect_s3_class(mirror, "ggplot")
    expect_s3_class(difference, "ggplot")
    expect_s3_class(ggplot2::ggplot_build(mirror), "ggplot_built")
    expect_s3_class(ggplot2::ggplot_build(difference), "ggplot_built")
})

test_that("plot.ISA() downsamples points but keeps every segment line", {
    isa <- seg_step()
    n_windows <- nrow(isa$out)

    full <- ggplot2::ggplot_build(plot(isa, plot.type = "diff"))
    expect_equal(nrow(full$data[[1L]]), n_windows)

    thinned <- ggplot2::ggplot_build(plot(
        isa,
        plot.type = "diff",
        max_points = 7L
    ))
    expect_equal(nrow(thinned$data[[1L]]), 7L)
    expect_equal(
        nrow(thinned$data[[3L]]),
        sum(isa$segments$putative_introgression)
    )

    mirror <- ggplot2::ggplot_build(plot(
        isa,
        plot.type = "mirror",
        max_points = 5L
    ))
    expect_equal(nrow(mirror$data[[1L]]), 5L)
    expect_equal(nrow(mirror$data[[2L]]), 5L)

    expect_equal(
        nrow(ggplot2::ggplot_build(plot(isa, max_points = Inf))$data[[1L]]),
        n_windows
    )
    expect_equal(
        nrow(ggplot2::ggplot_build(plot(isa, max_points = NULL))$data[[1L]]),
        n_windows
    )
    for (bad in list(-1, 0, NA_real_, "many")) {
        expect_error(plot(isa, max_points = bad), "positive number")
    }
})

test_that("segments() honours thresholds, widths and input validation", {
    permissive <- seg_step()
    strict <- seg_step(min_segment_median = 1e9)
    wide <- seg_step(min_width = 1000L)
    expect_true(any(permissive$segments$putative_introgression))
    expect_false(any(strict$segments$putative_introgression))
    expect_equal(nrow(strict$segments), nrow(permissive$segments))
    expect_equal(nrow(wide$segments), 0L)

    corrected <- step_rca()
    expect_error(segments(1), "not RCA")
    expect_error(segments(do.call(RCA, step_inputs()), target_pairs = tp), "corrected")
    for (bad in list(c(nope = "r1"), c(d1 = "nope"))) expect_error(segments(corrected, target_pairs = bad, verbose = FALSE), "not present")
    for (args in list(list(alpha = 0), list(alpha = 2), list(min_width = 1), list(undo_SD = 0), list(undo_SD = 11), list(seed = "x"), list(seed = c(1L, 2L)))) {
        expect_error(do.call(segments, c(list(x = corrected, target_pairs = tp, verbose = FALSE), args)))
    }
})

test_that("segments() is reproducible and leaves caller RNG untouched", {
    x <- step_rca()
    set.seed(999); first <- segments(x, target_pairs = tp, verbose = FALSE)
    set.seed(123); second <- segments(x, target_pairs = tp, verbose = FALSE)
    expect_identical(first$segments, second$segments)
    expect_identical(first$segments, segments(x, target_pairs = tp, verbose = FALSE)$segments)
    expect_equal(segments(x, target_pairs = tp, seed = 5L, verbose = FALSE)$param$seed, 5L)
    expect_equal(segments(x, target_pairs = tp, verbose = FALSE)$param$seed, 1L)
    set.seed(7); before <- stats::runif(3)
    set.seed(7); invisible(segments(x, target_pairs = tp, verbose = FALSE))
    expect_identical(before, stats::runif(3))
})

test_that("reverseWindows() re-segments affected pairs of an ISA", {
    isa <- seg_step(); reversed <- reverseWindows(isa, chr_ids = "d1")
    expect_s3_class(reversed, "ISA")
    expect_equal(reversed$param, isa$param)
    expect_setequal(reversed$meta$don_chr_ids, "d1"); expect_setequal(reversed$meta$rec_chr_ids, "r1")
    expect_true(nrow(reversed$segments) >= 1L)
    expect_false(identical(reversed$out$diff, isa$out$diff))
    expect_error(reverseWindows(isa, chr_ids = "nope"), "not present")
})

test_that("summary() and print() describe an ISA", {
    isa <- seg_step()
    expect_output(print(isa), "ISA object")
    expect_output(print(summary(isa, unit = "Mb")), "Putative introgressions found")
})
