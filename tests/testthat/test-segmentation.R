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
