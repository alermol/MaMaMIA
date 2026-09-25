test_that("RCA() builds a Read Count Array from windowed tables", {
    x <- small_rca()

    expect_s3_class(x, "RCA")
    expect_false(x$corrected)
    expect_equal(nrow(x$data), 4L * 40L)
    expect_setequal(unique(x$data$chr_id), c("d1", "d2", "r1", "r2"))
    expect_setequal(x$meta$don_chr_ids, c("d1", "d2"))
    expect_setequal(x$meta$rec_chr_ids, c("r1", "r2"))

    iid <- split(x$data$iid, x$data$chr_id)
    expect_true(all(vapply(
        iid,
        function(i) identical(as.integer(i), seq_len(40L)),
        logical(1L)
    )))

    expect_true(all(c("chr_id", "iid", "start", "stop", "cov", "gc") %in% names(x$data)))
    expect_true(all(c("chr_name", "subgenome") %in% names(x$data)))
})

test_that("RCA() accepts a matrix of windowed values", {
    inp <- small_inputs()
    inp$dcov <- as.matrix(inp$dcov)
    inp$dgc <- as.matrix(inp$dgc)

    expect_s3_class(do.call(RCA, inp), "RCA")
})

test_that("RCA() rejects malformed window tables", {
    inp <- small_inputs()
    inp$dcov <- inp$dcov[, 1:3]
    expect_error(do.call(RCA, inp), "4 columns")

    inp <- small_inputs()
    inp$meta <- inp$meta[, 1:2]
    expect_error(do.call(RCA, inp), "3 columns")

    inp <- small_inputs()
    inp$dcov$start[1L] <- inp$dcov$stop[1L]
    expect_error(do.call(RCA, inp), "start >= stop")

    inp <- small_inputs()
    inp$dcov$value[1L] <- NA_real_
    expect_error(do.call(RCA, inp), "must not contain NA")
})

test_that("RCA() requires contiguous, non-overlapping windows", {
    inp <- small_inputs()
    inp$dcov$start[3L] <- inp$dcov$start[3L] + 50
    expect_error(do.call(RCA, inp), "gaps")

    inp <- small_inputs()
    inp$dcov$start[3L] <- inp$dcov$start[3L] - 50
    expect_error(do.call(RCA, inp), "overlapping")
})

test_that("RCA() requires the same number of windows on every chromosome", {
    inp <- small_inputs()
    last_d1 <- max(inp$dcov$start[inp$dcov$chr_id == "d1"])
    inp$dcov <- inp$dcov[!(inp$dcov$chr_id == "d1" & inp$dcov$start == last_d1), ]
    inp$dgc <- inp$dgc[!(inp$dgc$chr_id == "d1" & inp$dgc$start == last_d1), ]

    expect_error(do.call(RCA, inp), "same number of windows")
})

test_that("RCA() requires disjoint donor and recipient subgenomes", {
    inp <- small_inputs()
    inp$meta$subgenome[inp$meta$chr_id == "r1"] <- "At"

    expect_error(do.call(RCA, inp), "disjoint")
})

test_that("RCA() requires coverage and GC tables to share coordinates", {
    inp <- small_inputs()
    shifted <- inp$dgc$chr_id == "d2"
    inp$dgc$start[shifted] <- inp$dgc$start[shifted] + 5
    inp$dgc$stop[shifted] <- inp$dgc$stop[shifted] + 5

    # Windows still tile, but no longer match the coverage windows.
    expect_error(do.call(RCA, inp))
})

test_that("RCA() requires chromosome metadata for every input chromosome", {
    inp <- small_inputs()
    inp$meta <- inp$meta[inp$meta$chr_id != "r2", ]

    expect_error(do.call(RCA, inp))
})

test_that("print() and summary() describe an RCA object", {
    x <- small_rca()

    expect_output(print(x), "RCA object")
    expect_output(print(x), "GC-corrected: FALSE")

    s <- summary(x)
    expect_s3_class(s, "summary.RCA")
    expect_equal(s$n_chromosomes, 4L)
    expect_false(s$corrected)
    expect_null(s$outliers)
    expect_output(print(s), "RCA Summary")
})
