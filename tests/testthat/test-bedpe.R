bedpe_columns <- c(
    "chrom1", "start1", "end1",
    "chrom2", "start2", "end2",
    "name", "score", "strand1", "strand2",
    "chrom_name1", "chrom_name2", "subgen_name1", "subgen_name2"
)

test_that("writeToBedpe() writes putative introgressions", {
    isa <- segments(step_rca(), target_pairs = c("d1" = "r1"), verbose = FALSE)
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    suppressMessages(writeToBedpe(isa, file = target))
    expect_true(file.exists(target))

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    n_putative <- sum(isa$segments$putative_introgression)
    expect_equal(nrow(out), n_putative)
    expect_true(n_putative >= 1L)
    expect_true(all(bedpe_columns %in% names(out)))
    expect_true(all(out$chrom1 == "d1"))
    expect_true(all(out$chrom2 == "r1"))
    expect_true(all(out$score >= 0))
    expect_true(all(out$subgen_name1 == "At"))
    expect_true(all(out$subgen_name2 == "A"))
})

test_that("writeToBedpe() writes a header only when nothing is flagged", {
    isa <- segments(
        step_rca(),
        target_pairs = c("d1" = "r1"),
        min_segment_median = 1e9,
        verbose = FALSE
    )
    expect_false(any(isa$segments$putative_introgression))

    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    expect_warning(
        suppressMessages(writeToBedpe(isa, file = target)),
        "header only"
    )

    header <- strsplit(readLines(target)[1L], "\t")[[1L]]
    expect_equal(header, bedpe_columns)
    expect_equal(nrow(utils::read.table(target, header = TRUE, sep = "\t")), 0L)
})

test_that("writeToBedpe() rejects non-ISA input", {
    expect_error(writeToBedpe(list(), tempfile(fileext = ".bedpe")), "ISA")
})
