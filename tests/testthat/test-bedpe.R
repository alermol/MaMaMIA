bedpe_columns <- c(
    "chrom1", "start1", "end1",
    "chrom2", "start2", "end2",
    "name", "score", "strand1", "strand2",
    "chrom_name1", "chrom_name2", "subgen_name1", "subgen_name2"
)

## A hand-built ISA with predictable coordinates and values, so the exact BEDPE
## bytes can be asserted without depending on CBS.
fake_bedpe_segments <- function(putative = c(TRUE, FALSE)) {
    data.frame(
        chr_id.don = rep("d1", length(putative)),
        chr_id.rec = rep("r1", length(putative)),
        loc.start = c(1L, 2L)[seq_along(putative)],
        loc.end = c(2L, 3L)[seq_along(putative)],
        num.mark = c(1L, 2L)[seq_along(putative)],
        seg.mean = c(0.5, 2.5)[seq_along(putative)],
        seg.sd = 0,
        seg.median = c(0.5, 2.5)[seq_along(putative)],
        seg.mad = 0,
        putative_introgression = putative,
        chr_name.don = rep("1D", length(putative)),
        subgenome.don = rep("D", length(putative)),
        chr_name.rec = rep("1R", length(putative)),
        subgenome.rec = rep("R", length(putative)),
        stringsAsFactors = FALSE
    )
}

fake_bedpe_isa <- function(segments = fake_bedpe_segments()) {
    out <- data.frame(
        iid = 1:3,
        chr_id.don = "d1",
        start.don = c(0L, 50L, 100L),
        stop.don = c(50L, 100L, 150L),
        chr_id.rec = "r1",
        start.rec = c(0L, 50L, 100L),
        stop.rec = c(50L, 100L, 150L),
        stringsAsFactors = FALSE
    )
    structure(
        list(
            out = out,
            segments = segments,
            meta = list(
                metrics = "seg.median",
                metrics_threshold = 0,
                value_col = "cor.gc"
            )
        ),
        class = "ISA"
    )
}

fake_bedpe_cohort <- function() {
    structure(
        list(
            isa = list(
                `sample-a` = fake_bedpe_isa(),
                `sample-b` = fake_bedpe_isa()
            )
        ),
        class = c("BamAnalysisResult", "list")
    )
}

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

test_that("the legacy layout is byte-for-byte unchanged when sample is NULL", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    suppressMessages(writeToBedpe(fake_bedpe_isa(), target))

    expected <- c(
        paste(bedpe_columns, collapse = "\t"),
        paste(
            c(
                "d1", "0", "100", "r1", "0", "100", "PIS1", "0.5",
                ".", ".", "1D", "1R", "D", "R"
            ),
            collapse = "\t"
        )
    )
    expect_identical(readLines(target), expected)
})

test_that("writeToBedpe() appends the sample column when requested", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    suppressMessages(writeToBedpe(fake_bedpe_isa(), target, sample = "sample-a"))

    header <- strsplit(readLines(target)[1L], "\t")[[1L]]
    expect_identical(header, c(bedpe_columns, "sample"))

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(out$sample, "sample-a")
})

test_that("a sample header is kept when no row is selected", {
    isa <- fake_bedpe_isa(segments = fake_bedpe_segments()[0, , drop = FALSE])
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    expect_warning(
        suppressMessages(writeToBedpe(isa, target, sample = "sample-a")),
        "header only"
    )
    header <- strsplit(readLines(target)[1L], "\t")[[1L]]
    expect_identical(header, c(bedpe_columns, "sample"))
    expect_equal(nrow(utils::read.table(target, header = TRUE, sep = "\t")), 0L)
})

test_that("writeToBedpe() validates its input", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    expect_error(writeToBedpe(list(), target), "ISA")
    expect_error(
        writeToBedpe(fake_bedpe_isa(), target, sample = c("a", "b")),
        "single non-empty string"
    )
    expect_error(
        writeToBedpe(fake_bedpe_isa(), target, sample = ""),
        "single non-empty string"
    )
})

test_that("writeCohortToBedpe() combines every putative call with sample identity", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    suppressMessages(writeCohortToBedpe(fake_bedpe_cohort(), target))

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(names(out), c(bedpe_columns, "sample"))
    expect_identical(nrow(out), 2L)
    expect_identical(out$sample, c("sample-a", "sample-b"))
    expect_true(all(out$name == "PIS1"))
})

test_that("writeCohortToBedpe() maps selected rows per sample", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    suppressMessages(writeCohortToBedpe(
        fake_bedpe_cohort(), target,
        selected = list(`sample-a` = 2L, `sample-b` = 1L)
    ))

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(nrow(out), 2L)
    expect_identical(out$sample, c("sample-a", "sample-b"))
    ## Row 2 of the fake ISA runs from window iid 2 to 3.
    expect_identical(out$start1, c(50L, 0L))
    expect_identical(out$end1, c(150L, 100L))
    expect_identical(out$score, c(2.5, 0.5))
})

test_that("writeCohortToBedpe() omits samples missing from selected", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    suppressMessages(writeCohortToBedpe(
        fake_bedpe_cohort(), target,
        selected = list(`sample-a` = 1L)
    ))

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(out$sample, "sample-a")
    expect_equal(nrow(out), 1L)
})

test_that("writeCohortToBedpe() keeps the header when nothing is selected", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    expect_warning(
        suppressMessages(writeCohortToBedpe(
            fake_bedpe_cohort(), target,
            selected = list(`sample-a` = integer(), `sample-b` = integer())
        )),
        "header only"
    )

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(names(out), c(bedpe_columns, "sample"))
    expect_equal(nrow(out), 0L)
})

test_that("writeCohortToBedpe() honours the samples subset", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    suppressMessages(writeCohortToBedpe(
        fake_bedpe_cohort(), target,
        samples = "sample-b"
    ))

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(unique(out$sample), "sample-b")
})

test_that("writeCohortToBedpe() accepts a named list of ISA objects", {
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    isas <- fake_bedpe_cohort()$isa
    suppressMessages(writeCohortToBedpe(isas, target))

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(out$sample, c("sample-a", "sample-b"))
})

test_that("writeCohortToBedpe() validates samples and selected rows", {
    cohort <- fake_bedpe_cohort()
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    expect_error(writeCohortToBedpe(list(), target), "non-empty named list")
    expect_error(writeCohortToBedpe(cohort, target, samples = "nope"), "Unknown sample")
    expect_error(
        writeCohortToBedpe(cohort, target, selected = list(nope = 1L)),
        "not exported"
    )
    expect_error(
        writeCohortToBedpe(cohort, target, selected = list(`sample-a` = 99L)),
        "between 1 and 2"
    )
    expect_error(
        writeCohortToBedpe(cohort, target, selected = list(`sample-a` = 1.5)),
        "whole numbers"
    )
    expect_error(
        writeCohortToBedpe(cohort, target, selected = list(1L)),
        "named list"
    )
})

test_that("writeToBedpe() rejects non-ISA input", {
    expect_error(writeToBedpe(list(), tempfile(fileext = ".bedpe")), "ISA")
})
