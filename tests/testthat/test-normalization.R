## T4.1: multi-sample dosage normalization.

dose_by_iid <- function(x, chr) {
    data <- x$data[x$data$chr_id == chr, , drop = FALSE]
    data$dose[order(data$iid)]
}

normalization_cohort <- function(donor, recipient = donor) {
    stopifnot(is.list(donor), is.list(recipient))
    stopifnot(identical(names(donor), names(recipient)))
    out <- Map(function(d, r) normalization_rca(d, r), donor, recipient)
    names(out) <- names(donor)
    out
}

test_that("median normalization produces the hand-computed dosage table", {
    donor <- list(
        A = c(2, 0, 10, 4),
        B = c(4, 0, 10, 4),
        C = c(6, 4, 10, 4)
    )
    rcas <- normalization_cohort(donor)
    normalized <- normalizeReadCounts(rcas)

    expect_named(normalized, c("A", "B", "C"))
    expect_equal(dose_by_iid(normalized$A, "d1"), c(2.5 / 4.5, 1, 1, 1),
        tolerance = 1e-12
    )
    expect_equal(dose_by_iid(normalized$B, "d1"), c(1, 1, 1, 1),
        tolerance = 1e-12
    )
    expect_equal(dose_by_iid(normalized$C, "d1"), c(6.5 / 4.5, 9, 1, 1),
        tolerance = 1e-12
    )
    ## Both sides use the same values, so their doses agree.
    expect_equal(dose_by_iid(normalized$A, "r1"), dose_by_iid(normalized$A, "d1"))

    for (sample in names(normalized)) {
        expect_identical(
            normalized[[sample]]$normalization,
            list(
                mode = "median", reference = NULL,
                by_side = TRUE, pseudocount = 0.5
            )
        )
        expect_identical(normalized[[sample]]$data$cor.gc,
            rcas[[sample]]$data$cor.gc
        )
        expect_true(all(normalized[[sample]]$data$ideal))
        expect_true(all(is.finite(normalized[[sample]]$data$dose)))
    }
})

test_that("reference normalization makes the reference sample exactly one", {
    donor <- list(
        A = c(2, 0, 10, 4),
        B = c(4, 0, 10, 4),
        C = c(6, 4, 10, 4)
    )
    rcas <- normalization_cohort(donor)
    normalized <- normalizeReadCounts(rcas, reference = "B")

    expect_true(all(normalized$B$data$dose == 1))
    expect_equal(dose_by_iid(normalized$A, "d1"), c(2.5 / 4.5, 1, 1, 1),
        tolerance = 1e-12
    )
    expect_equal(dose_by_iid(normalized$C, "d1"), c(6.5 / 4.5, 9, 1, 1),
        tolerance = 1e-12
    )
    expect_identical(
        normalized$A$normalization,
        list(
            mode = "reference", reference = "B",
            by_side = TRUE, pseudocount = 0.5
        )
    )
})

test_that("all-zero windows and zero reference windows get dose one", {
    donor <- list(
        A = c(2, 0, 10, 4),
        B = c(4, 0, 10, 4),
        C = c(6, 4, 10, 4)
    )
    zero_window <- lapply(donor, function(x) {
        x[[2L]] <- 0
        x
    })
    normalized <- normalizeReadCounts(normalization_cohort(zero_window))
    for (sample in names(normalized)) {
        expect_equal(dose_by_iid(normalized[[sample]], "d1")[2L], 1)
    }

    ## A zero reference window with positive samples uses the pseudocount ratio.
    donor$A[[2L]] <- 0
    donor$B[[2L]] <- 4
    donor$C[[2L]] <- 4
    normalized <- normalizeReadCounts(normalization_cohort(donor), reference = "A")
    expect_equal(dose_by_iid(normalized$A, "d1")[2L], 1)
    expect_equal(dose_by_iid(normalized$B, "d1")[2L], 4.5 / 0.5)
    expect_equal(dose_by_iid(normalized$C, "d1")[2L], 4.5 / 0.5)
})

test_that("missing values stay missing and are no longer ideal", {
    donor <- list(
        A = c(2, 0, 10, 4),
        B = c(4, 0, 10, 4),
        C = c(6, 4, 10, 4)
    )
    rcas <- normalization_cohort(donor)
    missing <- rcas
    idx <- which(missing$A$data$chr_id == "d1")[2L]
    missing$A$data$cor.gc[[idx]] <- NA_real_

    normalized <- normalizeReadCounts(missing)
    expect_true(is.na(dose_by_iid(normalized$A, "d1")[2L]))
    expect_false(normalized$A$data$ideal[[idx]])
    expect_equal(dose_by_iid(normalized$B, "d1")[2L], 0.5 / 2.5)
    expect_equal(dose_by_iid(normalized$C, "d1")[2L], 4.5 / 2.5)

    ## When every sample is missing, the window is not ideal anywhere.
    all_missing <- rcas
    for (sample in names(all_missing)) {
        idx <- which(all_missing[[sample]]$data$chr_id == "d1")[2L]
        all_missing[[sample]]$data$cor.gc[[idx]] <- NA_real_
    }
    normalized <- normalizeReadCounts(all_missing)
    for (sample in names(normalized)) {
        expect_true(is.na(dose_by_iid(normalized[[sample]], "d1")[2L]))
        idx <- which(normalized[[sample]]$data$chr_id == "d1")[2L]
        expect_false(normalized[[sample]]$data$ideal[[idx]])
    }
})

test_that("a single sample normalizes to one", {
    rcas <- normalization_cohort(list(A = c(2, 0, 10, 4)))
    normalized <- normalizeReadCounts(rcas)

    expect_true(all(normalized$A$data$dose == 1))
    expect_identical(normalized$A$normalization$mode, "median")
})

test_that("different window orders are aligned by (chr_id, iid)", {
    donor <- list(
        A = c(2, 0, 10, 4),
        B = c(4, 0, 10, 4),
        C = c(6, 4, 10, 4)
    )
    set.seed(42)
    shuffled <- normalization_cohort(donor)
    shuffled$B$data <- shuffled$B$data[
        sample(nrow(shuffled$B$data)),
        ,
        drop = FALSE
    ]
    normalized <- normalizeReadCounts(shuffled)
    expect_equal(dose_by_iid(normalized$B, "d1"), c(1, 1, 1, 1))
    expect_equal(dose_by_iid(normalized$A, "d1"), c(2.5 / 4.5, 1, 1, 1))

    ## Shuffling the first sample must not change any per-window result.
    shuffled_first <- normalization_cohort(donor)
    shuffled_first$A$data <- shuffled_first$A$data[
        sample(nrow(shuffled_first$A$data)),
        ,
        drop = FALSE
    ]
    normalized <- normalizeReadCounts(shuffled_first)
    expect_equal(dose_by_iid(normalized$A, "d1"), c(2.5 / 4.5, 1, 1, 1))
    expect_equal(dose_by_iid(normalized$C, "d1"), c(6.5 / 4.5, 9, 1, 1))
})

test_that("normalization is approximately scale invariant", {
    donor <- list(
        A = c(1000, 1000, 3000, 2000),
        B = c(2000, 1000, 3000, 2000),
        C = c(3000, 4000, 3000, 2000)
    )
    rcas <- normalization_cohort(donor)
    scaled <- lapply(rcas, function(x) {
        x$data$cor.gc <- x$data$cor.gc * 10
        x
    })
    plain <- normalizeReadCounts(rcas)
    multiplied <- normalizeReadCounts(scaled)
    expect_equal(
        dose_by_iid(plain$A, "d1"),
        dose_by_iid(multiplied$A, "d1"),
        tolerance = 1e-2
    )
})

test_that("by_side and joint normalization are different sensitivity modes", {
    rcas <- normalization_cohort(list(A = rep(10, 4)), list(A = rep(100, 4)))

    by_side <- normalizeReadCounts(rcas, by_side = TRUE)
    expect_equal(dose_by_iid(by_side$A, "d1")[1L], 1)
    expect_equal(dose_by_iid(by_side$A, "r1")[1L], 1)

    joint <- normalizeReadCounts(rcas, by_side = FALSE)
    expect_false(joint$A$normalization$by_side)
    expect_equal(dose_by_iid(joint$A, "d1")[1L], 10.5 / 55.5)
    expect_equal(dose_by_iid(joint$A, "r1")[1L], 100.5 / 55.5)
})

test_that("normalizeReadCounts() validates its input", {
    rcas <- normalization_cohort(list(
        A = c(2, 0, 10, 4),
        B = c(4, 0, 10, 4)
    ))

    expect_error(
        normalizeReadCounts(list()),
        "named list",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        normalizeReadCounts(list(rcas$A)),
        "one non-empty name",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        normalizeReadCounts(stats::setNames(list(rcas$A, rcas$B), c("A", "A"))),
        "duplicate",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        normalizeReadCounts(list(A = 1)),
        "not an RCA",
        class = "MaMaMIAAnalysisError"
    )

    uncorrected <- rcas
    uncorrected$A$corrected <- FALSE
    expect_error(
        normalizeReadCounts(uncorrected),
        "GC-corrected",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        normalizeReadCounts(rcas, reference = "nope"),
        "not one of",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        normalizeReadCounts(rcas, reference = 1),
        "single sample name",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        normalizeReadCounts(rcas, by_side = NA),
        "TRUE or FALSE",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        normalizeReadCounts(rcas, reference = "A", by_side = FALSE),
        "only available for median",
        class = "MaMaMIAAnalysisError"
    )
    for (bad in list(0, -1, NA_real_, Inf, "x", c(1, 2))) {
        expect_error(
            normalizeReadCounts(rcas, pseudocount = bad),
            "positive number",
            class = "MaMaMIAAnalysisError"
        )
    }

    different_chromosomes <- rcas
    different_chromosomes$B$data <- different_chromosomes$B$data[
        different_chromosomes$B$data$chr_id != "d2",
        ,
        drop = FALSE
    ]
    expect_error(
        normalizeReadCounts(different_chromosomes),
        "same chromosomes",
        class = "MaMaMIAAnalysisError"
    )

    different_windows <- rcas
    different_windows$B$data$iid[[1L]] <- 99L
    expect_error(
        normalizeReadCounts(different_windows),
        "missing window",
        class = "MaMaMIAAnalysisError"
    )

    different_coordinates <- rcas
    different_coordinates$B$data$start[[1L]] <-
        different_coordinates$B$data$start[[1L]] + 1
    expect_error(
        normalizeReadCounts(different_coordinates),
        "different coordinates",
        class = "MaMaMIAAnalysisError"
    )
})
