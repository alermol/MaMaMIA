## P6: length-stratified calibration and call classification.

calibration_calls <- function() {
    data.frame(
        sample = c(
            "C1", "C1", "C1", "C2", "C2", "C2",
            "P1", "P1", "P1"
        ),
        mode = "multi",
        group = "1",
        width_bp = c(6e6, 30e6, 500e6, 7e6, 40e6, 400e6, 6e6, 30e6, 500e6),
        score_mean = c(0.1, 0.2, 0.3, 0.15, 0.25, 0.35, 1.5, 1.6, 1.7),
        donor_mean = c(0.1, 0.2, 0.3, 0.15, 0.25, 0.35, 1.5, 1.6, 1.7),
        recipient_mean = 0,
        stringsAsFactors = FALSE
    )
}

test_that("calibration uses length bins and floors", {
    calls <- calibration_calls()
    calibrated <- calibrateSegments(
        calls, fdr = 0.01,
        length_bins_mb = c(5, 10, 25, 50, 100, Inf),
        floors = c(0.5, 0.4, 0.3, 0.25, 0.4)
    )
    expect_named(calibrated, c("calibration", "calls"))
    expect_true(all(c(
        "length_bin", "n_control", "q99", "threshold", "fdr", "version"
    ) %in% names(calibrated$calibration)))
    expect_identical(
        sort(unique(calibrated$calibration$length_bin)),
        sort(c("5-10", "10-25", "25-50", "50-100", ">100"))
    )
    expect_true(all(
        calibrated$calibration$threshold >= c(0.5, 0.4, 0.3, 0.25, 0.4)
    ))
    expect_identical(
        sort(unique(calibrated$calls$sample[calibrated$calls$pass])),
        "P1"
    )
    expect_true(all(calibrated$calls$pass[calibrated$calls$sample == "P1"]))
})

test_that("controls restrict the empirical null", {
    calls <- calibration_calls()
    calibrated <- calibrateSegments(
        calls, controls = "C1", fdr = 0.01,
        length_bins_mb = c(5, 10, 25, 50, 100, Inf),
        floors = c(0.5, 0.4, 0.3, 0.25, 0.4)
    )
    expect_true(all(calibrated$calibration$n_control <= 2L))
})

test_that("a sparse length bin borrows the mode-wide null", {
    calls <- calibration_calls()
    calls <- calls[calls$width_bp < 100e6, , drop = FALSE]
    calibrated <- calibrateSegments(
        calls, fdr = 0.01,
        length_bins_mb = c(5, 10, 25, 50, 100),
        floors = c(0, 0, 0, 0)
    )
    expect_true(all(is.finite(calibrated$calibration$threshold)))
    expect_identical(nrow(calibrated$calibration), 4L)
})

test_that("self-calibration does not let a target call define the null", {
    scores <- synthetic_score_table(event_effect = 2, seed = 21)
    calls <- segmentScores(scores, min_mb = 1)
    calibrated <- calibrateSegments(
        calls, fdr = 0.01,
        length_bins_mb = c(5, 10, 25, 50, 100, Inf),
        floors = c(0.3, 0.3, 0.3, 0.3, 0.3)
    )
    expect_gte(sum(calibrated$calls$pass), 1L)
})

test_that("calibration validates fdr and empty input", {
    calls <- calibration_calls()
    expect_error(calibrateSegments(calls, fdr = 0), "fdr")
    empty <- calibrateSegments(calls[0, , drop = FALSE])
    expect_identical(nrow(empty$calibration), 0L)
})

test_that("classifyCalls() assigns class and confidence", {
    calls <- data.frame(
        sample = c("S", "S", "S", "S"),
        group = "1",
        start_iid = 1:4,
        donor_mean = c(2, 1.2, 0.7, 0.1),
        recipient_mean = c(-1, -0.6, 0.1, 0.0),
        score_mean = c(2, 1.2, 0.7, 0.1),
        pass = c(TRUE, TRUE, FALSE, TRUE),
        stringsAsFactors = FALSE
    )
    out <- classifyCalls(calls)
    expect_identical(out$class, c(
        "substitution", "substitution", "addition", "other"
    ))
    expect_identical(out$confidence, c("high", "high", "moderate", "low"))
})

test_that("classifyCalls() handles empty input", {
    calls <- data.frame(
        sample = character(), group = character(),
        donor_mean = numeric(), recipient_mean = numeric(),
        score_mean = numeric(), stringsAsFactors = FALSE
    )
    out <- classifyCalls(calls)
    expect_identical(out$class, character())
    expect_identical(out$confidence, character())
})
