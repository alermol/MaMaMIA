## Candidate 7: one calibration object derived from the declared bins.

calibration_defaults <- function(bins = NULL, floors = NULL) {
    defaults <- .config_defaults()
    if (!is.null(bins)) {
        defaults$calibration$length_bins_mb <- bins
    }
    if (!is.null(floors)) {
        defaults$calibration$single_sample_floors <- floors
    }
    defaults
}

test_that("the bundled five-bin defaults resolve to the registry floors", {
    parameters <- resolveCalibrationParameters(.config_defaults())
    expect_s3_class(parameters, "MaMaMIACalibrationParameters")
    expect_identical(parameters$n_bins, 5L)
    expect_identical(
        parameters$labels,
        c("5-10", "10-25", "25-50", "50-100", ">100")
    )
    expect_equal(
        unname(parameters$floors), c(2.0, 1.5, 1.2, 1.0, 0.5)
    )
    expect_identical(parameters$fdr, 0.01)
    expect_null(parameters$controls)
    expect_identical(parameters$segmentation$alpha, 0.001)
    expect_identical(parameters$segmentation$min_width, 5L)
    expect_identical(parameters$segmentation$merge_gap, 100L)
    expect_identical(parameters$segmentation$min_mb, 5)
})

test_that("one, three, five and seven bins each get matching floors", {
    cases <- list(
        `1` = c(5, Inf),
        `3` = c(5, 25, 100, Inf),
        `5` = c(5, 10, 25, 50, 100, Inf),
        `7` = c(5, 10, 25, 50, 75, 100, 200, Inf)
    )
    expected <- c(1L, 3L, 5L, 7L)
    for (i in seq_along(cases)) {
        parameters <- resolveCalibrationParameters(
            calibration_defaults(bins = cases[[i]])
        )
        expect_identical(parameters$n_bins, expected[[i]])
        expect_identical(
            length(parameters$floors), expected[[i]]
        )
        expect_identical(length(parameters$labels), expected[[i]])
        expect_identical(
            names(parameters$floors), parameters$labels
        )
        expect_true(all(parameters$floors >= 0))
    }
})

test_that("runtime overrides replace defaults without name splicing", {
    parameters <- resolveCalibrationParameters(
        .config_defaults(),
        overrides = list(
            length_bins_mb = c(1, 10, Inf),
            fdr = 0.05,
            floors = c(0.7, 0.4),
            controls = c("C1", "C2"),
            alpha = 0.01,
            undo_SD = 2,
            min_width = 7L,
            merge_gap = 25L,
            merge_floor = -0.2,
            merge_min = 0.6,
            min_mb = 0.5
        )
    )
    expect_identical(parameters$length_bins_mb, c(1, 10, Inf))
    expect_identical(parameters$n_bins, 2L)
    expect_equal(unname(parameters$floors), c(0.7, 0.4))
    expect_identical(parameters$fdr, 0.05)
    expect_identical(parameters$controls, c("C1", "C2"))
    expect_identical(parameters$segmentation$alpha, 0.01)
    expect_identical(parameters$segmentation$undo_SD, 2)
    expect_identical(parameters$segmentation$min_width, 7L)
    expect_identical(parameters$segmentation$merge_gap, 25L)
    expect_equal(parameters$segmentation$merge_floor, -0.2)
    expect_equal(parameters$segmentation$merge_min, 0.6)
    expect_equal(parameters$segmentation$min_mb, 0.5)
})

test_that("single and multi modes choose different default floors", {
    single <- resolveCalibrationParameters(
        .config_defaults(), mode = "single"
    )
    multi <- resolveCalibrationParameters(
        .config_defaults(), mode = "multi"
    )
    expect_equal(
        unname(single$floors), c(2.0, 1.5, 1.2, 1.0, 0.5)
    )
    expect_equal(unname(multi$floors), rep(0.25, 5L))

    overridden <- resolveCalibrationParameters(
        .config_defaults(),
        overrides = list(floors = 0.4),
        mode = "multi"
    )
    expect_equal(unname(overridden$floors), rep(0.4, 5L))
})

test_that("configuration defaults can be passed as a full config", {
    config <- list(defaults = .config_defaults())
    parameters <- resolveCalibrationParameters(config)
    expect_identical(parameters$n_bins, 5L)

    null_defaults <- resolveCalibrationParameters(NULL)
    expect_identical(null_defaults$n_bins, 5L)
})

test_that("dashboard floor controls follow the declared bin count", {
    three <- list(
        defaults = calibration_defaults(bins = c(5, 25, 100, Inf))
    )
    labels <- bam_calibration_floor_labels(three)
    expect_identical(labels, c(
        "Floor 5-25 Mb", "Floor 25-100 Mb", "Floor >100 Mb"
    ))
    expect_identical(bam_calibration_floor_count(three), 3L)
    html <- as.character(bam_advanced_parameters_ui(
        three, list(calibration_floors = rep(0.25, 3L))
    ))
    expect_identical(
        lengths(regmatches(
            html,
            gregexpr('<input id="calibration_floor_[0-9]+"', html)
        ))[[1L]],
        3L
    )

    seven <- list(
        defaults = calibration_defaults(
            bins = c(5, 10, 25, 50, 75, 100, 200, Inf)
        )
    )
    expect_identical(bam_calibration_floor_count(seven), 7L)
    expect_length(bam_calibration_floor_labels(seven), 7L)
})

test_that("invalid calibration overrides are rejected clearly", {
    expect_error(
        resolveCalibrationParameters(
            .config_defaults(), overrides = list(fdr = 0)
        ),
        "`fdr`",
        class = "MaMaMIACalibrationError"
    )
    expect_error(
        resolveCalibrationParameters(
            .config_defaults(), overrides = list(floors = c(0.1, 0.2, 0.3))
        ),
        "one value per length bin",
        class = "MaMaMIACalibrationError"
    )
    expect_error(
        resolveCalibrationParameters(
            .config_defaults(),
            overrides = list(floors = c(0.1, -1, 0.2, 0.3, 0.4))
        ),
        "non-negative",
        class = "MaMaMIACalibrationError"
    )
    expect_error(
        resolveCalibrationParameters(
            .config_defaults(),
            overrides = list(length_bins_mb = c(5, 10))
        ),
        "end with Inf",
        class = "MaMaMIACalibrationError"
    )
})
