## P3: WindowMetrics construction.

window_metrics_fixture <- function() {
    skip_if_no_tiny_fixtures()
    config <- readMaMaMIAConfig(tiny_config_path())
    windows <- makeAnalysisWindows(
        tiny_assembly(config), tiny_contigs(), 100L
    )
    gc <- extractWindowGc(windows, tiny_assembly(config))
    header <- bamHeaderInfo(tiny_bam_path())
    measured <- countWindowMetrics(
        tiny_bam_path(), windows, "tiny-sample", header$read_groups
    )
    list(
        config = config, windows = windows, gc = gc, measured = measured
    )
}

test_that("buildWindowMetrics() implements the WindowMetrics contract", {
    fixture <- window_metrics_fixture()
    metrics <- buildWindowMetrics(
        fixture$measured, fixture$windows, fixture$gc
    )
    expect_named(metrics, c(
        "sample", "chr_id", "chr_name", "group", "role", "subgenome",
        "copy_number", "iid", "start", "stop", "length", "count", "covered",
        "breadth", "gc", "lambda_breadth", "lambda_count", "read_length"
    ))
    expect_identical(nrow(metrics), nrow(fixture$windows$table))
    expect_true(all(metrics$breadth >= 0 & metrics$breadth <= 1))
    expect_equal(
        metrics$breadth,
        metrics$covered / metrics$length,
        tolerance = 1e-12
    )
    expect_equal(
        metrics$lambda_breadth,
        -log(1 - metrics$breadth)
    )
    expect_equal(
        metrics$lambda_count,
        metrics$count * metrics$read_length / metrics$length
    )
    expect_identical(unique(metrics$read_length), 50)
    expect_identical(unique(metrics$group), "1")
    expect_identical(sort(unique(metrics$role)), c("donor", "recipient"))
})

test_that("buildWindowMetrics() rejects a mismatched GC table", {
    fixture <- window_metrics_fixture()
    gc <- fixture$gc
    gc$start[[1L]] <- gc$start[[1L]] + 1
    expect_error(
        buildWindowMetrics(fixture$measured, fixture$windows, gc),
        "GC is missing"
    )
})

test_that("buildWindowMetrics() requires a positive read length", {
    fixture <- window_metrics_fixture()
    measured <- fixture$measured
    measured$read_length <- NULL
    expect_error(
        buildWindowMetrics(measured, fixture$windows, fixture$gc),
        "read_length"
    )
    expect_silent(buildWindowMetrics(
        measured, fixture$windows, fixture$gc, read_length = 50
    ))
})

test_that("buildWindowMetrics() rejects windows outside the grid", {
    fixture <- window_metrics_fixture()
    measured <- fixture$measured
    measured$metrics$iid[[1L]] <- 0L
    measured$metrics$chr_id[[1L]] <- "unknown"
    expect_error(
        buildWindowMetrics(measured, fixture$windows, fixture$gc),
        "not in `windows`"
    )
})

test_that("buildWindowMetrics() keeps samples separate", {
    fixture <- window_metrics_fixture()
    measured <- fixture$measured
    second <- measured$metrics
    second$sample <- "second"
    measured$metrics$read_length <- 50
    second$read_length <- 50
    combined <- rbind(measured$metrics, second)
    metrics <- buildWindowMetrics(combined, fixture$windows, fixture$gc)
    expect_identical(sort(unique(metrics$sample)), c("second", "tiny-sample"))
    expect_identical(
        sum(metrics$count[metrics$sample == "second"]),
        sum(metrics$count[metrics$sample == "tiny-sample"])
    )
})
