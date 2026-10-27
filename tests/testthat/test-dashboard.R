## Shiny launcher and result helpers.

test_that("bam_run_parameters() resolves the run controls", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    params <- bam_run_parameters(
        list(resolution = 100), config, samples = "demo-sample"
    )
    expect_identical(params$resolution, 100L)
    expect_identical(params$metric, config$defaults$metric)
    expect_identical(params$score, config$defaults$score)
    expect_named(params, c(
        "resolution", "mapq", "exclude_flags", "cores", "metric", "score"
    ))
})

test_that("an unknown progress stage does not abort the analysis", {
    expect_silent(
        label <- bam_progress_label("metrics_done", list(n_windows = 10))
    )
    expect_identical(label, "metrics_done")
    expect_silent(bam_progress_label("brand_new_stage", list()))
})

test_that("progress labels cover the pipeline stages", {
    stages <- bam_progress_stages(2L)
    expect_true(all(c(
        "preflight", "index", "windows", "gc", "metrics", "score", "mask",
        "segmentation", "calibration"
    ) %in% stages))
    expect_identical(sum(stages == "metrics"), 2L)
    expect_match(bam_progress_label("metrics", list(total = 1000)), "reads")
    expect_match(
        bam_progress_label("score", list(mode = "single", metric = "breadth")),
        "single"
    )
    expect_match(
        bam_progress_label("mask", list(masked_fraction = 0.1)),
        "10%"
    )
    expect_match(
        bam_progress_label("calibration", list(n_calls = 3)),
        "3 calibrated"
    )
})

test_that("cohort summaries and warnings read WindowMetrics", {
    result <- list(
        samples = data.frame(
            sample = c("A", "B"),
            n_windows = c(2L, 2L),
            n_calls = c(1L, 0L),
            n_high_confidence = c(1L, 0L),
            masked_fraction = c(0.1, 0.1),
            mode = c("multi", "multi"),
            stringsAsFactors = FALSE
        ),
        metrics = data.frame(
            sample = c("A", "A", "B", "B"),
            role = c("donor", "recipient", "donor", "recipient"),
            count = c(10L, 20L, 0L, 0L),
            stringsAsFactors = FALSE
        )
    )
    summary <- bam_cohort_summary(result)
    expect_identical(summary$reads_total, c(30, 0))
    expect_identical(summary$n_calls, c(1L, 0L))
    warnings <- bam_coverage_warnings(result)
    expect_true(any(grepl("Sample 'B' has no reads", warnings)))
    expect_false(any(grepl("Sample 'A' has no reads", warnings)))
})

test_that("the launcher builds with the demo registry", {
    config <- system.file("extdata", "bam-demo", "assemblies.yml",
        package = "MaMaMIA")
    expect_true(file.exists(config))
    app <- bam_launcher_app(config = config)
    expect_s3_class(app, "shiny.appobj")
})

test_that("exploreIntrogressions() requires a BamAnalysisResult", {
    expect_error(
        exploreIntrogressions(list()),
        "BamAnalysisResult"
    )
})

test_that("the explorer score profile renders to a device", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    result <- bam_example_result(100)
    calls <- result$calls[result$calls$pass, , drop = FALSE]
    plot <- explorer_score_plot(result$scores, calls)
    expect_s3_class(plot, "ggplot")

    output <- tempfile(fileext = ".png")
    on.exit(unlink(output), add = TRUE)
    grDevices::png(output, width = 900, height = 500, type = "cairo")
    print(plot)
    grDevices::dev.off()
    expect_gt(file.info(output)$size, 0)
})
