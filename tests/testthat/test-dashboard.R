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
    expect_identical(
        params$donor_aggregation, config$defaults$donor_aggregation
    )
    expect_named(params, c(
        "resolution", "mapq", "exclude_flags", "cores", "metric", "score",
        "donor_aggregation"
    ))
})

test_that("the BAM chooser accepts exactly one file", {
    chooser <- as.character(bam_single_file_input())
    expect_length(chooser, 1L)
    expect_match(chooser, 'id="bams"', fixed = TRUE)
    expect_match(chooser, 'accept=".bam"', fixed = TRUE)
    expect_match(chooser, "Choose BAM file", fixed = TRUE)
    expect_match(chooser, "No BAM selected", fixed = TRUE)
    expect_false(grepl("multiple", chooser, fixed = TRUE))
})

test_that("the upload path helper rejects several files", {
    one <- data.frame(datapath = "a.bam", stringsAsFactors = FALSE)
    expect_identical(bam_single_upload_path(one), "a.bam")
    expect_identical(bam_single_upload_path(NULL), character())
    expect_identical(
        bam_single_upload_path(data.frame(datapath = character())),
        character()
    )
    two <- data.frame(
        datapath = c("a.bam", "b.bam"), stringsAsFactors = FALSE
    )
    expect_error(
        bam_single_upload_path(two),
        "one BAM per run",
        class = "MaMaMIAUploadError"
    )
})

test_that("donor aggregation defaults from the registry and is set under Advanced", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    expect_identical(config$defaults$donor_aggregation, "max")
    defaults <- bam_run_parameters(
        list(resolution = 100), config, samples = "demo-sample"
    )
    expect_identical(defaults$donor_aggregation, "max")
    meaned <- bam_run_parameters(
        list(resolution = 100, donor_aggregation = "mean"),
        config,
        samples = "demo-sample"
    )
    expect_identical(meaned$donor_aggregation, "mean")
    expect_error(
        bam_run_parameters(
            list(resolution = 100, donor_aggregation = "median"),
            config,
            samples = "demo-sample"
        ),
        "donor_aggregation"
    )
    html <- as.character(bam_sidebar_details_ui(config))
    expect_match(html, 'id="donor_aggregation"', fixed = TRUE)
    expect_match(html, "Advanced", fixed = TRUE)
})

test_that("the analysis key changes with the donor aggregation", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    preflight <- list(
        ok = TRUE,
        table = data.frame(
            path = "demo.bam", sample = "demo-sample",
            stringsAsFactors = FALSE
        ),
        assembly = "triticum_like__A__secale_like__R",
        config = config
    )
    max_key <- bam_analysis_key(preflight, 100, "max")
    mean_key <- bam_analysis_key(preflight, 100, "mean")
    expect_true(is.character(max_key) && nzchar(max_key))
    expect_false(identical(max_key, mean_key))
})

test_that("the species pair resolves every declared component", {
    fixture <- multi_component_fixture(samples = "sample1")
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    config <- readMaMaMIAConfig(fixture$config_path)
    resolved <- bam_resolve_species_pair(config, "recipient", "donor")
    expect_identical(resolved$recipient_components, c("A", "B"))
    expect_identical(resolved$donor_components, "R")
    expect_identical(
        sort(unique(resolved$assembly$contigs$recipient$component)),
        c("A", "B")
    )
    expect_identical(
        unique(resolved$assembly$contigs$donor$component), "R"
    )
})

test_that("the sidebar no longer offers component selectors", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    html <- as.character(bam_sidebar_species_ui(config))
    expect_false(grepl("recipient_component", html, fixed = TRUE))
    expect_false(grepl("donor_component", html, fixed = TRUE))
    expect_match(html, 'id="component_summary_ui"', fixed = TRUE)
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

test_that("the explorer has no chromosome group selector", {
    skip_if_not_installed("htmlwidgets")
    bam <- as.character(bam_results_ui())
    standalone <- as.character(explorer_controls_ui("explorer"))
    filters <- as.character(explorer_view_filters_ui("explorer"))
    for (html in list(bam, standalone, filters)) {
        expect_false(grepl("Chromosome groups", html, fixed = TRUE))
    }
    expect_match(filters, "Call class", fixed = TRUE)
    expect_match(filters, "Confidence", fixed = TRUE)
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
