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
    typed <- bam_run_parameters(
        list(resolution = 100), config,
        samples = "demo-sample", resolution = 1010L
    )
    expect_identical(typed$resolution, 1010L)
    expect_identical(params$metric, config$defaults$metric)
    expect_identical(params$score, config$defaults$score)
    expect_identical(
        params$donor_aggregation, config$defaults$donor_aggregation
    )
    expect_named(params, c(
        "resolution", "mapq", "exclude_flags", "cores", "metric", "score",
        "donor_aggregation", "control_samples", "calibration"
    ))
    expect_identical(params$control_samples, character(0))
    expect_null(params$calibration)
})

test_that("the resolution picker is a ten-window slider with three ticks", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    range <- bam_resolution_range(config$defaults)
    expect_identical(range$min, 100L)
    expect_identical(range$max, 5000L)
    expect_identical(range$ticks, c(100L, 1000L, 5000L))
    expect_identical(range$step, 10L)
    expect_identical(bam_default_resolution(config$defaults), 100L)
    expect_identical(bam_resolution_value(1010, config$defaults), 1010L)
    expect_error(bam_resolution_value(1005, config$defaults), "steps of 10")
    expect_error(bam_resolution_value(5001, config$defaults), "between")
    slider <- as.character(bam_resolution_slider_ui(config$defaults))
    expect_match(slider, 'id="resolution"', fixed = TRUE)
    expect_match(slider, 'data-min="100"', fixed = TRUE)
    expect_match(slider, 'data-max="5000"', fixed = TRUE)
    expect_match(slider, 'data-step="10"', fixed = TRUE)
    expect_match(slider, 'mamamia-resolution-ticks', fixed = TRUE)
    expect_match(slider, 'id="resolution_number"', fixed = TRUE)
    expect_match(slider, 'type="number"', fixed = TRUE)
    expect_match(slider, 'min="100"', fixed = TRUE)
    expect_match(slider, 'max="5000"', fixed = TRUE)
    expect_match(slider, 'step="10"', fixed = TRUE)
    expect_match(slider, "1,000", fixed = TRUE)
    expect_match(slider, "5,000", fixed = TRUE)
})

test_that("the example resolution slider is locked at its precomputed value", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    slider <- as.character(bam_example_resolution_slider_ui(config$defaults))
    expect_match(slider, 'id="resolution"', fixed = TRUE)
    expect_match(slider, 'data-min="100"', fixed = TRUE)
    expect_match(slider, 'data-max="5000"', fixed = TRUE)
    expect_match(slider, 'data-from="1000"', fixed = TRUE)
    expect_match(slider, 'resolution_locked', fixed = TRUE)
    expect_match(slider, 'id="resolution_number"', fixed = TRUE)
    expect_match(slider, 'disabled="disabled"', fixed = TRUE)
    expect_match(
        slider,
        "Fixed at 1,000 windows per chromosome.",
        fixed = TRUE
    )
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
    expect_match(html, 'id="donor_aggregation_ui"', fixed = TRUE)
    expect_match(html, "Advanced", fixed = TRUE)
})

test_that("control samples are selected under Advanced", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    samples <- c("S1", "S2", "S3")
    selected <- bam_run_parameters(
        list(resolution = 100, control_samples = c("S2", "S3")),
        config, samples = samples
    )
    expect_identical(selected$control_samples, c("S2", "S3"))
    expect_identical(selected$calibration$controls, c("S2", "S3"))
    expect_error(
        bam_run_parameters(
            list(resolution = 100, control_samples = "S9"),
            config, samples = samples
        ),
        "Unknown control sample"
    )
    single <- bam_run_parameters(
        list(resolution = 100, control_samples = "S1"),
        config, samples = "S1"
    )
    expect_identical(single$control_samples, character(0))
    expect_null(single$calibration)
    html <- as.character(bam_sidebar_details_ui(config))
    expect_match(html, 'id="control_samples_ui"', fixed = TRUE)
})

test_that("the analysis key changes with control samples", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    preflight <- list(
        ok = TRUE,
        table = data.frame(
            path = c("a.bam", "b.bam", "c.bam"),
            sample = c("S1", "S2", "S3"),
            stringsAsFactors = FALSE
        ),
        assembly = "triticum_like__A__secale_like__R",
        config = config
    )
    none <- bam_analysis_key(preflight, 100, "max", NULL)
    controls <- bam_analysis_key(preflight, 100, "max", c("S2"))
    all <- bam_analysis_key(preflight, 100, "max", c("S1", "S2", "S3"))
    expect_true(is.character(none) && nzchar(none))
    expect_false(identical(none, controls))
    expect_identical(none, all)
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

test_that("the species dropdowns use registry labels", {
    config <- bam_example_config()
    choices <- bam_species_choices(config)
    expect_identical(
        names(choices),
        c("Triticum-like (A)", "Secale-like (R)")
    )
    expect_identical(
        unname(choices[c("Triticum-like (A)", "Secale-like (R)")]),
        c("triticum_like", "secale_like")
    )
    expect_false(any(grepl(
        "triticum_like|secale_like", names(choices)
    )))
    html <- as.character(bam_sidebar_species_ui(config))
    expect_match(html, ">Triticum-like (A)</option>", fixed = TRUE)
    expect_match(html, ">Secale-like (R)</option>", fixed = TRUE)
})

test_that("the Assembly panel keeps its title", {
    config <- bam_example_config()
    html <- as.character(bam_sidebar_details_ui(config))
    expect_match(html, "Assembly", fixed = TRUE)
    expect_false(grepl('id="assembly_panel_title"', html, fixed = TRUE))
})

test_that("the preflight recovery card is calm and actionable", {
    assembly <- list(
        stage = "assembly",
        message = paste(
            "BAM '/tmp/session/001-0.bam' is missing target(s) required",
            "by assembly 'wheat__rye': CP155614.1, CP155613.1"
        ),
        suggestion = "Choose a matching species pair."
    )
    html <- as.character(bam_sidebar_recovery_ui(assembly))
    expect_match(html, "mamamia-recovery", fixed = TRUE)
    expect_match(html, "No matching assembly", fixed = TRUE)
    expect_match(html, "does not contain the contigs", fixed = TRUE)
    expect_match(html, "calculation log lists the missing targets", fixed = TRUE)
    expect_match(html, "Choose assembly", fixed = TRUE)
    expect_match(html, "Choose a matching species pair.", fixed = TRUE)
    expect_false(grepl("text-danger", html, fixed = TRUE))
    expect_false(grepl("/tmp/session", html, fixed = TRUE))
    expect_false(grepl("CP155614.1", html, fixed = TRUE))
    locked <- as.character(
        bam_sidebar_recovery_ui(assembly, show_action = FALSE)
    )
    expect_false(grepl("unlock_setup", locked, fixed = TRUE))
    generic <- as.character(bam_sidebar_recovery_ui(list(
        stage = "index",
        message = "The index could not be written.",
        suggestion = "Check the temporary directory."
    )))
    expect_match(generic, "Setup needs attention", fixed = TRUE)
    expect_match(generic, "Review setup", fixed = TRUE)
})

test_that("the sidebar no longer offers component selectors", {
    config <- readMaMaMIAConfig(
        system.file("extdata", "bam-demo", "assemblies.yml",
            package = "MaMaMIA")
    )
    html <- as.character(bam_sidebar_species_ui(config))
    expect_false(grepl("recipient_component", html, fixed = TRUE))
    expect_false(grepl("donor_component", html, fixed = TRUE))
    expect_match(html, 'id="recipient_species"', fixed = TRUE)
    expect_match(html, 'id="donor_species"', fixed = TRUE)
    expect_match(html, 'id="swap_species"', fixed = TRUE)
    expect_false(grepl('id="component_summary_ui"', html, fixed = TRUE))
})

test_that("example analysis controls render as read-only fields", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    html <- as.character(bam_sidebar_species_ui(
        bam_example_config(),
        selection = bam_example_species_selection(),
        locked = TRUE
    ))
    expect_false(grepl('id="recipient_species"', html, fixed = TRUE))
    expect_false(grepl('id="donor_species"', html, fixed = TRUE))
    expect_false(grepl('id="swap_species"', html, fixed = TRUE))
    expect_match(html, "Recipient species", fixed = TRUE)
    expect_match(html, "Donor species", fixed = TRUE)
    expect_match(html, "Fixed for the bundled example", fixed = TRUE)
    expect_false(grepl('id="component_summary_ui"', html, fixed = TRUE))
})

test_that("example analysis keys ignore stale run controls", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    example <- bam_example_preflight_state()
    expected <- bam_analysis_key(example, 1000, "max", NULL)
    stale <- bam_analysis_key(example, 5000, "mean", "example-sample")
    expect_identical(stale, expected)
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
})

test_that("progress labels stay short and ignore detail fields", {
    expect_identical(
        bam_progress_label("metrics", list(sample = "A", total = 1000)),
        "Counting reads"
    )
    expect_identical(
        bam_progress_label(
            "score", list(mode = "single", metric = "breadth")
        ),
        "Scoring windows"
    )
    expect_identical(
        bam_progress_label("mask", list(masked_fraction = 0.1)),
        "Masking windows"
    )
    expect_identical(
        bam_progress_label("calibration", list(n_calls = 3)),
        "Calibrating calls"
    )
    expect_identical(
        bam_progress_label("segmentation", list(n_segments = 18)),
        "Segmenting profiles"
    )
})

test_that("progress callbacks can update an inline status", {
    seen <- list()
    callback <- bam_progress_callback(
        1L,
        handler = function(stage, detail, progress) {
            seen[[length(seen) + 1L]] <<- list(
                stage = stage, progress = progress
            )
        }
    )
    callback("preflight", list())
    callback("index", list())
    expect_length(seen, 2L)
    expect_identical(seen[[1L]]$stage, "preflight")
    expect_equal(
        seen[[2L]]$progress,
        2 / length(bam_progress_stages(1L))
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
    expect_match(bam, "mamamia-sample-picker", fixed = TRUE)
    expect_match(bam, "mamamia-call-note", fixed = TRUE)
    expect_match(filters, "Call class", fixed = TRUE)
    expect_match(filters, "deletion", fixed = TRUE)
    expect_match(filters, "complex", fixed = TRUE)
    expect_match(filters, "Class status", fixed = TRUE)
    expect_match(filters, 'value="confident" selected', fixed = TRUE)
    expect_match(filters, "Copy state", fixed = TRUE)
    expect_match(filters, "balanced", fixed = TRUE)
    expect_match(filters, "Confidence", fixed = TRUE)
})

test_that("the export panel offers one results bundle", {
    html <- as.character(explorer_export_ui("explorer"))
    expect_match(html, "Results bundle", fixed = TRUE)
    expect_match(html, "download_bundle", fixed = TRUE)
    expect_match(html, "Download results bundle", fixed = TRUE)
    expect_match(html, "Selected calls", fixed = TRUE)
    expect_match(html, "download_bedpe", fixed = TRUE)
})

test_that("the empty sample picker is not hidden from Shiny", {
    ## Shiny suspends outputs whose wrapper has display:none, so hiding the
    ## empty picker would prevent it from ever rendering a multi-sample
    ## selector. The wrapper must stay visible and collapse its spacing.
    css <- bam_dashboard_css()
    expect_false(grepl(
        ".mamamia-sample-picker:empty", css, fixed = TRUE
    ))
    expect_match(
        css, ".mamamia-sample-picker .shiny-input-container", fixed = TRUE
    )
    expect_match(
        css, "margin-right: .55rem", fixed = TRUE
    )
})

test_that("notifications match bslib's specificity to stay compact", {
    ## bslib ships doubled class and id selectors for the toast panel and
    ## cards; a single-class rule is overridden and the popup grows again.
    css <- bam_dashboard_css()
    expect_match(
        css, "#shiny-notification-panel#shiny-notification-panel",
        fixed = TRUE
    )
    expect_match(
        css, ".shiny-notification.shiny-notification", fixed = TRUE
    )
    expect_match(
        css,
        ".shiny-notification.shiny-notification .shiny-progress-notification",
        fixed = TRUE
    )
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
    result <- bam_example_result(1000)
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
