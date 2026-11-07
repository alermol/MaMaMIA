## Shiny launcher and result helpers.

test_that("bam_run_parameters() resolves the run controls", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
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
    expect_identical(params$cores, mamamia_default_cores())
    expect_identical(
        params$donor_aggregation, config$defaults$donor_aggregation
    )
    expect_named(params, c(
        "resolution", "mapq", "exclude_flags", "cores", "metric", "score",
        "donor_aggregation", "control_samples", "calibration",
        "profile_name", "profile_base", "profile_parameters",
        "profile_changed"
    ))
    expect_identical(params$control_samples, character(0))
    expect_null(params$calibration)
    expect_identical(params$profile_name, "custom")
    expect_null(params$mask)
})

test_that("fractional Advanced integers are rejected", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
    )
    expect_error(
        bam_run_parameters(
            list(mapq = 20.9), config, samples = "demo-sample"
        ),
        "whole number"
    )
    expect_error(
        bam_run_parameters(
            list(calibration_merge_gap = 100.9), config,
            samples = "demo-sample"
        ),
        "whole number"
    )
    params <- bam_run_parameters(
        list(mapq = 20, calibration_merge_gap = 100), config,
        samples = "demo-sample"
    )
    expect_identical(params$mapq, 20L)
})

test_that("analysis profiles resolve complete parameter bundles", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
    )
    profiles <- mamamia_analysis_profiles()
    expect_false("component_panel" %in% names(profiles))
    expect_false(any(vapply(
        profiles,
        function(profile) "component_panel" %in% names(profile),
        logical(1)
    )))
    sensitive <- mamamia_analysis_profile("sensitive")
    expect_identical(sensitive$parameters$windows_per_chromosome, 100L)
    expect_identical(sensitive$parameters$metric, "combined")
    expect_identical(sensitive$parameters$score, "paired_log_ratio")
    expect_identical(sensitive$parameters$donor_aggregation, "max")
    expect_equal(
        sensitive$parameters$calibration$floors,
        c(1.5, 1.2, 1.0, 0.8, 0.4)
    )
    expect_false("mask" %in% names(sensitive$parameters))
    params <- bam_run_parameters(
        list(
            analysis_profile = "sensitive",
            control_samples = character(0)
        ),
        config, samples = "demo-sample"
    )
    expect_identical(params$resolution, 100L)
    expect_identical(params$metric, "combined")
    expect_identical(params$profile_name, "sensitive")
    expect_identical(params$profile_base, "sensitive")
    choices <- mamamia_analysis_profile_choices()
    expect_identical(
        unname(choices[c("Balanced", "Sensitive", "Custom")]),
        c("balanced", "sensitive", "custom")
    )
    expect_null(mamamia_analysis_profile("custom"))
    expect_error(mamamia_analysis_profile("not-a-profile"), "Unknown analysis profile")
    html <- as.character(bam_sidebar_details_ui(config))
    expect_match(html, 'id="analysis_profile_ui"', fixed = TRUE)
    expect_match(html, 'id="advanced_parameters_ui"', fixed = TRUE)

    changed <- bam_run_parameters(
        list(
            analysis_profile = "sensitive",
            metric = "breadth",
            control_samples = character(0)
        ),
        config, samples = "demo-sample",
        profile_base = "sensitive"
    )
    expect_identical(changed$profile_name, "custom")
    expect_identical(changed$profile_base, "sensitive")
    expect_identical(changed$metric, "breadth")
    expect_identical(changed$profile_changed, "metric")
    expect_identical(changed$resolution, 100L)
})

test_that("profile choices follow the configured resolutions", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
    )
    choices <- bam_analysis_profile_choices(config)
    expect_identical(
        unname(choices[c("Balanced", "Sensitive")]),
        c("balanced", "sensitive")
    )
    expect_identical(bam_default_profile(config), "balanced")

    narrow <- config
    narrow$defaults$allowed_windows <- 1000L
    narrow$defaults$windows_per_chromosome <- 1000L
    expect_identical(
        unname(bam_analysis_profile_choices(narrow)),
        c("balanced", "complex", "custom")
    )
    expect_error(
        bam_run_parameters(
            list(analysis_profile = "sensitive"),
            narrow, samples = "demo-sample"
        ),
        "this configuration allows"
    )

    only100 <- config
    only100$defaults$allowed_windows <- 100L
    only100$defaults$windows_per_chromosome <- 100L
    expect_identical(
        unname(bam_analysis_profile_choices(only100)),
        c("sensitive", "cnv_robust", "custom")
    )
    expect_identical(bam_default_profile(only100), "sensitive")

    none <- config
    none$defaults$allowed_windows <- 200L
    none$defaults$windows_per_chromosome <- 200L
    expect_identical(
        unname(bam_analysis_profile_choices(none)),
        "custom"
    )
    expect_identical(bam_default_profile(none), "custom")
    params <- bam_run_parameters(
        list(control_samples = character(0)), none, samples = "demo-sample"
    )
    expect_identical(params$resolution, 200L)
    expect_identical(params$profile_name, "custom")
})

test_that("Advanced exposes every profile parameter and marks defaults", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
    )
    defaults <- bam_profile_default_parameters(config, "sensitive")
    html <- as.character(bam_advanced_parameters_ui(config, defaults))
    ids <- c(
        "resolution_ui", "metric", "score", "donor_aggregation", "mapq",
        "calibration_fdr", paste0("calibration_floor_", 1:5),
        "calibration_min_mb", "calibration_merge_gap"
    )
    for (id in ids) {
        expect_match(html, sprintf('id="%s"', id), fixed = TRUE)
    }
    expect_match(
        html, "Cohort recurrence cutoff is derived automatically",
        fixed = TRUE
    )
    expect_false(grepl("mask_recurrence_threshold", html, fixed = TRUE))
    expect_match(html, "mamamia-param-floors", fixed = TRUE)
    expect_match(html, "mamamia-param-legend", fixed = TRUE)
    expect_match(html, "reset_profile_parameters", fixed = TRUE)
    expect_match(html, "mamamia-param-actions", fixed = TRUE)
    help_label <- as.character(bam_profile_parameter_label(
        "mapq", "Minimum MAPQ", FALSE,
        help = bam_parameter_help("mapq")
    ))
    expect_match(help_label, "bslib-tooltip", fixed = TRUE)
    expect_match(help_label, "mamamia-param-help", fixed = TRUE)
    help_ids <- c(
        "resolution", "metric", "score", "donor_aggregation", "mapq",
        "calibration_fdr", "calibration_min_mb",
        "calibration_merge_gap", "control_samples"
    )
    for (id in help_ids) {
        expect_true(nzchar(bam_parameter_help(id)))
    }
    expect_match(
        bam_parameter_help("calibration_floor_2", "Floor 10-25 Mb"),
        "10-25",
        fixed = TRUE
    )

    default_label <- as.character(
        bam_profile_parameter_label("mapq", "Minimum MAPQ", FALSE)
    )
    expect_match(default_label, "mamamia-param-default", fixed = TRUE)
    expect_false(grepl("mamamia-param-changed", default_label, fixed = TRUE))
    changed_label <- as.character(
        bam_profile_parameter_label("mapq", "Minimum MAPQ", TRUE)
    )
    expect_match(changed_label, "mamamia-param-changed", fixed = TRUE)

    current <- defaults
    current$mapq <- 30L
    current$calibration_floors[[2L]] <- 1.1
    expect_identical(
        bam_profile_changed_parameters(current, defaults),
        c("mapq", "calibration_floors")
    )
})

test_that("every profile parameter can be overridden into a Custom run", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
    )
    params <- bam_run_parameters(
        list(
            analysis_profile = "sensitive",
            resolution = 1000,
            metric = "count",
            score = "donor_only",
            donor_aggregation = "mean",
            mapq = 30,
            calibration_fdr = 0.02,
            calibration_floor_1 = 1.4,
            calibration_floor_2 = 1.1,
            calibration_floor_3 = 0.9,
            calibration_floor_4 = 0.7,
            calibration_floor_5 = 0.3,
            calibration_min_mb = 7,
            calibration_merge_gap = 50,
            control_samples = character(0)
        ),
        config, samples = "demo-sample",
        profile_base = "sensitive"
    )
    expect_identical(params$profile_name, "custom")
    expect_identical(params$profile_base, "sensitive")
    expect_identical(params$resolution, 1000L)
    expect_identical(params$metric, "count")
    expect_identical(params$score, "donor_only")
    expect_identical(params$donor_aggregation, "mean")
    expect_identical(params$mapq, 30L)
    expect_equal(params$calibration$fdr, 0.02)
    expect_equal(
        params$calibration$floors,
        c(1.4, 1.1, 0.9, 0.7, 0.3)
    )
    expect_equal(params$calibration$min_mb, 7)
    expect_identical(params$calibration$merge_gap, 50L)
    expect_setequal(
        params$profile_changed,
        c(
            "windows_per_chromosome", "metric", "score",
            "donor_aggregation", "mapq", "calibration_fdr",
            "calibration_floors", "calibration_min_mb",
            "calibration_merge_gap"
        )
    )
})

test_that("the resolution picker is a ten-window slider with three ticks", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
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
        mamamia_path("data", "bam-demo", "assemblies.yml")
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

test_that("component evidence stays in the bundle, not the Calls UI", {
    html <- as.character(bam_results_ui())
    expect_false(grepl("Component evidence", html, fixed = TRUE))
    expect_false(grepl('id="component_evidence"', html, fixed = TRUE))
    css <- bam_dashboard_css()
    expect_match(
        css,
        ".mamamia-main .card-body > .datatables.html-widget",
        fixed = TRUE
    )
    expect_match(css, "flex: 0 0 auto !important", fixed = TRUE)
    expect_false(grepl(".mamamia-component-evidence", css, fixed = TRUE))
})

test_that("the BAM chooser accepts exactly one file", {
    chooser <- as.character(bam_single_file_input())
    expect_length(chooser, 1L)
    expect_match(chooser, 'id="bams"', fixed = TRUE)
    expect_match(chooser, 'accept=".bam"', fixed = TRUE)
    expect_match(chooser, "Choose BAM file", fixed = TRUE)
    expect_match(chooser, "No BAM selected", fixed = TRUE)
    expect_false(grepl("multiple", chooser, fixed = TRUE))
    expect_false(grepl("data-max-bytes", chooser, fixed = TRUE))
    expect_false(grepl("bam_upload_help", chooser, fixed = TRUE))
    expect_match(chooser, "mamamia-upload-timer", fixed = TRUE)
    expect_match(chooser, "mamamia-upload-timer-value", fixed = TRUE)
    expect_match(chooser, "Time to finish", fixed = TRUE)
    expect_match(chooser, "--:--", fixed = TRUE)

    limited <- as.character(bam_single_file_input(5 * 1024^3))
    expect_match(
        limited, 'data-max-bytes="5368709120"', fixed = TRUE
    )
    expect_match(limited, "Max 5.0 GiB", fixed = TRUE)
    expect_match(limited, 'id="bam_upload_help"', fixed = TRUE)

    help <- as.character(bam_upload_help_ui(5 * 1024^3))
    expect_match(help, "5.0 GiB", fixed = TRUE)
    expect_match(help, "Before upload", fixed = TRUE)
    expect_match(help, "50–90%", fixed = TRUE)
    expect_match(help, "bamhusk", fixed = TRUE)
    expect_match(help, "husk --drop-seq", fixed = TRUE)
    expect_match(help, "mamamia-help-command", fixed = TRUE)
    expect_match(help, "QUAL", fixed = TRUE)
    help_modal <- as.character(bam_upload_help_modal(5 * 1024^3))
    expect_match(help_modal, "BAM upload limit", fixed = TRUE)
    expect_match(help_modal, "mamamia-tooltip-modal", fixed = TRUE)

    css <- bam_dashboard_css()
    expect_match(css, ".mamamia-upload-timer.is-visible", fixed = TRUE)
    expect_match(
        css, ".mamamia-upload-timer.is-running .mamamia-upload-timer-value",
        fixed = TRUE
    )
    expect_match(css, "font-variant-numeric: tabular-nums", fixed = TRUE)
})

test_that("upload caps apply only to the ShinyProxy deployment", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
    )
    withr::local_envvar(
        MAMAMIA_DEPLOYMENT_MODE = "direct",
        MAMAMIA_MAX_REQUEST_SIZE = NA
    )
    direct <- bam_launcher_state(config)
    expect_true(direct$ok)
    expect_true(is.infinite(direct$config$server$max_total_bytes))
    expect_identical(mamamia_max_request_size(), Inf)

    withr::local_envvar(
        MAMAMIA_DEPLOYMENT_MODE = "shinyproxy",
        MAMAMIA_MAX_REQUEST_SIZE = "1234"
    )
    proxy <- bam_launcher_state(config)
    expect_true(proxy$ok)
    expect_identical(
        proxy$config$server$max_total_bytes,
        config$server$max_total_bytes
    )
    expect_identical(mamamia_max_request_size(), 1234)
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

test_that("oversized uploads are rejected with a bamhusk suggestion", {
    uploads <- data.frame(
        name = "huge.bam", size = 11, datapath = "huge.bam",
        stringsAsFactors = FALSE
    )
    error <- tryCatch(
        bam_upload_limit(uploads, 10),
        error = function(e) e
    )
    expect_s3_class(error, "MaMaMIBamTooLargeError")
    expect_match(conditionMessage(error), "this server accepts at most")
    info <- bam_error_info(error)
    expect_identical(info$stage, "upload")
    expect_identical(info$size, 11)
    expect_identical(info$limit, 10)
    expect_match(info$suggestion, "bamhusk", fixed = TRUE)
    expect_silent(bam_upload_limit(uploads, 11))

    html <- as.character(bam_sidebar_recovery_ui(info))
    expect_match(html, "BAM is too large", fixed = TRUE)
    expect_match(
        html, "https://github.com/alermol/bamhusk", fixed = TRUE
    )
    expect_match(html, "husk --drop-seq", fixed = TRUE)
    expect_false(grepl("unlock_setup", html, fixed = TRUE))

    modal <- as.character(bam_upload_too_large_modal(info))
    expect_match(modal, "BAM is too large to upload", fixed = TRUE)
    expect_match(
        modal, "https://github.com/alermol/bamhusk", fixed = TRUE
    )
    expect_match(modal, "husk --drop-seq", fixed = TRUE)
    expect_match(modal, "clear_oversized_upload", fixed = TRUE)
    expect_identical(bam_format_size(5 * 1024^3), "5.0 GiB")

    log_lines <- bam_preflight_log(list(
        ok = FALSE, empty = FALSE, error = info, candidates = NULL
    ))
    expect_match(paste(log_lines, collapse = "\n"), "bamhusk", fixed = TRUE)

    payload_info <- bam_upload_payload_info(
        list(name = "huge.bam", size = 11), 10
    )
    expect_identical(payload_info$stage, "upload")
    expect_identical(payload_info$size, 11)
    expect_null(bam_upload_payload_info(
        list(name = "huge.bam", size = 11), 11
    ))
    expect_null(bam_upload_payload_info(
        list(name = "huge.bam", size = "not-a-number"), 10
    ))
})

test_that("donor aggregation defaults from the registry and is set under Advanced", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
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
    expect_match(html, 'id="advanced_parameters_ui"', fixed = TRUE)
    expect_match(html, "Advanced", fixed = TRUE)
})

test_that("control samples are selected under Advanced", {
    config <- readMaMaMIAConfig(
        mamamia_path("data", "bam-demo", "assemblies.yml")
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
        mamamia_path("data", "bam-demo", "assemblies.yml")
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
        mamamia_path("data", "bam-demo", "assemblies.yml")
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
        mamamia_path("data", "bam-demo", "assemblies.yml")
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

cohort_fixture <- function() {
    list(
        samples = data.frame(
            sample = c("A", "B"),
            n_windows = c(2L, 2L),
            n_calls = c(1L, 0L),
            n_high_confidence = c(1L, 0L),
            masked_fraction = c(0.1, 0.2),
            mode = c("multi", "multi"),
            stringsAsFactors = FALSE
        ),
        metrics = data.frame(
            sample = c("A", "A", "B", "B"),
            role = c("donor", "recipient", "donor", "recipient"),
            count = c(10L, 20L, 0L, 0L),
            stringsAsFactors = FALSE
        ),
        windows = list(table = data.frame(
            chr_id = c("r1", "r2", "d1", "d2"),
            role = c("recipient", "recipient", "donor", "donor"),
            group = c("1", "2", "1", "2"),
            stringsAsFactors = FALSE
        )),
        calls = data.frame(
            sample = c("A", "A", "B"),
            pass = c(TRUE, FALSE, FALSE),
            class = c("substitution", "other", "deletion"),
            class_status = c("confident", "needs_review", "needs_review"),
            confidence = c("high", "low", "low"),
            group = c("1", "1", "1"),
            width_bp = c(1e6, 2e6, 3e6),
            score_mean = c(1, 2, 3),
            stringsAsFactors = FALSE
        )
    )
}

test_that("cohort summaries and warnings read WindowMetrics", {
    result <- cohort_fixture()
    summary <- bam_cohort_summary(result)
    expect_identical(summary$reads_total, c(30, 0))
    expect_identical(summary$n_calls, c(1L, 0L))
    expect_equal(summary$reads_total_per_window, c(15, 0))
    expect_equal(summary$donor_read_share, c(1 / 3, NA))
    expect_equal(summary$masked_windows, c(0.2, 0.4))
    expect_identical(summary$n_all_calls, c(2L, 1L))
    expect_identical(summary$n_substitution, c(1L, 0L))
    expect_identical(summary$n_deletion, c(0L, 0L))
    expect_identical(summary$n_class_confident, c(1L, 0L))
    expect_identical(summary$n_class_needs_review, c(1L, 1L))
    expect_equal(summary$pass_rate, c(0.5, 0))
    expect_equal(summary$total_width_mb, c(1, NA))
    expect_equal(summary$mean_width_mb, c(1, NA))
    expect_equal(summary$max_width_mb, c(1, NA))
    expect_identical(summary$n_contigs, c(4L, 4L))
    expect_identical(summary$n_groups, c(2L, 2L))
    expect_identical(summary$n_groups_with_calls, c(1L, 0L))
    warnings <- bam_coverage_warnings(result)
    expect_true(any(grepl("Sample 'B' has no reads", warnings)))
    expect_false(any(grepl("Sample 'A' has no reads", warnings)))
})

test_that("every cohort column has a label, description and format", {
    summary <- bam_cohort_summary(cohort_fixture())
    metadata <- bam_cohort_column_metadata(names(summary))
    expect_identical(nrow(metadata), ncol(summary))
    expect_identical(anyDuplicated(metadata$column), 0L)
    expect_identical(anyDuplicated(metadata$label), 0L)
    expect_false(anyNA(metadata$label))
    expect_false(any(!nzchar(metadata$label)))
    expect_false(anyNA(metadata$description))
    expect_false(any(!nzchar(metadata$description)))
    expect_true(all(metadata$format %in% c(
        "text", "integer", "decimal1", "decimal2", "percent", "signif"
    )))
    expect_identical(
        metadata$description[metadata$column == "sample"],
        "Sample name from the BAM @RG SM tag."
    )
    expect_identical(
        metadata$description[metadata$column == "n_calls"],
        "Calls that pass calibration filters."
    )
    choices <- bam_cohort_extra_column_choices()
    expect_true(all(nzchar(names(choices))))
    expect_true("All calls" %in% names(choices))
    expect_true("n_all_calls" %in% unname(choices))
    expect_false(any(unname(choices) %in% bam_cohort_brief_columns()))
    expect_false("sample" %in% unname(choices))
    header <- as.character(bam_cohort_table_header(c("sample", "n_calls")))
    expect_match(
        header, "Sample name from the BAM @RG SM tag.", fixed = TRUE
    )
    expect_match(
        header, "Calls that pass calibration filters.", fixed = TRUE
    )
    expect_match(header, 'class="mamamia-call-label"', fixed = TRUE)
})

test_that("the cohort table supports tooltips, sorting and extra columns", {
    skip_if_not_installed("DT")
    summary <- bam_cohort_summary(cohort_fixture())
    table <- bam_cohort_table(
        summary, c("n_all_calls", "reads_total_per_window")
    )
    expect_s3_class(table, "datatables")
    expect_true(all(c(
        "n_all_calls", "reads_total_per_window"
    ) %in% names(table$x$data)))
    expect_true(isTRUE(table$x$options$orderMulti))
    definitions <- Filter(
        function(x) !is.null(x$orderSequence),
        table$x$options$columnDefs
    )
    expect_length(definitions, 1L)
    expect_identical(
        definitions[[1]]$orderSequence, c("asc", "desc", "")
    )
    expect_identical(definitions[[1]]$targets, "_all")
    expect_match(
        as.character(table$x$options$initComplete), "order([])",
        fixed = TRUE
    )
    header <- as.character(table$x$container)
    expect_match(header, "Passing calls", fixed = TRUE)
    expect_match(
        header, "Calls that pass calibration filters.", fixed = TRUE
    )
    expect_match(header, "Reads/window", fixed = TRUE)
})

test_that("the cohort columns picker is available in the Export tab", {
    ui <- as.character(bam_cohort_columns_ui())
    expect_match(ui, "Additional columns", fixed = TRUE)
    expect_match(ui, "Add columns", fixed = TRUE)
    expect_match(ui, 'id="cohort_columns"', fixed = TRUE)
    export <- as.character(bam_results_ui())
    expect_match(export, "cohort_ui", fixed = TRUE)
    expect_match(export, "mamamia-export-panel", fixed = TRUE)
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
    expect_false(grepl('value="confident" selected', filters, fixed = TRUE))
    expect_match(
        filters, 'id="explorer-pass_only".*checked'
    )
    expect_match(filters, "Copy state", fixed = TRUE)
    expect_match(filters, "balanced", fixed = TRUE)
    expect_match(filters, "Confidence", fixed = TRUE)
})

test_that("the explorer shows passing calls by default", {
    view <- as.character(explorer_view_filters_ui("explorer"))
    width <- as.character(explorer_call_filters_ui("explorer"))
    expect_false(grepl('value="confident" selected', view, fixed = TRUE))
    expect_match(view, 'id="explorer-pass_only".*checked')
    expect_match(width, 'value="0"', fixed = TRUE)
    expect_match(width, "Additional columns", fixed = TRUE)
    expect_match(width, "Add columns", fixed = TRUE)
})

test_that("the call table header carries column tooltips", {
    header <- as.character(bam_call_table_header(c(
        "group", "donor_chr_name", "recipient_mean_depletion"
    )))
    expect_match(
        header,
        "Chromosome homology group used for scoring",
        fixed = TRUE
    )
    expect_match(header, "Display name of the donor contig.", fixed = TRUE)
    expect_match(
        header,
        "Mean corrected log-depth depletion",
        fixed = TRUE
    )
    expect_match(header, 'class="mamamia-call-label"', fixed = TRUE)
})

test_that("the export panel offers one results bundle", {
    html <- as.character(explorer_export_ui("explorer"))
    expect_match(html, "mamamia-export-actions", fixed = TRUE)
    expect_match(html, "mamamia-export-card", fixed = TRUE)
    expect_match(html, "mamamia-download-row", fixed = TRUE)
    expect_match(html, "mamamia-download-spinner", fixed = TRUE)
    expect_match(html, "Preparing bundle...", fixed = TRUE)
    expect_match(html, "Results bundle", fixed = TRUE)
    expect_match(html, "download_bundle", fixed = TRUE)
    expect_match(html, "Download results bundle", fixed = TRUE)
    expect_match(html, "one folder per sample", fixed = TRUE)
    expect_match(html, "All calls", fixed = TRUE)
    expect_match(html, "download_all_bedpe", fixed = TRUE)
    expect_match(html, "Download all BEDPE", fixed = TRUE)
    expect_match(html, "Selected calls", fixed = TRUE)
    expect_match(html, "download_bedpe", fixed = TRUE)
    css <- bam_dashboard_css()
    expect_match(css, ".mamamia-download-row", fixed = TRUE)
    expect_match(css, ".mamamia-download-spinner", fixed = TRUE)
    expect_match(css, ".mamamia-download-status.is-visible", fixed = TRUE)
    expect_match(css, "@keyframes mamamia-download-spin", fixed = TRUE)
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
    config <- mamamia_path("data", "bam-demo", "assemblies.yml")
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

test_that("the result cache keeps a fixed session budget", {
    cache <- bam_result_cache(max_entries = 2L)
    cache$set("a", 1)
    cache$set("b", 2)
    cache$set("c", 3)
    expect_identical(cache$limit(), 2L)
    expect_identical(cache$size(), 2L)
    expect_null(cache$get("a"))
    expect_identical(cache$get("c"), 3)
})

test_that("species auto-detection matches the BAM header", {
    path <- test_path("fixtures", "tiny.bam")
    config <- list(species = list(
        tiny = list(
            label = "Tiny species",
            contigs = data.frame(
                id = c("d1", "r1"), length = c(20000, 20000),
                stringsAsFactors = FALSE
            )
        ),
        other = list(
            label = "Other species",
            contigs = data.frame(
                id = "absent", length = 1,
                stringsAsFactors = FALSE
            )
        )
    ))
    candidates <- bam_species_candidates(path, config)
    expect_true("tiny" %in% candidates$key)
    expect_true(candidates$exact[candidates$key == "tiny"])
    expect_false(candidates$exact[candidates$key == "other"])
    detected <- bam_detect_species(path, config)
    expect_identical(detected$species, "tiny")
    expect_identical(detected$label, "Tiny species")

    ambiguous <- config
    ambiguous$species$subset <- list(
        label = "Subset species",
        contigs = data.frame(
            id = "d1", length = 20000,
            stringsAsFactors = FALSE
        )
    )
    expect_identical(
        bam_detect_species(path, ambiguous)$species,
        "tiny"
    )
    ambiguous$species$same <- ambiguous$species$tiny
    expect_null(bam_detect_species(path, ambiguous))
})

test_that("species pair auto-detection prefers the current donor role", {
    config <- bam_example_config()
    recipient <- config$species$triticum_like$contigs
    donor <- config$species$secale_like$contigs
    targets <- data.frame(
        name = c(recipient$id, donor$id),
        length = c(recipient$length, donor$length),
        stringsAsFactors = FALSE
    )
    candidates <- bam_species_pair_candidate_table(targets, config)
    expect_true(any(candidates$exact))
    expect_identical(sum(candidates$exact), 2L)
    chosen <- bam_choose_species_pair(
        candidates[candidates$exact, , drop = FALSE],
        current_recipient = "triticum_like",
        current_donor = "secale_like"
    )
    expect_identical(as.character(chosen$recipient), "triticum_like")
    expect_identical(as.character(chosen$donor), "secale_like")

    swapped <- bam_choose_species_pair(
        candidates[candidates$exact, , drop = FALSE],
        current_recipient = "secale_like",
        current_donor = "triticum_like"
    )
    expect_identical(as.character(swapped$recipient), "secale_like")
    expect_identical(as.character(swapped$donor), "triticum_like")
})

test_that("single-sample runs get a cohort-mask warning", {
    warning_ui <- as.character(bam_cohort_mask_warning_ui(1L))
    expect_match(warning_ui, "Cohort recurrence unavailable", fixed = TRUE)
    expect_match(warning_ui, "at least two samples", fixed = TRUE)
    expect_null(bam_cohort_mask_warning_ui(2L))
    expect_null(bam_cohort_mask_warning_ui(0L))
    expect_null(bam_cohort_mask_warning_ui(NA_integer_))
})
