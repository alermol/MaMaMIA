## The bundled example must behave like a BAM run without touching a BAM.

skip_if_no_example <- function() {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("bslib")
    testthat::skip_if_not_installed("DT")
    testthat::skip_if_not_installed("gridlayout")
    testthat::skip_if_not_installed("shinyWidgets")
    if (!MaMaMIA:::bam_example_available()) {
        testthat::skip("bundled example data is not installed")
    }
}

demo_config_path <- function() {
    demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
    file.path(demo, "assemblies.yml")
}

output_html <- function(output) {
    paste(as.character(output), collapse = "")
}

example_inputs <- function(...) {
    defaults <- list(
        resolution = "100",
        reference_sample = ""
    )
    utils::modifyList(defaults, list(...))
}

test_that("the bundled example ships pre-computed results only", {
    skip_if_no_example()

    expect_true(MaMaMIA:::bam_example_available())
    expect_identical(MaMaMIA:::bam_example_resolutions(), c(100L, 1000L))
    expect_identical(MaMaMIA:::bam_example_default_resolution(), "1000")

    config <- MaMaMIA:::bam_example_config()
    expect_identical(
        sort(unique(as.integer(config$defaults$allowed_windows))),
        c(100L, 1000L)
    )
    expect_false("bam_root" %in% names(config$server))

    state <- MaMaMIA:::bam_example_preflight_state()
    expect_true(state$ok)
    expect_true(state$example)
    expect_identical(state$table$sample, "demo-sample")
    expect_identical(state$table$path, "demo.bam")
    expect_identical(state$table$index_action, "reuse")
    expect_null(state$staged)
    expect_false(file.exists(state$table$path))

    result <- MaMaMIA:::bam_example_result(100L)
    expect_s3_class(result, "BamAnalysisResult")
    expect_identical(result$windows$n_per_chromosome, 100L)
    expect_identical(nrow(result$windows$table), 200L)
    expect_identical(result$samples$sample, "demo-sample")
    expect_identical(result$samples$path, "demo.bam")
    expect_false(file.exists(result$samples$path))
})

test_that("the example rejects resolutions it does not ship", {
    skip_if_no_example()
    expect_error(
        MaMaMIA:::bam_example_result(10000L),
        "available at resolution",
        class = "MaMaMIAWindowError"
    )
})

test_that("the example replay uses the normal pipeline stage messages", {
    skip_if_no_example()
    result <- MaMaMIA:::bam_example_result(100L)
    seen <- list()
    MaMaMIA:::bam_example_progress(result, function(stage, detail) {
        seen[[length(seen) + 1L]] <<- list(stage = stage, detail = detail)
    })

    expect_identical(
        vapply(seen, `[[`, character(1), "stage"),
        MaMaMIA:::bam_progress_stages(1L, "none")
    )
    expect_match(
        MaMaMIA:::bam_progress_label("preflight", seen[[1L]]$detail),
        "Checking BAM headers \\(1 BAM file, 1 sample"
    )
    expect_match(
        MaMaMIA:::bam_progress_label("counts", seen[[5L]]$detail),
        "Counting reads for demo-sample"
    )
    expect_match(
        MaMaMIA:::bam_progress_label("correction", seen[[6L]]$detail),
        "Correcting GC bias for demo-sample"
    )
    expect_match(
        MaMaMIA:::bam_progress_label("segmentation", seen[[7L]]$detail),
        "Segmenting chromosome pairs for demo-sample"
    )
})

test_that("the server runs the bundled example without an upload", {
    skip_if_no_example()
    state <- MaMaMIA:::bam_launcher_state(demo_config_path())
    expect_true(state$ok)

    shiny::testServer(MaMaMIA:::bam_launcher_server(state, cache = FALSE), {
        do.call(session$setInputs, example_inputs(example_data = 1L))
        expect_true(example_mode())
        expect_true(sidebar_unlocked())
        pf <- preflight_state()
        expect_true(pf$ok)
        expect_true(pf$example)
        expect_identical(pf$table$sample, "demo-sample")

        session$setInputs(resolution = "100")
        expect_true(preflight_state()$ok)
        suppressWarnings(session$setInputs(run = 1L))

        result <- result()
        expect_s3_class(result, "BamAnalysisResult")
        expect_identical(result$windows$n_per_chromosome, 100L)
        expect_match(output$cohort_table, "demo-sample")
        expect_match(output$cohort_table, "none")
        log <- output$calculation_log
        expect_match(log, "Counting reads for demo-sample", fixed = TRUE)
        expect_match(log, "Analysis finished", fixed = TRUE)
    })
})

test_that("the Example data button unlocks the setup and writes the footer log", {
    skip_if_no_example()
    state <- MaMaMIA:::bam_launcher_state(demo_config_path())

    shiny::testServer(MaMaMIA:::bam_launcher_server(state, cache = FALSE), {
        expect_false(sidebar_unlocked())
        expect_null(output$sidebar_recovery)

        session$setInputs(example_data = 1L)
        expect_true(example_mode())
        expect_true(sidebar_unlocked())
        expect_true(preflight_state()$ok)

        sidebar <- output_html(output$sidebar_body)
        expect_match(sidebar, "Analysis")
        expect_match(sidebar, "Assembly")
        expect_false(grepl("Advanced", sidebar, fixed = TRUE))
        expect_match(sidebar, "Filters")
        expect_match(sidebar, "Export")
        species <- output_html(output$species_ui)
        expect_match(species, "Recipient species")
        expect_match(species, "Donor species")
        expect_match(species, "triticum_like")

        log <- output$calculation_log
        expect_match(log, "Preflight OK", fixed = TRUE)
        expect_match(log, "Detected assembly", fixed = TRUE)

        session$setInputs(example_data = 2L)
        expect_false(example_mode())
        expect_false(sidebar_unlocked())
    })
})

test_that("the example resolution guard disables unsupported runs", {
    skip_if_no_example()
    state <- MaMaMIA:::bam_launcher_state(demo_config_path())

    shiny::testServer(MaMaMIA:::bam_launcher_server(state, cache = FALSE), {
        do.call(session$setInputs, example_inputs(example_data = 1L))
        expect_true(preflight_state()$ok)
        session$setInputs(resolution = "10000")

        html <- paste(as.character(output$run_ui), collapse = "")
        expect_match(html, "supports resolution")
        expect_match(html, "disabled")
    })
})
