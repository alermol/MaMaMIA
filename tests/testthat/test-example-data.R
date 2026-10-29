## The bundled example is a pre-computed BamAnalysisResult.

test_that("the bundled example is available", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    config <- bam_example_config()
    expect_identical(config$version, 2L)
    expect_true("triticum_like" %in% names(config$species))
    expect_identical(bam_example_resolutions(), 100L)
    expect_identical(bam_example_default_resolution(), "100")
})

test_that("the bundled example replays a result and progress", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    result <- bam_example_result(100)
    expect_s3_class(result, "BamAnalysisResult")
    expect_named(result, c(
        "assembly", "config", "bams", "samples", "windows", "gc", "metrics",
        "scores", "masks", "calibration", "calls", "provenance"
    ))
    expect_gte(sum(result$calls$pass), 1L)
    expect_true(any(result$calls$recipient_assignment == "assigned"))
    expect_true(any(result$calls$class == "substitution"))
    stages <- character()
    bam_example_progress(
        result,
        function(stage, detail) stages <<- c(stages, stage),
        pause = 0
    )
    expect_true(all(c(
        "preflight", "index", "windows", "gc", "metrics", "score", "mask",
        "segmentation", "calibration"
    ) %in% stages))
})

test_that("the example preflight state resolves the assembly", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    state <- bam_example_preflight_state()
    expect_true(state$ok)
    expect_true(state$example)
    expect_identical(state$assembly, "triticum_like__A__secale_like__R")
    expect_true(nrow(state$table) >= 1L)
    expect_match(bam_example_note(), "Example dataset")
})

test_that("the bundled example satisfies the result contract", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    result <- bam_example_result(100)
    expect_silent(validateBamAnalysisResult(result))
    expect_true(all(c("start", "stop") %in% names(result$scores)))
    expect_true(all(result$scores$stop > result$scores$start))
    calls <- result$calls[result$calls$pass %in% TRUE, , drop = FALSE]
    expect_true(all(c(
        "recipient_observed_breadth", "recipient_copy_number",
        "affected_fraction", "affected_copies", "dosage_confidence"
    ) %in% names(calls)))
    expect_true(any(!is.na(calls$affected_copies)))
})
