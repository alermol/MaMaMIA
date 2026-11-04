## BamAnalysisResult contract validation.

test_that("validation requires the calibrated call columns", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    result <- bam_example_result(1000)

    calls <- result$calls
    calls$pass <- NULL
    broken <- result
    broken$calls <- calls
    expect_error(validateBamAnalysisResult(broken), "pass")

    metrics <- result$metrics
    metrics$copy_number <- NULL
    broken <- result
    broken$metrics <- metrics
    expect_error(validateBamAnalysisResult(broken), "copy_number")

    scores <- result$scores
    scores$score_type <- NULL
    broken <- result
    broken$scores <- scores
    expect_error(validateBamAnalysisResult(broken), "score_type")
})
