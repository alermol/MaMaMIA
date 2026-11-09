## Bundle QC schema validation.

bundle_qc_tables <- function(result, sample) {
    metrics <- result$metrics[
        as.character(result$metrics$sample) == sample, , drop = FALSE
    ]
    scores <- result$scores[
        as.character(result$scores$sample) == sample, , drop = FALSE
    ]
    list(
        metrics = bam_bundle_rename_columns(
            metrics, bam_bundle_metric_columns()
        ),
        scores = bam_bundle_rename_columns(
            scores, bam_bundle_score_columns()
        ),
        masks = bam_bundle_rename_columns(
            result$masks, bam_bundle_mask_columns()
        ),
        calibration = bam_bundle_rename_columns(
            result$calibration, bam_bundle_calibration_columns()
        )
    )
}

test_that("the bundle QC schema is derived from the result contract", {
    schema <- bam_bundle_qc_schema()
    contracts <- bam_result_contracts()
    expect_setequal(names(schema), c(
        "metrics", "scores", "masks", "calibration"
    ))
    for (name in names(schema)) {
        spec <- schema[[name]]
        expect_identical(
            spec$required,
            unname(spec$rename[contracts[[name]]])
        )
        expect_identical(anyDuplicated(spec$required), 0L)
        expect_identical(anyDuplicated(names(spec$rename)), 0L)
        expect_identical(
            spec$optional,
            setdiff(unname(spec$rename), spec$required)
        )
        expect_length(intersect(spec$required, spec$optional), 0L)
        expect_identical(
            bam_bundle_optional_columns(name), spec$optional
        )
    }
    expect_true(length(bam_bundle_optional_columns("scores")) > 0L)
})

test_that("optional QC columns may be absent without failing validation", {
    result <- bam_example_result(1000)
    sample <- unique(as.character(result$calls$sample))[[1L]]
    tables <- bundle_qc_tables(result, sample)
    optional <- bam_bundle_optional_columns("scores")
    tables$scores[optional] <- NULL
    expect_invisible(bam_bundle_validate_qc_tables(tables))
})

test_that("valid and genuinely empty QC tables pass validation", {
    result <- bam_example_result(1000)
    sample <- unique(as.character(result$calls$sample))[[1L]]
    tables <- bundle_qc_tables(result, sample)
    expect_invisible(bam_bundle_validate_qc_tables(tables))

    empty <- lapply(tables, function(table) table[0, , drop = FALSE])
    expect_invisible(bam_bundle_validate_qc_tables(empty))
    expect_invisible(bam_bundle_validate_qc_tables(
        empty[c("metrics", "scores")],
        include_masks = FALSE,
        include_calibration = FALSE
    ))
})

test_that("a missing QC column fails with the table and column name", {
    result <- bam_example_result(1000)
    sample <- unique(as.character(result$calls$sample))[[1L]]
    tables <- bundle_qc_tables(result, sample)
    tables$metrics$breadth <- NULL
    expect_error(
        bam_bundle_validate_qc_tables(tables),
        "metrics.*breadth"
    )

    tables <- bundle_qc_tables(result, sample)
    tables$scores <- NULL
    expect_error(
        bam_bundle_validate_qc_tables(tables),
        "scores.*missing"
    )
})

test_that("QC contract drift fails before the archive is created", {
    result <- bam_example_result(1000)
    sample <- unique(as.character(result$calls$sample))[[1L]]
    result$metrics$read_length <- NULL
    archive <- tempfile(fileext = ".zip")
    on.exit(unlink(archive), add = TRUE)
    expect_error(
        bam_write_result_bundle(result, archive, sample = sample),
        "metrics.*read_length_bp"
    )
    expect_false(file.exists(archive))
})
