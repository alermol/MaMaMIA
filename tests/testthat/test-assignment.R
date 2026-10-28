## Per-call assignment of donor segments to recipient contigs.

assignment_result <- function(r1_breadth = 0.0009, r2_breadth = 0.02) {
    result <- bam_example_result(100)
    extra <- result$assembly$contigs$recipient[1, , drop = FALSE]
    extra$id <- "r2"
    extra$name <- "1B"
    extra$subgenome <- "B"
    result$assembly$contigs$recipient <- rbind(
        result$assembly$contigs$recipient, extra
    )
    metric <- result$metrics[result$metrics$role == "recipient", ,
        drop = FALSE
    ][1, , drop = FALSE]
    metric$chr_id <- "r2"
    metric$chr_name <- "1B"
    metric$subgenome <- "B"
    metric$breadth <- r2_breadth
    metric$lambda_breadth <- -log(1 - r2_breadth)
    result$metrics <- rbind(result$metrics, metric)
    result$metrics$breadth[result$metrics$chr_id == "r1"] <- r1_breadth
    result$metrics$lambda_breadth[
        result$metrics$chr_id == "r1"
    ] <- -log(1 - r1_breadth)

    call <- result$calls[2, , drop = FALSE]
    call$start_iid <- 0L
    call$end_iid <- 99L
    call$start_bp <- 0
    call$end_bp <- 20000
    call$width_bp <- 20000
    call$n_windows <- 100L
    call$class <- "substitution"
    call$confidence <- "high"
    call$pass <- TRUE
    result$calls <- call
    result
}

test_that("the example calls are assigned to the only recipient chromosome", {
    result <- bam_example_result(100)
    calls <- assignCallPositions(
        result$calls, result$metrics, result$windows,
        result$scores, result$assembly
    )
    expect_true(all(!is.na(calls$recipient_chr_id)))
    expect_identical(unique(calls$recipient_subgenome), "A")
    expect_identical(unique(calls$recipient_chr_name), "1A")
    expect_identical(unique(calls$donor_subgenome), "R")
    expect_identical(unique(calls$donor_chr_name), "1R")
    expect_true(all(calls$recipient_start_bp >= 0))
    expect_true(all(calls$recipient_end_bp > calls$recipient_start_bp))
    expect_true(all(calls$recipient_assignment == "assigned"))
})

test_that("a lost recipient chromosome is the most probable host", {
    result <- assignment_result()
    calls <- assignCallPositions(
        result$calls, result$metrics, result$windows,
        result$scores, result$assembly
    )
    expect_identical(calls$recipient_subgenome, "A")
    expect_identical(calls$recipient_chr_name, "1A")
    expect_identical(calls$recipient_assignment, "assigned")
    expect_gt(calls$recipient_margin, 0)
    expect_identical(calls$donor_subgenome, "R")
})

test_that("a copy-neutral event is not given a confident host", {
    result <- assignment_result(r1_breadth = 0.01, r2_breadth = 0.01)
    calls <- assignCallPositions(
        result$calls, result$metrics, result$windows,
        result$scores, result$assembly
    )
    expect_identical(calls$recipient_assignment, "no_recipient_loss")
    expect_identical(calls$recipient_weight, 0.5, tolerance = 1e-6)
})
