## P5: cohort recurrence window masks.

multiply_donor <- function(metrics, group, from, to, multiplier,
                           samples = NULL) {
    idx <- metrics$role == "donor" & metrics$group == group &
        metrics$iid >= from & metrics$iid <= to
    if (!is.null(samples)) {
        idx <- idx & metrics$sample %in% samples
    }
    metrics$lambda_breadth[idx] <- metrics$lambda_breadth[idx] * multiplier
    metrics$lambda_count[idx] <- metrics$lambda_count[idx] * multiplier
    metrics$breadth[idx] <- 1 - exp(-metrics$lambda_breadth[idx])
    metrics$covered[idx] <- pmin(
        metrics$length[idx], round(metrics$breadth[idx] * metrics$length[idx])
    )
    metrics
}

test_that("a cohort mask removes recurrent windows", {
    metrics <- synthetic_window_metrics(
        samples = paste0("P", seq_len(12)), seed = 32
    )
    metrics <- multiply_donor(
        metrics, "1", 0L, 20L, 6, samples = paste0("P", 1:8)
    )
    metrics <- multiply_donor(
        metrics, "2", 100L, 130L, 6, samples = "P1"
    )
    scores <- scoreIntrogression(metrics, mode = "multi")
    mask <- buildWindowMask(
        metrics, scores = scores, mode = "multi"
    )
    expect_true(any(mask$mask[mask$group == "1" & mask$iid <= 20]))
    expect_false(any(mask$mask[mask$group == "2" & mask$iid >= 100 &
        mask$iid <= 130]))
    expect_identical(unique(mask$source), "cohort")
    expect_true(nzchar(unique(mask$version)))
    expect_identical(attr(mask, "threshold_source"), "automatic")
    expect_true(is.finite(attr(mask, "recurrence_threshold")))
})

test_that("the automatic cutoff follows the mixture", {
    metrics <- synthetic_window_metrics(
        samples = paste0("P", seq_len(12)), seed = 42
    )
    metrics <- multiply_donor(
        metrics, "1", 0L, 50L, 6, samples = paste0("P", 1:3)
    )
    scores <- scoreIntrogression(metrics, mode = "multi")
    automatic <- buildWindowMask(metrics, scores = scores, mode = "multi")
    expect_gt(sum(automatic$mask), 40L)
    expect_equal(attr(automatic, "recurrence_threshold"), 0.25)
    expect_identical(attr(automatic, "threshold_source"), "automatic")
})

test_that("a single sample returns an unmasked table", {
    metrics <- synthetic_window_metrics(samples = "P1")
    scores <- scoreIntrogression(metrics, mode = "single")
    mask <- buildWindowMask(metrics, scores = scores, mode = "single")
    expect_false(any(mask$mask))
    expect_identical(unique(mask$source), "none")
    expect_true(all(is.na(mask$f)))
    expect_true(all(is.na(mask$posterior)))
})

test_that("a multi mode with one sample is still unmasked", {
    metrics <- synthetic_window_metrics(samples = "P1")
    scores <- scoreIntrogression(metrics, mode = "single")
    mask <- buildWindowMask(metrics, scores = scores, mode = "multi")
    expect_false(any(mask$mask))
    expect_identical(unique(mask$source), "none")
})

test_that("applyWindowMask() sets masked scores to NA", {
    scores <- synthetic_score_table()
    mask <- data.frame(
        group = "1", iid = 0:10, start = 0:10 * 1e6,
        stop = (1:11) * 1e6, f = 1, posterior = 1, mask = TRUE,
        source = "cohort", version = "x", stringsAsFactors = FALSE
    )
    masked <- applyWindowMask(scores, mask)
    expect_true(all(is.na(masked$score[
        masked$group == "1" & masked$iid <= 10
    ])))
    expect_false(any(is.na(masked$score[
        masked$group == "1" & masked$iid > 10
    ])))
})

test_that("an all-NA mask leaves scores untouched", {
    scores <- data.frame(
        group = "1", iid = 0:2, score = c(1, 2, 3),
        stringsAsFactors = FALSE
    )
    mask <- data.frame(
        group = "1", iid = 0:2, mask = NA,
        stringsAsFactors = FALSE
    )
    expect_identical(applyWindowMask(scores, mask), scores)
})
