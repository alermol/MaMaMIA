## P5: window masks (panel, cohort, blacklist, physical fallback).

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

panel_fixture <- function() {
    metrics <- synthetic_window_metrics(
        samples = paste0("P", seq_len(12)), seed = 31
    )
    ## A recurrent terminal artifact in every panel sample.
    metrics <- multiply_donor(metrics, "1", 0L, 20L, 6)
    ## A unique event in one sample; it must not become a panel mask.
    metrics <- multiply_donor(metrics, "2", 100L, 130L, 6, samples = "P1")
    metrics
}

test_that("buildPanelMask() masks recurrent artifacts only", {
    metrics <- panel_fixture()
    panel <- buildPanelMask(metrics, z = 0.5, recurrence_threshold = 0.4)
    expect_identical(nrow(panel), 400L)
    expect_true(all(panel$mask[panel$group == "1" & panel$iid <= 20]))
    expect_false(any(panel$mask[panel$group == "2" & panel$iid >= 100 &
        panel$iid <= 130]))
    expect_true(mean(panel$mask) < 0.2)
    expect_identical(unique(panel$source), "panel")
    expect_true(nzchar(unique(panel$version)))
})

test_that("buildWindowMask() applies a panel mask and validates the grid", {
    metrics <- panel_fixture()
    panel <- buildPanelMask(metrics)
    target <- metrics[metrics$sample == "P1", , drop = FALSE]
    mask <- buildWindowMask(target, mode = "single", panel_mask = panel)
    expect_identical(mask$mask, panel$mask)
    expect_identical(unique(mask$source), "panel")
    expect_error(
        buildWindowMask(
            target, mode = "single",
            panel_mask = panel[-1L, , drop = FALSE]
        ),
        "does not cover"
    )
})

test_that("a small panel warns and uses a fixed fraction", {
    metrics <- synthetic_window_metrics(samples = paste0("P", 1:4))
    metrics <- multiply_donor(metrics, "1", 0L, 10L, 6)
    expect_warning(
        panel <- buildPanelMask(metrics),
        "only 4 sample"
    )
    expect_true(all(panel$mask[panel$group == "1" & panel$iid <= 10]))
})

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
        metrics, scores = scores, mode = "multi",
        physical_fallback_mb = 0
    )
    expect_true(any(mask$mask[mask$group == "1" & mask$iid <= 20]))
    expect_false(any(mask$mask[mask$group == "2" & mask$iid >= 100 &
        mask$iid <= 130]))
    expect_identical(unique(mask$source), "cohort")
})

test_that("a blacklist BED masks overlapping windows", {
    metrics <- synthetic_window_metrics()
    bed <- data.frame(
        chr_id = "d1S", start = 0, stop = 2e6, stringsAsFactors = FALSE
    )
    mask <- buildWindowMask(
        metrics, mode = "single", blacklist = bed,
        physical_fallback_mb = 0
    )
    expect_true(all(mask$mask[mask$group == "1" & mask$iid <= 1]))
    expect_false(any(mask$mask[mask$group == "1" & mask$iid > 1]))
    expect_identical(unique(mask$source[mask$mask]), "blacklist")
})

test_that("the physical fallback masks contig ends only", {
    metrics <- synthetic_window_metrics()
    mask <- buildWindowMask(
        metrics, mode = "single", physical_fallback_mb = 1
    )
    expect_true(mask$mask[mask$group == "1" & mask$iid == 0])
    expect_true(mask$mask[mask$group == "1" & mask$iid == 199])
    expect_false(mask$mask[mask$group == "1" & mask$iid == 100])
    expect_identical(unique(mask$source[mask$mask]), "physical")
})

test_that("applyWindowMask() sets masked scores to NA", {
    scores <- synthetic_score_table()
    mask <- data.frame(
        group = "1", iid = 0:10, start = 0:10 * 1e6,
        stop = (1:11) * 1e6, f = 1, posterior = 1, mask = TRUE,
        source = "panel", version = "x", stringsAsFactors = FALSE
    )
    masked <- applyWindowMask(scores, mask)
    expect_true(all(is.na(masked$score[
        masked$group == "1" & masked$iid <= 10
    ])))
    expect_false(any(is.na(masked$score[
        masked$group == "1" & masked$iid > 10
    ])))
})
