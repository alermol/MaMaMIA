## P6: CBS segmentation, gap merging and masking.

segment_positive <- function(calls) {
    calls[!is.na(calls$score_mean) & calls$score_mean > 0.5, , drop = FALSE]
}

test_that("merging three segments keeps the window-weighted mean", {
    calls <- data.frame(
        start_iid = c(0L, 10L, 20L),
        end_iid = c(4L, 14L, 24L),
        n_windows = c(5L, 5L, 5L),
        score_mean = c(1, 2, 10),
        stringsAsFactors = FALSE
    )
    profile <- stats::setNames(rep(0, 30L), as.character(0:29))
    merged <- MaMaMIA:::.segment_merge(
        calls, profile, merge_gap = 10L, merge_floor = -1, min_mb = 0
    )
    expect_equal(nrow(merged), 1L)
    expect_equal(merged$score_mean, (1 * 5 + 2 * 5 + 10 * 5) / 15)
    expect_equal(merged$n_windows, 25L)
    expect_equal(merged$start_iid, 0L)
    expect_equal(merged$end_iid, 24L)
})

test_that("masked windows do not contribute to segment evidence", {
    scores <- synthetic_score_table(
        event_from = 40L, event_to = 60L, event_effect = 2
    )
    inside <- scores$group == "1" & scores$iid == 50L
    scores$donor_residual[inside] <- 1e6
    mask <- data.frame(
        group = "1", iid = 50L, mask = TRUE,
        stringsAsFactors = FALSE
    )
    calls <- segmentScores(scores, mask = mask, min_mb = 1)
    hit <- calls[
        calls$group == "1" & calls$start_iid <= 50L &
            calls$end_iid >= 50L, ,
        drop = FALSE
    ]
    expect_gt(nrow(hit), 0L)
    expect_lt(max(abs(hit$donor_mean)), 100)
})

test_that("CBS finds a single event within a few windows", {
    scores <- synthetic_score_table(event_effect = 2)
    calls <- segmentScores(scores, min_mb = 1)
    hit <- segment_positive(calls)
    expect_gte(nrow(hit), 1L)
    event <- hit[which.max(hit$score_mean), ]
    expect_lte(abs(event$start_iid - 50), 3L)
    expect_lte(abs(event$end_iid - 90), 3L)
})

test_that("a whole-chromosome event is one merged call", {
    scores <- synthetic_score_table(
        event_from = 0L, event_to = 199L, event_effect = 1.5
    )
    calls <- segmentScores(scores, min_mb = 1)
    hit <- segment_positive(calls)
    expect_equal(nrow(hit[hit$group == "1", ]), 1L)
    expect_equal(hit$start_iid[hit$group == "1"], 0L)
    expect_equal(hit$end_iid[hit$group == "1"], 199L)
})

test_that("two events separated by a true null are not merged", {
    scores <- synthetic_score_table(
        event_from = 20L, event_to = 50L, event_effect = 1.5
    )
    idx <- scores$group == "1" & scores$iid >= 150 & scores$iid <= 180
    scores$score[idx] <- scores$score[idx] + 1.5
    calls <- segmentScores(scores, min_mb = 1)
    hit <- segment_positive(calls)
    hit <- hit[hit$group == "1", ]
    expect_equal(nrow(hit), 2L)
    expect_true(min(hit$end_iid) < 150)
})

test_that("masked windows are NA before segmentation", {
    scores <- synthetic_score_table(event_effect = 2)
    mask <- data.frame(
        group = "1", iid = 0:60, start = 0:60 * 1e6,
        stop = (1:61) * 1e6, f = 1, posterior = 1, mask = TRUE,
        source = "panel", version = "test", stringsAsFactors = FALSE
    )
    calls <- segmentScores(scores, mask = mask, min_mb = 1)
    hit <- segment_positive(calls)
    hit <- hit[hit$group == "1", ]
    expect_gte(nrow(hit), 1L)
    expect_gte(min(hit$start_iid), 60L)
})

test_that("calls shorter than min_mb are dropped", {
    scores <- synthetic_score_table(
        event_from = 50L, event_to = 55L, event_effect = 2
    )
    calls <- segmentScores(scores, min_mb = 20)
    expect_equal(sum(calls$score_mean > 0.5), 0L)
})

test_that("pure noise produces no calibrated calls", {
    scores <- synthetic_score_table(
        event_group = NULL, noise = 0.05, seed = 11
    )
    calls <- segmentScores(scores, min_mb = 1)
    calibrated <- calibrateSegments(
        calls, fdr = 0.01, floors = c(0.3, 0.3, 0.3, 0.3, 0.3)
    )
    expect_equal(sum(calibrated$calls$pass), 0L)
})

test_that("segmentation validates its input", {
    expect_error(segmentScores(data.frame()), "ScoreTable")
    scores <- synthetic_score_table()
    expect_error(
        segmentScores(scores, min_width = 0),
        "Invalid segmentation"
    )
})
