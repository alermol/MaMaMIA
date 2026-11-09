## Pure explorer filtering and module seams.

explorer_filter_fixture <- function() {
    data.frame(
        group = c("G1", "G2", "G2", "G3", "G4"),
        class = c("substitution", "addition", "deletion", NA, "complex"),
        class_status = c(
            "confident", "needs_review", "ambiguous", "confident", NA
        ),
        copy_state = c("balanced", "gain", NA, "loss", "neutral"),
        confidence = c("high", "low", "moderate", "high", NA),
        width_bp = c(2e6, 0.5e6, 4e6, 3e6, NA_real_),
        pass = c(TRUE, FALSE, TRUE, NA, TRUE),
        recipient_group = c(NA, "G1", NA, NA, "G2"),
        recipient_cross_group = c(FALSE, TRUE, FALSE, NA, TRUE),
        recipient_assignment = c(
            "assigned", "no_recipient_loss", "assigned", NA,
            "no_recipient_loss"
        ),
        stringsAsFactors = FALSE
    )
}

test_that("explorer_calls() passes NULL and empty calls through", {
    expect_identical(
        explorer_calls(NULL),
        list(calls = NULL, display_group = character(), unplaced = 0L)
    )
    empty <- explorer_filter_fixture()[0, , drop = FALSE]
    out <- explorer_calls(empty, list(pass_only = TRUE))
    expect_identical(out$calls, empty)
    expect_identical(out$display_group, character())
    expect_identical(out$unplaced, 0L)
})

test_that("explorer_calls() counts and optionally hides unplaced calls", {
    calls <- explorer_filter_fixture()
    out <- explorer_calls(calls)
    expect_identical(out$unplaced, 2L)
    expect_identical(
        out$calls$group, c("G1", "G2", "G2", "G3", "G4")
    )

    hidden <- explorer_calls(calls, list(hide_unplaced = TRUE))
    expect_identical(hidden$calls$group, c("G1", "G2", "G3"))
    expect_identical(hidden$unplaced, 2L)
    expect_identical(hidden$display_group, c("G1", "G2", "G3"))

    passed <- explorer_calls(
        calls, list(pass_only = TRUE, hide_unplaced = TRUE)
    )
    expect_identical(passed$calls$group, c("G1", "G2"))
    expect_identical(passed$unplaced, 1L)
})

test_that("explorer_calls() applies pass-only without treating NA as pass", {
    out <- explorer_calls(explorer_filter_fixture(), list(pass_only = TRUE))
    expect_identical(out$calls$group, c("G1", "G2", "G4"))
    expect_true(all(isTRUE_vec(out$calls$pass)))
})

test_that("group selection includes calls displayed under the recipient", {
    calls <- explorer_filter_fixture()
    out <- explorer_calls(calls, list(groups = "G1"))
    expect_identical(out$calls$group, c("G1", "G2"))
    expect_identical(out$display_group, c("G1", "G1"))

    out <- explorer_calls(calls, list(groups = c("G2", "G3")))
    expect_identical(
        out$calls$group, c("G2", "G2", "G3", "G4")
    )
    expect_identical(
        out$display_group, c("G1", "G2", "G3", "G2")
    )
})

test_that("explorer_display_group() keeps cross-group semantics", {
    calls <- explorer_filter_fixture()
    expect_identical(
        explorer_display_group(calls),
        c("G1", "G1", "G2", "G3", "G2")
    )
    no_recipient <- calls[, c("group", "pass"), drop = FALSE]
    expect_identical(
        explorer_display_group(no_recipient),
        as.character(no_recipient$group)
    )
})

test_that("explorer_calls() filters class, status, copy and confidence", {
    calls <- explorer_filter_fixture()
    out <- explorer_calls(calls, list(class = c("addition", "complex")))
    expect_identical(out$calls$group, c("G2", "G4"))

    out <- explorer_calls(
        calls, list(class_status = c("confident", "ambiguous"))
    )
    expect_identical(out$calls$group, c("G1", "G2", "G3"))

    out <- explorer_calls(calls, list(copy_state = c("gain", "loss")))
    expect_identical(out$calls$group, c("G2", "G3"))

    out <- explorer_calls(calls, list(confidence = "high"))
    expect_identical(out$calls$group, c("G1", "G3"))
})

test_that("explorer_calls() applies the minimum width in base pairs", {
    calls <- explorer_filter_fixture()
    out <- explorer_calls(calls, list(min_width_bp = 2e6))
    ## A missing width stays in the index as an NA row, matching the
    ## dashboard's historical subsetting semantics.
    expect_identical(out$calls$group, c("G1", "G2", "G3", NA))
    expect_equal(out$calls$width_bp, c(2e6, 4e6, 3e6, NA_real_))

    expect_identical(
        explorer_calls(calls, list(min_width_bp = 0))$calls,
        calls
    )
    expect_identical(
        explorer_calls(calls, list(min_width_bp = NA_real_))$calls,
        calls
    )
    expect_identical(
        explorer_calls(calls, list(min_width_bp = c(1, 2)))$calls,
        calls
    )
})

test_that("explorer_calls() keeps filter precedence", {
    calls <- explorer_filter_fixture()
    out <- explorer_calls(calls, list(
        pass_only = TRUE,
        groups = c("G1", "G2"),
        class = c("substitution", "deletion"),
        class_status = c("confident", "ambiguous"),
        copy_state = c("balanced", "loss"),
        confidence = "high",
        min_width_bp = 1e6
    ))
    expect_identical(out$calls$group, "G1")
    expect_identical(out$display_group, "G1")
})

test_that("dashboard theme and CSS helpers render from the theme module", {
    css <- bam_dashboard_css()
    expect_type(css, "character")
    expect_length(css, 1L)
    expect_match(css, "--mm-accent", fixed = TRUE)
    expect_match(css, ".mamamia-header", fixed = TRUE)
    theme <- bam_dashboard_theme()
    expect_false(is.null(theme))
    expect_gt(length(as.character(mamamia_favicon())), 0L)
})

test_that("cohort export card renders with and without a result", {
    before <- as.character(bam_cohort_export_ui(NULL))
    expect_match(before, "Run the analysis to enable the exports.",
        fixed = TRUE
    )
    result <- structure(
        list(samples = data.frame(sample = c("A", "B"))),
        class = "BamAnalysisResult"
    )
    card <- as.character(bam_cohort_export_ui(result))
    expect_match(card, "mamamia-export-cohort", fixed = TRUE)
    expect_match(card, "cohort_table", fixed = TRUE)
    expect_match(card, "Download all samples", fixed = TRUE)
})
