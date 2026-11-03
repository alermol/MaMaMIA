## Component evidence diagnostics.

component_evidence_fixture <- function(n = 5L) {
    iid <- seq_len(n)
    recipient <- data.frame(
        sample = "S1", group = "1", role = "recipient", component = "A",
        iid = iid, lambda_breadth = 0.01, lambda_count = 0.01,
        stringsAsFactors = FALSE
    )
    donor <- rbind(
        data.frame(
            sample = "S1", group = "1", role = "donor", component = "St",
            iid = iid, lambda_breadth = 0.02, lambda_count = 0.02,
            stringsAsFactors = FALSE
        ),
        data.frame(
            sample = "S1", group = "1", role = "donor", component = "V",
            iid = iid, lambda_breadth = 0.005, lambda_count = 0.005,
            stringsAsFactors = FALSE
        )
    )
    metrics <- rbind(recipient, donor)
    scores <- data.frame(
        sample = "S1", group = "1", iid = iid,
        score = 2, donor_mean = 0.02, donor_component = "St",
        metric = "breadth", score_type = "paired_log_ratio",
        stringsAsFactors = FALSE
    )
    calls <- data.frame(
        sample = "S1", group = "1", start_iid = 1L, end_iid = n,
        class = "addition", pass = TRUE, donor_component = "St",
        stringsAsFactors = FALSE
    )
    structure(
        list(
            metrics = metrics,
            scores = scores,
            calls = calls,
            assembly = list(
                contigs = list(
                    donor = data.frame(
                        group = "1", component = c("St", "V"),
                        stringsAsFactors = FALSE
                    )
                )
            )
        ),
        class = c("BamAnalysisResult", "list")
    )
}

test_that("component evidence reports every donor component per call", {
    evidence <- componentEvidence(component_evidence_fixture(), sample = "S1")
    expect_identical(nrow(evidence), 2L)
    expect_identical(
        evidence$donor_component[order(evidence$component_rank)],
        c("St", "V")
    )
    selected <- evidence[evidence$component_selected, , drop = FALSE]
    expect_identical(selected$donor_component, "St")
    expect_identical(selected$component_rank, 1L)
    expect_true(selected$component_supported)
    expect_gt(
        evidence$mean_component_score[evidence$donor_component == "St"],
        evidence$mean_component_score[evidence$donor_component == "V"]
    )
})

test_that("component evidence handles an empty call table", {
    fixture <- component_evidence_fixture()
    fixture$calls <- fixture$calls[0, , drop = FALSE]
    evidence <- componentEvidence(fixture, sample = "S1")
    expect_identical(nrow(evidence), 0L)
    expect_true("component_rank" %in% names(evidence))
})

test_that("component evidence keeps the group-level call separate", {
    evidence <- componentEvidence(component_evidence_fixture(), sample = "S1")
    expect_identical(unique(evidence$call_class), "addition")
    expect_true(all(evidence$call_pass))
    expect_identical(unique(evidence$chromosome_group), "1")
    expect_identical(unique(evidence$n_component_contigs), 1L)
})
