## Call classification from donor and recipient evidence.

test_that("classification uses both donor and recipient evidence", {
    calls <- data.frame(
        sample = "S1",
        group = "1",
        start_iid = 1:4,
        donor_mean = c(1.2, 1.2, 0.1, 0.1),
        recipient_mean = c(-1.2, 0, -1.2, 0),
        score_mean = c(2, 1, 1, 0),
        pass = TRUE,
        stringsAsFactors = FALSE
    )
    out <- classifyCalls(calls)
    expect_identical(
        out$class,
        c("substitution", "addition", "deletion", "other")
    )
    expect_identical(
        out$class_evidence,
        c("donor+recipient", "donor", "recipient", "none")
    )
    expect_identical(
        out$class_confidence,
        c("high", "high", "high", "low")
    )
    expect_identical(
        out$class_status,
        c("confident", "confident", "confident", "needs_review")
    )
    expect_identical(out$confidence, c("high", "high", "low", "low"))
})

test_that("class cutoffs are calibrated from control segments", {
    controls <- data.frame(
        sample = "C",
        mode = "multi",
        length_bin = ">100",
        donor_mean = seq(0, 1.9, length.out = 20),
        recipient_mean = -seq(0, 1.9, length.out = 20),
        stringsAsFactors = FALSE
    )
    thresholds <- .call_class_thresholds(controls, controls = "C")
    threshold <- thresholds[["multi|>100"]]
    expect_gt(threshold$donor_mod, 0.5)
    expect_lt(threshold$loss_mod, -0.5)
    expect_gte(threshold$donor_high, threshold$donor_mod)
    expect_lte(threshold$loss_high, threshold$loss_mod)
})

test_that("total depth separates balanced, gain, loss and compound events", {
    calls <- data.frame(
        sample = "S1",
        group = "1",
        start_iid = 1:5,
        donor_mean = c(1.2, 1.2, 0.1, 1.2, 0.1),
        recipient_mean = c(-1.2, 0, -1.2, -1.2, 0),
        total_mean = c(0, 1.2, -1.2, 1.2, 1.2),
        score_mean = c(2, 1, 1, 2, 1),
        pass = TRUE,
        stringsAsFactors = FALSE
    )
    out <- classifyCalls(calls)
    expect_identical(
        out$class,
        c("substitution", "addition", "deletion", "complex", "other")
    )
    expect_identical(
        out$copy_state,
        c("balanced", "gain", "loss", "compound", "gain")
    )
    expect_identical(
        out$total_evidence,
        c("neutral", "gain", "loss", "gain", "gain")
    )
    expect_identical(out$class_status[[4L]], "needs_review")
    expect_identical(out$class_status[[5L]], "ambiguous")
})

test_that("component heterogeneity downgrades a confident class", {
    calls <- data.frame(
        sample = "S1",
        group = "1",
        start_iid = 1:2,
        donor_mean = c(1.2, 1.2),
        recipient_mean = c(-1.2, -1.2),
        total_mean = c(0, 0),
        score_mean = c(2, 2),
        pass = TRUE,
        donor_component_spread = c(0, 3),
        recipient_component_spread = c(0, 0),
        n_donor_supported = c(1, 1),
        stringsAsFactors = FALSE
    )
    out <- classifyCalls(calls)
    expect_identical(out$component_status, c("consistent", "heterogeneous"))
    expect_identical(out$class_confidence, c("high", "moderate"))
    expect_identical(out$class_status, c("confident", "ambiguous"))
})

test_that("low mappability downgrades a confident class", {
    calls <- data.frame(
        sample = "S1",
        group = "1",
        start_iid = 1:2,
        donor_mean = c(1.2, 1.2),
        recipient_mean = c(-1.2, -1.2),
        total_mean = c(0, 0),
        score_mean = c(2, 2),
        pass = TRUE,
        mapq_fraction = c(0.9, 0.2),
        stringsAsFactors = FALSE
    )
    out <- classifyCalls(calls)
    expect_identical(out$mappability_status, c("ok", "low"))
    expect_identical(out$class_confidence, c("high", "moderate"))
    expect_identical(out$class_status, c("confident", "ambiguous"))
})
