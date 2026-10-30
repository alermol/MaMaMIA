## P4: median polish and ploidy-agnostic scoring.

score_metrics <- function(n = 200L, groups = c("1", "2"),
                          donor_components = "S",
                          recipient_components = c("A"),
                          donor_lambda = 0.1,
                          recipient_lambda = 0.1,
                          gc_effect = 0,
                          sample = "S1") {
    rows <- list()
    for (group in groups) {
        iid <- seq.int(0L, n - 1L)
        gc <- 0.4 + 0.1 * sin(iid / 20)
        for (component in donor_components) {
            lambda <- donor_lambda * exp(
                gc_effect * ((gc - 0.5) / 0.1)^2
            ) * exp(stats::rnorm(n, 0, 0.05))
            rows[[length(rows) + 1L]] <- data.frame(
                sample = sample, chr_id = paste0("d", group, component),
                chr_name = paste0(group, component), group = group,
                role = "donor", component = component, copy_number = 2,
                iid = iid, start = iid * 1e6, stop = (iid + 1) * 1e6,
                length = 1e6, count = round(lambda * 1e6 / 100),
                covered = pmin(1e6, round((1 - exp(-lambda)) * 1e6)),
                breadth = 1 - exp(-lambda), gc = gc,
                lambda_breadth = lambda, lambda_count = lambda,
                read_length = 100, stringsAsFactors = FALSE
            )
        }
        for (component in recipient_components) {
            lambda <- recipient_lambda * exp(
                gc_effect * ((gc - 0.5) / 0.1)^2
            ) * exp(stats::rnorm(n, 0, 0.05))
            rows[[length(rows) + 1L]] <- data.frame(
                sample = sample, chr_id = paste0("r", group, component),
                chr_name = paste0(group, component), group = group,
                role = "recipient", component = component, copy_number = 2,
                iid = iid, start = iid * 1e6, stop = (iid + 1) * 1e6,
                length = 1e6, count = round(lambda * 1e6 / 100),
                covered = pmin(1e6, round((1 - exp(-lambda)) * 1e6)),
                breadth = 1 - exp(-lambda), gc = gc,
                lambda_breadth = lambda, lambda_count = lambda,
                read_length = 100, stringsAsFactors = FALSE
            )
        }
    }
    do.call(rbind, rows)
}

add_donor_event <- function(metrics, group = "1", from = 50L, to = 90L,
                            multiplier = 2) {
    idx <- metrics$role == "donor" & metrics$group == group &
        metrics$iid >= from & metrics$iid <= to
    metrics$lambda_breadth[idx] <- metrics$lambda_breadth[idx] * multiplier
    metrics$lambda_count[idx] <- metrics$lambda_count[idx] * multiplier
    metrics$breadth[idx] <- 1 - exp(-metrics$lambda_breadth[idx])
    metrics$covered[idx] <- pmin(
        metrics$length[idx],
        round(metrics$breadth[idx] * metrics$length[idx])
    )
    metrics
}

add_recipient_event <- function(metrics, group = "1", from = 50L, to = 90L,
                                multiplier = 0.5) {
    idx <- metrics$role == "recipient" & metrics$group == group &
        metrics$iid >= from & metrics$iid <= to
    metrics$lambda_breadth[idx] <- metrics$lambda_breadth[idx] * multiplier
    metrics$lambda_count[idx] <- metrics$lambda_count[idx] * multiplier
    metrics$breadth[idx] <- 1 - exp(-metrics$lambda_breadth[idx])
    metrics$covered[idx] <- pmin(
        metrics$length[idx],
        round(metrics$breadth[idx] * metrics$length[idx])
    )
    metrics
}

test_that("medianPolish() recovers an additive two-way structure", {
    set.seed(1)
    a <- c(-1, 0, 2)
    b <- seq(-0.5, 0.5, length.out = 20)
    x <- outer(b, a, "+") + matrix(rnorm(60, 0, 0.01), nrow = 20)
    polished <- medianPolish(x)
    expect_true(polished$converged)
    expect_lt(max(abs(polished$residuals)), 0.05)
    expect_equal(polished$column, a, tolerance = 0.05)
})

test_that("medianPolish() handles missing values", {
    x <- matrix(c(1, 2, NA, 4, 5, 6), nrow = 2)
    polished <- medianPolish(x)
    expect_identical(dim(polished$residuals), dim(x))
    expect_true(is.na(polished$residuals[1, 2]))
})

test_that("a diploid substitution gives a positive score", {
    set.seed(2)
    metrics <- score_metrics(recipient_components = "A")
    metrics <- add_donor_event(metrics, multiplier = 2)
    scored <- scoreIntrogression(metrics, mode = "single")
    event <- scored$score[scored$group == "1" & scored$iid >= 50 &
        scored$iid <= 90]
    control <- scored$score[scored$group == "1" & scored$iid < 40]
    expect_gt(mean(event), 0.4)
    expect_lt(abs(mean(control)), 0.1)
})

test_that("tetraploid and hexaploid substitution signals follow log ploidy", {
    set.seed(3)
    tetra <- add_donor_event(
        score_metrics(recipient_components = c("A", "B")),
        multiplier = 2
    )
    tetra_score <- mean(
        scoreIntrogression(tetra, mode = "single")$score[
            scoreIntrogression(tetra, mode = "single")$group == "1" &
                scoreIntrogression(tetra, mode = "single")$iid >= 50 &
                scoreIntrogression(tetra, mode = "single")$iid <= 90
        ]
    )
    expect_gt(tetra_score, 0.5)

    set.seed(4)
    hexa <- add_donor_event(
        score_metrics(recipient_components = c("A", "B", "D")),
        multiplier = 1.5
    )
    hexa_scored <- scoreIntrogression(hexa, mode = "single")
    hexa_score <- mean(
        hexa_scored$score[
            hexa_scored$group == "1" & hexa_scored$iid >= 50 &
                hexa_scored$iid <= 90
        ]
    )
    expect_gt(hexa_score, 0.2)
    expect_lt(hexa_score, 0.8)
})

test_that("multi-donor max selects the source component and sum doubles", {
    set.seed(5)
    metrics <- score_metrics(
        donor_components = c("S1", "S2"),
        recipient_components = c("A", "B")
    )
    idx <- metrics$role == "donor" & metrics$component == "S2" &
        metrics$group == "1" & metrics$iid >= 50 & metrics$iid <= 90
    metrics$lambda_breadth[idx] <- metrics$lambda_breadth[idx] * 3
    metrics$lambda_count[idx] <- metrics$lambda_count[idx] * 3
    metrics$breadth[idx] <- 1 - exp(-metrics$lambda_breadth[idx])
    metrics$covered[idx] <- pmin(
        metrics$length[idx], round(metrics$breadth[idx] * metrics$length[idx])
    )
    maxed <- scoreIntrogression(
        metrics, mode = "single", donor_aggregation = "max"
    )
    event <- maxed[maxed$group == "1" & maxed$iid >= 50 &
        maxed$iid <= 90, ]
    expect_true(all(event$donor_component == "S2"))
    expect_gt(mean(event$score), 0.5)

    summed <- scoreIntrogression(
        metrics, mode = "single", donor_aggregation = "sum"
    )
    sum_event <- summed[summed$group == "1" & summed$iid >= 50 &
        summed$iid <= 90, ]
    expect_gt(mean(sum_event$donor_mean), mean(event$donor_mean))
})

test_that("donor aggregation defaults to max", {
    set.seed(15)
    metrics <- score_metrics(
        donor_components = c("S1", "S2"),
        recipient_components = "A"
    )
    idx <- metrics$role == "donor" & metrics$component == "S2" &
        metrics$group == "1" & metrics$iid >= 50 & metrics$iid <= 90
    metrics$lambda_breadth[idx] <- metrics$lambda_breadth[idx] * 3
    metrics$lambda_count[idx] <- metrics$lambda_count[idx] * 3
    metrics$breadth[idx] <- 1 - exp(-metrics$lambda_breadth[idx])
    metrics$covered[idx] <- pmin(
        metrics$length[idx], round(metrics$breadth[idx] * metrics$length[idx])
    )
    default <- scoreIntrogression(metrics, mode = "single")
    meaned <- scoreIntrogression(
        metrics, mode = "single", donor_aggregation = "mean"
    )
    event <- default$group == "1" & default$iid >= 50 & default$iid <= 90
    expect_true(all(default$donor_component[event] == "S2"))
    expect_gt(
        mean(default$score[event]),
        mean(meaned$score[event])
    )
})

test_that("median polish removes a constant donor ploidy offset", {
    set.seed(6)
    m1 <- score_metrics(sample = "A", donor_lambda = 0.2)
    m2 <- score_metrics(sample = "B", donor_lambda = 0.2)
    metrics <- rbind(m1, m2)
    scored <- scoreIntrogression(metrics, mode = "multi")
    expect_lt(abs(mean(scored$score)), 0.15)
    expect_lt(stats::sd(scored$score), 0.2)
})

test_that("a cohort keeps the event in the affected sample", {
    set.seed(7)
    positive <- add_donor_event(score_metrics(sample = "P"), multiplier = 2)
    controls <- lapply(c("C1", "C2", "C3", "C4"), function(sample) {
        score_metrics(sample = sample)
    })
    scored <- scoreIntrogression(
        do.call(rbind, c(list(positive), controls)), mode = "multi"
    )
    p <- scored$score[scored$sample == "P" & scored$group == "1" &
        scored$iid >= 50 & scored$iid <= 90]
    c <- scored$score[scored$sample != "P" & scored$group == "1" &
        scored$iid >= 50 & scored$iid <= 90]
    expect_gt(mean(p), 0.4)
    expect_lt(abs(mean(c)), 0.3)
})

test_that("the GC spline removes a shared GC gradient", {
    set.seed(8)
    metrics <- score_metrics(gc_effect = 0.8)
    raw <- scoreIntrogression(metrics, mode = "single", gc_correction = FALSE)
    corrected <- scoreIntrogression(metrics, mode = "single")
    expect_lt(stats::sd(corrected$score), stats::sd(raw$score))
    expect_lt(abs(mean(corrected$score)), 0.1)
})

test_that("zero windows do not produce NaN or Inf", {
    metrics <- score_metrics(n = 50L)
    metrics$lambda_breadth <- 0
    metrics$lambda_count <- 0
    metrics$breadth <- 0
    metrics$covered <- 0
    metrics$count <- 0
    scored <- scoreIntrogression(metrics, mode = "single")
    expect_true(all(is.finite(scored$score)))
    expect_true(all(is.finite(scored$donor_mean)))
    expect_true(all(is.finite(scored$recipient_mean)))
})

test_that("the copy-number model reports donor copies and fraction", {
    metrics <- score_metrics(
        donor_lambda = 0.2, recipient_lambda = 0.1,
        recipient_components = c("A", "B", "D")
    )
    scored <- scoreIntrogression(
        metrics, mode = "single", sample_composition = c(A = 2, B = 2, D = 2)
    )
    expect_true(all(is.finite(scored$donor_copies)))
    expect_true(all(scored$donor_copies > 0))
    expect_true(all(scored$donor_fraction > 0 & scored$donor_fraction < 1))
})

test_that("donor-only scoring is available as a fallback", {
    set.seed(9)
    metrics <- add_donor_event(score_metrics(), multiplier = 3)
    scored <- scoreIntrogression(
        metrics, mode = "single", score_type = "donor_only"
    )
    event <- scored$score[scored$group == "1" & scored$iid >= 50 &
        scored$iid <= 90]
    expect_gt(mean(event), 0.5)
    expect_identical(unique(scored$score_type), "donor_only")
})

test_that("the total-depth residual separates balanced and unbalanced events", {
    set.seed(21)
    base <- score_metrics(recipient_components = "A")
    substitution <- base
    sub_donor <- substitution$role == "donor" & substitution$group == "1" &
        substitution$iid >= 50 & substitution$iid <= 90
    sub_recipient <- substitution$role == "recipient" &
        substitution$group == "1" & substitution$iid >= 50 &
        substitution$iid <= 90
    substitution$lambda_breadth[sub_donor] <- 0.2
    substitution$lambda_breadth[sub_recipient] <- 0
    substitution$breadth[sub_donor] <- 1 - exp(-0.2)
    substitution$breadth[sub_recipient] <- 0
    substitution$covered[sub_donor] <- round(
        substitution$breadth[sub_donor] * substitution$length[sub_donor]
    )
    substitution$covered[sub_recipient] <- 0
    scored <- scoreIntrogression(substitution, mode = "single")
    event <- scored$group == "1" & scored$iid >= 50 & scored$iid <= 90
    expect_lt(abs(mean(scored$total_residual[event])), 0.2)

    addition <- add_donor_event(base, multiplier = 2)
    scored <- scoreIntrogression(addition, mode = "single")
    expect_gt(mean(scored$total_residual[event]), 0.2)

    deletion <- add_recipient_event(base, multiplier = 0.5)
    scored <- scoreIntrogression(deletion, mode = "single")
    expect_lt(mean(scored$total_residual[event]), -0.2)
})

test_that("component evidence is centred before measuring spread", {
    set.seed(22)
    metrics <- score_metrics(recipient_components = c("A", "B", "D"))
    scored <- scoreIntrogression(
        metrics, mode = "single",
        sample_composition = c(A = 2, B = 2, D = 2)
    )
    expect_true(all(is.finite(scored$donor_component_spread)))
    expect_true(all(is.finite(scored$recipient_component_spread)))
    expect_true(all(is.finite(scored$n_donor_supported)))
    ## Component baselines are removed, so normal windows have a small spread.
    expect_lt(stats::median(scored$recipient_component_spread), 1)
})
