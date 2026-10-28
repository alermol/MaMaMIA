#!/usr/bin/env Rscript
## A compact 12-scenario synthetic benchmark for the introgression score.
##
## The simulator follows the study design (depth, GC amplitude, spatial
## overdispersion, background, recipient CNV and event size) but uses the
## package's WindowMetrics/ScoreTable/Calls functions, so it measures the
## implementation rather than a reimplementation of the method.

suppressPackageStartupMessages(library(MaMaMIA))

`%||%` <- function(x, y) if (is.null(x)) y else x

set.seed(42)
n_windows <- 200L
groups <- as.character(1:5)
read_length <- 100L
length_bp <- 1e6
n_rep <- 6L

scenarios <- list(
    "baseline=base" = list(),
    "depth=low" = list(depth = 0.05),
    "depth=high" = list(depth = 0.2),
    "gc=none" = list(gc_amp = 0),
    "gc=strong" = list(gc_amp = 0.8),
    "overdispersion=low" = list(rho = 0.2, sd = 0.02),
    "overdispersion=high" = list(rho = 0.9, sd = 0.3),
    "background=low" = list(bg = 0.05),
    "background=high" = list(bg = 0.15),
    "recipient_cnv=deletion" = list(cnv = TRUE),
    "size=small" = list(event = 20L),
    "size=large" = list(event = 80L)
)

ar1 <- function(n, rho, sd) {
    e <- numeric(n)
    e[[1L]] <- stats::rnorm(1, 0, sd)
    for (i in seq.int(2L, n)) {
        e[[i]] <- rho * e[[i - 1L]] + sqrt(1 - rho^2) * stats::rnorm(1, 0, sd)
    }
    e
}

simulate_one <- function(options, rep) {
    depth <- options$depth %||% 0.1
    gc_amp <- options$gc_amp %||% 0.3
    rho <- options$rho %||% 0.7
    sd <- options$sd %||% 0.1
    bg <- options$bg %||% 0.1
    event_windows <- options$event %||% 40L
    samples <- c("P1", "P2", "C1", "C2", "C3", "C4")
    rows <- list()
    for (sample in samples) {
        sample_depth <- depth * exp(stats::rnorm(1, 0, 0.2))
        for (group in groups) {
            iid <- seq.int(0L, n_windows - 1L)
            gc <- 0.4 + 0.1 * sin(iid / 15 + as.integer(group))
            donor <- sample_depth * bg * exp(
                gc_amp * ((gc - 0.5) / 0.1)^2 + ar1(n_windows, rho, sd)
            )
            recipient <- sample_depth * bg * exp(
                gc_amp * ((gc - 0.5) / 0.1)^2 + ar1(n_windows, rho, sd)
            )
            if (identical(group, "1") && sample %in% c("P1", "P2")) {
                idx <- iid >= 50L & iid < 50L + event_windows
                donor[idx] <- donor[idx] * 2.5
            }
            if (isTRUE(options$cnv) && identical(group, "3")) {
                recipient <- recipient * 0.5
            }
            for (role in c("donor", "recipient")) {
                lambda <- if (role == "donor") donor else recipient
                count <- stats::rpois(n_windows, lambda * length_bp / read_length)
                covered <- pmin(
                    length_bp,
                    stats::rbinom(
                        n_windows, length_bp, 1 - exp(-lambda)
                    )
                )
                rows[[length(rows) + 1L]] <- data.frame(
                    sample = sample, chr_id = paste0(role, group),
                    chr_name = paste0(group, if (role == "donor") "S" else "A"),
                    group = group, role = role,
                    component = if (role == "donor") "S" else "A",
                    copy_number = 2, iid = iid, start = iid * length_bp,
                    stop = (iid + 1L) * length_bp, length = length_bp,
                    count = as.integer(count), covered = as.integer(covered),
                    breadth = covered / length_bp, gc = gc,
                    lambda_breadth = -log(pmax(1 - covered / length_bp, 1e-12)),
                    lambda_count = count * read_length / length_bp,
                    read_length = read_length, stringsAsFactors = FALSE
                )
            }
        }
    }
    metrics <- do.call(rbind, rows)
    truth <- data.frame(
        sample = c("P1", "P2"), group = "1",
        start = 50L, end = 50L + event_windows - 1L,
        stringsAsFactors = FALSE
    )
    list(metrics = metrics, truth = truth)
}

evaluate <- function(scored, truth, mode) {
    event <- scored$score[scored$group == "1" &
        scored$iid >= truth$start[[1L]] &
        scored$iid <= truth$end[[1L]]]
    control <- scored$score[scored$group != "1"]
    auc <- if (length(event) > 0L && length(control) > 0L) {
        as.numeric(wilcox.test(event, control, exact = FALSE)$statistic) /
            (length(event) * length(control))
    } else {
        NA_real_
    }
    calls <- segmentScores(scored, min_mb = 0.5)
    calibrated <- calibrateSegments(
        calls,
        controls = if (identical(mode, "multi")) c("C1", "C2", "C3", "C4") else NULL,
        mode = mode,
        length_bins_mb = c(0.5, 5, 10, 25, 50, 100, Inf),
        fdr = 0.01,
        floors = if (identical(mode, "single")) {
            c(0.3, 0.3, 0.3, 0.3, 0.3, 0.3)
        } else {
            rep(0.25, 6)
        }
    )
    pass <- calibrated$calls[calibrated$calls$pass, , drop = FALSE]
    detected <- any(
        pass$group == truth$group[[1L]] &
            pass$start_iid <= truth$end[[1L]] &
            pass$end_iid >= truth$start[[1L]]
    )
    list(
        auc = auc, detection = as.numeric(detected),
        n_calls = nrow(pass),
        control_calls = sum(
            !pass$group %in% truth$group &
                pass$width_bp >= 20e6
        )
    )
}

results <- list()
for (scenario in names(scenarios)) {
    for (rep in seq_len(n_rep)) {
        simulated <- simulate_one(scenarios[[scenario]], rep)
        # Single-sample: each positive sample is scored and calibrated alone.
        for (sample in c("P1", "P2")) {
            scored <- scoreIntrogression(
                simulated$metrics[simulated$metrics$sample == sample, ],
                mode = "single"
            )
            truth <- simulated$truth[simulated$truth$sample == sample, ]
            values <- evaluate(scored, truth, "single")
            results[[length(results) + 1L]] <- data.frame(
                scenario = scenario, rep = rep, mode = "single",
                sample = sample, auc = values$auc,
                detection = values$detection, n_calls = values$n_calls,
                control_calls = values$control_calls,
                stringsAsFactors = FALSE
            )
        }
        # Multi-sample: one cohort score, one call set.
        scored <- scoreIntrogression(simulated$metrics, mode = "multi")
        for (sample in c("P1", "P2")) {
            truth <- simulated$truth[simulated$truth$sample == sample, ]
            values <- evaluate(
                scored[scored$sample == sample, ], truth, "multi"
            )
            results[[length(results) + 1L]] <- data.frame(
                scenario = scenario, rep = rep, mode = "multi",
                sample = sample, auc = values$auc,
                detection = values$detection, n_calls = values$n_calls,
                control_calls = values$control_calls,
                stringsAsFactors = FALSE
            )
        }
    }
}
out <- do.call(rbind, results)
path <- file.path(
    "test_data", "introgression-study", "results",
    "synthetic_benchmark_pipeline.tsv"
)
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)

summary_table <- aggregate(
    cbind(detection, auc, control_calls) ~ scenario + mode,
    data = out, FUN = mean
)
print(summary_table, row.names = FALSE)
cat("\nOverall detection by mode:\n")
print(tapply(out$detection, out$mode, mean))
cat("Wrote", path, "\n")
