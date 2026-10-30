## Assignment calibration helpers.
##
## Assignment is a second detection problem: a call is compared with every
## recipient contig and the best host must beat what a width-matched random
## interval would produce. These helpers build the per-contig recipient
## residuals, interval statistics, empirical null thresholds and non-overlapping
## interval schedule used by assignCallPositions().

#' Per-contig corrected recipient residuals
#'
#' A cohort is corrected with median polish per contig and window; a single
#' sample is corrected with a GC spline per contig. `depletion` is the negative
#' residual, so a deletion is positive. `expected_breadth` is the residual's
#' baseline back-transformed to breadth units for dosage.
#' @noRd
.assignment_recipient_residuals <- function(recipient, mode, masks = NULL,
                                            eps = 1e-6) {
    if (!is.data.frame(recipient) || nrow(recipient) == 0L) {
        return(recipient)
    }
    recipient$log_lambda <- log(pmax(as.numeric(recipient$lambda_breadth), 0) +
        eps)
    multi <- identical(mode, "multi") &&
        length(unique(as.character(recipient$sample))) > 1L
    if (multi) {
        recipient$residual <- .side_correction_multi(
            recipient$log_lambda,
            as.character(recipient$sample),
            as.character(recipient$chr_id),
            as.integer(recipient$iid)
        )
    } else {
        recipient$residual <- NA_real_
        for (chr_id in unique(as.character(recipient$chr_id))) {
            idx <- as.character(recipient$chr_id) == chr_id
            recipient$residual[idx] <- .side_correction_single(
                recipient$log_lambda[idx],
                recipient$gc[idx]
            )
        }
    }
    recipient$depletion <- -as.numeric(recipient$residual)
    recipient$masked <- FALSE
    if (!is.null(masks) && is.data.frame(masks) &&
        all(c("group", "iid", "mask") %in% names(masks))) {
        key <- paste(recipient$group, recipient$iid, sep = ":")
        mask_key <- paste(masks$group, masks$iid, sep = ":")
        recipient$masked <- as.logical(
            masks$mask[match(key, mask_key)]
        ) %in% TRUE
        recipient$masked[is.na(recipient$masked)] <- FALSE
    }
    ## Expected breadth must not be derived as observed/exp(residual): a zero
    ## breadth gives zero expected breadth and hides a complete deletion. Use
    ## a cross-sample median (cohort) or a contig median (single sample)
    ## instead.
    if (multi) {
        key <- paste(
            as.character(recipient$chr_id), as.integer(recipient$iid),
            sep = ":"
        )
        values <- as.numeric(recipient$lambda_breadth)
        values[recipient$masked] <- NA_real_
        expected_lambda <- tapply(
            values, key,
            function(x) {
                x <- x[is.finite(x)]
                if (length(x) == 0L) NA_real_ else stats::median(x)
            }
        )
        recipient$expected_lambda <- as.numeric(expected_lambda[key])
    } else {
        expected_lambda <- tapply(
            as.numeric(recipient$lambda_breadth),
            as.character(recipient$chr_id),
            function(values) {
                values <- values[is.finite(values)]
                if (length(values) == 0L) NA_real_ else stats::median(values)
            }
        )
        recipient$expected_lambda <- as.numeric(
            expected_lambda[as.character(recipient$chr_id)]
        )
    }
    expected <- 1 - exp(-pmax(recipient$expected_lambda, 0))
    expected[!is.finite(expected)] <- NA_real_
    recipient$expected_breadth <- pmin(pmax(expected, 0), 1)
    contig_length <- tapply(
        as.numeric(recipient$stop), as.character(recipient$chr_id), max,
        na.rm = TRUE
    )
    recipient$contig_length <- as.numeric(
        contig_length[as.character(recipient$chr_id)]
    )
    recipient
}

#' Mean, standard error and lower bound of an interval
#' @noRd
.assignment_interval_stats <- function(values, inside) {
    x <- as.numeric(values[inside])
    x <- x[is.finite(x)]
    n <- length(x)
    if (n == 0L) {
        return(list(
            n = 0L, mean = NA_real_, se = NA_real_,
            lower = NA_real_, z = NA_real_
        ))
    }
    estimate <- mean(x)
    if (n < 2L) {
        return(list(
            n = n, mean = estimate, se = 0,
            lower = estimate, z = Inf
        ))
    }
    se <- stats::sd(x) / sqrt(.assignment_effective_n(x))
    if (!is.finite(se)) {
        se <- 0
    }
    lower <- estimate - 1.645 * se
    z <- if (se > 0) estimate / se else Inf
    list(n = n, mean = estimate, se = se, lower = lower, z = z)
}

#' Effective sample size for a window profile
#'
#' Windows are autocorrelated, so the number of rows overstates the
#' information in a run. A small autocorrelation correction keeps the standard
#' error from being anticonservative without requiring a full block bootstrap.
#' @noRd
.assignment_effective_n <- function(values) {
    n <- length(values)
    if (n < 4L || stats::sd(values) == 0) {
        return(as.numeric(n))
    }
    acf <- tryCatch(
        stats::acf(values, lag.max = min(10L, n - 1L), plot = FALSE,
                   na.action = stats::na.pass)$acf,
        error = function(e) NULL
    )
    if (is.null(acf) || length(acf) < 2L) {
        return(as.numeric(n))
    }
    rho <- sum(acf[-1L], na.rm = TRUE)
    n_eff <- n / max(1, 1 + 2 * rho)
    min(max(n_eff, 1), n)
}

#' Empirical null thresholds for assignment
#'
#' For each mode and length bin, random intervals of the same width are drawn
#' on the null samples' recipient contigs. The 95th/99th percentiles of the
#' best width-matched lower-bound score replace the fixed depletion thresholds;
#' the null best-over-second difference sets the assignment margin, and the
#' null cross-group-over-same-group difference sets the cross-group margin.
#' @noRd
.assignment_null_thresholds <- function(
    calls, recipient, contigs, controls = NULL,
    permutations = 64L, min_calls = 5L, min_null = 20L,
    same_group_bonus = 1
) {
    if (!is.data.frame(calls) || nrow(calls) < min_calls ||
        !is.data.frame(recipient) || nrow(recipient) == 0L) {
        return(NULL)
    }
    permutations <- max(8L, as.integer(permutations)[[1L]])
    samples <- unique(as.character(recipient$sample))
    null_samples <- if (!is.null(controls) && length(controls) > 0L) {
        intersect(as.character(controls), samples)
    } else {
        samples
    }
    if (length(null_samples) == 0L) {
        return(NULL)
    }
    n_per <- max(as.integer(recipient$iid), na.rm = TRUE) + 1L
    if (!is.finite(n_per) || n_per < 1L) {
        return(NULL)
    }
    contig_ids <- unique(as.character(recipient$chr_id))
    contig_group <- vapply(contig_ids, function(chr_id) {
        as.character(recipient$group[as.character(recipient$chr_id) ==
            chr_id][[1L]])
    }, character(1))

    ## A per-window median across the null samples is a robust cohort null: a
    ## real event in one or two samples cannot raise the threshold for every
    ## call. The pooled within-window spread supplies the standard error.
    profiles <- list()
    for (chr_id in contig_ids) {
        rows <- as.character(recipient$chr_id) == chr_id &
            as.character(recipient$sample) %in% null_samples
        x <- recipient[rows, , drop = FALSE]
        if (nrow(x) == 0L) {
            next
        }
        iid <- as.integer(x$iid) + 1L
        samples_here <- unique(as.character(x$sample))
        matrix_dep <- matrix(
            NA_real_, nrow = length(samples_here), ncol = n_per
        )
        matrix_mask <- matrix(
            TRUE, nrow = length(samples_here), ncol = n_per
        )
        for (k in seq_along(samples_here)) {
            idx <- as.character(x$sample) == samples_here[[k]]
            matrix_dep[k, iid[idx]] <- as.numeric(x$depletion[idx])
            matrix_mask[k, iid[idx]] <- as.logical(x$masked[idx]) %in% TRUE
        }
        valid_count <- rowSums(is.finite(matrix_dep) & !matrix_mask)
        profile <- suppressWarnings(
            apply(matrix_dep, 2L, stats::median, na.rm = TRUE)
        )
        profile[valid_count < 1L] <- NA_real_
        spread <- suppressWarnings(
            apply(matrix_dep, 2L, stats::sd, na.rm = TRUE)
        )
        spread[!is.finite(spread)] <- 0
        profiles[[chr_id]] <- list(
            dep = profile,
            spread = spread,
            valid = is.finite(profile)
        )
    }
    if (length(profiles) == 0L) {
        return(NULL)
    }
    window_depletion <- unlist(lapply(profiles, function(p) {
        p$dep[p$valid]
    }), use.names = FALSE)
    window_depletion <- window_depletion[is.finite(window_depletion)]
    run_threshold <- if (length(window_depletion) >= min_null) {
        max(0.5, as.numeric(stats::quantile(
            window_depletion, 0.95, na.rm = TRUE, names = FALSE
        )))
    } else {
        NULL
    }

    keys <- if ("length_bin" %in% names(calls)) {
        as.character(calls$length_bin)
    } else {
        rep("all", nrow(calls))
    }
    keys[is.na(keys) | !nzchar(keys)] <- "all"
    modes <- if ("mode" %in% names(calls)) {
        as.character(calls$mode)
    } else {
        rep("all", nrow(calls))
    }
    modes[is.na(modes) | !nzchar(modes)] <- "all"
    key <- paste(modes, keys, sep = "|")

    null_rows <- vector("list", nrow(calls))
    for (i in seq_len(nrow(calls))) {
        width <- as.integer(calls$end_iid[[i]]) -
            as.integer(calls$start_iid[[i]]) + 1L
        if (!is.finite(width) || width < 1L || width > n_per) {
            next
        }
        group <- as.character(calls$group[[i]])
        best <- numeric(0)
        second <- numeric(0)
        cross <- numeric(0)
        best_same <- numeric(0)
        for (b in seq_len(permutations)) {
            values <- numeric(length(contig_ids))
            for (j in seq_along(contig_ids)) {
                profile <- profiles[[contig_ids[[j]]]]
                if (is.null(profile)) {
                    values[[j]] <- NA_real_
                    next
                }
                starts <- which(profile$valid)
                starts <- starts[starts + width - 1L <= n_per]
                if (length(starts) == 0L) {
                    values[[j]] <- NA_real_
                    next
                }
                index <- ((b - 1L) * 37L + (j - 1L) * 101L) %%
                    length(starts) + 1L
                lo <- starts[[index]]
                hi <- lo + width - 1L
                active <- profile$valid[lo:hi]
                n <- sum(active)
                if (n < 1L) {
                    values[[j]] <- NA_real_
                    next
                }
                dep <- profile$dep[lo:hi][active]
                est <- mean(dep)
                se <- mean(profile$spread[lo:hi][active], na.rm = TRUE) /
                    sqrt(n)
                if (!is.finite(se)) {
                    se <- 0
                }
                values[[j]] <- est - 1.645 * se +
                    if (identical(contig_group[[j]], group)) {
                        same_group_bonus
                    } else {
                        0
                    }
            }
            values <- values[is.finite(values)]
            if (length(values) == 0L) {
                next
            }
            ordered <- sort(values, decreasing = TRUE)
            best <- c(best, ordered[[1L]])
            if (length(ordered) > 1L) {
                second <- c(second, ordered[[2L]])
            }
            same <- values[contig_group == group]
            other <- values[contig_group != group]
            if (length(same) > 0L) {
                best_same <- c(best_same, max(same))
            }
            if (length(other) > 0L) {
                cross <- c(cross, max(other))
            }
        }
        if (length(best) < min_null) {
            next
        }
        margin <- best - second[seq_len(min(length(best), length(second)))]
        cross_margin <- if (length(cross) > 0L && length(best_same) > 0L) {
            cross[seq_len(min(length(cross), length(best_same)))] -
                best_same[seq_len(min(length(cross), length(best_same)))]
        } else {
            numeric(0)
        }
        null_rows[[i]] <- data.frame(
            key = key[[i]],
            mode = modes[[i]],
            min_depletion = as.numeric(stats::quantile(
                best, 0.95, na.rm = TRUE, names = FALSE
            )),
            strong_depletion = as.numeric(stats::quantile(
                best, 0.99, na.rm = TRUE, names = FALSE
            )),
            assignment_margin = if (length(margin) > 0L) {
                max(0, as.numeric(stats::quantile(
                    margin, 0.95, na.rm = TRUE, names = FALSE
                )))
            } else {
                NA_real_
            },
            cross_group_margin = if (length(cross_margin) > 0L) {
                max(0, as.numeric(stats::quantile(
                    cross_margin, 0.99, na.rm = TRUE, names = FALSE
                )))
            } else {
                NA_real_
            },
            n_null = length(best),
            run_threshold = if (is.null(run_threshold)) NA_real_ else run_threshold,
            stringsAsFactors = FALSE
        )
    }
    null_rows <- null_rows[!vapply(null_rows, is.null, logical(1))]
    if (length(null_rows) == 0L) {
        return(NULL)
    }
    null <- do.call(rbind, null_rows)
    rownames(null) <- NULL

    ## Shrink each bin toward its mode-wide estimate. A handful of calls in one
    ## length bin should not set an extreme threshold.
    out <- list()
    for (current_mode in unique(null$mode)) {
        mode_rows <- null[null$mode == current_mode, , drop = FALSE]
        columns <- c(
            "min_depletion", "strong_depletion", "assignment_margin",
            "cross_group_margin", "run_threshold"
        )
        global <- lapply(columns, function(column) {
            values <- mode_rows[[column]]
            values <- values[is.finite(values)]
            if (length(values) == 0L) NA_real_ else stats::median(values)
        })
        names(global) <- columns
        for (current_key in unique(mode_rows$key)) {
            rows <- mode_rows[mode_rows$key == current_key, , drop = FALSE]
            weight <- nrow(rows) / (nrow(rows) + 20)
            shrunk <- lapply(columns, function(column) {
                values <- rows[[column]]
                values <- values[is.finite(values)]
                local_value <- if (length(values) == 0L) {
                    NA_real_
                } else {
                    stats::median(values)
                }
                if (!is.finite(local_value)) {
                    return(global[[column]])
                }
                if (!is.finite(global[[column]])) {
                    return(local_value)
                }
                weight * local_value + (1 - weight) * global[[column]]
            })
            names(shrunk) <- columns
            out[[current_key]] <- c(
                shrunk,
                list(
                    n_null = sum(rows$n_null, na.rm = TRUE),
                    mode = current_mode
                )
            )
        }
    }
    out
}

#' Per-call assignment thresholds
#' @noRd
.assignment_thresholds <- function(calls, null, class_thresholds,
                                   defaults) {
    n <- nrow(calls)
    mode <- if ("mode" %in% names(calls)) {
        as.character(calls$mode)
    } else {
        rep("all", n)
    }
    bin <- if ("length_bin" %in% names(calls)) {
        as.character(calls$length_bin)
    } else {
        rep("all", n)
    }
    mode[is.na(mode) | !nzchar(mode)] <- "all"
    bin[is.na(bin) | !nzchar(bin)] <- "all"
    key <- paste(mode, bin, sep = "|")
    result <- data.frame(
        min_depletion = rep(defaults$min_depletion, n),
        strong_depletion = rep(defaults$strong_depletion, n),
        assignment_margin = rep(defaults$assignment_margin, n),
        cross_group_margin = rep(defaults$cross_group_margin, n),
        run_threshold = rep(defaults$run_threshold, n),
        stringsAsFactors = FALSE
    )
    if (!is.null(null)) {
        for (i in seq_len(n)) {
            entry <- null[[key[[i]]]]
            if (is.null(entry)) {
                next
            }
            for (column in names(result)) {
                value <- entry[[column]]
                if (!is.null(value) && length(value) == 1L &&
                    is.finite(value)) {
                    result[[column]][[i]] <- as.numeric(value)
                }
            }
        }
    }
    result$donor_mod <- defaults$donor_mod
    result$loss_depletion <- defaults$loss_depletion
    result$total_mod <- defaults$total_mod
    result$total_high <- defaults$total_high
    if (!is.null(class_thresholds)) {
        for (i in seq_len(n)) {
            entry <- class_thresholds[[key[[i]]]]
            if (is.null(entry)) {
                next
            }
            if (!is.null(entry$donor_mod) && is.finite(entry$donor_mod)) {
                result$donor_mod[[i]] <- as.numeric(entry$donor_mod)
            }
            if (!is.null(entry$loss_mod) && is.finite(entry$loss_mod)) {
                result$loss_depletion[[i]] <- max(
                    0.5, as.numeric(-entry$loss_mod)
                )
            }
            if (!is.null(entry$total_mod) && is.finite(entry$total_mod)) {
                result$total_mod[[i]] <- as.numeric(entry$total_mod)
            }
            if (!is.null(entry$total_high) && is.finite(entry$total_high)) {
                result$total_high[[i]] <- as.numeric(entry$total_high)
            }
        }
    }
    result
}

#' Weighted interval scheduling
#'
#' Chooses a non-overlapping subset of candidate runs on one contig that
#' maximises the total score. Returns the row indices of the chosen runs.
#' @noRd
.assignment_interval_schedule <- function(start, end, score) {
    n <- length(start)
    if (n == 0L) {
        return(integer())
    }
    keep <- is.finite(start) & is.finite(end) & is.finite(score) & score > 0
    if (!any(keep)) {
        return(integer())
    }
    idx <- which(keep)
    start <- as.integer(start[idx])
    end <- as.integer(end[idx])
    score <- as.numeric(score[idx])
    order_end <- order(end, start)
    idx <- idx[order_end]
    start <- start[order_end]
    end <- end[order_end]
    score <- score[order_end]
    predecessor <- integer(length(idx))
    for (i in seq_along(idx)) {
        candidates <- which(end < start[[i]])
        predecessor[[i]] <- if (length(candidates) == 0L) {
            0L
        } else {
            max(candidates)
        }
    }
    best <- numeric(length(idx) + 1L)
    choice <- integer(length(idx) + 1L)
    for (i in seq_along(idx)) {
        take <- score[[i]] + best[[predecessor[[i]] + 1L]]
        skip <- best[[i]]
        if (take >= skip) {
            best[[i + 1L]] <- take
            choice[[i + 1L]] <- i
        } else {
            best[[i + 1L]] <- skip
            choice[[i + 1L]] <- 0L
        }
    }
    chosen <- integer()
    i <- length(idx)
    while (i > 0L) {
        if (choice[[i + 1L]] == i) {
            chosen <- c(chosen, idx[[i]])
            i <- predecessor[[i]]
        } else {
            i <- i - 1L
        }
    }
    sort(chosen)
}
