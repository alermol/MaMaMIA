## Assign a group-level donor segment to its most probable recipient position.
##
## Detection and calibration stay group-level: a call says that a chromosome
## group carries donor-like content. This layer asks which recipient contig is
## the most likely host. It compares every call with every recipient contig
## (all-vs-all), using a cohort breadth baseline when several samples are
## available and a same-contig outside baseline for a single sample. Candidate
## scores combine depth depletion, a deletion-run width match and a same-group
## homology bonus. A per-sample assignment then chooses at most one call per
## recipient contig (with a second pass for extra non-overlapping segments) and
## a dummy "no host" option, so an addition is not forced onto a normal
## chromosome. A cross-group choice is reported as a translocation candidate
## and is only marked assigned when the call's own group also has a reliable
## candidate.

.assignment_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAAssignmentError", fmt, ...)
}

#' Pick the donor component recorded by the score
#' @noRd
.assignment_donor_component <- function(scores, sample, group, start_iid,
                                         end_iid) {
    rows <- scores[
        as.character(scores$sample) == sample &
            as.character(scores$group) == group &
            scores$iid >= start_iid & scores$iid <= end_iid,
        ,
        drop = FALSE
    ]
    donors <- as.character(rows$donor_component)
    donors <- donors[!is.na(donors) & nzchar(donors)]
    if (length(donors) == 0L) {
        return(NA_character_)
    }
    counts <- sort(table(donors), decreasing = TRUE)
    names(counts)[[1L]]
}

#' Resolve the donor contig and its coordinates for one call
#' @noRd
.assignment_donor_contig <- function(call, scores, assembly, table) {
    sample <- as.character(call$sample[[1L]])
    group <- as.character(call$group[[1L]])
    start_iid <- as.integer(call$start_iid[[1L]])
    end_iid <- as.integer(call$end_iid[[1L]])
    donor_component <- .assignment_donor_component(
        scores, sample, group, start_iid, end_iid
    )
    donor_contigs <- assembly$contigs$donor
    if (is.null(donor_contigs) || nrow(donor_contigs) == 0L) {
        return(list(
            component = donor_component,
            chr_id = NA_character_,
            chr_name = NA_character_,
            length = NA_real_,
            start_bp = as.numeric(call$start_bp[[1L]]),
            end_bp = as.numeric(call$end_bp[[1L]])
        ))
    }
    rows <- donor_contigs[
        as.character(donor_contigs$group) == group &
            (!is.na(donor_component) &
                as.character(donor_contigs$component) == donor_component),
        ,
        drop = FALSE
    ]
    if (nrow(rows) == 0L) {
        rows <- donor_contigs[
            as.character(donor_contigs$group) == group, ,
            drop = FALSE
        ]
    }
    if (nrow(rows) == 0L) {
        return(list(
            component = donor_component,
            chr_id = NA_character_,
            chr_name = NA_character_,
            length = NA_real_,
            start_bp = as.numeric(call$start_bp[[1L]]),
            end_bp = as.numeric(call$end_bp[[1L]])
        ))
    }
    chr_id <- as.character(rows$id[[1L]])
    chr_name <- as.character(rows$name[[1L]])
    component <- as.character(rows$component[[1L]])
    if (is.na(donor_component)) {
        donor_component <- component
    }
    windows <- table[
        as.character(table$chr_id) == chr_id &
            table$iid >= start_iid & table$iid <= end_iid,
        ,
        drop = FALSE
    ]
    length_bp <- suppressWarnings(max(
        as.numeric(rows$length[[1L]]),
        as.numeric(table$stop[as.character(table$chr_id) == chr_id]),
        na.rm = TRUE
    ))
    if (!is.finite(length_bp) || length_bp <= 0) {
        length_bp <- as.numeric(call$width_bp[[1L]])
    }
    list(
        component = donor_component,
        chr_id = chr_id,
        chr_name = chr_name,
        length = length_bp,
        start_bp = if (nrow(windows) > 0L) {
            min(as.numeric(windows$start), na.rm = TRUE)
        } else {
            as.numeric(call$start_bp[[1L]])
        },
        end_bp = if (nrow(windows) > 0L) {
            max(as.numeric(windows$stop), na.rm = TRUE)
        } else {
            as.numeric(call$end_bp[[1L]])
        }
    )
}

#' Hungarian assignment for a rectangular score matrix
#'
#' Maximises the total score when every row is assigned to a distinct column.
#' The caller pads with dummy columns so that the matrix has at least as many
#' columns as rows. Returns one column index per row.
#' @noRd
.assignment_hungarian <- function(score) {
    score <- as.matrix(score)
    n <- nrow(score)
    m <- ncol(score)
    if (n == 0L) {
        return(integer())
    }
    if (n > m) {
        .assignment_error("The assignment matrix needs at least as many columns as rows")
    }
    max_score <- max(score, na.rm = TRUE)
    if (!is.finite(max_score)) {
        max_score <- 0
    }
    cost <- max_score - score
    cost[!is.finite(cost)] <- 1e9
    u <- rep(0, n + 1)
    v <- rep(0, m + 1)
    p <- integer(m + 1)
    way <- integer(m + 1)
    for (i in seq_len(n)) {
        p[[1L]] <- i
        j0 <- 1L
        minv <- rep(Inf, m + 1)
        used <- rep(FALSE, m + 1)
        repeat {
            used[[j0]] <- TRUE
            i0 <- p[[j0]]
            delta <- Inf
            j1 <- 0L
            for (j in seq_len(m)) {
                if (!used[[j + 1L]]) {
                    cur <- cost[i0, j] - u[[i0 + 1L]] - v[[j + 1L]]
                    if (cur < minv[[j + 1L]]) {
                        minv[[j + 1L]] <- cur
                        way[[j + 1L]] <- j0
                    }
                    if (minv[[j + 1L]] < delta) {
                        delta <- minv[[j + 1L]]
                        j1 <- j + 1L
                    }
                }
            }
            for (j in seq_len(m + 1L)) {
                if (used[[j]]) {
                    u[[p[[j]] + 1L]] <- u[[p[[j]] + 1L]] + delta
                    v[[j]] <- v[[j]] - delta
                } else {
                    minv[[j]] <- minv[[j]] - delta
                }
            }
            j0 <- j1
            if (p[[j0]] == 0L) {
                break
            }
        }
        repeat {
            j1 <- way[[j0]]
            p[[j0]] <- p[[j1]]
            j0 <- j1
            if (j0 == 1L) {
                break
            }
        }
    }
    assignment <- integer(n)
    for (j in seq_len(m)) {
        if (p[[j + 1L]] != 0L) {
            assignment[[p[[j + 1L]]]] <- j
        }
    }
    assignment
}

#' Cohort breadth baseline for one sample
#'
#' For every recipient contig and window, takes a high quantile of breadth
#' across the other samples. The current sample is excluded so an event in one
#' sample cannot raise its own expected baseline. A quantile rather than the
#' maximum keeps one outlier sample from inflating the expectation.
#' @noRd
.assignment_cohort_baseline <- function(recipient, sample, quantile = 0.75,
                                        min_others = 2L) {
    samples <- unique(as.character(recipient$sample))
    others <- setdiff(samples, sample)
    if (length(others) < min_others) {
        return(NULL)
    }
    rows <- recipient[as.character(recipient$sample) %in% others, ,
        drop = FALSE
    ]
    split_rows <- split(rows, as.character(rows$chr_id))
    out <- lapply(split_rows, function(x) {
        iids <- sort(unique(x$iid))
        values <- vapply(others, function(other) {
            y <- x[as.character(x$sample) == other, , drop = FALSE]
            if (nrow(y) == 0L) {
                return(rep(NA_real_, length(iids)))
            }
            y$breadth[match(iids, y$iid)]
        }, numeric(length(iids)))
        if (is.null(dim(values))) {
            values <- matrix(values, ncol = 1L)
        }
        expected <- suppressWarnings(apply(
            values, 1L, stats::quantile,
            probs = quantile, na.rm = TRUE, names = FALSE
        ))
        list(iid = iids, expected = as.numeric(expected))
    })
    out
}

#' Sample median recipient breadth
#' @noRd
.assignment_sample_median <- function(recipient, sample) {
    x <- as.numeric(recipient$breadth[
        as.character(recipient$sample) == sample
    ])
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0L) {
        return(NA_real_)
    }
    stats::median(x, na.rm = TRUE)
}

#' Candidate scores for one call against every recipient contig
#' @noRd
.assignment_call_candidates <- function(
    call, recipient, baseline, sample_median, contigs, donor,
    min_expected_breadth, relative_floor, run_threshold, min_depletion,
    strong_depletion, width_lambda, whole_mismatch_penalty,
    same_group_bonus, eps
) {
    sample <- as.character(call$sample[[1L]])
    group <- as.character(call$group[[1L]])
    start_iid <- as.integer(call$start_iid[[1L]])
    end_iid <- as.integer(call$end_iid[[1L]])
    if (is.na(start_iid) || is.na(end_iid)) {
        return(NULL)
    }
    donor_width_bp <- as.numeric(call$width_bp[[1L]])
    if (!is.finite(donor_width_bp) || donor_width_bp <= 0) {
        donor_width_bp <- as.numeric(donor$end_bp) - as.numeric(donor$start_bp)
    }
    donor_fraction <- if (is.finite(donor$length) && donor$length > 0) {
        donor_width_bp / donor$length
    } else {
        1
    }
    floor <- max(
        as.numeric(min_expected_breadth),
        as.numeric(relative_floor) * as.numeric(sample_median),
        na.rm = TRUE
    )
    if (!is.finite(floor) || floor <= 0) {
        floor <- as.numeric(min_expected_breadth)
    }

    rows <- lapply(seq_len(nrow(contigs)), function(j) {
        chr_id <- as.character(contigs$chr_id[[j]])
        x <- recipient[
            as.character(recipient$sample) == sample &
                as.character(recipient$chr_id) == chr_id,
            ,
            drop = FALSE
        ]
        if (nrow(x) == 0L) {
            return(NULL)
        }
        x <- x[order(x$iid), , drop = FALSE]
        inside <- x$iid >= start_iid & x$iid <= end_iid
        if (!any(inside)) {
            return(NULL)
        }
        base <- baseline[[chr_id]]
        if (!is.null(base)) {
            expected <- base$expected[match(x$iid, base$iid)]
            observed <- as.numeric(x$breadth)
            depletion <- log(pmax(expected, 0) + eps) -
                log(pmax(observed, 0) + eps)
            informative <- is.finite(expected) &
                expected >= min_expected_breadth
            depletion[!informative] <- 0
            inside_depletion <- mean(depletion[inside], na.rm = TRUE)
            expected_breadth <- mean(expected[inside], na.rm = TRUE)
            observed_breadth <- mean(observed[inside], na.rm = TRUE)
            active <- informative & is.finite(depletion) &
                depletion > run_threshold
            idx <- which(inside)
            lo <- min(idx)
            hi <- max(idx)
            run_found <- any(active[inside], na.rm = TRUE)
            if (run_found) {
                while (lo > 1L && isTRUE(active[[lo - 1L]])) {
                    lo <- lo - 1L
                }
                while (hi < nrow(x) && isTRUE(active[[hi + 1L]])) {
                    hi <- hi + 1L
                }
            }
            run_start_iid <- as.integer(x$iid[[lo]])
            run_end_iid <- as.integer(x$iid[[hi]])
            run_width_bp <- as.numeric(x$stop[[hi]]) - as.numeric(x$start[[lo]])
        } else {
            fraction <- mean(inside)
            outside <- !inside
            expected_breadth <- if (fraction >= 0.95 || sum(outside) < 10L) {
                sample_median
            } else {
                mean(as.numeric(x$breadth[outside]), na.rm = TRUE)
            }
            observed_breadth <- mean(
                as.numeric(x$breadth[inside]), na.rm = TRUE
            )
            inside_depletion <- log(
                pmax(expected_breadth, 0) + eps
            ) - log(pmax(observed_breadth, 0) + eps)
            run_found <- is.finite(inside_depletion) &&
                inside_depletion > run_threshold
            idx <- which(inside)
            run_start_iid <- as.integer(x$iid[[min(idx)]])
            run_end_iid <- as.integer(x$iid[[max(idx)]])
            run_width_bp <- as.numeric(x$stop[[max(idx)]]) -
                as.numeric(x$start[[min(idx)]])
        }
        if (!is.finite(inside_depletion)) {
            inside_depletion <- -Inf
        }
        recipient_length <- suppressWarnings(max(
            as.numeric(x$stop), na.rm = TRUE
        ))
        run_fraction <- if (is.finite(recipient_length) &&
            recipient_length > 0) {
            run_width_bp / recipient_length
        } else {
            1
        }
        both_whole <- donor_fraction >= 0.95 && run_fraction >= 0.95
        whole_mismatch <- xor(donor_fraction >= 0.95, run_fraction >= 0.95)
        width_penalty <- if (both_whole) {
            0
        } else if (whole_mismatch) {
            -as.numeric(whole_mismatch_penalty)
        } else {
            -as.numeric(width_lambda) *
                abs(log((donor_width_bp + eps) / (run_width_bp + eps)))
        }
        reliable <- is.finite(expected_breadth) &&
            expected_breadth >= floor &&
            is.finite(inside_depletion) &&
            inside_depletion >= min_depletion &&
            !(donor_fraction >= 0.95 && !both_whole)
        strong <- reliable && run_found &&
            expected_breadth >= 2 * min_expected_breadth &&
            inside_depletion >= strong_depletion
        same_group <- as.character(contigs$group[[j]]) == group
        score <- if (reliable) {
            inside_depletion + if (same_group) same_group_bonus else 0
        } else {
            -1e6
        }
        score <- score + if (reliable) width_penalty else 0
        data.frame(
            chr_id = chr_id,
            chr_name = as.character(contigs$chr_name[[j]]),
            component = as.character(contigs$component[[j]]),
            group = as.character(contigs$group[[j]]),
            copy_number = as.numeric(contigs$copy_number[[j]]),
            inside_depletion = as.numeric(inside_depletion),
            expected_breadth = as.numeric(expected_breadth),
            observed_breadth = as.numeric(observed_breadth),
            run_start_iid = as.integer(run_start_iid),
            run_end_iid = as.integer(run_end_iid),
            run_width_bp = as.numeric(run_width_bp),
            run_fraction = as.numeric(run_fraction),
            donor_width_bp = as.numeric(donor_width_bp),
            donor_fraction = as.numeric(donor_fraction),
            width_penalty = as.numeric(width_penalty),
            run_found = isTRUE(run_found),
            reliable = isTRUE(reliable),
            strong = isTRUE(strong),
            pass = if ("pass" %in% names(call)) {
                isTRUE(as.logical(call$pass[[1L]]))
            } else {
                TRUE
            },
            cross_group = !same_group,
            score = as.numeric(score),
            stringsAsFactors = FALSE
        )
    })
    rows <- do.call(rbind, rows)
    if (is.null(rows) || nrow(rows) == 0L) {
        return(NULL)
    }
    rows
}

#' Do two closed integer intervals overlap?
#' @noRd
.assignment_overlap <- function(a, b) {
    if (is.null(a) || is.null(b) || length(a) < 2L || length(b) < 2L) {
        return(FALSE)
    }
    max(a[[1L]], b[[1L]]) <= min(a[[2L]], b[[2L]])
}

#' Match all calls of one sample to recipient contigs
#' @noRd
.assignment_match_sample <- function(candidates, n_calls, contigs) {
    n_contigs <- nrow(contigs)
    candidate_index <- matrix(
        NA_integer_, nrow = n_calls, ncol = n_contigs
    )
    if (!is.null(candidates) && nrow(candidates) > 0L) {
        for (r in seq_len(nrow(candidates))) {
            j <- match(
                as.character(candidates$chr_id[[r]]),
                as.character(contigs$chr_id)
            )
            if (!is.na(j)) {
                candidate_index[candidates$call[[r]], j] <- r
            }
        }
    }
    priority <- if (!is.null(candidates) && "pass" %in% names(candidates)) {
        ifelse(isTRUE_vec(candidates$pass), 100, 0)
    } else if (!is.null(candidates)) {
        rep(0, nrow(candidates))
    } else {
        numeric()
    }
    call_priority <- rep(0, n_calls)
    if (length(priority) > 0L) {
        call_priority[candidates$call] <- priority
    }
    ## Passing calls are optimised lexicographically before non-passing calls:
    ## one point of a passing call must outweigh every non-passing candidate
    ## combined, otherwise a non-passing call can take the best host and push
    ## the passing call onto a worse one. Non-passing calls still fill the
    ## hosts that remain.
    pass_weight <- 1e6
    weight <- ifelse(call_priority > 0, pass_weight, 1)
    score <- matrix(-1e6, nrow = n_calls, ncol = n_contigs)
    for (i in seq_len(n_calls)) {
        ## An impossible cell must be as unattractive as an unreliable
        ## candidate on the same scale.
        score[i, ] <- -1e6 * weight[[i]]
        for (j in seq_len(n_contigs)) {
            r <- candidate_index[i, j]
            if (!is.na(r)) {
                score[i, j] <- candidates$score[[r]] * weight[[i]]
            }
        }
    }
    ## Each call has its own dummy column. A positive candidate score always
    ## beats "no host"; a negative one does not.
    dummy <- matrix(-1e6, nrow = n_calls, ncol = n_calls)
    diag(dummy) <- 0
    full <- cbind(score, dummy)
    assignment <- .assignment_hungarian(full)
    chosen <- rep(NA_integer_, n_calls)
    assigned_runs <- vector("list", n_contigs)
    for (i in seq_len(n_calls)) {
        j <- assignment[[i]]
        if (j <= n_contigs) {
            r <- candidate_index[i, j]
            if (!is.na(r)) {
                chosen[[i]] <- r
                assigned_runs[[j]] <- c(
                    assigned_runs[[j]],
                    list(c(
                        as.integer(candidates$run_start_iid[[r]]),
                        as.integer(candidates$run_end_iid[[r]])
                    ))
                )
            }
        }
    }
    ## A chromosome can host more than one segment at different positions.
    ## Give calls that lost the one-per-contig matching a second chance on a
    ## non-overlapping interval.
    dummy_rows <- which(is.na(chosen))
    if (length(dummy_rows) > 0L) {
        best_score <- vapply(dummy_rows, function(i) {
            rows <- which(candidate_index[i, ] > 0L)
            if (length(rows) == 0L) {
                return(-Inf)
            }
            r <- candidate_index[i, rows]
            max(candidates$score[r] + priority[r], na.rm = TRUE)
        }, numeric(1))
        dummy_rows <- dummy_rows[order(-best_score)]
        for (i in dummy_rows) {
            order_rows <- which(candidate_index[i, ] > 0L)
            if (length(order_rows) == 0L) {
                next
            }
            order_rows <- order_rows[order(
                -(candidates$score[candidate_index[i, order_rows]] +
                    priority[candidate_index[i, order_rows]])
            )]
            for (j in order_rows) {
                r <- candidate_index[i, j]
                if (!isTRUE(candidates$reliable[[r]]) ||
                    candidates$score[[r]] <= 0) {
                    next
                }
                run <- c(
                    as.integer(candidates$run_start_iid[[r]]),
                    as.integer(candidates$run_end_iid[[r]])
                )
                conflict <- any(vapply(
                    assigned_runs[[j]],
                    function(other) .assignment_overlap(run, other),
                    logical(1)
                ))
                if (!conflict) {
                    chosen[[i]] <- r
                    assigned_runs[[j]] <- c(assigned_runs[[j]], list(run))
                    break
                }
            }
        }
    }
    chosen
}

#' Estimate how many recipient copies an assigned donor segment replaced
#'
#' The recipient breadth is converted to a per-copy depth scale with
#' `lambda = -log(1 - breadth)`. Losing `k` of `n` copies gives
#' `lambda_observed / lambda_expected = (n - k) / n`. The estimate is only
#' reported for an assigned call with donor support; a recipient deletion
#' without donor signal and a copy-neutral addition have no dosage.
#' @noRd
.assignment_dosage <- function(assignment, donor_mean, expected_breadth,
                               observed_breadth, copy_number,
                               donor_cutoff = 0.5,
                               expected_floor = 0.005,
                               eps = 1e-6) {
    unavailable <- list(
        affected_fraction = NA_real_,
        affected_copies = NA_integer_,
        confidence = "unavailable"
    )
    if (!identical(assignment, "assigned") ||
        !is.finite(donor_mean) || donor_mean <= donor_cutoff ||
        !is.finite(expected_breadth) || expected_breadth <= 0 ||
        !is.finite(observed_breadth) || observed_breadth < 0 ||
        !is.finite(copy_number) || copy_number < 1) {
        return(unavailable)
    }
    copies_total <- as.integer(round(copy_number))
    if (copies_total < 1L) {
        return(unavailable)
    }
    lambda <- function(breadth) {
        -log(pmax(1 - pmin(as.numeric(breadth), 1 - eps), eps))
    }
    expected_lambda <- lambda(expected_breadth)
    observed_lambda <- lambda(observed_breadth)
    if (!is.finite(expected_lambda) || expected_lambda <= 0 ||
        !is.finite(observed_lambda)) {
        return(unavailable)
    }
    fraction <- 1 - observed_lambda / expected_lambda
    fraction <- min(max(fraction, 0), 1)
    affected <- copies_total * fraction
    copies <- as.integer(round(affected))
    distance <- abs(affected - copies)
    if (fraction < 0.25) {
        return(list(
            affected_fraction = fraction,
            affected_copies = NA_integer_,
            confidence = "unavailable"
        ))
    }
    if (distance > 0.35) {
        return(list(
            affected_fraction = fraction,
            affected_copies = NA_integer_,
            confidence = "low"
        ))
    }
    confidence <- if (donor_mean >= 1 && expected_breadth >= expected_floor &&
        distance <= 0.2) {
        "high"
    } else {
        "moderate"
    }
    list(
        affected_fraction = fraction,
        affected_copies = copies,
        confidence = confidence
    )
}

#' Add assignment columns without clearing existing values
#' @noRd
.assignment_empty <- function(calls) {
    n <- nrow(calls)
    if (!"class_group" %in% names(calls)) {
        calls$class_group <- if ("class" %in% names(calls)) {
            as.character(calls$class)
        } else {
            rep(NA_character_, n)
        }
    } else if ("class" %in% names(calls)) {
        ## Re-running the assignment must start from the group-level class,
        ## not from a previous assignment-aware class.
        calls$class <- calls$class_group
    }
    if (!"class_confidence" %in% names(calls)) {
        calls$class_confidence <- if ("confidence" %in% names(calls)) {
            as.character(calls$confidence)
        } else {
            rep(NA_character_, n)
        }
    }
    if (!"class_evidence" %in% names(calls)) {
        calls$class_evidence <- rep(NA_character_, n)
    }
    if (!"class_status" %in% names(calls)) {
        calls$class_status <- rep("needs_review", n)
    }
    calls$recipient_component <- rep(NA_character_, n)
    calls$recipient_chr_id <- rep(NA_character_, n)
    calls$recipient_chr_name <- rep(NA_character_, n)
    calls$recipient_group <- rep(NA_character_, n)
    calls$recipient_start_bp <- rep(NA_real_, n)
    calls$recipient_end_bp <- rep(NA_real_, n)
    calls$recipient_run_start_bp <- rep(NA_real_, n)
    calls$recipient_run_end_bp <- rep(NA_real_, n)
    calls$recipient_deletion_width_bp <- rep(NA_real_, n)
    calls$recipient_expected_breadth <- rep(NA_real_, n)
    calls$recipient_observed_breadth <- rep(NA_real_, n)
    calls$recipient_copy_number <- rep(NA_real_, n)
    calls$recipient_mean_depletion <- rep(NA_real_, n)
    calls$affected_fraction <- rep(NA_real_, n)
    calls$affected_copies <- rep(NA_integer_, n)
    calls$dosage_confidence <- rep("unavailable", n)
    calls$recipient_weight <- rep(NA_real_, n)
    calls$recipient_margin <- rep(NA_real_, n)
    calls$recipient_assignment <- rep(NA_character_, n)
    calls$recipient_cross_group <- rep(FALSE, n)
    calls$recipient_reliable <- rep(FALSE, n)
    calls$recipient_strong <- rep(FALSE, n)
    calls$recipient_alternative <- rep(NA_character_, n)
    calls$recipient_alternative_score <- rep(NA_real_, n)
    calls$donor_component <- rep(NA_character_, n)
    calls$donor_chr_id <- rep(NA_character_, n)
    calls$donor_chr_name <- rep(NA_character_, n)
    calls$donor_start_bp <- rep(NA_real_, n)
    calls$donor_end_bp <- rep(NA_real_, n)
    calls
}

#' Assign calls to recipient contigs
#'
#' Adds the assigned recipient contig and coordinates, the assignment weight,
#' margin and flag, the donor contig and coordinates, and the diagnostics used
#' by the karyotype. The original group-level class is kept in `class_group`;
#' `class` may change to `substitution` when the assigned recipient chromosome
#' carries a strong deletion.
#'
#' @param calls A Calls data frame.
#' @param metrics A WindowMetrics data frame.
#' @param windows An AnalysisWindows list.
#' @param scores A ScoreTable data frame.
#' @param assembly The resolved assembly list with `contigs`.
#' @param same_group_bonus Score bonus for a recipient contig in the call group.
#' @param min_expected_breadth Absolute baseline floor.
#' @param relative_floor Baseline floor relative to the sample median breadth.
#' @param run_threshold Depletion threshold for a deletion run.
#' @param min_depletion Minimum depletion for a candidate host.
#' @param strong_depletion Depletion required for a confident assignment.
#' @param width_lambda Penalty per log width mismatch.
#' @param whole_mismatch_penalty Penalty when only one side is full-length.
#' @param assignment_margin Minimum margin over the next candidate.
#' @param cross_group_margin Minimum margin over the best same-group candidate.
#' @param substitution_cutoff Depletion above which an assigned call is a
#'   substitution.
#' @param addition_cutoff Donor mean above which a call is an addition.
#' @param eps Depth offset.
#'
#' @return The Calls table with assignment columns.
#' @noRd
assignCallPositions <- function(
    calls, metrics, windows, scores, assembly,
    same_group_bonus = 1, min_expected_breadth = 0.001,
    relative_floor = 0.1, run_threshold = 1, min_depletion = 0.5,
    strong_depletion = 1.5, width_lambda = 2, whole_mismatch_penalty = 8,
    assignment_margin = 0.5, cross_group_margin = 2,
    substitution_cutoff = 1.5, addition_cutoff = 0.5, eps = 1e-6
) {
    if (!is.data.frame(calls)) {
        .assignment_error("`calls` must be a data frame")
    }
    calls <- .assignment_empty(calls)
    if (nrow(calls) == 0L) {
        return(calls)
    }
    if (!is.data.frame(metrics) ||
        !all(c("sample", "chr_id", "chr_name", "group", "role",
               "component", "iid", "lambda_breadth", "gc", "breadth") %in%
            names(metrics))) {
        .assignment_error("`metrics` must be a WindowMetrics table")
    }
    if (!is.list(windows) || is.null(windows$table)) {
        .assignment_error("`windows` must be an AnalysisWindows list")
    }
    if (!is.data.frame(scores) ||
        !all(c("sample", "group", "iid", "donor_component") %in%
            names(scores))) {
        .assignment_error("`scores` must be a ScoreTable with donor_component")
    }
    recipient <- metrics[
        as.character(metrics$role) == "recipient", , drop = FALSE
    ]
    if (nrow(recipient) == 0L) {
        return(calls)
    }
    table <- windows$table
    mode <- unique(as.character(scores$mode))
    if (length(mode) != 1L || !mode %in% c("single", "multi")) {
        mode <- if (length(unique(recipient$sample)) > 1L) "multi" else "single"
    }
    contigs <- assembly$contigs$recipient
    if (is.null(contigs) || nrow(contigs) == 0L) {
        contigs <- unique(recipient[, c(
            "chr_id", "chr_name", "component", "group"
        )])
    } else {
        contigs <- data.frame(
            chr_id = as.character(contigs$id),
            chr_name = as.character(contigs$name),
            component = as.character(contigs$component),
            group = as.character(contigs$group),
            copy_number = if ("copy_number" %in% names(contigs)) {
                as.numeric(contigs$copy_number)
            } else {
                rep(2, nrow(contigs))
            },
            stringsAsFactors = FALSE
        )
    }
    rownames(contigs) <- NULL

    for (sample in unique(as.character(calls$sample))) {
        idx <- which(as.character(calls$sample) == sample)
        if (length(idx) == 0L) {
            next
        }
        baseline <- if (identical(mode, "multi")) {
            .assignment_cohort_baseline(recipient, sample)
        } else {
            NULL
        }
        sample_median <- .assignment_sample_median(recipient, sample)
        donor <- lapply(idx, function(i) {
            .assignment_donor_contig(
                calls[i, , drop = FALSE], scores, assembly, table
            )
        })
        candidates <- lapply(seq_along(idx), function(k) {
            i <- idx[[k]]
            rows <- .assignment_call_candidates(
                calls[i, , drop = FALSE], recipient, baseline,
                sample_median, contigs, donor[[k]],
                min_expected_breadth, relative_floor, run_threshold,
                min_depletion, strong_depletion, width_lambda,
                whole_mismatch_penalty, same_group_bonus, eps
            )
            if (!is.null(rows) && nrow(rows) > 0L) {
                rows$call <- i
            }
            rows
        })
        candidates <- do.call(rbind, candidates)
        if (!is.null(candidates) && nrow(candidates) > 0L) {
            candidates$call <- match(candidates$call, idx)
        }
        chosen <- .assignment_match_sample(
            candidates, length(idx), contigs
        )
        for (k in seq_along(idx)) {
            i <- idx[[k]]
            call <- calls[i, , drop = FALSE]
            donor_info <- donor[[k]]
            calls$donor_component[[i]] <- donor_info$component
            calls$donor_chr_id[[i]] <- donor_info$chr_id
            calls$donor_chr_name[[i]] <- donor_info$chr_name
            calls$donor_start_bp[[i]] <- donor_info$start_bp
            calls$donor_end_bp[[i]] <- donor_info$end_bp

            rows <- if (!is.null(candidates) && nrow(candidates) > 0L) {
                which(candidates$call == k)
            } else {
                integer()
            }
            best_row <- chosen[[k]]
            best <- if (!is.na(best_row) && length(rows) > 0L) {
                candidates[best_row, , drop = FALSE]
            } else {
                NULL
            }
            real_rows <- rows[candidates$score[rows] > -1e5]
            other_rows <- setdiff(real_rows, best_row)
            alternative <- if (length(other_rows) > 0L) {
                other_rows[which.max(candidates$score[other_rows])]
            } else {
                NA_integer_
            }
            same_group_rows <- real_rows[
                as.character(candidates$group[real_rows]) ==
                    as.character(call$group[[1L]])
            ]
            same_group_best <- if (length(same_group_rows) > 0L) {
                max(candidates$score[same_group_rows], na.rm = TRUE)
            } else {
                -Inf
            }

            if (is.null(best)) {
                calls$recipient_assignment[[i]] <- "no_recipient_loss"
                calls$recipient_weight[[i]] <- 0
                calls$recipient_margin[[i]] <- NA_real_
                calls$recipient_cross_group[[i]] <- FALSE
                calls$recipient_reliable[[i]] <- FALSE
                calls$recipient_strong[[i]] <- FALSE
                if (!is.na(alternative)) {
                    calls$recipient_alternative[[i]] <-
                        as.character(candidates$chr_name[[alternative]])
                    calls$recipient_alternative_score[[i]] <-
                        as.numeric(candidates$score[[alternative]])
                }
                calls$class[[i]] <- calls$class_group[[i]]
                next
            }

            call_interval <- c(
                as.integer(call$start_iid[[1L]]),
                as.integer(call$end_iid[[1L]])
            )
            recipient_rows <- table[
                as.character(table$chr_id) ==
                    as.character(best$chr_id[[1L]]) &
                    table$iid >= call_interval[[1L]] &
                    table$iid <= call_interval[[2L]],
                ,
                drop = FALSE
            ]
            run_rows <- table[
                as.character(table$chr_id) ==
                    as.character(best$chr_id[[1L]]) &
                    table$iid >= as.integer(best$run_start_iid[[1L]]) &
                    table$iid <= as.integer(best$run_end_iid[[1L]]),
                ,
                drop = FALSE
            ]
            weights <- exp(
                candidates$score[real_rows] -
                    max(candidates$score[real_rows], na.rm = TRUE)
            )
            weights <- weights / sum(weights, na.rm = TRUE)
            weight <- weights[match(best_row, real_rows)]
            if (length(weight) == 0L || !is.finite(weight)) {
                weight <- 0
            }
            alternatives <- c(
                candidates$score[other_rows],
                0
            )
            margin <- if (length(alternatives) > 0L) {
                as.numeric(best$score[[1L]]) - max(alternatives, na.rm = TRUE)
            } else {
                NA_real_
            }
            is_cross <- isTRUE(best$cross_group[[1L]])
            assignment <- if (!isTRUE(best$reliable[[1L]])) {
                "ambiguous"
            } else if (is_cross && !is.finite(same_group_best)) {
                ## A cross-group host without any reliable candidate in the
                ## call's own group is a translocation candidate, not a
                ## confirmed position. This keeps a copy-neutral event from
                ## being "assigned" to an unrelated deletion.
                "ambiguous"
            } else if (is_cross &&
                (as.numeric(best$score[[1L]]) - same_group_best) <
                    cross_group_margin) {
                "ambiguous"
            } else if (!isTRUE(best$strong[[1L]])) {
                "ambiguous"
            } else if (is.finite(margin) && margin < assignment_margin) {
                "ambiguous"
            } else {
                "assigned"
            }

            calls$recipient_component[[i]] <- best$component[[1L]]
            calls$recipient_chr_id[[i]] <- best$chr_id[[1L]]
            calls$recipient_chr_name[[i]] <- best$chr_name[[1L]]
            calls$recipient_group[[i]] <- best$group[[1L]]
            calls$recipient_start_bp[[i]] <- if (nrow(recipient_rows) > 0L) {
                min(as.numeric(recipient_rows$start), na.rm = TRUE)
            } else {
                NA_real_
            }
            calls$recipient_end_bp[[i]] <- if (nrow(recipient_rows) > 0L) {
                max(as.numeric(recipient_rows$stop), na.rm = TRUE)
            } else {
                NA_real_
            }
            calls$recipient_run_start_bp[[i]] <- if (nrow(run_rows) > 0L) {
                min(as.numeric(run_rows$start), na.rm = TRUE)
            } else {
                NA_real_
            }
            calls$recipient_run_end_bp[[i]] <- if (nrow(run_rows) > 0L) {
                max(as.numeric(run_rows$stop), na.rm = TRUE)
            } else {
                NA_real_
            }
            calls$recipient_deletion_width_bp[[i]] <-
                if (isTRUE(best$run_found[[1L]])) {
                    as.numeric(best$run_width_bp[[1L]])
                } else {
                    NA_real_
                }
            calls$recipient_expected_breadth[[i]] <-
                as.numeric(best$expected_breadth[[1L]])
            calls$recipient_mean_depletion[[i]] <-
                as.numeric(best$inside_depletion[[1L]])
            calls$recipient_weight[[i]] <- as.numeric(weight)
            calls$recipient_margin[[i]] <- margin
            calls$recipient_assignment[[i]] <- assignment
            calls$recipient_cross_group[[i]] <- is_cross
            calls$recipient_reliable[[i]] <- isTRUE(best$reliable[[1L]])
            calls$recipient_strong[[i]] <- isTRUE(best$strong[[1L]])
            if (!is.na(alternative)) {
                calls$recipient_alternative[[i]] <-
                    as.character(candidates$chr_name[[alternative]])
                calls$recipient_alternative_score[[i]] <-
                    as.numeric(candidates$score[[alternative]])
            }
            dosage <- .assignment_dosage(
                assignment = assignment,
                donor_mean = as.numeric(call$donor_mean[[1L]]),
                expected_breadth = as.numeric(best$expected_breadth[[1L]]),
                observed_breadth = as.numeric(best$observed_breadth[[1L]]),
                copy_number = as.numeric(best$copy_number[[1L]])
            )
            calls$recipient_observed_breadth[[i]] <-
                as.numeric(best$observed_breadth[[1L]])
            calls$recipient_copy_number[[i]] <-
                as.numeric(best$copy_number[[1L]])
            calls$affected_fraction[[i]] <- dosage$affected_fraction
            calls$affected_copies[[i]] <- dosage$affected_copies
            calls$dosage_confidence[[i]] <- dosage$confidence
            assigned_class <- calls$class_group[[i]]
            assigned_evidence <- calls$class_evidence[[i]]
            assigned_confidence <- calls$class_confidence[[i]]
            assigned_status <- calls$class_status[[i]]
            donor_mean <- as.numeric(call$donor_mean[[1L]])
            donor_supported <- is.finite(donor_mean) &&
                donor_mean > addition_cutoff
            local_loss <- isTRUE(best$strong[[1L]]) &&
                as.numeric(best$inside_depletion[[1L]]) >= substitution_cutoff
            width_ratio <- if (isTRUE(best$run_found[[1L]]) &&
                is.finite(best$run_width_bp[[1L]]) &&
                best$run_width_bp[[1L]] > 0 &&
                is.finite(best$donor_width_bp[[1L]]) &&
                best$donor_width_bp[[1L]] > 0) {
                as.numeric(best$donor_width_bp[[1L]]) /
                    as.numeric(best$run_width_bp[[1L]])
            } else {
                NA_real_
            }
            if (identical(assignment, "assigned")) {
                if (donor_supported && local_loss) {
                    if (is.finite(width_ratio) &&
                        (width_ratio < 0.5 || width_ratio > 2)) {
                        assigned_class <- "complex"
                        assigned_evidence <- "donor+recipient; width mismatch"
                        assigned_confidence <- "low"
                        assigned_status <- "needs_review"
                    } else {
                        assigned_class <- "substitution"
                        assigned_evidence <- "donor+recipient"
                        assigned_confidence <- if (isTRUE(best$strong[[1L]]) &&
                            donor_mean >= 1) {
                            "high"
                        } else {
                            "moderate"
                        }
                        assigned_status <- if (assigned_confidence == "high") {
                            "confident"
                        } else {
                            "ambiguous"
                        }
                    }
                } else if (donor_supported) {
                    assigned_class <- "addition"
                    assigned_evidence <- "donor"
                    assigned_confidence <- if (donor_mean >= 1) {
                        "high"
                    } else {
                        "moderate"
                    }
                    assigned_status <- if (assigned_confidence == "high") {
                        "confident"
                    } else {
                        "ambiguous"
                    }
                } else if (local_loss) {
                    assigned_class <- "deletion"
                    assigned_evidence <- "recipient"
                    assigned_confidence <- if (isTRUE(best$strong[[1L]])) {
                        "moderate"
                    } else {
                        "low"
                    }
                    assigned_status <- "ambiguous"
                }
            }
            calls$class[[i]] <- assigned_class
            calls$class_evidence[[i]] <- assigned_evidence
            calls$class_confidence[[i]] <- assigned_confidence
            calls$class_status[[i]] <- assigned_status
        }
    }
    calls
}
