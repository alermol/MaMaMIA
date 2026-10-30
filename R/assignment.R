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

#' Mean of a score column over one call interval
#' @noRd
.assignment_segment_mean <- function(scores, sample, group, start_iid,
                                     end_iid, column) {
    if (!column %in% names(scores)) {
        return(NA_real_)
    }
    rows <- as.character(scores$sample) == as.character(sample) &
        as.character(scores$group) == as.character(group) &
        as.integer(scores$iid) >= as.integer(start_iid) &
        as.integer(scores$iid) <= as.integer(end_iid)
    values <- as.numeric(scores[[column]][rows])
    values <- values[is.finite(values)]
    if (length(values) == 0L) {
        return(NA_real_)
    }
    mean(values)
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

#' Candidate scores for one call against every recipient contig
#' @noRd
.assignment_call_candidates <- function(
    call, recipient, contigs, donor,
    min_expected_breadth, relative_floor, run_threshold, min_depletion,
    strong_depletion, width_lambda, whole_mismatch_penalty,
    same_group_bonus, eps, sample_median = NULL
) {
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
    if (is.null(sample_median)) {
        sample_breadth <- if (is.list(recipient) &&
            !is.data.frame(recipient)) {
            unlist(lapply(recipient, function(x) {
                x$breadth[!x$masked & is.finite(x$breadth)]
            }), use.names = FALSE)
        } else {
            recipient$breadth[
                !recipient$masked & is.finite(recipient$breadth)
            ]
        }
        sample_median <- if (length(sample_breadth) > 0L) {
            stats::median(sample_breadth, na.rm = TRUE)
        } else {
            NA_real_
        }
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
        x <- if (is.list(recipient) && !is.data.frame(recipient)) {
            recipient[[chr_id]]
        } else {
            recipient[
                as.character(recipient$chr_id) == chr_id,
                ,
                drop = FALSE
            ]
        }
        if (is.null(x) || nrow(x) == 0L) {
            return(NULL)
        }
        x <- x[order(x$iid), , drop = FALSE]
        inside <- x$iid >= start_iid & x$iid <= end_iid
        usable <- inside & !x$masked & is.finite(x$depletion)
        if (!any(usable)) {
            return(NULL)
        }
        stats <- .assignment_interval_stats(x$depletion, usable)
        if (!is.finite(stats$mean)) {
            return(NULL)
        }
        expected_breadth <- mean(
            x$expected_breadth[usable], na.rm = TRUE
        )
        observed_breadth <- mean(x$breadth[usable], na.rm = TRUE)
        mapq_fraction <- if ("mapq_high_fraction" %in% names(x)) {
            values <- x$mapq_high_fraction[usable]
            values <- values[is.finite(values)]
            if (length(values) == 0L) NA_real_ else mean(values)
        } else {
            NA_real_
        }
        active <- !x$masked & is.finite(x$depletion) &
            x$depletion > run_threshold
        idx <- which(inside & !x$masked)
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
        contig_length <- suppressWarnings(max(
            as.numeric(x$contig_length), as.numeric(x$stop),
            na.rm = TRUE
        ))
        run_fraction <- if (is.finite(contig_length) && contig_length > 0) {
            run_width_bp / contig_length
        } else {
            1
        }
        both_whole <- donor_fraction >= 0.95 && run_fraction >= 0.95
        whole_mismatch <- xor(donor_fraction >= 0.95, run_fraction >= 0.95)
        n_per <- max(as.integer(x$iid), na.rm = TRUE) + 1L
        tolerance <- if (is.finite(n_per) && n_per > 0) 1 / n_per else 0
        donor_start <- donor$start_bp / donor$length
        donor_end <- donor$end_bp / donor$length
        if (!is.finite(donor_start) || !is.finite(donor_end) ||
            donor_end <= donor_start || !is.finite(donor$length) ||
            donor$length <= 0) {
            donor_start <- 0
            donor_end <- donor_fraction
        }
        run_start <- as.numeric(x$start[[lo]]) / contig_length
        run_end <- as.numeric(x$stop[[hi]]) / contig_length
        if (!is.finite(run_start) || !is.finite(run_end) ||
            run_end <= run_start) {
            run_start <- 0
            run_end <- run_fraction
        }
        d0 <- max(0, donor_start - tolerance)
        d1 <- min(1, donor_end + tolerance)
        r0 <- max(0, run_start - tolerance)
        r1 <- min(1, run_end + tolerance)
        denominator <- max(d1 - d0, r1 - r0)
        reciprocal_overlap <- if (is.finite(denominator) &&
            denominator > 0) {
            max(0, min(d1, r1) - max(d0, r0)) / denominator
        } else {
            0
        }
        fraction_ratio <- if (is.finite(run_fraction) &&
            run_fraction > 0) {
            donor_fraction / run_fraction
        } else {
            NA_real_
        }
        width_penalty <- if (both_whole) {
            0
        } else if (whole_mismatch) {
            -as.numeric(whole_mismatch_penalty)
        } else if (is.finite(fraction_ratio) && fraction_ratio > 0) {
            -as.numeric(width_lambda) * abs(log(fraction_ratio))
        } else {
            -as.numeric(whole_mismatch_penalty)
        }
        reliable <- is.finite(expected_breadth) &&
            expected_breadth >= floor &&
            is.finite(stats$lower) &&
            stats$lower >= min_depletion &&
            !(donor_fraction >= 0.95 && !both_whole)
        strong <- reliable && run_found &&
            expected_breadth >= 2 * min_expected_breadth &&
            stats$lower >= strong_depletion
        same_group <- as.character(contigs$group[[j]]) == group
        score <- if (reliable) {
            stats$lower + if (same_group) same_group_bonus else 0
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
            inside_depletion = as.numeric(stats$mean),
            inside_se = as.numeric(stats$se),
            inside_lower = as.numeric(stats$lower),
            inside_z = as.numeric(stats$z),
            inside_n = as.integer(stats$n),
            expected_breadth = as.numeric(expected_breadth),
            observed_breadth = as.numeric(observed_breadth),
            mapq_fraction = as.numeric(mapq_fraction),
            run_start_iid = as.integer(run_start_iid),
            run_end_iid = as.integer(run_end_iid),
            run_width_bp = as.numeric(run_width_bp),
            run_fraction = as.numeric(run_fraction),
            donor_width_bp = as.numeric(donor_width_bp),
            donor_fraction = as.numeric(donor_fraction),
            reciprocal_overlap = as.numeric(reciprocal_overlap),
            fraction_ratio = as.numeric(fraction_ratio),
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
    ## For each contig, choose the non-overlapping subset with the highest
    ## total score, then let the next contig place the calls that remain.
    dummy_rows <- which(is.na(chosen))
    if (length(dummy_rows) > 0L) {
        contig_order <- order(vapply(seq_len(n_contigs), function(j) {
            cells <- candidate_index[dummy_rows, j]
            rows <- dummy_rows[!is.na(cells) & cells > 0L]
            if (length(rows) == 0L) {
                return(-Inf)
            }
            r <- candidate_index[rows, j]
            max(
                candidates$score[r] +
                    ifelse(priority[r] > 0, 1e6, 0),
                na.rm = TRUE
            )
        }, numeric(1)), decreasing = TRUE)
        for (j in contig_order) {
            cells <- candidate_index[dummy_rows, j]
            rows <- dummy_rows[!is.na(cells) & cells > 0L]
            rows <- rows[is.na(chosen[rows])]
            if (length(rows) == 0L) {
                next
            }
            r <- candidate_index[rows, j]
            keep <- isTRUE_vec(candidates$reliable[r]) &
                candidates$score[r] > 0
            if (length(assigned_runs[[j]]) > 0L) {
                conflict <- vapply(seq_along(r), function(k) {
                    run <- c(
                        as.integer(candidates$run_start_iid[[r[[k]]]]),
                        as.integer(candidates$run_end_iid[[r[[k]]]])
                    )
                    any(vapply(
                        assigned_runs[[j]],
                        function(other) .assignment_overlap(run, other),
                        logical(1)
                    ))
                }, logical(1))
                keep <- keep & !conflict
            }
            rows <- rows[keep]
            r <- r[keep]
            if (length(rows) == 0L) {
                next
            }
            schedule <- .assignment_interval_schedule(
                candidates$run_start_iid[r],
                candidates$run_end_iid[r],
                candidates$score[r] +
                    ifelse(priority[r] > 0, 1e6, 0)
            )
            for (s in schedule) {
                i <- rows[[s]]
                if (!is.na(chosen[[i]])) {
                    next
                }
                rr <- r[[s]]
                chosen[[i]] <- rr
                assigned_runs[[j]] <- c(
                    assigned_runs[[j]],
                    list(c(
                        as.integer(candidates$run_start_iid[[rr]]),
                        as.integer(candidates$run_end_iid[[rr]])
                    ))
                )
            }
        }
    }
    chosen
}

#' Estimate how many recipient copies an assigned donor segment replaced
#'
#' The recipient breadth is converted to a per-copy depth scale with
#' `lambda = -log(1 - breadth)`. Losing `k` of `n` copies gives
#' `lambda_observed / lambda_expected = (n - k) / n`. A donor-only gain is
#' reported from `donor_copies`, so additions are no longer always
#' `unavailable`; a recipient loss without donor support is reported from the
#' loss fraction, at lower confidence when no total-depth support is present.
#' @noRd
.assignment_dosage <- function(assignment, donor_mean, expected_breadth,
                               observed_breadth, copy_number,
                               donor_copies = NA_real_,
                               total_mean = NA_real_,
                               donor_cutoff = 0.5,
                               total_cutoff = Inf,
                               expected_floor = 0.005,
                               eps = 1e-6) {
    unavailable <- list(
        affected_fraction = NA_real_,
        affected_copies = NA_integer_,
        confidence = "unavailable"
    )
    if (!identical(assignment, "assigned") ||
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
    loss_fraction <- NA_real_
    if (is.finite(expected_breadth) && expected_breadth > 0 &&
        is.finite(observed_breadth) && observed_breadth >= 0) {
        expected_lambda <- lambda(expected_breadth)
        observed_lambda <- lambda(observed_breadth)
        if (is.finite(expected_lambda) && expected_lambda > 0 &&
            is.finite(observed_lambda)) {
            loss_fraction <- 1 - observed_lambda / expected_lambda
            loss_fraction <- min(max(loss_fraction, 0), 1)
        }
    }
    loss_supported <- is.finite(loss_fraction) && loss_fraction >= 0.25
    gain_supported <- is.finite(donor_mean) && donor_mean > donor_cutoff &&
        is.finite(donor_copies) && donor_copies >= 0.25
    total_loss <- is.finite(total_mean) && is.finite(total_cutoff) &&
        total_mean <= -total_cutoff

    if (loss_supported && gain_supported) {
        affected <- copies_total * loss_fraction
        copies <- as.integer(round(affected))
        distance <- abs(affected - copies)
        if (distance > 0.35) {
            return(list(
                affected_fraction = loss_fraction,
                affected_copies = NA_integer_,
                confidence = "low"
            ))
        }
        agreement <- abs(donor_copies - affected)
        confidence <- if (donor_mean >= 1 && agreement <= 0.5 &&
            expected_breadth >= expected_floor) {
            "high"
        } else {
            "moderate"
        }
        return(list(
            affected_fraction = loss_fraction,
            affected_copies = copies,
            confidence = confidence
        ))
    }
    if (gain_supported) {
        copies <- as.integer(round(donor_copies))
        if (copies < 1L) {
            return(unavailable)
        }
        confidence <- if (donor_mean >= 1 && donor_copies >= 1) {
            "high"
        } else {
            "moderate"
        }
        return(list(
            affected_fraction = donor_copies / copies_total,
            affected_copies = copies,
            confidence = confidence
        ))
    }
    if (loss_supported) {
        affected <- copies_total * loss_fraction
        copies <- as.integer(round(affected))
        distance <- abs(affected - copies)
        if (distance > 0.35) {
            return(list(
                affected_fraction = loss_fraction,
                affected_copies = NA_integer_,
                confidence = "low"
            ))
        }
        confidence <- if (distance <= 0.2 &&
            (donor_mean >= 1 || total_loss)) {
            "high"
        } else {
            "moderate"
        }
        return(list(
            affected_fraction = loss_fraction,
            affected_copies = copies,
            confidence = confidence
        ))
    }
    unavailable
}

#' Copy state implied by a class
#' @noRd
.assignment_copy_state_from_class <- function(class) {
    class <- as.character(class)
    ifelse(
        class == "substitution", "balanced",
        ifelse(
            class == "addition", "gain",
            ifelse(
                class == "deletion", "loss",
                ifelse(class == "complex", "compound", "neutral")
            )
        )
    )
}

#' Component heterogeneity status
#'
#' A group-level call can hide a one-component event. A large donor or
#' recipient component spread, or more than one supported donor component,
#' marks the call as heterogeneous rather than failing it.
#' @noRd
.assignment_component_status <- function(call) {
    donor_spread <- if ("donor_component_spread" %in% names(call)) {
        as.numeric(call$donor_component_spread[[1L]])
    } else {
        0
    }
    recipient_spread <- if ("recipient_component_spread" %in% names(call)) {
        as.numeric(call$recipient_component_spread[[1L]])
    } else {
        0
    }
    supported <- if ("n_donor_supported" %in% names(call)) {
        as.numeric(call$n_donor_supported[[1L]])
    } else {
        0
    }
    if ((is.finite(donor_spread) && donor_spread > 2) ||
        (is.finite(recipient_spread) && recipient_spread > 2) ||
        (is.finite(supported) && supported > 1)) {
        "heterogeneous"
    } else {
        "consistent"
    }
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
    copy_state_default <- .assignment_copy_state_from_class(calls$class_group)
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
    calls$recipient_mapq_fraction <- rep(NA_real_, n)
    calls$recipient_copy_number <- rep(NA_real_, n)
    calls$recipient_mean_depletion <- rep(NA_real_, n)
    calls$recipient_depletion_se <- rep(NA_real_, n)
    calls$recipient_depletion_lower <- rep(NA_real_, n)
    calls$recipient_depletion_z <- rep(NA_real_, n)
    calls$recipient_width_overlap <- rep(NA_real_, n)
    calls$recipient_fraction_ratio <- rep(NA_real_, n)
    calls$donor_copies <- rep(NA_real_, n)
    if (!"total_mean" %in% names(calls)) {
        calls$total_mean <- rep(NA_real_, n)
    }
    if (!"copy_state_group" %in% names(calls)) {
        calls$copy_state_group <- if ("copy_state" %in% names(calls)) {
            as.character(calls$copy_state)
        } else {
            copy_state_default
        }
    }
    calls$copy_state <- calls$copy_state_group
    if (!"total_evidence" %in% names(calls)) {
        calls$total_evidence <- rep(NA_character_, n)
    }
    calls$component_status <- rep(NA_character_, n)
    if (!"mapq_fraction" %in% names(calls)) {
        calls$mapq_fraction <- rep(NA_real_, n)
    }
    calls$mappability_status <- ifelse(
        is.finite(as.numeric(calls$mapq_fraction)) &
            as.numeric(calls$mapq_fraction) < 0.5,
        "low", "ok"
    )
    calls$mappability_status[
        !is.finite(as.numeric(calls$mapq_fraction))
    ] <- NA_character_
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
    masks = NULL, class_thresholds = NULL, controls = NULL,
    assignment_calibration = TRUE, null_permutations = 64L,
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

    recipient_residuals <- .assignment_recipient_residuals(
        recipient, mode, masks = masks, eps = eps
    )
    null <- if (isTRUE(assignment_calibration)) {
        .assignment_null_thresholds(
            calls, recipient_residuals, contigs,
            controls = controls,
            permutations = null_permutations,
            same_group_bonus = same_group_bonus
        )
    } else {
        NULL
    }
    threshold_table <- .assignment_thresholds(
        calls, null, class_thresholds,
        defaults = list(
            min_depletion = min_depletion,
            strong_depletion = strong_depletion,
            assignment_margin = assignment_margin,
            cross_group_margin = cross_group_margin,
            run_threshold = run_threshold,
            donor_mod = addition_cutoff,
            loss_depletion = substitution_cutoff,
            total_mod = Inf,
            total_high = Inf
        )
    )

    for (sample in unique(as.character(calls$sample))) {
        idx <- which(as.character(calls$sample) == sample)
        if (length(idx) == 0L) {
            next
        }
        recipient_sample <- recipient_residuals[
            as.character(recipient_residuals$sample) == sample, ,
            drop = FALSE
        ]
        sample_median <- stats::median(
            recipient_sample$breadth[
                !recipient_sample$masked &
                    is.finite(recipient_sample$breadth)
            ],
            na.rm = TRUE
        )
        recipient_sample <- split(
            recipient_sample, as.character(recipient_sample$chr_id)
        )
        donor <- lapply(idx, function(i) {
            .assignment_donor_contig(
                calls[i, , drop = FALSE], scores, assembly, table
            )
        })
        candidates <- lapply(seq_along(idx), function(k) {
            i <- idx[[k]]
            rows <- .assignment_call_candidates(
                calls[i, , drop = FALSE], recipient_sample, contigs,
                donor[[k]],
                min_expected_breadth, relative_floor,
                threshold_table$run_threshold[[i]],
                threshold_table$min_depletion[[i]],
                threshold_table$strong_depletion[[i]],
                width_lambda, whole_mismatch_penalty,
                same_group_bonus, eps, sample_median = sample_median
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
            calls$component_status[[i]] <- .assignment_component_status(call)
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
            min_depletion <- threshold_table$min_depletion[[i]]
            strong_depletion <- threshold_table$strong_depletion[[i]]
            assignment_margin <- threshold_table$assignment_margin[[i]]
            cross_group_margin <- threshold_table$cross_group_margin[[i]]
            assignment <- if (!isTRUE(best$reliable[[1L]]) ||
                !isTRUE(best$strong[[1L]])) {
                "ambiguous"
            } else if (is_cross && !is.finite(same_group_best)) {
                ## A cross-group host without any reliable candidate in the
                ## call's own group is an explicit translocation candidate.
                "cross_group"
            } else if (is_cross &&
                (as.numeric(best$score[[1L]]) - same_group_best) <
                    cross_group_margin) {
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
            calls$recipient_depletion_se[[i]] <-
                as.numeric(best$inside_se[[1L]])
            calls$recipient_depletion_lower[[i]] <-
                as.numeric(best$inside_lower[[1L]])
            calls$recipient_depletion_z[[i]] <-
                as.numeric(best$inside_z[[1L]])
            calls$recipient_width_overlap[[i]] <-
                as.numeric(best$reciprocal_overlap[[1L]])
            calls$recipient_fraction_ratio[[i]] <-
                as.numeric(best$fraction_ratio[[1L]])
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
            donor_copies <- .assignment_segment_mean(
                scores, sample, call$group[[1L]],
                call$start_iid[[1L]], call$end_iid[[1L]], "donor_copies"
            )
            total_mean <- if ("total_mean" %in% names(call)) {
                as.numeric(call$total_mean[[1L]])
            } else {
                NA_real_
            }
            calls$donor_copies[[i]] <- donor_copies
            calls$total_mean[[i]] <- total_mean
            total_mod <- threshold_table$total_mod[[i]]
            total_gain <- is.finite(total_mean) && is.finite(total_mod) &&
                total_mean > total_mod
            total_loss <- is.finite(total_mean) && is.finite(total_mod) &&
                total_mean < -total_mod
            calls$total_evidence[[i]] <- if (total_gain) {
                "gain"
            } else if (total_loss) {
                "loss"
            } else if (is.finite(total_mean) && is.finite(total_mod)) {
                "neutral"
            } else {
                NA_character_
            }
            dosage <- .assignment_dosage(
                assignment = assignment,
                donor_mean = as.numeric(call$donor_mean[[1L]]),
                expected_breadth = as.numeric(best$expected_breadth[[1L]]),
                observed_breadth = as.numeric(best$observed_breadth[[1L]]),
                copy_number = as.numeric(best$copy_number[[1L]]),
                donor_copies = donor_copies,
                total_mean = total_mean,
                donor_cutoff = threshold_table$donor_mod[[i]],
                total_cutoff = total_mod
            )
            calls$recipient_observed_breadth[[i]] <-
                as.numeric(best$observed_breadth[[1L]])
            calls$recipient_mapq_fraction[[i]] <-
                as.numeric(best$mapq_fraction[[1L]])
            if (is.finite(calls$recipient_mapq_fraction[[i]]) &&
                calls$recipient_mapq_fraction[[i]] < 0.5) {
                calls$mappability_status[[i]] <- "low"
            }
            calls$recipient_copy_number[[i]] <-
                as.numeric(best$copy_number[[1L]])
            calls$affected_fraction[[i]] <- dosage$affected_fraction
            calls$affected_copies[[i]] <- dosage$affected_copies
            calls$dosage_confidence[[i]] <- dosage$confidence
            assigned_class <- calls$class_group[[i]]
            assigned_evidence <- calls$class_evidence[[i]]
            assigned_confidence <- calls$class_confidence[[i]]
            assigned_status <- calls$class_status[[i]]
            assigned_copy_state <- calls$copy_state_group[[i]]
            donor_mean <- as.numeric(call$donor_mean[[1L]])
            donor_supported <- is.finite(donor_mean) &&
                donor_mean > threshold_table$donor_mod[[i]]
            local_loss <- isTRUE(best$strong[[1L]]) &&
                as.numeric(best$inside_lower[[1L]]) >=
                    threshold_table$loss_depletion[[i]]
            width_ok <- is.finite(best$reciprocal_overlap[[1L]]) &&
                best$reciprocal_overlap[[1L]] >= 0.5 &&
                is.finite(best$fraction_ratio[[1L]]) &&
                best$fraction_ratio[[1L]] >= 0.5 &&
                best$fraction_ratio[[1L]] <= 2
            if (identical(assignment, "assigned")) {
                if (donor_supported && local_loss) {
                    if (!width_ok) {
                        assigned_class <- "complex"
                        assigned_evidence <- "donor+recipient; width mismatch"
                        assigned_confidence <- "low"
                        assigned_status <- "needs_review"
                        assigned_copy_state <- "compound"
                    } else if (total_gain || total_loss) {
                        assigned_class <- "complex"
                        assigned_evidence <- "donor+recipient; total imbalance"
                        assigned_confidence <- "low"
                        assigned_status <- "needs_review"
                        assigned_copy_state <- "compound"
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
                        assigned_copy_state <- "balanced"
                    }
                } else if (donor_supported) {
                    if (total_loss) {
                        assigned_class <- "complex"
                        assigned_evidence <- "donor; total loss"
                        assigned_confidence <- "low"
                        assigned_status <- "needs_review"
                        assigned_copy_state <- "compound"
                    } else {
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
                        assigned_copy_state <- "gain"
                    }
                } else if (local_loss) {
                    if (total_gain) {
                        assigned_class <- "complex"
                        assigned_evidence <- "recipient; total gain"
                        assigned_confidence <- "low"
                        assigned_status <- "needs_review"
                        assigned_copy_state <- "compound"
                    } else {
                        assigned_class <- "deletion"
                        assigned_evidence <- "recipient"
                        assigned_confidence <- if (isTRUE(best$strong[[1L]]) &&
                            total_loss) {
                            "high"
                        } else if (isTRUE(best$strong[[1L]])) {
                            "moderate"
                        } else {
                            "low"
                        }
                        assigned_status <- if (assigned_confidence == "high") {
                            "confident"
                        } else {
                            "ambiguous"
                        }
                        assigned_copy_state <- "loss"
                    }
                }
            } else if (identical(assignment, "cross_group")) {
                assigned_status <- "needs_review"
            }
            calls$copy_state[[i]] <- assigned_copy_state
            if ((identical(calls$component_status[[i]], "heterogeneous") ||
                identical(calls$mappability_status[[i]], "low")) &&
                identical(assigned_confidence, "high")) {
                assigned_confidence <- "moderate"
                assigned_status <- "ambiguous"
            }
            calls$class[[i]] <- assigned_class
            calls$class_evidence[[i]] <- assigned_evidence
            calls$class_confidence[[i]] <- assigned_confidence
            calls$class_status[[i]] <- assigned_status
        }
    }
    calls
}
