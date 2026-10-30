## Call classification and confidence.

.calls_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIACallError", fmt, ...)
}

#' @noRd
.call_class_evidence <- function(donor_mean, recipient_mean,
                                 total_mean = NULL,
                                 donor_cutoff = 0.5,
                                 loss_cutoff = -0.5,
                                 total_cutoff = Inf) {
    donor_gain <- is.finite(donor_mean) & donor_mean > donor_cutoff
    recipient_loss <- is.finite(recipient_mean) & recipient_mean < loss_cutoff
    recipient_gain <- is.finite(recipient_mean) &
        recipient_mean > abs(loss_cutoff)
    if (is.null(total_mean)) {
        total_mean <- rep(NA_real_, length(donor_mean))
    }
    total_gain <- is.finite(total_mean) & is.finite(total_cutoff) &
        total_mean > total_cutoff
    total_loss <- is.finite(total_mean) & is.finite(total_cutoff) &
        total_mean < -total_cutoff
    total_neutral <- !total_gain & !total_loss
    class <- ifelse(
        donor_gain & recipient_loss, "substitution",
        ifelse(
            donor_gain & recipient_gain, "complex",
            ifelse(
                donor_gain, "addition",
                ifelse(
                    recipient_loss, "deletion",
                    ifelse(recipient_gain, "other", "other")
                )
            )
        )
    )
    evidence <- ifelse(
        donor_gain & (recipient_loss | recipient_gain), "donor+recipient",
        ifelse(
            donor_gain, "donor",
            ifelse(recipient_loss | recipient_gain, "recipient", "none")
        )
    )
    list(
        class = as.character(class),
        evidence = as.character(evidence),
        donor_gain = donor_gain,
        recipient_loss = recipient_loss,
        recipient_gain = recipient_gain,
        total_gain = total_gain,
        total_loss = total_loss,
        total_neutral = total_neutral
    )
}

#' Empirical donor/recipient/total cutoffs for call classes
#' @noRd
.call_class_thresholds <- function(calls, controls,
                                   addition_cutoff = 0.5,
                                   substitution_cutoff = -0.5,
                                   high = 1.0,
                                   min_control = 5L) {
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
    sample <- as.character(calls$sample)
    has_controls <- !is.null(controls) && length(controls) > 0L &&
        any(sample %in% as.character(controls))
    pool <- if (has_controls) {
        sample %in% as.character(controls)
    } else {
        rep(FALSE, n)
    }
    donor <- as.numeric(calls$donor_mean)
    loss <- as.numeric(calls$recipient_mean)
    total <- if ("total_mean" %in% names(calls)) {
        as.numeric(calls$total_mean)
    } else {
        rep(NA_real_, n)
    }
    defaults <- list(
        donor_mod = addition_cutoff,
        donor_high = high,
        loss_mod = substitution_cutoff,
        loss_high = -high,
        total_mod = 0.5,
        total_high = 1.0
    )
    make <- function(idx) {
        d <- donor[idx & is.finite(donor)]
        r <- loss[idx & is.finite(loss)]
        t <- abs(total[idx & is.finite(total)])
        if (length(d) < min_control || length(r) < min_control) {
            return(NULL)
        }
        donor_mod <- max(
            addition_cutoff,
            as.numeric(stats::quantile(d, 0.95, na.rm = TRUE))
        )
        donor_high <- max(
            high, donor_mod,
            as.numeric(stats::quantile(d, 0.99, na.rm = TRUE))
        )
        loss_mod <- min(
            substitution_cutoff,
            as.numeric(stats::quantile(r, 0.05, na.rm = TRUE))
        )
        loss_high <- min(
            -high, loss_mod,
            as.numeric(stats::quantile(r, 0.01, na.rm = TRUE))
        )
        total_mod <- if (length(t) >= min_control) {
            max(0.5, as.numeric(stats::quantile(t, 0.95, na.rm = TRUE)))
        } else {
            defaults$total_mod
        }
        total_high <- if (length(t) >= min_control) {
            max(1.0, total_mod,
                as.numeric(stats::quantile(t, 0.99, na.rm = TRUE)))
        } else {
            defaults$total_high
        }
        list(
            donor_mod = donor_mod,
            donor_high = donor_high,
            loss_mod = loss_mod,
            loss_high = loss_high,
            total_mod = total_mod,
            total_high = total_high,
            n = length(d)
        )
    }
    out <- list()
    if (!has_controls) {
        for (m in unique(mode)) {
            for (b in unique(bin[mode == m])) {
                out[[paste(m, b, sep = "|")]] <- defaults
            }
        }
        return(out)
    }
    for (m in unique(mode)) {
        global <- make(pool & mode == m)
        if (is.null(global)) {
            global <- make(mode == m)
        }
        if (is.null(global)) {
            global <- defaults
        }
        for (b in unique(bin[mode == m])) {
            local <- make(pool & mode == m & bin == b)
            if (is.null(local)) {
                local <- global
            } else if (!is.null(global) && local$n < 100L) {
                ## Partial pooling: a small control bin is shrunk toward the
                ## mode-wide cutoff instead of setting an extreme threshold.
                weight <- local$n / (local$n + 20)
                shrink <- function(value, reference) {
                    if (!is.finite(value)) {
                        return(reference)
                    }
                    if (!is.finite(reference)) {
                        return(value)
                    }
                    weight * value + (1 - weight) * reference
                }
                local$donor_mod <- max(
                    addition_cutoff,
                    shrink(local$donor_mod, global$donor_mod)
                )
                local$donor_high <- max(
                    high, local$donor_mod,
                    shrink(local$donor_high, global$donor_high)
                )
                local$loss_mod <- min(
                    substitution_cutoff,
                    shrink(local$loss_mod, global$loss_mod)
                )
                local$loss_high <- min(
                    -high, local$loss_mod,
                    shrink(local$loss_high, global$loss_high)
                )
                local$total_mod <- max(
                    0.5, shrink(local$total_mod, global$total_mod)
                )
                local$total_high <- max(
                    1.0, local$total_mod,
                    shrink(local$total_high, global$total_high)
                )
            }
            out[[paste(m, b, sep = "|")]] <- local
        }
    }
    out
}

#' @noRd
.call_class_status <- function(class_confidence, pass) {
    ifelse(
        !isTRUE_vec(pass), "needs_review",
        ifelse(
            class_confidence == "high", "confident",
            ifelse(class_confidence == "moderate", "ambiguous", "needs_review")
        )
    )
}

#' @noRd
.call_component_status <- function(calls, spread_cutoff = 2.0) {
    n <- nrow(calls)
    spread <- function(name) {
        if (!name %in% names(calls)) {
            return(rep(0, n))
        }
        values <- as.numeric(calls[[name]])
        values[!is.finite(values)] <- 0
        values
    }
    donor_spread <- spread("donor_component_spread")
    recipient_spread <- spread("recipient_component_spread")
    supported <- if ("n_donor_supported" %in% names(calls)) {
        values <- as.numeric(calls$n_donor_supported)
        values[!is.finite(values)] <- 0
        values
    } else {
        rep(0, n)
    }
    ifelse(
        donor_spread > spread_cutoff |
            recipient_spread > spread_cutoff |
            supported > 1,
        "heterogeneous", "consistent"
    )
}

#' @noRd
.call_mappability_status <- function(calls, cutoff = 0.5) {
    n <- nrow(calls)
    if (!"mapq_fraction" %in% names(calls)) {
        return(rep(NA_character_, n))
    }
    values <- as.numeric(calls$mapq_fraction)
    ifelse(
        is.finite(values) & values < cutoff, "low", "ok"
    )
}

#' @noRd
.call_copy_state <- function(evidence, class) {
    n <- length(class)
    out <- rep("neutral", n)
    out[evidence$donor_gain & evidence$recipient_loss] <- "balanced"
    out[evidence$donor_gain & !evidence$recipient_loss] <- "gain"
    out[!evidence$donor_gain & evidence$recipient_loss] <- "loss"
    out[evidence$total_gain | evidence$total_loss] <- "compound"
    out[class == "complex"] <- "compound"
    out[class == "other" & evidence$total_gain] <- "gain"
    out[class == "other" & evidence$total_loss] <- "loss"
    out
}

#' Classify introgression calls from donor, recipient and total evidence
#'
#' The classifier uses both sides of the paired score and, when available, a
#' relative total-depth residual. A call is a `substitution` when donor
#' coverage is gained and recipient coverage is lost with a neutral total, an
#' `addition` when donor is gained with a total gain, a `deletion` when
#' recipient is lost with a total loss, and `other` when no side is strong.
#' A donor/recipient event with a total imbalance is `complex`; so is a
#' donor/recipient width mismatch refined during assignment.
#'
#' `class_evidence` records which side(s) supported the call, `class_confidence`
#' is a class-level confidence, and `class_status` is `confident`, `ambiguous`
#' or `needs_review`. `copy_state` is a separate axis: `balanced`, `gain`,
#' `loss`, `compound` or `neutral`. `total_evidence` records
#' `gain`, `loss` or `neutral`. The older `confidence` column remains the
#' donor-effect confidence. Cutoffs are calibrated from control segments when
#' enough controls are available.
#'
#' @param calls A `Calls` data frame from [calibrateSegments()].
#' @param controls Optional control sample names used to calibrate the class
#'   cutoffs. When `NULL`, or when none of the requested controls are present,
#'   the fixed floors are used.
#' @param thresholds Optional precomputed `.call_class_thresholds()` result.
#' @param substitution_cutoff Minimum recipient-loss cutoff.
#' @param addition_cutoff Minimum donor-gain cutoff.
#' @param high,moderate Effect floors for high and moderate class confidence.
#' @param min_control Minimum number of control segments per mode/length bin
#'   before empirical cutoffs are used.
#' @param component_spread_cutoff Log-depth spread above which a group-level
#'   call is flagged as component-heterogeneous.
#'
#' @return The calls with `class`, `copy_state`, `total_evidence`,
#'   `class_evidence`, `class_confidence`, `class_status`,
#'   `component_status`, `mappability_status` and `confidence`.
#' @export
classifyCalls <- function(calls, controls = NULL, thresholds = NULL,
                          substitution_cutoff = -0.5,
                          addition_cutoff = 0.5, high = 1.0,
                          moderate = 0.5, min_control = 5L,
                          component_spread_cutoff = 2.0) {
    if (!is.data.frame(calls)) {
        .calls_error("`calls` must be a data frame")
    }
    if (nrow(calls) == 0L) {
        calls$class <- character(0)
        calls$copy_state <- character(0)
        calls$total_evidence <- character(0)
        calls$class_evidence <- character(0)
        calls$class_confidence <- character(0)
        calls$class_status <- character(0)
        calls$component_status <- character(0)
        calls$mappability_status <- character(0)
        calls$confidence <- character(0)
        return(calls)
    }
    required <- c("donor_mean", "recipient_mean", "score_mean")
    missing <- setdiff(required, names(calls))
    if (length(missing) > 0L) {
        .calls_error(
            "`calls` is missing column(s): %s", paste(missing, collapse = ", ")
        )
    }
    pass <- if ("pass" %in% names(calls)) {
        as.logical(calls$pass) %in% TRUE
    } else {
        rep(TRUE, nrow(calls))
    }
    if (is.null(thresholds)) {
        thresholds <- .call_class_thresholds(
            calls, controls,
            addition_cutoff = addition_cutoff,
            substitution_cutoff = substitution_cutoff,
            high = high,
            min_control = min_control
        )
    }
    mode <- if ("mode" %in% names(calls)) {
        as.character(calls$mode)
    } else {
        rep("all", nrow(calls))
    }
    bin <- if ("length_bin" %in% names(calls)) {
        as.character(calls$length_bin)
    } else {
        rep("all", nrow(calls))
    }
    mode[is.na(mode) | !nzchar(mode)] <- "all"
    bin[is.na(bin) | !nzchar(bin)] <- "all"
    key <- paste(mode, bin, sep = "|")
    thresholds <- thresholds[key]
    donor_mod <- vapply(thresholds, `[[`, numeric(1), "donor_mod")
    donor_high <- vapply(thresholds, `[[`, numeric(1), "donor_high")
    loss_mod <- vapply(thresholds, `[[`, numeric(1), "loss_mod")
    loss_high <- vapply(thresholds, `[[`, numeric(1), "loss_high")
    total_mod <- vapply(thresholds, `[[`, numeric(1), "total_mod")
    total_high <- vapply(thresholds, `[[`, numeric(1), "total_high")
    total_mean <- if ("total_mean" %in% names(calls)) {
        as.numeric(calls$total_mean)
    } else {
        rep(NA_real_, nrow(calls))
    }

    evidence <- .call_class_evidence(
        calls$donor_mean, calls$recipient_mean, total_mean,
        donor_cutoff = donor_mod,
        loss_cutoff = loss_mod,
        total_cutoff = total_mod
    )
    calls$class <- evidence$class
    calls$class_evidence <- evidence$evidence
    calls$copy_state <- .call_copy_state(evidence, calls$class)
    calls$total_evidence <- ifelse(
        evidence$total_gain, "gain",
        ifelse(evidence$total_loss, "loss", "neutral")
    )
    calls$total_evidence[!is.finite(total_mean)] <- NA_character_

    both <- evidence$donor_gain & evidence$recipient_loss
    donor_only <- evidence$donor_gain & !evidence$recipient_loss
    loss_only <- !evidence$donor_gain & evidence$recipient_loss
    total_known <- is.finite(total_mean) & is.finite(total_mod)
    total_imbalance <- evidence$total_gain | evidence$total_loss
    total_neutral <- evidence$total_neutral | !total_known

    ## A total imbalance turns a two-sided event into a compound event, and a
    ## donor gain with a total loss or a recipient loss with a total gain is
    ## internally inconsistent.
    compound <- (both & total_imbalance) |
        (donor_only & evidence$total_loss) |
        (loss_only & evidence$total_gain) |
        (evidence$donor_gain & evidence$recipient_gain)
    calls$class[compound] <- "complex"
    calls$class_evidence[compound] <- ifelse(
        both[compound],
        "donor+recipient; total imbalance",
        "total imbalance"
    )
    calls$copy_state[compound] <- "compound"
    calls$copy_state[evidence$donor_gain & !evidence$recipient_loss &
        !compound] <- "gain"
    calls$copy_state[!evidence$donor_gain & evidence$recipient_loss &
        !compound] <- "loss"

    high_ok <- (both & !compound & calls$donor_mean >= donor_high &
        calls$recipient_mean <= loss_mod & total_neutral) |
        (donor_only & !compound & calls$donor_mean >= donor_high &
            (total_mean >= total_high | !total_known)) |
        (loss_only & !compound & calls$recipient_mean <= loss_high &
            (total_mean <= -total_high | !total_known))
    moderate_ok <- (both & !compound & calls$donor_mean >= donor_mod &
        calls$recipient_mean <= loss_mod) |
        (donor_only & !compound & calls$donor_mean >= donor_mod) |
        (loss_only & !compound & calls$recipient_mean <= loss_mod) |
        compound
    calls$class_confidence <- ifelse(
        !is.finite(calls$donor_mean) | !is.finite(calls$recipient_mean),
        "low",
        ifelse(
            pass & high_ok, "high",
            ifelse(moderate_ok, "moderate", "low")
        )
    )
    calls$class_status <- .call_class_status(calls$class_confidence, pass)
    calls$class_status[compound] <- "needs_review"
    calls$class_status[calls$class == "other" &
        (evidence$total_gain | evidence$total_loss)] <- "ambiguous"

    calls$component_status <- .call_component_status(
        calls, spread_cutoff = component_spread_cutoff
    )
    calls$mappability_status <- .call_mappability_status(calls)
    downgrade <- (calls$component_status == "heterogeneous" |
        calls$mappability_status %in% "low") &
        calls$class_confidence == "high"
    calls$class_confidence[downgrade] <- "moderate"
    calls$class_status[downgrade & calls$class_status == "confident"] <-
        "ambiguous"

    calls$confidence <- ifelse(
        pass & calls$donor_mean >= high,
        "high",
        ifelse(calls$donor_mean >= moderate, "moderate", "low")
    )
    calls <- calls[order(calls$sample, calls$group, calls$start_iid), ,
        drop = FALSE]
    rownames(calls) <- NULL
    calls
}
