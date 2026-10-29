## Call classification and confidence.

.calls_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIACallError", fmt, ...)
}

#' @noRd
.call_class_evidence <- function(donor_mean, recipient_mean,
                                 donor_cutoff = 0.5,
                                 loss_cutoff = -0.5) {
    donor_gain <- is.finite(donor_mean) & donor_mean > donor_cutoff
    recipient_loss <- is.finite(recipient_mean) & recipient_mean < loss_cutoff
    class <- ifelse(
        donor_gain & recipient_loss, "substitution",
        ifelse(
            donor_gain, "addition",
            ifelse(recipient_loss, "deletion", "other")
        )
    )
    evidence <- ifelse(
        donor_gain & recipient_loss, "donor+recipient",
        ifelse(donor_gain, "donor",
            ifelse(recipient_loss, "recipient", "none")
        )
    )
    list(
        class = as.character(class),
        evidence = as.character(evidence),
        donor_gain = donor_gain,
        recipient_loss = recipient_loss
    )
}

#' Empirical donor/recipient cutoffs for call classes
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
    pool <- rep(TRUE, n)
    if (!is.null(controls) && length(controls) > 0L) {
        pool <- sample %in% as.character(controls)
        if (!any(pool)) {
            pool <- rep(TRUE, n)
        }
    }
    donor <- as.numeric(calls$donor_mean)
    loss <- as.numeric(calls$recipient_mean)
    defaults <- list(
        donor_mod = addition_cutoff,
        donor_high = high,
        loss_mod = substitution_cutoff,
        loss_high = -high
    )
    make <- function(idx) {
        d <- donor[idx & is.finite(donor)]
        r <- loss[idx & is.finite(loss)]
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
        list(
            donor_mod = donor_mod,
            donor_high = donor_high,
            loss_mod = loss_mod,
            loss_high = loss_high
        )
    }
    out <- list()
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

#' Classify introgression calls from donor and recipient evidence
#'
#' The classifier uses both sides of the paired score. A call is a
#' `substitution` when donor coverage is gained and recipient coverage is lost,
#' an `addition` when only the donor side is gained, a `deletion` when only the
#' recipient side is lost, and `other` when neither side has a strong effect.
#' `class_evidence` records which side(s) supported the call,
#' `class_confidence` is a class-level confidence, and `class_status` is
#' `confident`, `ambiguous` or `needs_review`. The older `confidence` column
#' remains the donor-effect confidence. Donor and recipient cutoffs are
#' calibrated from control segments when enough controls are available.
#'
#' @param calls A `Calls` data frame from [calibrateSegments()].
#' @param controls Optional control sample names used to calibrate the class
#'   cutoffs. When `NULL`, every sample contributes.
#' @param substitution_cutoff Minimum recipient-loss cutoff.
#' @param addition_cutoff Minimum donor-gain cutoff.
#' @param high,moderate Effect floors for high and moderate class confidence.
#' @param min_control Minimum number of control segments per mode/length bin
#'   before empirical cutoffs are used.
#'
#' @return The calls with `class`, `class_evidence`, `class_confidence`,
#'   `class_status` and `confidence`.
#' @export
classifyCalls <- function(calls, controls = NULL,
                          substitution_cutoff = -0.5,
                          addition_cutoff = 0.5, high = 1.0,
                          moderate = 0.5, min_control = 5L) {
    if (!is.data.frame(calls)) {
        .calls_error("`calls` must be a data frame")
    }
    if (nrow(calls) == 0L) {
        calls$class <- character(0)
        calls$class_evidence <- character(0)
        calls$class_confidence <- character(0)
        calls$class_status <- character(0)
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
    thresholds <- .call_class_thresholds(
        calls, controls,
        addition_cutoff = addition_cutoff,
        substitution_cutoff = substitution_cutoff,
        high = high,
        min_control = min_control
    )
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

    evidence <- .call_class_evidence(
        calls$donor_mean, calls$recipient_mean,
        donor_cutoff = donor_mod,
        loss_cutoff = loss_mod
    )
    calls$class <- evidence$class
    calls$class_evidence <- evidence$evidence
    both <- evidence$donor_gain & evidence$recipient_loss
    donor_only <- evidence$donor_gain & !evidence$recipient_loss
    loss_only <- !evidence$donor_gain & evidence$recipient_loss
    high_ok <- (both & calls$donor_mean >= donor_high &
        calls$recipient_mean <= loss_high) |
        (donor_only & calls$donor_mean >= donor_high) |
        (loss_only & calls$recipient_mean <= loss_high)
    moderate_ok <- (both & calls$donor_mean >= donor_mod &
        calls$recipient_mean <= loss_mod) |
        (donor_only & calls$donor_mean >= donor_mod) |
        (loss_only & calls$recipient_mean <= loss_mod)
    calls$class_confidence <- ifelse(
        !is.finite(calls$donor_mean) | !is.finite(calls$recipient_mean),
        "low",
        ifelse(
            pass & high_ok, "high",
            ifelse(moderate_ok, "moderate", "low")
        )
    )
    calls$class_status <- .call_class_status(calls$class_confidence, pass)
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
