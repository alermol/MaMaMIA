## Call classification and confidence.

.calls_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIACallError", fmt, ...)
}

#' Classify introgression calls
#'
#' A call is a `substitution` when the recipient side lost coverage
#' (`recipient_mean < -0.5`), an `addition` when the donor side gained coverage
#' (`donor_mean > 0.5`), and `other` otherwise. Confidence combines the donor
#' effect with the calibration pass flag.
#'
#' @param calls A `Calls` data frame from [calibrateSegments()].
#' @param substitution_cutoff Recipient residual below which a call is a
#'   substitution.
#' @param addition_cutoff Donor residual above which a call is an addition.
#' @param high,moderate Donor-residual cutoffs for `high` and `moderate`
#'   confidence.
#'
#' @return The calls with `class` and `confidence`.
#' @export
classifyCalls <- function(calls, substitution_cutoff = -0.5,
                          addition_cutoff = 0.5, high = 1.0,
                          moderate = 0.5) {
    if (!is.data.frame(calls)) {
        .calls_error("`calls` must be a data frame")
    }
    if (nrow(calls) == 0L) {
        calls$class <- character(0)
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
    calls$class <- ifelse(
        calls$recipient_mean < substitution_cutoff,
        "substitution",
        ifelse(calls$donor_mean > addition_cutoff, "addition", "other")
    )
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
