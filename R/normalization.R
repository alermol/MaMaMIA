## Multi-sample dosage normalization (T4.1).
##
## `normalizeReadCounts()` converts the GC-corrected coverage of several
## samples into a dosage ratio. It runs after `correctReadCounts()` and before
## CBS (D6), writes `dose` into every RCA and leaves `cor.gc` untouched.
##
## Reference mode uses one sample as the denominator:
##   dose_s = (C_s + eps) / (C_ref + eps)
## Median mode uses the window-wise median across all samples:
##   dose_s = (C_s + eps) / (median_s(C_s) + eps)
## `by_side = TRUE` (the default) computes the denominator separately for the
## donor and recipient windows. `by_side = FALSE` pools both sides per window
## position (`iid`) and is only a sensitivity option.

.normalization_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAAnalysisError", fmt, ...)
}

.normalization_flag <- function(x, what) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        .normalization_error("`%s` must be TRUE or FALSE", what)
    }
    x
}

.normalization_pseudocount <- function(x) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
        !is.finite(x) || x <= 0) {
        .normalization_error("`pseudocount` must be a single positive number")
    }
    as.numeric(x)
}

.normalization_inputs <- function(rcas) {
    if (!is.list(rcas) || inherits(rcas, "RCA") || length(rcas) == 0L) {
        .normalization_error(
            "`rcas` must be a non-empty named list of RCA objects"
        )
    }
    sample_names <- names(rcas)
    if (is.null(sample_names) || anyNA(sample_names) ||
        any(!nzchar(sample_names))) {
        .normalization_error("`rcas` must have one non-empty name per sample")
    }
    if (anyDuplicated(sample_names)) {
        .normalization_error(
            "`rcas` has duplicate sample name(s): %s",
            paste(unique(sample_names[duplicated(sample_names)]), collapse = ", ")
        )
    }

    required <- c("chr_id", "iid", "start", "stop", "cor.gc", "ideal")
    for (i in seq_along(rcas)) {
        sample <- sample_names[[i]]
        x <- rcas[[i]]
        if (!inherits(x, "RCA")) {
            .normalization_error("Sample '%s' is not an RCA object", sample)
        }
        if (!isTRUE(x$corrected)) {
            .normalization_error(
                "Sample '%s' must be GC-corrected before normalization",
                sample
            )
        }
        if (!is.data.frame(x$data) || nrow(x$data) == 0L) {
            .normalization_error("Sample '%s' has no window data", sample)
        }
        missing <- setdiff(required, names(x$data))
        if (length(missing) > 0L) {
            .normalization_error(
                "Sample '%s' is missing columns: %s",
                sample, paste(missing, collapse = ", ")
            )
        }
        if (!is.numeric(x$data$cor.gc)) {
            .normalization_error(
                "Sample '%s': `cor.gc` must be numeric", sample
            )
        }
        if (!is.logical(x$data$ideal)) {
            .normalization_error(
                "Sample '%s': `ideal` must be logical", sample
            )
        }
        if (anyNA(x$data$chr_id) || anyNA(x$data$iid)) {
            .normalization_error(
                "Sample '%s' has missing `chr_id` or `iid` values", sample
            )
        }
        if (anyDuplicated(.normalization_keys(x$data))) {
            .normalization_error(
                "Sample '%s' has duplicated (chr_id, iid) windows", sample
            )
        }
        for (side in c("don_chr_ids", "rec_chr_ids")) {
            ids <- x$meta[[side]]
            if (is.null(ids) || length(ids) == 0L) {
                .normalization_error(
                    "Sample '%s' has no `%s`", sample, side
                )
            }
        }
    }
    invisible(sample_names)
}

.normalization_keys <- function(data) {
    paste(as.character(data$chr_id), as.character(data$iid), sep = "\r")
}

## Align every sample to the first sample's window order. All samples must
## describe exactly the same (chr_id, iid) windows with the same coordinates.
.normalization_align <- function(rcas, sample_names) {
    base <- rcas[[1L]]$data
    base_keys <- .normalization_keys(base)
    base_chr <- unique(as.character(base$chr_id))
    base_don <- as.character(rcas[[1L]]$meta$don_chr_ids)
    base_rec <- as.character(rcas[[1L]]$meta$rec_chr_ids)

    aligned <- vector("list", length(rcas))
    positions <- vector("list", length(rcas))
    names(aligned) <- sample_names
    aligned[[1L]] <- base
    positions[[1L]] <- seq_len(nrow(base))

    for (i in seq_along(rcas)[-1L]) {
        sample <- sample_names[[i]]
        data <- rcas[[i]]$data
        chr <- unique(as.character(data$chr_id))

        missing_chr <- setdiff(base_chr, chr)
        extra_chr <- setdiff(chr, base_chr)
        if (length(missing_chr) > 0L || length(extra_chr) > 0L) {
            .normalization_error(
                "Sample '%s' does not have the same chromosomes as sample '%s'%s%s",
                sample, sample_names[[1L]],
                if (length(missing_chr) > 0L) {
                    paste0("; missing: ", paste(missing_chr, collapse = ", "))
                } else {
                    ""
                },
                if (length(extra_chr) > 0L) {
                    paste0("; extra: ", paste(extra_chr, collapse = ", "))
                } else {
                    ""
                }
            )
        }

        don <- as.character(rcas[[i]]$meta$don_chr_ids)
        rec <- as.character(rcas[[i]]$meta$rec_chr_ids)
        if (!setequal(don, base_don) || !setequal(rec, base_rec)) {
            .normalization_error(
                "Sample '%s' has a different donor/recipient split than sample '%s'",
                sample, sample_names[[1L]]
            )
        }

        keys <- .normalization_keys(data)
        pos <- match(base_keys, keys)
        if (anyNA(pos)) {
            missing_window <- which(is.na(pos))[[1L]]
            .normalization_error(
                "Sample '%s' is missing window %s:%s",
                sample, base$chr_id[[missing_window]], base$iid[[missing_window]]
            )
        }
        if (length(keys) != length(base_keys)) {
            .normalization_error(
                "Sample '%s' has %d windows; sample '%s' has %d",
                sample, length(keys), sample_names[[1L]], length(base_keys)
            )
        }

        ordered <- data[pos, , drop = FALSE]
        same_coords <- as.numeric(ordered$start) == as.numeric(base$start) &
            as.numeric(ordered$stop) == as.numeric(base$stop)
        if (anyNA(same_coords) || any(!same_coords)) {
            bad <- which(is.na(same_coords) | !same_coords)[[1L]]
            .normalization_error(
                "Sample '%s' window %s:%s has different coordinates than sample '%s'",
                sample, base$chr_id[[bad]], base$iid[[bad]], sample_names[[1L]]
            )
        }

        aligned[[i]] <- ordered
        positions[[i]] <- pos
    }

    list(aligned = aligned, positions = positions)
}

.normalization_matrix <- function(aligned, index) {
    out <- matrix(
        NA_real_,
        nrow = length(index),
        ncol = length(aligned),
        dimnames = list(NULL, names(aligned))
    )
    for (j in seq_along(aligned)) {
        out[, j] <- as.numeric(aligned[[j]]$cor.gc[index])
    }
    ## Infinite corrected coverage can only come from a broken fit; treat it
    ## as missing so that the resulting dose is NA and the window is not ideal.
    out[!is.finite(out)] <- NA_real_
    out
}

.normalization_row_median <- function(values) {
    if (ncol(values) == 1L) {
        return(values[, 1L])
    }
    apply(values, 1L, stats::median, na.rm = TRUE)
}

## A joint baseline pools donor and recipient windows that share the same
## position index (`iid`). It is only used when `by_side = FALSE`.
.normalization_joint_baseline <- function(donor, recipient, iid_donor, iid_recipient) {
    combined <- rbind(donor, recipient)
    iid <- c(iid_donor, iid_recipient)
    values <- as.vector(combined)
    labels <- rep(iid, times = ncol(combined))
    baseline <- vapply(
        split(values, labels),
        stats::median,
        numeric(1L),
        na.rm = TRUE
    )
    list(
        donor = as.numeric(baseline[as.character(iid_donor)]),
        recipient = as.numeric(baseline[as.character(iid_recipient)])
    )
}

.normalization_dosage <- function(values, denominator, pseudocount) {
    dose <- (values + pseudocount) / (denominator + pseudocount)
    ## Both sides of the ratio are zero: the pseudocount formula already gives
    ## 1, but state the rule explicitly so it survives future formula changes.
    dose[!is.na(values) & values == 0 &
        !is.na(denominator) & denominator == 0] <- 1
    dose[is.na(denominator)] <- NA_real_
    dose[!is.finite(dose)] <- NA_real_
    dose
}

#' Normalize corrected coverage to dosage across samples
#'
#' @description
#' Computes a per-window dosage ratio from the GC-corrected coverage of a named
#' list of [RCA()] objects. Normalization runs after [correctReadCounts()] and
#' before [segments()] (D6): it reads `cor.gc`, writes `dose`, and leaves
#' `cor.gc` unchanged.
#'
#' With `reference`, the reference sample's corrected coverage is the
#' denominator for every other sample (`dose = (C + eps) / (C_ref + eps)`).
#' Without it, the denominator is the window-wise median across samples
#' (`dose = (C + eps) / (median(C) + eps)`). `eps` is `pseudocount`, so the
#' ratio stays finite when one side is zero; a window that is zero in every
#' sample (or zero in the reference and the sample) gets dose `1`.
#'
#' By default the donor and recipient windows are normalized independently
#' (`by_side = TRUE`), because they come from different reference genomes and
#' mappability regimes. `by_side = FALSE` is an explicit sensitivity option: it
#' pools donor and recipient windows with the same position index (`iid`) into a
#' single per-position baseline.
#'
#' Windows with a non-finite dose are marked `ideal = FALSE` so [segments()]
#' and the plots skip them. Missing (`NA`) corrected coverage stays missing in
#' the dose; infinite input coverage is treated as missing.
#'
#' @param rcas Named list of GC-corrected `RCA` objects. Every sample must
#'   describe exactly the same windows (`chr_id`, `iid`) and chromosome sets.
#' @param reference Optional sample name. When supplied, reference mode is used
#'   and `reference` must be one of `names(rcas)`.
#' @param by_side Logical; if `TRUE` (default), normalize donor and recipient
#'   windows independently. If `FALSE`, use one joint baseline per window
#'   position across both sides. Joint normalization is only available in
#'   median mode.
#' @param pseudocount Positive additive constant `eps` in the dosage ratio
#'   (default `0.5`).
#'
#' @return The input list, with `dose` added to every `RCA$data` table and a
#'   `normalization` list (`mode`, `reference`, `by_side`, `pseudocount`)
#'   stored in every `RCA`.
#' @seealso [correctReadCounts()], [segments()], [analyzeBamInput()]
#' @export
normalizeReadCounts <- function(rcas, reference = NULL, by_side = TRUE,
                                pseudocount = 0.5) {
    sample_names <- .normalization_inputs(rcas)
    by_side <- .normalization_flag(by_side, "by_side")
    pseudocount <- .normalization_pseudocount(pseudocount)

    reference_index <- NULL
    if (!is.null(reference)) {
        if (!is.character(reference) || length(reference) != 1L ||
            is.na(reference) || !nzchar(reference)) {
            .normalization_error(
                "`reference` must be NULL or a single sample name"
            )
        }
        reference_index <- match(reference, sample_names)
        if (is.na(reference_index)) {
            .normalization_error(
                "Reference sample '%s' is not one of the input samples: %s",
                reference, paste(sample_names, collapse = ", ")
            )
        }
    }
    mode <- if (is.null(reference_index)) "median" else "reference"
    if (identical(mode, "reference") && !by_side) {
        .normalization_error(
            "`by_side = FALSE` is only available for median normalization"
        )
    }

    aligned_info <- .normalization_align(rcas, sample_names)
    aligned <- aligned_info$aligned
    base <- aligned[[1L]]

    don_ids <- as.character(rcas[[1L]]$meta$don_chr_ids)
    rec_ids <- as.character(rcas[[1L]]$meta$rec_chr_ids)
    don_index <- which(base$chr_id %in% don_ids)
    rec_index <- which(base$chr_id %in% rec_ids)
    if (length(don_index) == 0L || length(rec_index) == 0L) {
        .normalization_error(
            "The RCA data must contain donor and recipient windows"
        )
    }
    if (length(don_index) + length(rec_index) != nrow(base)) {
        .normalization_error(
            "Every window must belong to exactly one donor or recipient chromosome"
        )
    }

    donor <- .normalization_matrix(aligned, don_index)
    recipient <- .normalization_matrix(aligned, rec_index)

    if (by_side) {
        donor_denominator <- if (mode == "reference") {
            donor[, reference_index]
        } else {
            .normalization_row_median(donor)
        }
        recipient_denominator <- if (mode == "reference") {
            recipient[, reference_index]
        } else {
            .normalization_row_median(recipient)
        }
    } else {
        baseline <- .normalization_joint_baseline(
            donor, recipient,
            base$iid[don_index], base$iid[rec_index]
        )
        donor_denominator <- baseline$donor
        recipient_denominator <- baseline$recipient
    }

    dose_donor <- .normalization_dosage(
        donor, donor_denominator, pseudocount
    )
    dose_recipient <- .normalization_dosage(
        recipient, recipient_denominator, pseudocount
    )

    dose <- matrix(
        NA_real_,
        nrow = nrow(base),
        ncol = length(rcas),
        dimnames = list(NULL, sample_names)
    )
    dose[don_index, ] <- dose_donor
    dose[rec_index, ] <- dose_recipient

    normalization <- list(
        mode = mode,
        reference = if (mode == "reference") reference else NULL,
        by_side = by_side,
        pseudocount = pseudocount
    )

    for (i in seq_along(rcas)) {
        out <- rcas[[i]]
        dose_original <- rep(NA_real_, nrow(out$data))
        dose_original[aligned_info$positions[[i]]] <- dose[, i]
        out$data$dose <- dose_original
        out$data$ideal <- out$data$ideal & is.finite(dose_original)
        out$normalization <- normalization
        rcas[[i]] <- out
    }

    rcas
}
