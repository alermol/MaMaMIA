## Introgression scoring for any donor/recipient ploidy.
##
## The primary score is the paired donor/recipient log-ratio of per-subgenome
## depth:
##
##   y = log(lambda_donor + eps) - log(lambda_recipient + eps)
##
## `lambda` is breadth-derived (-log(1 - breadth)) or count-derived
## (count * read_length / length). Recipient homoeologs are averaged within a
## chromosome group; donor homoeologs are averaged, summed or max-selected.
## A single sample is GC-corrected with a natural spline and centred; a cohort
## is corrected with median polish, which removes the shared window effect and
## each sample's library/ploidy offset.

.scoring_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAScoringError", fmt, ...)
}

.scoring_eps <- 1e-4

.scoring_metric <- function(metric) {
    match.arg(metric, c("breadth", "count", "combined"))
}

.scoring_lambda <- function(metrics, metric, eps = .scoring_eps) {
    if (identical(metric, "breadth")) {
        return(as.numeric(metrics$lambda_breadth))
    }
    if (identical(metric, "count")) {
        return(as.numeric(metrics$lambda_count))
    }
    0.5 * (log(as.numeric(metrics$lambda_breadth) + eps) +
        log(as.numeric(metrics$lambda_count) + eps))
}

.validate_score_metrics <- function(metrics) {
    required <- c(
        "sample", "chr_id", "chr_name", "group", "role", "subgenome", "iid",
        "start", "stop", "length", "count", "covered", "breadth", "gc",
        "lambda_breadth", "lambda_count", "read_length"
    )
    missing <- setdiff(required, names(metrics))
    if (length(missing) > 0L) {
        .scoring_error(
            "`metrics` is missing column(s): %s", paste(missing, collapse = ", ")
        )
    }
    if (!is.data.frame(metrics) || nrow(metrics) == 0L) {
        .scoring_error("`metrics` must be a non-empty WindowMetrics data frame")
    }
    metrics$sample <- as.character(metrics$sample)
    metrics$group <- as.character(metrics$group)
    metrics$role <- as.character(metrics$role)
    metrics$subgenome <- as.character(metrics$subgenome)
    metrics$iid <- as.integer(metrics$iid)
    if (anyNA(metrics$group) || any(!nzchar(metrics$group))) {
        .scoring_error(
            "Every window must have a chromosome group; add `group` to the species contigs"
        )
    }
    if (anyNA(metrics$iid)) {
        .scoring_error("`metrics$iid` must not contain NA")
    }
    if (!all(c("donor", "recipient") %in% unique(metrics$role))) {
        .scoring_error("`metrics` must contain donor and recipient windows")
    }
    metrics
}

.scoring_composition <- function(metrics, sample_composition = NULL) {
    recipient <- unique(metrics$subgenome[metrics$role == "recipient"])
    if (is.null(sample_composition)) {
        if ("copy_number" %in% names(metrics)) {
            values <- vapply(recipient, function(s) {
                x <- metrics$copy_number[
                    metrics$role == "recipient" & metrics$subgenome == s
                ]
                if (length(x) == 0L) 2 else as.numeric(x[[1L]])
            }, numeric(1))
        } else {
            values <- stats::setNames(rep(2, length(recipient)), recipient)
        }
        return(values)
    }
    values <- unlist(sample_composition, use.names = TRUE)
    if (is.null(values) || is.null(names(values)) || anyNA(names(values)) ||
        any(!nzchar(names(values)))) {
        .scoring_error(
            "`sample_composition` must be a named numeric vector, for example c(A = 2, B = 2, D = 2)"
        )
    }
    unknown <- setdiff(names(values), recipient)
    if (length(unknown) > 0L) {
        .scoring_error(
            "`sample_composition` names subgenome(s) not present in the recipient data: %s",
            paste(unknown, collapse = ", ")
        )
    }
    values <- vapply(seq_along(values), function(i) {
        x <- values[[i]]
        if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
            !is.finite(x) || x < 0 || x != floor(x)) {
            .scoring_error(
                "`sample_composition` values must be non-negative integers"
            )
        }
        as.numeric(x)
    }, numeric(1))
    names(values) <- names(unlist(sample_composition, use.names = TRUE))
    missing <- setdiff(recipient, names(values))
    if (length(missing) > 0L) {
        values <- c(values, stats::setNames(rep(2, length(missing)), missing))
    }
    values[recipient]
}

.side_correction_single <- function(log_lambda, gc) {
    y <- as.numeric(log_lambda)
    frame <- data.frame(y = y, gc = as.numeric(gc))
    valid <- is.finite(frame$y) & is.finite(frame$gc)
    if (sum(valid) < 20L || length(unique(frame$gc[valid])) < 3L) {
        return(y - stats::median(y, na.rm = TRUE))
    }
    keep <- valid & frame$y <= stats::quantile(
        frame$y[valid], 0.90, na.rm = TRUE
    )
    fit <- tryCatch(
        stats::lm(y ~ splines::ns(gc, df = 3), data = frame[keep, , drop = FALSE]),
        error = function(e) NULL
    )
    if (is.null(fit)) {
        return(y - stats::median(y, na.rm = TRUE))
    }
    predicted <- rep(NA_real_, length(y))
    predicted[valid] <- as.numeric(stats::predict(fit, newdata = frame[valid, , drop = FALSE]))
    out <- y - predicted
    out - stats::median(out, na.rm = TRUE)
}

.side_correction_multi <- function(y, sample, group, iid) {
    key <- paste(group, iid, sep = ":")
    samples <- unique(sample)
    keys <- unique(key)
    matrix_y <- matrix(
        NA_real_, nrow = length(keys), ncol = length(samples),
        dimnames = list(keys, samples)
    )
    row <- match(key, keys)
    column <- match(sample, samples)
    matrix_y[cbind(row, column)] <- y
    polished <- medianPolish(matrix_y)
    residual <- polished$residuals[cbind(row, column)]
    as.numeric(residual)
}

.score_one_group <- function(donor, recipient, metric, donor_aggregation,
                             composition, eps = .scoring_eps,
                             score_type = "paired_log_ratio") {
    donor_lambda <- .scoring_lambda(donor, metric, eps)
    recipient_lambda <- .scoring_lambda(recipient, metric, eps)
    iid <- sort(unique(donor$iid))
    donor_matrix <- matrix(
        NA_real_, nrow = length(iid), ncol = length(unique(donor$subgenome)),
        dimnames = list(iid, sort(unique(donor$subgenome)))
    )
    donor_matrix[cbind(match(donor$iid, iid), match(donor$subgenome,
        colnames(donor_matrix)))] <- donor_lambda
    recipient_matrix <- matrix(
        NA_real_, nrow = length(iid), ncol = length(unique(recipient$subgenome)),
        dimnames = list(iid, sort(unique(recipient$subgenome)))
    )
    recipient_matrix[cbind(match(recipient$iid, iid),
        match(recipient$subgenome, colnames(recipient_matrix)))] <-
        recipient_lambda

    active <- names(composition)[composition > 0]
    active <- intersect(active, colnames(recipient_matrix))
    if (length(active) == 0L) {
        .scoring_error(
            "No recipient subgenome has a positive copy number in `sample_composition`"
        )
    }
    recipient_mean <- rowMeans(
        recipient_matrix[, active, drop = FALSE], na.rm = TRUE
    )
    recipient_mean[!is.finite(recipient_mean)] <- 0

    donor_source <- rep(NA_character_, length(iid))
    if (identical(donor_aggregation, "max")) {
        donor_mean <- rep(NA_real_, length(iid))
        for (i in seq_along(iid)) {
            values <- donor_matrix[i, , drop = TRUE]
            scores <- log(values + eps) - log(recipient_mean[[i]] + eps)
            best <- which.max(scores)
            if (length(best) == 0L || !is.finite(scores[[best]])) {
                next
            }
            donor_mean[[i]] <- values[[best]]
            donor_source[[i]] <- colnames(donor_matrix)[[best]]
        }
    } else if (identical(donor_aggregation, "sum")) {
        donor_mean <- rowSums(donor_matrix, na.rm = TRUE)
        donor_source[] <- if (ncol(donor_matrix) == 1L) {
            colnames(donor_matrix)
        } else {
            NA_character_
        }
    } else {
        donor_mean <- rowMeans(donor_matrix, na.rm = TRUE)
        donor_source[] <- if (ncol(donor_matrix) == 1L) {
            colnames(donor_matrix)
        } else {
            NA_character_
        }
    }
    donor_mean[!is.finite(donor_mean)] <- 0
    score <- if (identical(score_type, "donor_only")) {
        log(donor_mean + eps)
    } else {
        log(donor_mean + eps) - log(recipient_mean + eps)
    }

    gc_donor <- vapply(iid, function(i) {
        mean(donor$gc[donor$iid == i], na.rm = TRUE)
    }, numeric(1))
    gc_recipient <- vapply(iid, function(i) {
        mean(recipient$gc[recipient$iid == i], na.rm = TRUE)
    }, numeric(1))

    group_length <- max(recipient$stop, na.rm = TRUE)
    n_per <- length(iid)
    list(
        iid = iid,
        start = floor(iid * group_length / n_per),
        stop = floor((iid + 1L) * group_length / n_per),
        score = score,
        donor_mean = donor_mean,
        recipient_mean = recipient_mean,
        donor_subgenome = donor_source,
        n_donor = rep(ncol(donor_matrix), length(iid)),
        n_recipient = rep(length(active), length(iid)),
        gc_donor = gc_donor,
        gc_recipient = gc_recipient,
        log_donor = log(donor_mean + eps),
        log_recipient = log(recipient_mean + eps),
        length = group_length,
        chr_name = recipient$chr_name[[1L]],
        n_per = n_per
    )
}

#' Score introgression windows
#'
#' Computes the paired donor/recipient log-ratio (or the donor-only fallback)
#' for every sample, chromosome group and window. Supports any number of donor
#' and recipient subgenomes. Single samples are GC-corrected with a natural
#' spline; cohorts use median polish. The optional copy-number model reports
#' donor copies and donor fraction.
#'
#' @param metrics A `WindowMetrics` data frame from [buildWindowMetrics()].
#' @param mode `"single"` or `"multi"`.
#' @param metric `"breadth"`, `"count"` or `"combined"`.
#' @param donor_aggregation `"mean"`, `"max"` or `"sum"`.
#' @param recipient_aggregation Only `"mean"` is supported.
#' @param sample_composition Optional named copy-number vector for recipient
#'   subgenomes, for example `c(A = 2, B = 2, D = 2)`.
#' @param gc_correction Correct the single-sample score for GC.
#' @param score_type `"paired_log_ratio"` or `"donor_only"`.
#'
#' @return A `ScoreTable` data frame.
#' @export
scoreIntrogression <- function(
    metrics,
    mode = c("single", "multi"),
    metric = c("breadth", "count", "combined"),
    donor_aggregation = c("mean", "max", "sum"),
    recipient_aggregation = "mean",
    sample_composition = NULL,
    gc_correction = TRUE,
    score_type = c("paired_log_ratio", "donor_only")
) {
    metrics <- .validate_score_metrics(metrics)
    mode <- match.arg(mode, c("single", "multi"))
    metric <- .scoring_metric(metric)
    donor_aggregation <- match.arg(donor_aggregation, c("mean", "max", "sum"))
    score_type <- match.arg(score_type, c("paired_log_ratio", "donor_only"))
    if (!identical(recipient_aggregation, "mean")) {
        .scoring_error("Only `recipient_aggregation = 'mean'` is supported")
    }
    if (length(gc_correction) != 1L || is.na(gc_correction) ||
        !is.logical(gc_correction)) {
        .scoring_error("`gc_correction` must be true or false")
    }
    samples <- unique(metrics$sample)
    if (identical(mode, "multi") && length(samples) == 1L) {
        mode <- "single"
    }
    composition <- .scoring_composition(metrics, sample_composition)
    eps <- .scoring_eps

    pieces <- list()
    for (sample in samples) {
        x <- metrics[metrics$sample == sample, , drop = FALSE]
        for (group in sort(unique(x$group))) {
            donor <- x[x$role == "donor" & x$group == group, , drop = FALSE]
            recipient <- x[
                x$role == "recipient" & x$group == group, ,
                drop = FALSE
            ]
            if (nrow(donor) == 0L || nrow(recipient) == 0L) {
                next
            }
            donor <- donor[order(donor$subgenome, donor$iid), , drop = FALSE]
            recipient <- recipient[
                order(recipient$subgenome, recipient$iid), ,
                drop = FALSE
            ]
            scored <- .score_one_group(
                donor, recipient, metric, donor_aggregation, composition, eps,
                score_type
            )
            pieces[[length(pieces) + 1L]] <- data.frame(
                sample = sample,
                group = group,
                iid = scored$iid,
                start = scored$start,
                stop = scored$stop,
                mode = mode,
                metric = metric,
                score_type = score_type,
                score = scored$score,
                donor_mean = scored$donor_mean,
                recipient_mean = scored$recipient_mean,
                donor_subgenome = scored$donor_subgenome,
                n_donor = scored$n_donor,
                n_recipient = scored$n_recipient,
                chr_name = scored$chr_name,
                length = scored$length,
                n_per_chromosome = scored$n_per,
                gc_donor = scored$gc_donor,
                gc_recipient = scored$gc_recipient,
                log_donor = scored$log_donor,
                log_recipient = scored$log_recipient,
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(pieces) == 0L) {
        .scoring_error("No donor/recipient group could be scored")
    }
    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    out$score_raw <- out$score

    if (identical(mode, "single")) {
        score <- numeric(nrow(out))
        donor_residual <- numeric(nrow(out))
        recipient_residual <- numeric(nrow(out))
        for (sample in unique(out$sample)) {
            idx <- which(out$sample == sample)
            score[idx] <- if (isTRUE(gc_correction)) {
                if (identical(score_type, "donor_only")) {
                    .side_correction_single(
                        out$score[idx], out$gc_donor[idx]
                    )
                } else {
                    .score_correct_single_pair(
                        out$score[idx], out$gc_donor[idx],
                        out$gc_recipient[idx]
                    )
                }
            } else {
                out$score[idx] - stats::median(out$score[idx], na.rm = TRUE)
            }
            donor_residual[idx] <- .side_correction_single(
                out$log_donor[idx], out$gc_donor[idx]
            )
            recipient_residual[idx] <- .side_correction_single(
                out$log_recipient[idx], out$gc_recipient[idx]
            )
        }
        out$score <- score
        out$donor_residual <- donor_residual
        out$recipient_residual <- recipient_residual
    } else {
        out$score <- .side_correction_multi(
            out$score, out$sample, out$group, out$iid
        )
        out$donor_residual <- .side_correction_multi(
            out$log_donor, out$sample, out$group, out$iid
        )
        out$recipient_residual <- .side_correction_multi(
            out$log_recipient, out$sample, out$group, out$iid
        )
    }

    kappa_total <- sum(composition[composition > 0])
    active <- names(composition)[composition > 0]
    out$lambda_copy <- out$recipient_mean * length(active) / kappa_total
    out$donor_copies <- ifelse(
        is.finite(out$lambda_copy) & out$lambda_copy > 0,
        out$donor_mean / out$lambda_copy,
        NA_real_
    )
    out$donor_fraction <- ifelse(
        is.finite(out$donor_copies),
        out$donor_copies / (out$donor_copies + kappa_total),
        NA_real_
    )
    out <- out[order(
        out$sample, out$group, out$iid
    ), , drop = FALSE]
    rownames(out) <- NULL
    out[, c(
        "sample", "group", "iid", "start", "stop", "mode", "metric",
        "score_type", "score", "donor_mean",
        "recipient_mean", "donor_subgenome", "n_donor", "n_recipient",
        "donor_copies", "donor_fraction", "chr_name", "length",
        "n_per_chromosome", "gc_donor", "gc_recipient", "log_donor",
        "log_recipient", "donor_residual", "recipient_residual", "lambda_copy",
        "score_raw"
    )]
}

.score_correct_single_pair <- function(score, gc_donor, gc_recipient) {
    y <- as.numeric(score)
    frame <- data.frame(
        y = y,
        gc_donor = as.numeric(gc_donor),
        gc_recipient = as.numeric(gc_recipient)
    )
    valid <- is.finite(frame$y) & is.finite(frame$gc_donor) &
        is.finite(frame$gc_recipient)
    if (sum(valid) < 20L || length(unique(frame$gc_donor[valid])) < 3L ||
        length(unique(frame$gc_recipient[valid])) < 3L) {
        return(y - stats::median(y, na.rm = TRUE))
    }
    keep <- valid & frame$y <= stats::quantile(
        frame$y[valid], 0.90, na.rm = TRUE
    )
    fit <- tryCatch(
        stats::lm(
            y ~ splines::ns(gc_donor, df = 3) +
                splines::ns(gc_recipient, df = 3),
            data = frame[keep, , drop = FALSE]
        ),
        error = function(e) NULL
    )
    if (is.null(fit)) {
        return(y - stats::median(y, na.rm = TRUE))
    }
    predicted <- rep(NA_real_, length(y))
    predicted[valid] <- as.numeric(stats::predict(
        fit, newdata = frame[valid, , drop = FALSE]
    ))
    out <- y - predicted
    out - stats::median(out, na.rm = TRUE)
}
