## Window-level masks for recurrent artifacts.
##
## The only mask source is a cohort recurrence mask computed from the current
## samples (`source = "cohort"`). A single-sample run has no cohort and
## therefore returns an unmasked table with `source = "none"`.
##
## Masked windows are set to NA before CBS, not removed after segmentation, so
## segment boundaries and gap merging still see the full window grid.

.mask_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAMaskError", fmt, ...)
}

.mask_em2 <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2L) {
        return(list(
            w = c(1, 0), mu = c(stats::median(x), stats::median(x)),
            sd = c(1e-6, 1e-6)
        ))
    }
    mu <- c(0.02, 0.75)
    sd <- c(0.05, 0.15)
    w <- c(0.8, 0.2)
    for (k in seq_len(200L)) {
        d1 <- w[[1L]] * stats::dnorm(x, mu[[1L]], max(sd[[1L]], 1e-6))
        d2 <- w[[2L]] * stats::dnorm(x, mu[[2L]], max(sd[[2L]], 1e-6))
        denominator <- d1 + d2 + 1e-300
        p2 <- d2 / denominator
        p1 <- 1 - p2
        w_new <- c(mean(p1), mean(p2))
        mu_new <- c(
            sum(p1 * x) / max(sum(p1), 1e-12),
            sum(p2 * x) / max(sum(p2), 1e-12)
        )
        sd_new <- c(
            sqrt(sum(p1 * (x - mu_new[[1L]])^2) / max(sum(p1), 1e-12)),
            sqrt(sum(p2 * (x - mu_new[[2L]])^2) / max(sum(p2), 1e-12))
        )
        if (!all(is.finite(c(w_new, mu_new, sd_new)))) {
            break
        }
        converged <- max(abs(c(
            w_new - w, mu_new - mu, sd_new - sd
        ))) < 1e-8
        w <- w_new
        mu <- mu_new
        sd <- pmax(sd_new, 1e-6)
        if (converged) {
            break
        }
    }
    order <- order(mu)
    list(w = w[order], mu = mu[order], sd = sd[order])
}

.mask_posterior <- function(x, fit) {
    d1 <- fit$w[[1L]] * stats::dnorm(x, fit$mu[[1L]], max(fit$sd[[1L]], 1e-6))
    d2 <- fit$w[[2L]] * stats::dnorm(x, fit$mu[[2L]], max(fit$sd[[2L]], 1e-6))
    d2 / (d1 + d2 + 1e-300)
}

.mask_metrics_frame <- function(metrics) {
    required <- c("group", "iid", "start", "stop", "chr_id")
    missing <- setdiff(required, names(metrics))
    if (length(missing) > 0L) {
        .mask_error(
            "`metrics` is missing column(s): %s", paste(missing, collapse = ", ")
        )
    }
    ## A mask is indexed by chromosome group and relative window. Use the
    ## donor contig as the coordinate frame so every group/iid appears once.
    if ("role" %in% names(metrics) && any(metrics$role == "donor")) {
        metrics <- metrics[metrics$role == "donor", , drop = FALSE]
    }
    out <- unique(metrics[, required, drop = FALSE])
    out <- out[order(out$group, out$iid), , drop = FALSE]
    rownames(out) <- NULL
    out
}

.mask_from_fraction <- function(fraction, threshold) {
    fit <- .mask_em2(fraction$f)
    fraction$posterior <- .mask_posterior(fraction$f, fit)
    fraction$mask <- fraction$posterior > 0.5 & fraction$f > threshold
    attr(fraction, "fit") <- fit
    fraction
}

#' Build a cohort recurrence mask
#'
#' The recurrence mask is the only supported masking method. It requires at
#' least two samples in `scores`; a single-sample run returns an all-false
#' `MaskTable` with `source = "none"`. The cutoff is derived from the fitted
#' two-component mixture.
#'
#' @param metrics A `WindowMetrics` data frame.
#' @param scores Optional `ScoreTable` used for the cohort recurrence mask.
#' @param mode `"single"` or `"multi"`. Only `"multi"` can compute a mask.
#' @param z Score threshold that defines a high residual.
#'
#' @return A `MaskTable` data frame with one row per window. The derived
#'   recurrence threshold is stored in the `"recurrence_threshold"` attribute.
buildWindowMask <- function(metrics, scores = NULL,
                            mode = c("single", "multi"), z = 0.5) {
    if (!is.data.frame(metrics) || nrow(metrics) == 0L) {
        .mask_error("`metrics` must be a non-empty WindowMetrics table")
    }
    mode <- match.arg(mode, c("single", "multi"))
    frame <- .mask_metrics_frame(metrics)
    out <- data.frame(
        group = frame$group,
        iid = frame$iid,
        start = as.numeric(frame$start),
        stop = as.numeric(frame$stop),
        f = NA_real_,
        posterior = NA_real_,
        mask = FALSE,
        source = "none",
        version = NA_character_,
        stringsAsFactors = FALSE
    )
    attr(out, "recurrence_threshold") <- NA_real_
    attr(out, "threshold_source") <- "none"
    if (!identical(mode, "multi") || is.null(scores)) {
        return(out)
    }
    if (!is.data.frame(scores) ||
        !all(c("group", "iid", "score") %in% names(scores)) ||
        !"sample" %in% names(scores)) {
        .mask_error(
            "`scores` must be a ScoreTable with sample, group, iid and score"
        )
    }
    if (length(unique(as.character(scores$sample))) < 2L) {
        return(out)
    }
    keys <- unique(scores[, c("group", "iid")])
    values <- if ("score_raw" %in% names(scores)) {
        scores$score_raw
    } else {
        scores$score
    }
    key <- paste(scores$group, scores$iid, sep = "\r")
    recurrence <- tapply(values, key, function(x) {
        mean(x > z, na.rm = TRUE)
    })
    f <- as.numeric(recurrence[
        paste(keys$group, keys$iid, sep = "\r")
    ])
    fraction <- data.frame(
        group = keys$group, iid = keys$iid, f = f,
        stringsAsFactors = FALSE
    )
    n_samples <- length(unique(as.character(scores$sample)))
    threshold_floor <- if (n_samples < 10L) 0.5 + 1e-9 else -Inf
    fraction <- .mask_from_fraction(fraction, threshold_floor)
    fraction$mask[is.na(fraction$mask)] <- FALSE
    masked <- fraction$mask %in% TRUE
    derived_threshold <- if (any(masked)) {
        min(fraction$f[masked], na.rm = TRUE)
    } else {
        NA_real_
    }
    position <- match(
        paste(out$group, out$iid, sep = ":"),
        paste(fraction$group, fraction$iid, sep = ":")
    )
    if (anyNA(position)) {
        .mask_error("The cohort score and the window grid do not match")
    }
    out$f <- fraction$f[position]
    out$posterior <- fraction$posterior[position]
    out$mask <- fraction$mask[position]
    out$source <- "cohort"
    out$version <- hash_object(list(
        f = round(fraction$f, 8),
        z = z,
        threshold = if (is.finite(derived_threshold)) {
            derived_threshold
        } else {
            NA_real_
        },
        samples = sort(unique(as.character(scores$sample)))
    ))
    attr(out, "recurrence_threshold") <- derived_threshold
    attr(out, "threshold_source") <- "automatic"
    out
}

#' Apply a mask to a score profile
#'
#' @param scores A `ScoreTable`.
#' @param mask A `MaskTable`.
#'
#' @return `scores` with `score` set to `NA` in masked windows.
applyWindowMask <- function(scores, mask) {
    if (is.null(mask) || nrow(mask) == 0L ||
        !any(mask$mask %in% TRUE)) {
        return(scores)
    }
    key <- paste(scores$group, scores$iid, sep = ":")
    masked <- mask$mask %in% TRUE
    mask_key <- paste(mask$group[masked], mask$iid[masked], sep = ":")
    scores$score[key %in% mask_key] <- NA_real_
    scores
}
