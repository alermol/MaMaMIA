## Window-level masks for recurrent artifacts.
##
## Three sources are supported, in this order:
##   1. a fixed panel mask built from a reference cohort (`source = panel`);
##   2. a cohort recurrence mask computed from the current samples
##      (`source = cohort`);
##   3. an external BED blacklist (`source = blacklist`).
## When none of those applies, a physical terminal fallback masks a fixed
## number of megabases at both contig ends (`source = physical`).
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

#' Build a panel recurrence mask
#'
#' Computes the single-sample score for every panel sample, the fraction of
#' samples above `z` in each window, and fits a two-component mixture to those
#' fractions. Windows in the recurrent (upper) component with a fraction above
#' `recurrence_threshold` are masked.
#'
#' @param panel_metrics A `WindowMetrics` data frame for the panel.
#' @param z Score threshold that defines a high residual in one sample.
#' @param recurrence_threshold Minimum fraction of panel samples.
#' @param output Optional TSV path for the mask.
#' @param metric,donor_aggregation Scoring options passed to
#'   [scoreIntrogression()].
#'
#' @return A `MaskTable` data frame.
#' @export
buildPanelMask <- function(panel_metrics, z = 0.5,
                           recurrence_threshold = 0.4, output = NULL,
                           metric = "breadth",
                           donor_aggregation = "mean") {
    if (!is.data.frame(panel_metrics) || nrow(panel_metrics) == 0L) {
        .mask_error("`panel_metrics` must be a non-empty WindowMetrics table")
    }
    n_panel <- length(unique(panel_metrics$sample))
    if (n_panel < 10L) {
        warning(
            sprintf(
                "Panel has only %d sample(s); using a fixed fraction > 0.5 instead of the two-component posterior",
                n_panel
            ),
            call. = FALSE
        )
    }
    scored <- scoreIntrogression(
        panel_metrics,
        mode = "single",
        metric = metric,
        donor_aggregation = donor_aggregation
    )
    keys <- unique(scored[, c("group", "iid")])
    f <- vapply(seq_len(nrow(keys)), function(i) {
        values <- scored$score[
            scored$group == keys$group[[i]] & scored$iid == keys$iid[[i]]
        ]
        mean(values > z, na.rm = TRUE)
    }, numeric(1))
    fraction <- data.frame(
        group = keys$group,
        iid = keys$iid,
        f = f,
        stringsAsFactors = FALSE
    )
    if (n_panel < 10L) {
        fraction$posterior <- NA_real_
        fraction$mask <- fraction$f > 0.5
    } else {
        fraction <- .mask_from_fraction(fraction, recurrence_threshold)
    }
    frame <- .mask_metrics_frame(panel_metrics)
    position <- match(
        paste(fraction$group, fraction$iid, sep = ":"),
        paste(frame$group, frame$iid, sep = ":")
    )
    if (anyNA(position)) {
        .mask_error("The panel score and the window grid do not match")
    }
    out <- data.frame(
        group = frame$group,
        iid = frame$iid,
        start = as.numeric(frame$start),
        stop = as.numeric(frame$stop),
        f = NA_real_,
        posterior = NA_real_,
        mask = FALSE,
        source = "panel",
        version = hash_object(list(
            f = round(fraction$f, 8),
            z = z,
            threshold = recurrence_threshold,
            n = n_panel
        )),
        stringsAsFactors = FALSE
    )
    out$f[position] <- fraction$f
    out$posterior[position] <- fraction$posterior
    out$mask[position] <- fraction$mask
    if (!is.null(output)) {
        utils::write.table(
            out, output, sep = "\t", quote = FALSE, row.names = FALSE
        )
    }
    out
}

.read_mask_table <- function(panel_mask) {
    if (is.character(panel_mask)) {
        if (length(panel_mask) != 1L || !file.exists(panel_mask)) {
            .mask_error("Panel mask file does not exist: %s", panel_mask)
        }
        panel_mask <- utils::read.delim(
            panel_mask,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }
    if (!is.data.frame(panel_mask) ||
        !all(c("group", "iid", "mask") %in% names(panel_mask))) {
        .mask_error(
            "`panel_mask` must be a data frame or TSV with group, iid and mask"
        )
    }
    panel_mask$group <- as.character(panel_mask$group)
    panel_mask$iid <- as.integer(panel_mask$iid)
    panel_mask$mask <- as.logical(panel_mask$mask)
    panel_mask
}

.read_blacklist <- function(blacklist) {
    if (is.null(blacklist)) {
        return(NULL)
    }
    if (is.character(blacklist)) {
        if (length(blacklist) != 1L || !file.exists(blacklist)) {
            .mask_error("Blacklist file does not exist: %s", blacklist)
        }
        blacklist <- utils::read.delim(
            blacklist,
            header = FALSE,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        if (ncol(blacklist) < 3L) {
            .mask_error(
                "Blacklist BED must have at least three columns (chr, start, end)"
            )
        }
        names(blacklist)[seq_len(3L)] <- c("chr_id", "start", "stop")
    }
    if (!is.data.frame(blacklist) ||
        !all(c("chr_id", "start", "stop") %in% names(blacklist))) {
        .mask_error(
            "`blacklist` must be a BED path or a data frame with chr_id, start and stop"
        )
    }
    blacklist$chr_id <- as.character(blacklist$chr_id)
    blacklist$start <- as.numeric(blacklist$start)
    blacklist$stop <- as.numeric(blacklist$stop)
    blacklist
}

#' Build a window mask
#'
#' @param metrics A `WindowMetrics` data frame.
#' @param scores Optional `ScoreTable` used for the cohort recurrence mask.
#' @param mode `"single"` or `"multi"`.
#' @param panel_mask Optional panel `MaskTable` or TSV path.
#' @param blacklist Optional BED path or data frame.
#' @param physical_fallback_mb Width in Mb masked at both contig ends when no
#'   other source applies. `0` disables the fallback.
#' @param recurrence_threshold Minimum recurrence fraction for a cohort mask.
#' @param z Score threshold that defines a high residual.
#' @param config Optional validated configuration.
#'
#' @return A `MaskTable` data frame with one row per window.
#' @export
buildWindowMask <- function(metrics, scores = NULL,
                            mode = c("single", "multi"),
                            panel_mask = NULL, blacklist = NULL,
                            physical_fallback_mb = 1,
                            recurrence_threshold = 0.4, z = 0.5,
                            config = NULL) {
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
    used <- FALSE

    if (!is.null(panel_mask)) {
        panel <- .read_mask_table(panel_mask)
        position <- match(
            paste(out$group, out$iid, sep = ":"),
            paste(panel$group, panel$iid, sep = ":")
        )
        if (anyNA(position)) {
            .mask_error(
                "The panel mask does not cover every analysis window; rebuild it with the same window grid"
            )
        }
        out$f <- panel$f[position]
        out$posterior <- panel$posterior[position]
        out$mask <- panel$mask[position]
        out$source <- "panel"
        out$version <- panel$version[position] %||% NA_character_
        used <- TRUE
    } else if (identical(mode, "multi") && !is.null(scores)) {
        if (!is.data.frame(scores) ||
            !all(c("group", "iid", "score") %in% names(scores))) {
            .mask_error("`scores` must be a ScoreTable with group, iid and score")
        }
        keys <- unique(scores[, c("group", "iid")])
        f <- vapply(seq_len(nrow(keys)), function(i) {
            values <- if ("score_raw" %in% names(scores)) {
                scores$score_raw
            } else {
                scores$score
            }
            values <- values[
                scores$group == keys$group[[i]] &
                    scores$iid == keys$iid[[i]]
            ]
            mean(values > z, na.rm = TRUE)
        }, numeric(1))
        fraction <- data.frame(
            group = keys$group, iid = keys$iid, f = f,
            stringsAsFactors = FALSE
        )
        n_samples <- length(unique(scores$sample))
        effective_threshold <- recurrence_threshold
        if (n_samples < 10L) {
            ## In a tiny cohort a unique event has f = 1/n; do not mask it. A
            ## shared artifact in every sample still has f = 1.
            effective_threshold <- max(effective_threshold, 0.5 + 1e-9)
        }
        fraction <- .mask_from_fraction(fraction, effective_threshold)
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
            threshold = effective_threshold,
            samples = sort(unique(scores$sample))
        ))
        used <- TRUE
    }

    if (!is.null(blacklist)) {
        bed <- .read_blacklist(blacklist)
        ids <- as.character(frame$chr_id)
        starts <- as.numeric(frame$start)
        stops <- as.numeric(frame$stop)
        mask <- rep(FALSE, length(ids))
        bed_key <- paste(bed$chr_id, bed$start, bed$stop, sep = ":")
        for (i in seq_len(nrow(bed))) {
            mask <- mask | (
                ids == bed$chr_id[[i]] &
                    stops > bed$start[[i]] &
                    starts < bed$stop[[i]]
            )
        }
        if (!used) {
            out$source <- "blacklist"
        }
        out$mask <- out$mask | mask
        out$version <- ifelse(
            is.na(out$version), hash_object(bed_key), out$version
        )
        used <- TRUE
    }

    if (!used && is.numeric(physical_fallback_mb) &&
        length(physical_fallback_mb) == 1L &&
        is.finite(physical_fallback_mb) && physical_fallback_mb > 0) {
        fallback_bp <- physical_fallback_mb * 1e6
        lengths <- tapply(
            as.numeric(frame$stop), as.character(frame$group), max,
            na.rm = TRUE
        )
        group_length <- as.numeric(lengths[out$group])
        ## A fixed megabase fallback must not swallow a short contig (for
        ## example the 20 kb test fixture); cap it at 10% of the contig.
        end_mask <- pmin(fallback_bp, 0.1 * group_length)
        out$mask <- out$start < end_mask | out$stop > (group_length - end_mask)
        out$source <- "physical"
        out$version <- hash_object(list(
            fallback_mb = physical_fallback_mb,
            lengths = lengths
        ))
    }
    out
}

#' Apply a mask to a score profile
#'
#' @param scores A `ScoreTable`.
#' @param mask A `MaskTable`.
#'
#' @return `scores` with `score` set to `NA` in masked windows.
#' @export
applyWindowMask <- function(scores, mask) {
    if (is.null(mask) || nrow(mask) == 0L || !any(mask$mask)) {
        return(scores)
    }
    key <- paste(scores$group, scores$iid, sep = ":")
    mask_key <- paste(mask$group[mask$mask], mask$iid[mask$mask], sep = ":")
    scores$score[key %in% mask_key] <- NA_real_
    scores
}
