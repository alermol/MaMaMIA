## CBS segmentation of the score profile, with gap merging.

.segmentation_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIASegmentationError", fmt, ...)
}

.segment_one <- function(x, alpha, undo_SD, min_width) {
    x <- x[order(x$iid), , drop = FALSE]
    y <- as.numeric(x$score)
    keep <- is.finite(y)
    if (sum(keep) < max(20L, min_width)) {
        return(NULL)
    }
    cna <- DNAcopy::CNA(
        genomdat = y[keep],
        chrom = rep("1", sum(keep)),
        maploc = as.integer(x$iid[keep]),
        data.type = "logratio"
    )
    smoothed <- DNAcopy::smooth.CNA(cna)
    segmented <- tryCatch(
        DNAcopy::segment(
            smoothed,
            verbose = 0,
            undo.splits = "sdundo",
            alpha = alpha,
            min.width = min_width,
            undo.SD = undo_SD
        ),
        error = function(e) NULL
    )
    if (is.null(segmented) || is.null(segmented$output) ||
        nrow(segmented$output) == 0L) {
        return(NULL)
    }
    output <- segmented$output
    data.frame(
        start_iid = as.integer(output$loc.start),
        end_iid = as.integer(output$loc.end),
        n_windows = as.integer(output$num.mark),
        score_mean = as.numeric(output$seg.mean),
        stringsAsFactors = FALSE
    )
}

.segment_merge <- function(calls, profile, merge_gap, merge_floor,
                           min_mb, merge_min = 0.25) {
    if (nrow(calls) == 0L) {
        return(calls)
    }
    calls <- calls[order(calls$start_iid), , drop = FALSE]
    positive <- calls$score_mean > merge_min
    merged <- list()
    current <- NULL
    for (i in seq_len(nrow(calls))) {
        row <- calls[i, , drop = FALSE]
        if (!positive[[i]]) {
            if (!is.null(current)) {
                merged[[length(merged) + 1L]] <- current
                current <- NULL
            }
            merged[[length(merged) + 1L]] <- row
            next
        }
        if (is.null(current)) {
            current <- row
            next
        }
        gap_start <- current$end_iid + 1L
        gap_end <- row$start_iid - 1L
        gap_mean <- if (gap_end >= gap_start) {
            mean(profile[as.character(gap_start:gap_end)], na.rm = TRUE)
        } else {
            0
        }
        if (!is.finite(gap_mean)) {
            gap_mean <- 0
        }
        if (row$start_iid - current$end_iid - 1L <= merge_gap &&
            gap_mean > merge_floor) {
            weight_current <- current$n_windows
            weight_row <- row$n_windows
            current$score_mean <- (
                current$score_mean * weight_current +
                    row$score_mean * weight_row
            ) / (weight_current + weight_row)
            current$end_iid <- row$end_iid
            current$n_windows <- current$end_iid - current$start_iid + 1L
        } else {
            merged[[length(merged) + 1L]] <- current
            current <- row
        }
    }
    if (!is.null(current)) {
        merged[[length(merged) + 1L]] <- current
    }
    out <- do.call(rbind, merged)
    rownames(out) <- NULL
    out[out$end_iid >= out$start_iid, , drop = FALSE]
}

#' Segment score profiles with CBS
#'
#' Runs circular binary segmentation on each sample × chromosome-group score
#' profile. Masked windows are set to `NA` and excluded from the CBS fit.
#' Adjacent positive segments separated by a short, near-null gap are merged.
#'
#' @param scores A `ScoreTable` from [scoreIntrogression()].
#' @param mask Optional `MaskTable` from [buildWindowMask()].
#' @param alpha CBS significance level.
#' @param undo_SD `undo.splits = "sdundo"` parameter.
#' @param min_width Minimum CBS segment width in windows.
#' @param merge_gap Maximum gap in windows that can be merged.
#' @param merge_floor Minimum mean score in the gap to allow merging.
#' @param merge_min Minimum segment mean for a segment to be merged; null
#'   segments below this value stay separate and are removed by calibration.
#' @param min_mb Minimum call width in Mb.
#' @param lengths Optional named vector of contig lengths by group.
#'
#' @return A `Calls` data frame without calibration thresholds.
#' @export
segmentScores <- function(scores, mask = NULL, alpha = 0.001, undo_SD = 1,
                          min_width = 5L, merge_gap = 100L,
                          merge_floor = -0.5, min_mb = 5,
                          merge_min = 0.25, lengths = NULL) {
    if (!is.data.frame(scores) || nrow(scores) == 0L) {
        .segmentation_error("`scores` must be a non-empty ScoreTable")
    }
    required <- c(
        "sample", "group", "iid", "score", "donor_residual",
        "recipient_residual", "length", "n_per_chromosome", "mode", "metric"
    )
    missing <- setdiff(required, names(scores))
    if (length(missing) > 0L) {
        .segmentation_error(
            "`scores` is missing column(s): %s", paste(missing, collapse = ", ")
        )
    }
    alpha <- as.numeric(alpha)[[1L]]
    undo_SD <- as.numeric(undo_SD)[[1L]]
    min_width <- as.integer(min_width)[[1L]]
    merge_gap <- as.integer(merge_gap)[[1L]]
    if (anyNA(c(alpha, undo_SD, min_width, merge_gap)) ||
        min_width < 1L || merge_gap < 0L) {
        .segmentation_error("Invalid segmentation parameters")
    }
    if (!is.null(mask)) {
        scores <- applyWindowMask(scores, mask)
    }
    if (is.null(lengths)) {
        lengths <- tapply(scores$length, scores$group, max, na.rm = TRUE)
    }
    if (is.null(names(lengths))) {
        .segmentation_error("`lengths` must be named by chromosome group")
    }

    pieces <- list()
    for (sample in unique(scores$sample)) {
        for (group in sort(unique(scores$group))) {
            x <- scores[
                scores$sample == sample & scores$group == group, ,
                drop = FALSE
            ]
            if (nrow(x) < max(20L, min_width)) {
                next
            }
            segments <- .segment_one(x, alpha, undo_SD, min_width)
            if (is.null(segments)) {
                next
            }
            profile <- stats::setNames(x$score, x$iid)
            segments <- .segment_merge(
                segments, profile, merge_gap, merge_floor, min_mb, merge_min
            )
            if (nrow(segments) == 0L) {
                next
            }
            group_length <- as.numeric(lengths[[group]])
            n_per <- as.integer(x$n_per_chromosome[[1L]])
            if (!is.finite(group_length) || group_length <= 0 || n_per < 1L) {
                .segmentation_error(
                    "No valid length for chromosome group '%s'", group
                )
            }
            segments$sample <- sample
            segments$group <- group
            segments$mode <- x$mode[[1L]]
            segments$metric <- x$metric[[1L]]
            segments$method <- if ("score_type" %in% names(x)) {
                if (x$score_type[[1L]] %in% c("donor_only", "donor-only")) {
                    "donor_only"
                } else {
                    "paired_lr"
                }
            } else {
                "paired_lr"
            }
            segments$start_bp <- floor(segments$start_iid * group_length / n_per)
            segments$end_bp <- floor(
                (segments$end_iid + 1L) * group_length / n_per
            )
            segments$width_bp <- segments$end_bp - segments$start_bp
            segments$donor_mean <- vapply(seq_len(nrow(segments)), function(i) {
                idx <- x$iid >= segments$start_iid[[i]] &
                    x$iid <= segments$end_iid[[i]]
                mean(x$donor_residual[idx], na.rm = TRUE)
            }, numeric(1))
            segments$recipient_mean <- vapply(
                seq_len(nrow(segments)), function(i) {
                    idx <- x$iid >= segments$start_iid[[i]] &
                        x$iid <= segments$end_iid[[i]]
                    mean(x$recipient_residual[idx], na.rm = TRUE)
                }, numeric(1)
            )
            if ("total_residual" %in% names(x)) {
                segments$total_mean <- vapply(
                    seq_len(nrow(segments)), function(i) {
                        idx <- x$iid >= segments$start_iid[[i]] &
                            x$iid <= segments$end_iid[[i]]
                        mean(x$total_residual[idx], na.rm = TRUE)
                    }, numeric(1)
                )
            } else {
                segments$total_mean <- NA_real_
            }
            for (column in c(
                "donor_component_spread", "recipient_component_spread",
                "n_donor_supported", "mapq_fraction"
            )) {
                if (column %in% names(x)) {
                    segments[[column]] <- vapply(
                        seq_len(nrow(segments)), function(i) {
                            idx <- x$iid >= segments$start_iid[[i]] &
                                x$iid <= segments$end_iid[[i]]
                            mean(x[[column]][idx], na.rm = TRUE)
                        }, numeric(1)
                    )
                } else {
                    segments[[column]] <- NA_real_
                }
            }
            pieces[[length(pieces) + 1L]] <- segments
        }
    }
    if (length(pieces) == 0L) {
        return(data.frame(
            sample = character(0), group = character(0),
            start_iid = integer(0), end_iid = integer(0),
            start_bp = numeric(0), end_bp = numeric(0),
            width_bp = numeric(0), score_mean = numeric(0),
            donor_mean = numeric(0), recipient_mean = numeric(0),
            total_mean = numeric(0), donor_component_spread = numeric(0),
            recipient_component_spread = numeric(0),
            n_donor_supported = numeric(0), mapq_fraction = numeric(0),
            n_windows = integer(0), method = character(0),
            mode = character(0), metric = character(0),
            stringsAsFactors = FALSE
        ))
    }
    out <- do.call(rbind, pieces)
    out <- out[out$width_bp >= min_mb * 1e6, , drop = FALSE]
    out <- out[, c(
        "sample", "group", "start_iid", "end_iid", "start_bp", "end_bp",
        "width_bp", "score_mean", "donor_mean", "recipient_mean",
        "total_mean", "donor_component_spread",
        "recipient_component_spread", "n_donor_supported",
        "mapq_fraction", "n_windows", "method", "mode", "metric"
    )]
    rownames(out) <- NULL
    out
}
