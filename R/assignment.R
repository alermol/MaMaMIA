## Assign a group-level donor segment to its most probable recipient position.
##
## The detection and calibration stay group-level: a call says that a
## chromosome group carries donor-like content. This layer asks which
## recipient chromosome in that group is the most likely host. It compares the
## per-recipient-subgenome depth over the call interval after removing each
## sample's GC trend and median depth; the most depleted recipient chromosome
## is the best candidate for a substitution or deletion. A copy-neutral event
## has no depth loss, so the best candidate is returned but flagged as
## unresolved rather than presented as a confident position.

.assignment_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAAssignmentError", fmt, ...)
}

#' Normalise recipient depth for positional assignment
#'
#' The group-level score uses the recipient mean. For assignment we need the
#' individual recipient contigs, so each sample's recipient log-depth is
#' GC-corrected and centred with the same single-sample smoother used by the
#' score. Centring the whole sample preserves differences between A, B and D.
#' @noRd
.assignment_recipient_depth <- function(metrics, eps = 1e-6) {
    out <- metrics
    out$log_depth <- log(pmax(as.numeric(out$lambda_breadth), eps))
    out$normalised_depth <- NA_real_
    for (sample in unique(out$sample)) {
        rows <- out$sample == sample
        out$normalised_depth[rows] <- .side_correction_single(
            out$log_depth[rows], out$gc[rows]
        )
    }
    out
}

#' Pick the donor subgenome recorded by the score
#' @noRd
.assignment_donor_subgenome <- function(scores, sample, group, start_iid,
                                         end_iid) {
    rows <- scores[
        as.character(scores$sample) == sample &
            as.character(scores$group) == group &
            scores$iid >= start_iid & scores$iid <= end_iid,
        ,
        drop = FALSE
    ]
    donors <- as.character(rows$donor_subgenome)
    donors <- donors[!is.na(donors) & nzchar(donors)]
    if (length(donors) == 0L) {
        return(NA_character_)
    }
    counts <- sort(table(donors), decreasing = TRUE)
    names(counts)[[1L]]
}

#' Assign calls to recipient contigs
#'
#' Adds `recipient_subgenome`, `recipient_chr_id`, `recipient_chr_name`,
#' `recipient_start_bp`, `recipient_end_bp`, `recipient_weight`,
#' `recipient_margin`, `recipient_assignment`, `donor_subgenome`,
#' `donor_chr_id`, `donor_chr_name`, `donor_start_bp` and `donor_end_bp` to
#' the Calls table.
#'
#' @param calls A Calls data frame.
#' @param metrics A WindowMetrics data frame.
#' @param windows An AnalysisWindows list.
#' @param scores A ScoreTable data frame.
#' @param assembly The resolved assembly list with `contigs`.
#' @param eps Depth offset.
#'
#' @return The Calls table with assignment columns.
#' @noRd
assignCallPositions <- function(calls, metrics, windows, scores, assembly,
                                eps = 1e-6) {
    empty <- function() {
        calls$recipient_subgenome <- character(nrow(calls))
        calls$recipient_chr_id <- character(nrow(calls))
        calls$recipient_chr_name <- character(nrow(calls))
        calls$recipient_start_bp <- numeric(nrow(calls))
        calls$recipient_end_bp <- numeric(nrow(calls))
        calls$recipient_weight <- numeric(nrow(calls))
        calls$recipient_margin <- numeric(nrow(calls))
        calls$recipient_assignment <- character(nrow(calls))
        calls$donor_subgenome <- character(nrow(calls))
        calls$donor_chr_id <- character(nrow(calls))
        calls$donor_chr_name <- character(nrow(calls))
        calls$donor_start_bp <- numeric(nrow(calls))
        calls$donor_end_bp <- numeric(nrow(calls))
        calls
    }
    if (!is.data.frame(calls)) {
        .assignment_error("`calls` must be a data frame")
    }
    if (nrow(calls) == 0L) {
        return(empty())
    }
    calls <- empty()
    if (!is.data.frame(metrics) ||
        !all(c("sample", "chr_id", "chr_name", "group", "role",
               "subgenome", "iid", "lambda_breadth", "gc") %in%
            names(metrics))) {
        .assignment_error("`metrics` must be a WindowMetrics table")
    }
    if (!is.list(windows) || is.null(windows$table)) {
        .assignment_error("`windows` must be an AnalysisWindows list")
    }
    if (!is.data.frame(scores) ||
        !all(c("sample", "group", "iid", "donor_subgenome") %in%
            names(scores))) {
        .assignment_error("`scores` must be a ScoreTable with donor_subgenome")
    }
    recipient_metrics <- metrics[
        as.character(metrics$role) == "recipient", , drop = FALSE
    ]
    if (nrow(recipient_metrics) == 0L) {
        return(empty())
    }
    recipient_metrics <- .assignment_recipient_depth(recipient_metrics, eps)
    table <- windows$table
    donor_contigs <- assembly$contigs$donor
    if (is.null(donor_contigs)) {
        donor_contigs <- data.frame()
    }

    for (i in seq_len(nrow(calls))) {
        call <- calls[i, , drop = FALSE]
        sample <- as.character(call$sample[[1L]])
        group <- as.character(call$group[[1L]])
        start_iid <- as.integer(call$start_iid[[1L]])
        end_iid <- as.integer(call$end_iid[[1L]])
        contig_rows <- recipient_metrics[
            as.character(recipient_metrics$sample) == sample &
                as.character(recipient_metrics$group) == group,
            ,
            drop = FALSE
        ]
        contig_rows <- contig_rows[
            is.finite(contig_rows$normalised_depth), , drop = FALSE
        ]
        if (nrow(contig_rows) == 0L) {
            next
        }
        split_key <- paste(
            contig_rows$subgenome, contig_rows$chr_id, sep = "\r"
        )
        candidates <- lapply(split(contig_rows, split_key), function(rows) {
            inside <- rows[
                rows$iid >= start_iid & rows$iid <= end_iid, ,
                drop = FALSE
            ]
            if (nrow(inside) == 0L) {
                return(NULL)
            }
            outside <- rows[
                rows$iid < start_iid | rows$iid > end_iid, ,
                drop = FALSE
            ]
            fraction <- nrow(inside) / nrow(rows)
            baseline <- if (fraction >= 0.95 || nrow(outside) < 10L) {
                0
            } else {
                mean(outside$normalised_depth, na.rm = TRUE)
            }
            data.frame(
                subgenome = as.character(rows$subgenome[[1L]]),
                chr_id = as.character(rows$chr_id[[1L]]),
                chr_name = as.character(rows$chr_name[[1L]]),
                normalised_depth = mean(inside$normalised_depth, na.rm = TRUE),
                baseline = baseline,
                depletion = baseline -
                    mean(inside$normalised_depth, na.rm = TRUE),
                stringsAsFactors = FALSE
            )
        })
        candidates <- do.call(rbind, candidates)
        if (is.null(candidates) || nrow(candidates) == 0L) {
            next
        }
        candidates <- candidates[order(-candidates$depletion), , drop = FALSE]
        best <- candidates[1L, , drop = FALSE]
        margin <- if (nrow(candidates) > 1L) {
            best$depletion[[1L]] - candidates$depletion[[2L]]
        } else {
            NA_real_
        }
        weights <- exp(candidates$depletion - max(candidates$depletion))
        weights <- weights / sum(weights)
        best_weight <- weights[[1L]]
        assignment <- if (nrow(candidates) == 1L) {
            "assigned"
        } else if (margin >= 0.5) {
            "assigned"
        } else if (margin >= 0.15) {
            "ambiguous"
        } else {
            "no_recipient_loss"
        }

        recipient_rows <- table[
            as.character(table$chr_id) == best$chr_id[[1L]] &
                table$iid >= start_iid & table$iid <= end_iid,
            ,
            drop = FALSE
        ]
        calls$recipient_subgenome[[i]] <- best$subgenome[[1L]]
        calls$recipient_chr_id[[i]] <- best$chr_id[[1L]]
        calls$recipient_chr_name[[i]] <- best$chr_name[[1L]]
        calls$recipient_start_bp[[i]] <- if (nrow(recipient_rows) > 0L) {
            min(recipient_rows$start, na.rm = TRUE)
        } else {
            NA_real_
        }
        calls$recipient_end_bp[[i]] <- if (nrow(recipient_rows) > 0L) {
            max(recipient_rows$stop, na.rm = TRUE)
        } else {
            NA_real_
        }
        calls$recipient_weight[[i]] <- best_weight
        calls$recipient_margin[[i]] <- margin
        calls$recipient_assignment[[i]] <- assignment

        donor_subgenome <- .assignment_donor_subgenome(
            scores, sample, group, start_iid, end_iid
        )
        donor_rows <- donor_contigs[
            as.character(donor_contigs$group) == group &
                (!is.na(donor_subgenome) &
                    as.character(donor_contigs$subgenome) == donor_subgenome),
            ,
            drop = FALSE
        ]
        if (nrow(donor_rows) == 0L) {
            donor_rows <- donor_contigs[
                as.character(donor_contigs$group) == group, ,
                drop = FALSE
            ]
        }
        if (nrow(donor_rows) > 0L) {
            donor_chr <- as.character(donor_rows$id[[1L]])
            donor_windows <- table[
                as.character(table$chr_id) == donor_chr &
                    table$iid >= start_iid & table$iid <= end_iid,
                ,
                drop = FALSE
            ]
            calls$donor_subgenome[[i]] <- if (is.na(donor_subgenome)) {
                NA_character_
            } else {
                donor_subgenome
            }
            calls$donor_chr_id[[i]] <- donor_chr
            calls$donor_chr_name[[i]] <- as.character(donor_rows$name[[1L]])
            calls$donor_start_bp[[i]] <- if (nrow(donor_windows) > 0L) {
                min(donor_windows$start, na.rm = TRUE)
            } else {
                NA_real_
            }
            calls$donor_end_bp[[i]] <- if (nrow(donor_windows) > 0L) {
                max(donor_windows$stop, na.rm = TRUE)
            } else {
                NA_real_
            }
        }
    }
    calls
}
