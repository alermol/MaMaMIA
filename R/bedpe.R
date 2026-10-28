## BEDPE export for introgression calls.
##
## A call is group-level: `start_iid`/`end_iid` are relative windows shared by
## the homoeologs in one chromosome group. For a BEDPE we materialise one donor
## interval and one recipient interval. When assignCallPositions() has added
## the assigned contigs and coordinates, those are used. Otherwise the donor
## component is the one recorded by a `max` aggregation (or the first donor
## contig in the group) and the recipient interval uses the first recipient
## component in the group; the `group` column and the score fields retain the
## group-level meaning.

.bedpe_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIABedpeError", fmt, ...)
}

.annotate_call_intervals <- function(calls, windows) {
    if (nrow(calls) == 0L) {
        calls$chr_id.don <- character(0)
        calls$start.don <- numeric(0)
        calls$stop.don <- numeric(0)
        calls$chr_id.rec <- character(0)
        calls$start.rec <- numeric(0)
        calls$stop.rec <- numeric(0)
        return(calls)
    }
    if (!is.list(windows) || is.null(windows$table)) {
        .bedpe_error("`windows` must be an AnalysisWindows list")
    }
    table <- windows$table
    n_per <- windows$n_per_chromosome
    donor <- table[table$role == "donor", , drop = FALSE]
    recipient <- table[table$role == "recipient", , drop = FALSE]
    pick <- function(role_table, group, start_iid, end_iid,
                     component = NULL) {
        rows <- role_table[role_table$group == group, , drop = FALSE]
        if (nrow(rows) == 0L) {
            return(list(chr = NA_character_, start = NA_real_, stop = NA_real_))
        }
        if (is.null(component) || is.na(component) ||
            !component %in% rows$component) {
            component <- rows$component[[1L]]
        }
        rows <- rows[rows$component == component, , drop = FALSE]
        start_row <- rows[rows$iid == start_iid, , drop = FALSE]
        if (nrow(start_row) == 0L) {
            return(list(chr = NA_character_, start = NA_real_, stop = NA_real_))
        }
        length_bp <- max(rows$stop, na.rm = TRUE)
        list(
            chr = as.character(start_row$chr_id[[1L]]),
            start = as.numeric(start_row$start[[1L]]),
            stop = floor((end_iid + 1) * length_bp / n_per)
        )
    }
    donor_intervals <- lapply(seq_len(nrow(calls)), function(i) {
        pick(
            donor, calls$group[[i]], calls$start_iid[[i]], calls$end_iid[[i]],
            if ("donor_component" %in% names(calls)) {
                calls$donor_component[[i]]
            } else {
                NULL
            }
        )
    })
    recipient_intervals <- lapply(seq_len(nrow(calls)), function(i) {
        pick(
            recipient, calls$group[[i]], calls$start_iid[[i]],
            calls$end_iid[[i]]
        )
    })
    calls$chr_id.don <- vapply(donor_intervals, `[[`, character(1), "chr")
    calls$start.don <- vapply(donor_intervals, `[[`, numeric(1), "start")
    calls$stop.don <- vapply(donor_intervals, `[[`, numeric(1), "stop")
    calls$chr_id.rec <- vapply(recipient_intervals, `[[`, character(1), "chr")
    calls$start.rec <- vapply(recipient_intervals, `[[`, numeric(1), "start")
    calls$stop.rec <- vapply(recipient_intervals, `[[`, numeric(1), "stop")

    ## Prefer the per-call assignment produced by assignCallPositions(); the
    ## group fallback above is only for older calls without those columns.
    if (all(c("recipient_chr_id", "recipient_start_bp",
              "recipient_end_bp") %in% names(calls))) {
        ok <- !is.na(calls$recipient_chr_id)
        calls$chr_id.rec[ok] <- as.character(calls$recipient_chr_id[ok])
        calls$start.rec[ok] <- as.numeric(calls$recipient_start_bp[ok])
        calls$stop.rec[ok] <- as.numeric(calls$recipient_end_bp[ok])
    }
    if (all(c("donor_chr_id", "donor_start_bp", "donor_end_bp") %in%
        names(calls))) {
        ok <- !is.na(calls$donor_chr_id)
        calls$chr_id.don[ok] <- as.character(calls$donor_chr_id[ok])
        calls$start.don[ok] <- as.numeric(calls$donor_start_bp[ok])
        calls$stop.don[ok] <- as.numeric(calls$donor_end_bp[ok])
    }
    if ("recipient_assignment" %in% names(calls)) {
        no_host <- !is.na(calls$recipient_assignment) &
            calls$recipient_assignment == "no_recipient_loss"
        calls$chr_id.rec[no_host] <- NA_character_
        calls$start.rec[no_host] <- NA_real_
        calls$stop.rec[no_host] <- NA_real_
    }
    calls
}

#' Write calls to a BEDPE file
#'
#' Writes one row per call with donor and recipient intervals plus the score,
#' class, confidence, method and mode. When the calls carry a positional
#' assignment, the recipient interval is the assigned contig and coordinates;
#' otherwise the first contig in the chromosome group is used. Calls with
#' `recipient_assignment = "no_recipient_loss"` are kept with `NA` recipient
#' coordinates because they have no host interval. Masked and failing calls are
#' excluded by default.
#'
#' @param x A `BamAnalysisResult` or a `Calls` data frame with `windows`
#'   supplied through `windows`.
#' @param file Output path.
#' @param samples Optional sample names to export.
#' @param selected Optional row indices (in `x$calls`) to export.
#' @param pass_only Export only calls with `pass == TRUE`.
#' @param windows Optional `AnalysisWindows` when `x` is a Calls data frame.
#'
#' @return The output path, invisibly.
#' @export
writeCallsToBedpe <- function(x, file, samples = NULL, selected = NULL,
                              pass_only = TRUE, windows = NULL) {
    if (inherits(x, "BamAnalysisResult")) {
        calls <- x$calls
        if (is.null(windows)) {
            windows <- x$windows
        }
    } else if (is.data.frame(x)) {
        calls <- x
    } else {
        .bedpe_error("`x` must be a BamAnalysisResult or a Calls data frame")
    }
    if (is.null(windows)) {
        .bedpe_error("`windows` is required when `x` is a Calls data frame")
    }
    if (!is.character(file) || length(file) != 1L || is.na(file) ||
        !nzchar(file)) {
        .bedpe_error("`file` must be a single output path")
    }
    if (nrow(calls) > 0L && !"chr_id.don" %in% names(calls)) {
        calls <- .annotate_call_intervals(calls, windows)
    }
    if (!is.null(selected)) {
        selected <- suppressWarnings(as.integer(selected))
        calls <- calls[seq_len(nrow(calls)) %in% selected, , drop = FALSE]
    }
    if (!is.null(samples)) {
        calls <- calls[calls$sample %in% samples, , drop = FALSE]
    }
    if (isTRUE(pass_only) && "pass" %in% names(calls)) {
        calls <- calls[isTRUE_vec(calls$pass), , drop = FALSE]
    }
    if (nrow(calls) == 0L) {
        utils::write.table(
            data.frame(
                sample = character(0), chr_id.don = character(0),
                start.don = numeric(0), stop.don = numeric(0),
                chr_id.rec = character(0), start.rec = numeric(0),
                stop.rec = numeric(0), width_bp = numeric(0),
                score_mean = numeric(0), donor_mean = numeric(0),
                recipient_mean = numeric(0), class = character(0),
                confidence = character(0), method = character(0),
                mode = character(0), group = character(0),
                class_group = character(0),
                recipient_component = character(0),
                recipient_chr_name = character(0),
                recipient_group = character(0),
                recipient_assignment = character(0),
                recipient_cross_group = logical(0),
                recipient_weight = numeric(0),
                recipient_expected_breadth = numeric(0),
                recipient_deletion_width_bp = numeric(0),
                recipient_alternative = character(0),
                recipient_mean_depletion = numeric(0),
                donor_component = character(0),
                donor_chr_name = character(0),
                stringsAsFactors = FALSE
            ),
            file, sep = "\t", quote = FALSE, row.names = FALSE
        )
        return(invisible(file))
    }
    columns <- c(
        "sample", "chr_id.don", "start.don", "stop.don",
        "chr_id.rec", "start.rec", "stop.rec", "width_bp",
        "score_mean", "donor_mean", "recipient_mean", "class",
        "class_group", "confidence", "method", "mode", "group",
        "start_iid", "end_iid", "recipient_component", "recipient_chr_name",
        "recipient_group", "recipient_assignment", "recipient_cross_group",
        "recipient_weight", "recipient_expected_breadth",
        "recipient_deletion_width_bp", "recipient_alternative",
        "recipient_mean_depletion", "donor_component", "donor_chr_name"
    )
    out <- calls[, intersect(columns, names(calls)), drop = FALSE]
    utils::write.table(
        out, file, sep = "\t", quote = FALSE, row.names = FALSE
    )
    invisible(file)
}

isTRUE_vec <- function(x) {
    !is.na(x) & as.logical(x)
}
