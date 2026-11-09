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

#' Canonical BEDPE column order
#'
#' The same list declares the empty and non-empty BEDPE schema. Call fields
#' come from the canonical call schema; the six interval columns are added by
#' `.annotate_call_intervals()`.
#' @noRd
bam_bedpe_columns <- function() {
    c(
        "sample", "chr_id.don", "start.don", "stop.don",
        "chr_id.rec", "start.rec", "stop.rec", "width_bp",
        "score_mean", "donor_mean", "recipient_mean", "class",
        "class_group", "class_evidence", "class_confidence", "class_status",
        "copy_state", "total_evidence", "component_status",
        "mappability_status", "mapq_fraction",
        "confidence", "method", "mode", "group",
        "start_iid", "end_iid", "pass", "recipient_component",
        "recipient_chr_name", "recipient_group", "recipient_assignment",
        "recipient_cross_group", "recipient_weight",
        "recipient_expected_breadth", "recipient_observed_breadth",
        "recipient_deletion_width_bp",
        "recipient_alternative", "recipient_mean_depletion",
        "recipient_depletion_se", "recipient_depletion_lower",
        "recipient_depletion_z", "recipient_width_overlap",
        "recipient_fraction_ratio", "recipient_mapq_fraction",
        "total_mean",
        "donor_component", "donor_chr_name"
    )
}

.bam_bedpe_types <- function() {
    types <- call_schema_types()
    types[c("chr_id.don", "chr_id.rec")] <- "character"
    types[c("start.don", "stop.don", "start.rec", "stop.rec")] <- "numeric"
    types
}

#' Project calls onto the canonical BEDPE schema
#'
#' Missing columns become typed `NA` columns so an empty export and a populated
#' export always describe the same fields.
#' @noRd
.bam_bedpe_table_schema <- function(calls) {
    columns <- bam_bedpe_columns()
    types <- .bam_bedpe_types()
    n <- if (is.data.frame(calls)) nrow(calls) else 0L
    out <- lapply(columns, function(name) {
        if (is.data.frame(calls) && name %in% names(calls)) {
            calls[[name]]
        } else {
            .call_schema_vector(unname(types[[name]]), n)
        }
    })
    names(out) <- columns
    as.data.frame(
        out, stringsAsFactors = FALSE, check.names = FALSE, optional = TRUE
    )
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
#' class, group class, class evidence/confidence/status, confidence, method,
#' mode, calibration pass flag and the recipient assignment fields. When the
#' calls carry a positional
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
.bedpe_table <- function(x, samples = NULL, selected = NULL,
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
    calls <- explorer_calls(calls, list(pass_only = pass_only))$calls
    if (nrow(calls) == 0L) {
        return(.bam_bedpe_table_schema(NULL))
    }
    .bam_bedpe_table_schema(calls)
}

writeCallsToBedpe <- function(x, file, samples = NULL, selected = NULL,
                              pass_only = TRUE, windows = NULL) {
    if (!is.character(file) || length(file) != 1L || is.na(file) ||
        !nzchar(file)) {
        .bedpe_error("`file` must be a single output path")
    }
    out <- .bedpe_table(
        x,
        samples = samples,
        selected = selected,
        pass_only = pass_only,
        windows = windows
    )
    utils::write.table(
        out, file, sep = "\t", quote = FALSE, row.names = FALSE
    )
    invisible(file)
}

isTRUE_vec <- function(x) {
    !is.na(x) & as.logical(x)
}
