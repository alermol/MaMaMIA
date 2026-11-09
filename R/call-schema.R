## Canonical schema for one introgression call.
##
## The call record is used in several shapes: the validated `Calls` data frame,
## the BEDPE export, the descriptive bundle dictionary and the karyotype
## payload. This module is the single declaration of the call fields and their
## type, default, display label, descriptive export name and description.
##
## Consumers derive their names and metadata from `call_schema()` instead of
## repeating the field list. Keep new fields here first, then migrate any
## consumer that needs them.

.call_schema_fields <- local({
    name <- c(
        "sample", "group", "start_iid", "end_iid", "start_bp", "end_bp",
        "width_bp", "score_mean", "donor_mean", "recipient_mean",
        "total_mean", "donor_component_spread",
        "recipient_component_spread", "n_donor_supported",
        "mapq_fraction", "n_windows", "method", "mode", "metric",
        "length_bin", "threshold", "pass", "class", "class_evidence",
        "copy_state", "total_evidence", "class_confidence",
        "class_status", "component_status", "mappability_status",
        "confidence", "class_group", "recipient_component",
        "recipient_chr_id", "recipient_chr_name", "recipient_group",
        "recipient_start_bp", "recipient_end_bp",
        "recipient_run_start_bp", "recipient_run_end_bp",
        "recipient_deletion_width_bp", "recipient_expected_breadth",
        "recipient_observed_breadth", "recipient_mapq_fraction",
        "recipient_mean_depletion",
        "recipient_depletion_se", "recipient_depletion_lower",
        "recipient_depletion_z", "recipient_width_overlap",
        "recipient_fraction_ratio", "copy_state_group",
        "recipient_weight", "recipient_margin",
        "recipient_assignment", "recipient_cross_group",
        "recipient_reliable", "recipient_strong",
        "recipient_alternative", "recipient_alternative_score",
        "donor_component", "donor_chr_id", "donor_chr_name",
        "donor_start_bp", "donor_end_bp"
    )
    type <- c(
        sample = "character",
        group = "character",
        start_iid = "integer",
        end_iid = "integer",
        start_bp = "numeric",
        end_bp = "numeric",
        width_bp = "numeric",
        score_mean = "numeric",
        donor_mean = "numeric",
        recipient_mean = "numeric",
        total_mean = "numeric",
        donor_component_spread = "numeric",
        recipient_component_spread = "numeric",
        n_donor_supported = "numeric",
        mapq_fraction = "numeric",
        n_windows = "integer",
        method = "character",
        mode = "character",
        metric = "character",
        length_bin = "character",
        threshold = "numeric",
        pass = "logical",
        class = "character",
        class_evidence = "character",
        copy_state = "character",
        total_evidence = "character",
        class_confidence = "character",
        class_status = "character",
        component_status = "character",
        mappability_status = "character",
        confidence = "character",
        class_group = "character",
        recipient_component = "character",
        recipient_chr_id = "character",
        recipient_chr_name = "character",
        recipient_group = "character",
        recipient_start_bp = "numeric",
        recipient_end_bp = "numeric",
        recipient_run_start_bp = "numeric",
        recipient_run_end_bp = "numeric",
        recipient_deletion_width_bp = "numeric",
        recipient_expected_breadth = "numeric",
        recipient_observed_breadth = "numeric",
        recipient_mapq_fraction = "numeric",
        recipient_mean_depletion = "numeric",
        recipient_depletion_se = "numeric",
        recipient_depletion_lower = "numeric",
        recipient_depletion_z = "numeric",
        recipient_width_overlap = "numeric",
        recipient_fraction_ratio = "numeric",
        copy_state_group = "character",
        recipient_weight = "numeric",
        recipient_margin = "numeric",
        recipient_assignment = "character",
        recipient_cross_group = "logical",
        recipient_reliable = "logical",
        recipient_strong = "logical",
        recipient_alternative = "character",
        recipient_alternative_score = "numeric",
        donor_component = "character",
        donor_chr_id = "character",
        donor_chr_name = "character",
        donor_start_bp = "numeric",
        donor_end_bp = "numeric"
    )
    if (!identical(names(type), name)) {
        stop("The call schema type names do not match the field names",
            call. = FALSE
        )
    }
    type <- unname(type)
    label <- c(
        "Sample", "Group", "Start window", "End window",
        "Group start (bp)", "Group end (bp)", "Width (bp)", "Score",
        "Donor mean", "Recipient mean", "Total mean", "Donor spread",
        "Recipient spread", "Donor sets", "MAPQ fraction", "Windows",
        "Method", "Mode", "Metric", "Length bin", "Threshold", "Pass",
        "Class", "Class evidence", "Copy state", "Total evidence",
        "Class confidence", "Status", "Component status", "Mappability",
        "Confidence", "Group class", "Recipient set", "Recipient ID",
        "Recipient", "Recipient group", "Recipient start (bp)",
        "Recipient end (bp)", "Depletion start (bp)",
        "Depletion end (bp)", "Depletion width (bp)",
        "Expected breadth", "Observed breadth", "Recipient MAPQ",
        "Mean depletion", "Depletion SE", "Depletion lower",
        "Depletion Z", "Width overlap", "Width ratio", "Group copy state",
        "Recipient weight", "Recipient margin", "Placement",
        "Cross-group", "Recipient reliable", "Recipient strong",
        "Alternative recipient", "Alternative score", "Donor set",
        "Donor ID", "Donor", "Donor start (bp)", "Donor end (bp)"
    )
    export <- c(
        "sample_name", "chromosome_group", "start_window_index",
        "end_window_index", "group_start_bp", "group_end_bp",
        "call_width_bp", "mean_introgression_score",
        "mean_donor_residual", "mean_recipient_residual",
        "mean_total_depth_residual", "donor_component_spread",
        "recipient_component_spread", "supporting_donor_components",
        "high_mapq_read_fraction", "windows_in_call",
        "scoring_method", "analysis_mode", "coverage_metric",
        "call_length_bin", "calibration_threshold",
        "passes_calibration", "call_class", "call_class_evidence",
        "copy_number_state", "total_depth_evidence", "call_class_confidence",
        "review_status", "component_consistency", "mappability_quality",
        "donor_confidence", "group_level_call_class",
        "recipient_copy_set", "recipient_contig_id",
        "recipient_chromosome", "recipient_chromosome_group",
        "recipient_start_bp", "recipient_end_bp",
        "recipient_deletion_run_start_bp",
        "recipient_deletion_run_end_bp",
        "recipient_deletion_run_width_bp",
        "recipient_expected_breadth", "recipient_observed_breadth",
        "recipient_high_mapq_read_fraction",
        "recipient_mean_depletion",
        "recipient_depletion_standard_error",
        "recipient_depletion_lower_bound",
        "recipient_depletion_z_score",
        "recipient_width_reciprocal_overlap",
        "donor_to_recipient_width_ratio",
        "group_level_copy_state",
        "recipient_assignment_weight", "recipient_assignment_margin",
        "recipient_assignment_status", "recipient_is_cross_group",
        "recipient_candidate_reliable", "recipient_candidate_strong",
        "alternative_recipient_chromosome",
        "alternative_recipient_score", "donor_copy_set",
        "donor_contig_id", "donor_chromosome",
        "donor_start_bp", "donor_end_bp"
    )
    description <- c(
        "Sample name read from the BAM @RG SM tag.",
        "Chromosome homology group used for scoring and segmentation, for example 7.",
        "First analysis window index in the call, zero-based within the group.",
        "Last analysis window index in the call, zero-based and inclusive.",
        "Call start in group-level coordinates, in base pairs, scaled to the donor group length.",
        "Call end in group-level coordinates, in base pairs.",
        "Call width in base pairs, calculated as end_bp minus start_bp.",
        "Mean introgression score across the call; higher values mean stronger donor evidence.",
        "Mean donor residual across the call; positive values mean donor enrichment.",
        "Mean recipient residual across the call; negative values mean recipient loss.",
        "Mean total-depth residual across the call; positive values mean a total gain and negative values a total loss.",
        "Spread of donor residuals across donor copy sets; large values flag component heterogeneity.",
        "Spread of recipient residuals across recipient copy sets; large values flag component heterogeneity.",
        "Number of donor copy sets that support the call.",
        "Mean fraction of high-MAPQ reads across the call; values below 0.5 flag low mappability.",
        "Number of analysis windows covered by the call.",
        "Scoring method: paired_lr for a paired donor/recipient log ratio, or donor_only.",
        "Analysis mode: single for one sample, or multi for a cohort.",
        "Coverage metric used for scoring: breadth, count or combined.",
        "Segment-length bin used to select the calibration threshold.",
        "Calibrated score threshold for this call's mode and length bin.",
        "TRUE when score_mean exceeds the calibrated threshold.",
        "Final call class after recipient assignment: substitution, addition, deletion, complex or other.",
        "Evidence that supported the call: donor, recipient, donor+recipient, total imbalance or none.",
        "Copy-number interpretation: balanced, gain, loss, compound or neutral.",
        "Total-depth evidence: gain, loss, neutral, or blank when unavailable.",
        "Class-level confidence: high, moderate or low.",
        "Review status: confident, ambiguous or needs_review.",
        "Component consistency: consistent or heterogeneous.",
        "Mappability flag: ok or low; low means the high-MAPQ read fraction is below 0.5.",
        "Legacy donor-effect confidence: high, moderate or low.",
        "Group-level call class before recipient assignment.",
        "Copy-set label of the assigned recipient contig.",
        "Assembly contig identifier of the assigned recipient.",
        "Display name of the assigned recipient contig.",
        "Homology group of the assigned recipient contig; differs from chromosome_group for cross-group placements.",
        "Start of the corresponding interval on the assigned recipient contig, in base pairs.",
        "End of that interval on the assigned recipient contig, in base pairs.",
        "Start of the contiguous recipient depletion run used for assignment, in base pairs.",
        "End of that depletion run, in base pairs.",
        "Width of the recipient depletion run, in base pairs; blank when no run was found.",
        "Baseline breadth expected for the recipient contig.",
        "Mean breadth observed inside the call interval on the recipient contig.",
        "Mean high-MAPQ read fraction in the recipient interval.",
        "Mean corrected log-depth depletion in the recipient interval; positive values mean loss.",
        "Standard error of the mean recipient depletion.",
        "Lower confidence bound of the recipient depletion estimate.",
        "Z score of the recipient depletion estimate.",
        "Reciprocal overlap between the donor segment and the recipient depletion run, from 0 to 1.",
        "Donor width fraction divided by recipient run width fraction; values near 1 mean matched widths.",
        "Copy state before recipient assignment.",
        "Weight of the chosen recipient among plausible candidates, from 0 to 1.",
        "Score difference between the chosen recipient and the best alternative; larger values are more decisive.",
        "Recipient placement outcome: assigned, ambiguous, cross_group or no_recipient_loss.",
        "TRUE when the assigned recipient contig is outside the call's homology group.",
        "TRUE when the recipient candidate passed the calibrated reliability screen.",
        "TRUE when the recipient candidate passed the stronger depletion screen.",
        "Display name of the next-best recipient contig, when available.",
        "Assignment score of the next-best recipient contig.",
        "Copy-set label of the donor contig that contributed the call.",
        "Assembly contig identifier of the donor contig.",
        "Display name of the donor contig.",
        "Donor segment start in base pairs.",
        "Donor segment end in base pairs."
    )
    if (!identical(
        length(name), length(type)
    ) || !identical(length(name), length(label)) ||
        !identical(length(name), length(export)) ||
        !identical(length(name), length(description))) {
        stop("The call schema columns have unequal lengths", call. = FALSE)
    }
    if (anyDuplicated(name) > 0L) {
        stop("The call schema has duplicate field names", call. = FALSE)
    }
    if (anyDuplicated(export) > 0L) {
        stop("The call schema has duplicate export names", call. = FALSE)
    }
    default <- lapply(type, function(value) {
        switch(value,
            character = NA_character_,
            numeric = NA_real_,
            integer = NA_integer_,
            logical = NA,
            stop("Unknown call schema type: ", value, call. = FALSE)
        )
    })
    names(default) <- name
    default[["pass"]] <- TRUE
    default[["recipient_cross_group"]] <- FALSE
    default[["recipient_reliable"]] <- FALSE
    default[["recipient_strong"]] <- FALSE
    data.frame(
        name = name,
        type = type,
        default = I(default),
        label = label,
        export = export,
        description = description,
        stringsAsFactors = FALSE
    )
})

#' Canonical call schema
#'
#' One row per call field with its type, default, display label, descriptive
#' export name and description.
#' @noRd
call_schema <- function() {
    .call_schema_fields
}

#' Canonical call field names
#' @noRd
call_schema_names <- function() {
    .call_schema_fields$name
}

#' Canonical call field types
#' @noRd
call_schema_types <- function() {
    stats::setNames(.call_schema_fields$type, .call_schema_fields$name)
}

#' Canonical call field defaults
#' @noRd
call_schema_defaults <- function() {
    stats::setNames(.call_schema_fields$default, .call_schema_fields$name)
}

#' Canonical call display labels
#' @noRd
call_schema_labels <- function() {
    stats::setNames(.call_schema_fields$label, .call_schema_fields$name)
}

#' Canonical call export names
#' @noRd
call_schema_export <- function() {
    stats::setNames(.call_schema_fields$export, .call_schema_fields$name)
}

.call_schema_vector <- function(type, n = 0L) {
    switch(type,
        character = character(n),
        numeric = numeric(n),
        integer = integer(n),
        logical = logical(n),
        stop("Unknown call schema type: ", type, call. = FALSE)
    )
}

#' Empty calls data frame with canonical columns
#'
#' @param columns Optional subset of canonical field names.
#' @param n Number of rows, zero by default.
#' @noRd
call_schema_empty <- function(columns = NULL, n = 0L) {
    if (is.null(columns)) {
        columns <- call_schema_names()
    }
    columns <- as.character(columns)
    types <- call_schema_types()
    absent <- setdiff(columns, names(types))
    if (length(absent) > 0L) {
        stop(
            "Unknown call schema column(s): ",
            paste(absent, collapse = ", "),
            call. = FALSE
        )
    }
    out <- lapply(columns, function(name) {
        .call_schema_vector(unname(types[[name]]), n)
    })
    names(out) <- columns
    as.data.frame(
        out, stringsAsFactors = FALSE, check.names = FALSE, optional = TRUE
    )
}

.call_schema_cast <- function(value, type) {
    switch(type,
        character = as.character(value),
        numeric = as.numeric(value),
        integer = as.integer(value),
        logical = as.logical(value),
        stop("Unknown call schema type: ", type, call. = FALSE)
    )
}

#' Build a canonical call record
#'
#' Coerces values to the canonical field types and fills missing fields with
#' their schema defaults. Extra fields are kept after the canonical fields so
#' callers such as the karyotype payload can add presentation-only values.
#'
#' @param values Named list of values.
#' @param include_missing Fill absent canonical fields with their defaults.
#' @noRd
call_record <- function(values = list(), include_missing = TRUE) {
    if (!is.list(values)) {
        stop("`values` must be a list", call. = FALSE)
    }
    schema <- call_schema()
    supplied <- names(values)
    if (is.null(supplied)) {
        supplied <- character()
    }
    canonical <- schema$name[schema$name %in% supplied]
    if (isTRUE(include_missing)) {
        canonical <- schema$name
    }
    out <- lapply(canonical, function(name) {
        row <- match(name, schema$name)
        value <- if (name %in% supplied) {
            values[[name]]
        } else {
            schema$default[[row]]
        }
        .call_schema_cast(value, schema$type[[row]])
    })
    names(out) <- canonical
    for (name in setdiff(supplied, schema$name)) {
        out[[name]] <- values[[name]]
    }
    out
}

#' Normalize a call record to the canonical schema
#' @noRd
call_record_normalize <- function(x) {
    call_record(x, include_missing = TRUE)
}
