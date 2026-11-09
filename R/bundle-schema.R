## Declared schema for the QC tables in a result bundle.
##
## The bundle renames raw result columns to descriptive export names. This
## module declares, for every QC table, the rename map and the columns that
## must exist after renaming. `bam_bundle_validate_qc_tables()` fails before
## any file is written when the result no longer satisfies the contract.

## Descriptive names for the per-window metric table.
bam_bundle_metric_columns <- function() {
    c(
        sample = "sample_name",
        chr_id = "contig_id",
        chr_name = "contig_name",
        group = "chromosome_group",
        role = "role",
        component = "component",
        copy_number = "copy_number",
        iid = "window_index",
        start = "start_bp",
        stop = "end_bp",
        length = "window_length_bp",
        count = "read_count",
        covered = "covered_bases",
        mapq_high = "high_mapq_read_count",
        mapq_high_fraction = "high_mapq_read_fraction",
        breadth = "breadth",
        gc = "gc_fraction",
        lambda_breadth = "expected_breadth",
        lambda_count = "expected_read_count",
        read_length = "read_length_bp"
    )
}

## Descriptive names for the per-window score table.
bam_bundle_score_columns <- function() {
    c(
        sample = "sample_name",
        group = "chromosome_group",
        iid = "window_index",
        start = "start_bp",
        stop = "end_bp",
        mode = "analysis_mode",
        metric = "coverage_metric",
        score_type = "score_type",
        score = "introgression_score",
        donor_mean = "mean_donor_coverage",
        recipient_mean = "mean_recipient_coverage",
        donor_component = "donor_component",
        n_donor = "n_donor_contigs",
        n_recipient = "n_recipient_contigs",
        donor_component_spread = "donor_component_spread",
        recipient_component_spread = "recipient_component_spread",
        n_donor_supported = "supporting_donor_contigs",
        mapq_fraction = "high_mapq_read_fraction",
        chr_name = "group_name",
        length = "group_length_bp",
        n_per_chromosome = "windows_per_chromosome",
        gc_donor = "donor_gc_fraction",
        gc_recipient = "recipient_gc_fraction",
        log_donor = "log_donor_coverage",
        log_recipient = "log_recipient_coverage",
        log_total = "log_total_coverage",
        donor_residual = "donor_residual",
        recipient_residual = "recipient_residual",
        total_residual = "total_depth_residual",
        score_raw = "raw_introgression_score"
    )
}

## Descriptive names for the mask table.
bam_bundle_mask_columns <- function() {
    c(
        group = "chromosome_group",
        iid = "window_index",
        start = "start_bp",
        stop = "end_bp",
        f = "mask_signal",
        posterior = "mask_posterior",
        mask = "masked",
        source = "mask_source",
        version = "mask_version"
    )
}

## Descriptive names for the calibration table.
bam_bundle_calibration_columns <- function() {
    c(
        mode = "analysis_mode",
        length_bin = "call_length_bin",
        n_control = "control_calls",
        q99 = "null_q99",
        threshold = "calibration_threshold",
        fdr = "target_fdr",
        version = "calibration_version"
    )
}

#' Declared QC tables for the result bundle
#'
#' Each table declares its descriptive rename map, the columns the validation
#' contract requires after renaming, and the file it is written to.
#' @noRd
bam_bundle_qc_schema <- function() {
    contracts <- bam_result_contracts()
    spec <- function(table, rename, file, required_source) {
        unmapped <- setdiff(required_source, names(rename))
        if (length(unmapped) > 0L) {
            stop(
                "Bundle QC schema drift in `", table,
                "`: no export name for contract column(s): ",
                paste(unmapped, collapse = ", "),
                call. = FALSE
            )
        }
        required <- unname(rename[required_source])
        if (anyDuplicated(required) > 0L) {
            stop(
                "Bundle QC schema drift in `", table,
                "`: duplicate export name(s): ",
                paste(unique(required[duplicated(required)]), collapse = ", "),
                call. = FALSE
            )
        }
        list(
            table = table,
            file = file,
            rename = rename,
            required = required,
            optional = setdiff(unname(rename), required)
        )
    }
    list(
        metrics = spec(
            "metrics", bam_bundle_metric_columns(),
            "qc/metrics-windows.csv.gz", contracts$metrics
        ),
        scores = spec(
            "scores", bam_bundle_score_columns(),
            "qc/scores-windows.csv.gz", contracts$scores
        ),
        masks = spec(
            "masks", bam_bundle_mask_columns(),
            "qc/masks-windows.csv.gz", contracts$masks
        ),
        calibration = spec(
            "calibration", bam_bundle_calibration_columns(),
            "qc/calibration.csv", contracts$calibration
        )
    )
}

#' Optional descriptive columns for a bundle QC table
#'
#' Optional columns are part of the declared export vocabulary but are not
#' guaranteed by the result contract. They may be absent in a valid result.
#' @param table One of `"metrics"`, `"scores"`, `"masks"` or `"calibration"`.
#' @return A character vector of optional descriptive column names.
#' @noRd
bam_bundle_optional_columns <- function(table) {
    schema <- bam_bundle_qc_schema()
    spec <- schema[[as.character(table)[[1L]]]]
    if (is.null(spec)) {
        stop("Unknown bundle QC table: ", table, call. = FALSE)
    }
    spec$optional
}

#' Validate renamed QC tables before writing a bundle
#'
#' A zero-row table with the required columns is a genuine empty result and
#' passes. A missing table, a non-data-frame, or a missing required column is
#' contract drift and fails with the table and column names in the message.
#'
#' @param tables Named list of renamed QC tables: `metrics`, `scores`, and
#'   optionally `masks` and `calibration`.
#' @param include_masks Validate the mask table.
#' @param include_calibration Validate the calibration table.
#' @param required_tables Tables that must always be present.
#' @return `tables`, invisibly.
#' @noRd
bam_bundle_validate_qc_tables <- function(tables, include_masks = TRUE,
                                          include_calibration = TRUE,
                                          required_tables = c(
                                              "metrics", "scores"
                                          )) {
    schema <- bam_bundle_qc_schema()
    wanted <- as.character(required_tables)
    if (isTRUE(include_masks)) {
        wanted <- c(wanted, "masks")
    }
    if (isTRUE(include_calibration)) {
        wanted <- c(wanted, "calibration")
    }
    wanted <- unique(wanted)
    if (!is.list(tables)) {
        stop("Bundle QC tables must be a named list", call. = FALSE)
    }
    for (table in wanted) {
        spec <- schema[[table]]
        if (!table %in% names(tables)) {
            stop(
                "Bundle QC table `", table, "` is missing (expected ",
                spec$file, ")",
                call. = FALSE
            )
        }
        value <- tables[[table]]
        if (is.null(value)) {
            stop(
                "Bundle QC table `", table, "` is missing (expected ",
                spec$file, ")",
                call. = FALSE
            )
        }
        if (!is.data.frame(value)) {
            stop(
                "Bundle QC table `", table, "` must be a data frame",
                call. = FALSE
            )
        }
        missing <- setdiff(spec$required, names(value))
        if (length(missing) > 0L) {
            stop(
                "Bundle QC table `", table, "` is missing required column(s): ",
                paste(missing, collapse = ", "),
                call. = FALSE
            )
        }
    }
    invisible(tables)
}
