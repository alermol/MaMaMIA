## Validation of the BamAnalysisResult contract.

#' Validate a BamAnalysisResult
#'
#' Checks the named fields and the WindowMetrics/ScoreTable/MaskTable/
#' Calibration/Calls column contracts.
#'
#' @param x A `BamAnalysisResult`.
#'
#' @return `x`, invisibly.
#' @export
validateBamAnalysisResult <- function(x) {
    if (!inherits(x, "BamAnalysisResult") || !is.list(x)) {
        stop("`x` must be a BamAnalysisResult", call. = FALSE)
    }
    required <- c(
        "assembly", "config", "bams", "samples", "windows", "gc", "metrics",
        "scores", "masks", "calibration", "calls", "provenance"
    )
    missing <- setdiff(required, names(x))
    if (length(missing) > 0L) {
        stop(
            "BamAnalysisResult is missing field(s): ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    contracts <- list(
        metrics = c(
            "sample", "chr_id", "chr_name", "group", "role", "component",
            "iid", "start", "stop", "length", "count", "covered", "breadth",
            "gc", "lambda_breadth", "lambda_count", "read_length"
        ),
        scores = c(
            "sample", "group", "iid", "start", "stop", "mode", "metric",
            "score", "donor_mean", "recipient_mean", "donor_component",
            "n_donor", "n_recipient", "donor_copies", "donor_fraction"
        ),
        masks = c(
            "group", "iid", "start", "stop", "f", "posterior", "mask",
            "source", "version"
        ),
        calibration = c(
            "mode", "length_bin", "n_control", "q99", "threshold", "fdr",
            "version"
        ),
        calls = c(
            "sample", "group", "start_iid", "end_iid", "start_bp", "end_bp",
            "width_bp", "score_mean", "donor_mean", "recipient_mean", "class",
            "class_group", "confidence", "method", "mode",
            "recipient_component", "recipient_chr_id", "recipient_chr_name",
            "recipient_group", "recipient_start_bp", "recipient_end_bp",
            "recipient_run_start_bp", "recipient_run_end_bp",
            "recipient_deletion_width_bp", "recipient_expected_breadth",
            "recipient_observed_breadth", "recipient_copy_number",
            "recipient_mean_depletion", "affected_fraction",
            "affected_copies", "dosage_confidence", "recipient_weight",
            "recipient_margin", "recipient_assignment",
            "recipient_cross_group", "recipient_reliable",
            "recipient_strong", "recipient_alternative",
            "recipient_alternative_score", "donor_component", "donor_chr_id",
            "donor_chr_name", "donor_start_bp", "donor_end_bp"
        )
    )
    for (field in names(contracts)) {
        table <- x[[field]]
        if (!is.data.frame(table)) {
            stop("`", field, "` must be a data frame", call. = FALSE)
        }
        absent <- setdiff(contracts[[field]], names(table))
        if (length(absent) > 0L) {
            stop(
                "`", field, "` is missing column(s): ",
                paste(absent, collapse = ", "),
                call. = FALSE
            )
        }
    }
    invisible(x)
}
