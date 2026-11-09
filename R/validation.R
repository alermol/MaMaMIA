## Validation of the BamAnalysisResult contract.

#' Column contracts for a BamAnalysisResult
#' @noRd
bam_result_contracts <- function() {
    list(
        metrics = c(
            "sample", "chr_id", "chr_name", "group", "role", "component",
            "copy_number", "iid", "start", "stop", "length", "count",
            "covered", "mapq_high", "mapq_high_fraction", "breadth", "gc",
            "lambda_breadth", "lambda_count", "read_length"
        ),
        scores = c(
            "sample", "group", "iid", "start", "stop", "mode", "metric",
            "score_type", "score", "donor_mean", "recipient_mean",
            "donor_component", "n_donor", "n_recipient",
            "donor_component_spread",
            "recipient_component_spread", "n_donor_supported",
            "mapq_fraction", "total_residual"
        ),
        masks = c(
            "group", "iid", "start", "stop", "f", "posterior", "mask",
            "source", "version"
        ),
        calibration = c(
            "mode", "length_bin", "n_control", "q99", "threshold", "fdr",
            "version"
        ),
        ## The call record is declared once in R/call-schema.R.
        calls = call_schema_names()
    )
}

#' Validate a BamAnalysisResult
#'
#' Checks the named fields and the WindowMetrics/ScoreTable/MaskTable/
#' Calibration/Calls column contracts.
#'
#' @param x A `BamAnalysisResult`.
#'
#' @return `x`, invisibly.
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
    contracts <- bam_result_contracts()
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
