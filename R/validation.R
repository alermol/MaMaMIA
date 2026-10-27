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
            "sample", "chr_id", "chr_name", "group", "role", "subgenome",
            "iid", "start", "stop", "length", "count", "covered", "breadth",
            "gc", "lambda_breadth", "lambda_count", "read_length"
        ),
        scores = c(
            "sample", "group", "iid", "start", "stop", "mode", "metric",
            "score", "donor_mean", "recipient_mean", "donor_subgenome",
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
            "confidence", "method", "mode"
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
