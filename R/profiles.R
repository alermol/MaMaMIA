## Named analysis profiles.
##
## Profiles bundle every parameter that changes the candidate set: assembly
## components are selected at assembly resolution, while the run profile
## controls window resolution, metric, score, donor aggregation, masking,
## calibration and segmentation. Species names never enter these contracts.

.mamamia_profile_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAProfileError", fmt, ...)
}

#' Named analysis profiles
#'
#' @return A named list of profile contracts. Each contract has a `label`, a
#'   `description`, a `component_panel` flag and a `parameters` list suitable
#'   for [analyzeBamInput()].
#' @export
mamamia_analysis_profiles <- function() {
    balanced <- list(
        windows_per_chromosome = 1000L,
        metric = "breadth",
        score = "paired_log_ratio",
        donor_aggregation = "max",
        mapq = 20L,
        mask = list(
            physical_fallback_mb = 1,
            recurrence_threshold = 0.4
        ),
        calibration = list(
            fdr = 0.01,
            floors = c(2.0, 1.5, 1.2, 1.0, 0.5),
            min_mb = 5,
            merge_gap = 100
        )
    )
    sensitive <- balanced
    sensitive$windows_per_chromosome <- 100L
    sensitive$metric <- "combined"
    sensitive$calibration$floors <- c(1.5, 1.2, 1.0, 0.8, 0.4)
    cnv_robust <- balanced
    cnv_robust$windows_per_chromosome <- 100L
    cnv_robust$score <- "donor_only"
    cnv_robust$calibration$floors <- c(1.5, 1.2, 1.0, 0.8, 0.4)
    complex <- balanced
    complex$metric <- "combined"
    component_panel <- balanced

    list(
        balanced = list(
            label = "Balanced",
            description = paste(
                "Default candidate discovery for normal or high coverage;",
                "breadth, paired donor/recipient log-ratio and max donor",
                "aggregation."
            ),
            component_panel = FALSE,
            parameters = balanced
        ),
        sensitive = list(
            label = "Sensitive",
            description = paste(
                "Lower-resolution, combined breadth/count scoring with",
                "relaxed floors for low-coverage or subtle candidates."
            ),
            component_panel = FALSE,
            parameters = sensitive
        ),
        cnv_robust = list(
            label = "CNV robust",
            description = paste(
                "Donor-only scoring at low resolution for samples with",
                "suspected recipient copy-number changes."
            ),
            component_panel = FALSE,
            parameters = cnv_robust
        ),
        complex = list(
            label = "Complex",
            description = paste(
                "Combined breadth/count scoring for multi-component donors;",
                "pair with a component-isolation assembly when J/St/V",
                "resolution matters."
            ),
            component_panel = FALSE,
            parameters = complex
        ),
        component_panel = list(
            label = "Component panel",
            description = paste(
                "Balanced scoring plus the per-call component evidence",
                "panel. Use the component-isolation assembly to test one",
                "donor component at a time."
            ),
            component_panel = TRUE,
            parameters = component_panel
        )
    )
}

#' Resolve one named analysis profile
#'
#' @param name Profile name, or `NULL`/`"custom"` for the existing controls.
#' @return A profile contract, or `NULL` for a custom run.
#' @export
mamamia_analysis_profile <- function(name = NULL) {
    if (is.null(name) || length(name) == 0L) {
        return(NULL)
    }
    name <- as.character(name)[[1L]]
    if (is.na(name) || !nzchar(name) || identical(name, "custom")) {
        return(NULL)
    }
    profiles <- mamamia_analysis_profiles()
    profile <- profiles[[name]]
    if (is.null(profile)) {
        .mamamia_profile_error(
            "Unknown analysis profile '%s'; available profiles: %s",
            name, paste(names(profiles), collapse = ", ")
        )
    }
    profile
}

#' Choices for the dashboard analysis-profile selector
#' @noRd
mamamia_analysis_profile_choices <- function() {
    profiles <- mamamia_analysis_profiles()
    labels <- vapply(profiles, `[[`, character(1), "label")
    c(stats::setNames(labels, names(profiles)), "Custom" = "custom")
}
