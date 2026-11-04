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
#'   `short_description`, a `description` and a `parameters` list suitable for
#'   [analyzeBamInput()].
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

    list(
        balanced = list(
            label = "Balanced",
            short_description = paste(
                "Default candidate discovery for normal or high coverage."
            ),
            description = paste(
                "Default candidate discovery for normal or high coverage;",
                "breadth, paired donor/recipient log-ratio and max donor",
                "aggregation."
            ),
            parameters = balanced
        ),
        sensitive = list(
            label = "Sensitive",
            short_description = paste(
                "Relaxed, low-resolution discovery for low-coverage samples."
            ),
            description = paste(
                "Lower-resolution, combined breadth/count scoring with",
                "relaxed floors for low-coverage or subtle candidates."
            ),
            parameters = sensitive
        ),
        cnv_robust = list(
            label = "CNV robust",
            short_description = paste(
                "Donor-only scoring when recipient copy number varies."
            ),
            description = paste(
                "Donor-only scoring at low resolution for samples with",
                "suspected recipient copy-number changes."
            ),
            parameters = cnv_robust
        ),
        complex = list(
            label = "Complex",
            short_description = paste(
                "Combined scoring for multi-component donors."
            ),
            description = paste(
                "Combined breadth/count scoring for multi-component donors;",
                "pair with a component-isolation assembly when J/St/V",
                "resolution matters."
            ),
            parameters = complex
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

#' Flatten one profile's parameters for forms, provenance and bundles
#'
#' Keeps every profile parameter in one scalar-or-vector list so the dashboard
#' can edit each value and the bundle can record exactly what ran.
#'
#' @param parameters A `parameters` list from [mamamia_analysis_profiles()].
#' @return A named list with scalar parameters and a five-value
#'   `calibration_floors` vector.
#' @noRd
mamamia_profile_flat_parameters <- function(parameters) {
    if (!is.list(parameters)) {
        .mamamia_profile_error("Profile parameters must be a list")
    }
    mask <- parameters$mask
    calibration <- parameters$calibration
    if (!is.list(mask)) {
        mask <- list()
    }
    if (!is.list(calibration)) {
        calibration <- list()
    }
    list(
        windows_per_chromosome = as.integer(
            parameters$windows_per_chromosome
        ),
        metric = as.character(parameters$metric),
        score = as.character(parameters$score),
        donor_aggregation = as.character(parameters$donor_aggregation),
        mapq = as.integer(parameters$mapq),
        mask_physical_fallback_mb = as.numeric(
            mask$physical_fallback_mb
        ),
        mask_recurrence_threshold = as.numeric(
            mask$recurrence_threshold
        ),
        calibration_fdr = as.numeric(calibration$fdr),
        calibration_floors = as.numeric(calibration$floors),
        calibration_min_mb = as.numeric(calibration$min_mb),
        calibration_merge_gap = as.integer(calibration$merge_gap)
    )
}

#' Choices for the dashboard analysis-profile selector
#' @noRd
mamamia_analysis_profile_choices <- function() {
    profiles <- mamamia_analysis_profiles()
    labels <- vapply(profiles, `[[`, character(1), "label")
    c(stats::setNames(names(profiles), labels), "Custom" = "custom")
}
