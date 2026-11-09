## Calibration and segmentation parameters as one object.
##
## The registry declares length bins, an FDR and single-sample floors; runtime
## calibration supplies floors, bins, an FDR and segmentation controls. Both
## used to be spliced together by name in several places. This module resolves
## the registry defaults and the runtime overrides once, so the pipeline, the
## dashboard and the result provenance all read the same values.

.calibration_defaults <- function(defaults) {
    if (is.null(defaults)) {
        return(.config_defaults())
    }
    if (!is.list(defaults)) {
        .mamamia_error(
            "MaMaMIACalibrationError", "`defaults` must be a mapping"
        )
    }
    if (!is.null(defaults$defaults) && is.list(defaults$defaults)) {
        defaults <- defaults$defaults
    }
    defaults
}

.calibration_number <- function(x, what, min = 0, max = Inf) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) || !is.finite(x) ||
        x < min || x > max) {
        .mamamia_error(
            "MaMaMIACalibrationError",
            "`%s` must be a single number in [%s, %s]",
            what,
            format(min, scientific = FALSE),
            format(max, scientific = FALSE)
        )
    }
    as.numeric(x)
}

## Stretch a registered floor vector to a different number of bins. The
## bundled five-bin curve is interpolated; other lengths are stretched the same
## way so a runtime bin override cannot leave floors and labels out of step.
.calibration_stretch_floors <- function(floors, n) {
    floors <- as.numeric(floors)
    if (length(floors) == n) {
        return(floors)
    }
    if (length(floors) == 1L) {
        return(rep(floors, n))
    }
    if (length(floors) == 0L || anyNA(floors) || any(!is.finite(floors))) {
        .mamamia_error(
            "MaMaMIACalibrationError",
            "Calibration floors must be finite numbers"
        )
    }
    positions <- if (length(floors) == 1L) {
        0
    } else {
        seq(0, 1, length.out = length(floors))
    }
    stats::approx(positions, floors, xout = seq(0, 1, length.out = n),
        rule = 2
    )$y
}

#' Resolve calibration and segmentation parameters
#'
#' Merges the registry calibration defaults with runtime overrides and returns
#' one object. `bins` holds the `breaks` and `labels` used by calibration,
#' `floors` is the selected per-bin floor (single-sample floors in single mode,
#' the cohort floor in multi mode), and `segmentation` holds every CBS and
#' merge argument.
#'
#' @param defaults Registry defaults (`config$defaults`), a full configuration,
#'   or `NULL` for the built-in defaults.
#' @param overrides Optional runtime calibration/segmentation list. Recognized
#'   names are `length_bins_mb`, `fdr`, `floors`, `single_sample_floors`,
#'   `controls`, `alpha`, `undo_SD`, `min_width`, `merge_gap`, `merge_floor`,
#'   `merge_min` and `min_mb`.
#' @param mode `"single"` or `"multi"`; selects the default floor vector.
#'
#' @return A `MaMaMIACalibrationParameters` list with `bins`, `length_bins_mb`,
#'   `labels`, `n_bins`, `floors`, `single_sample_floors`, `controls`, `fdr`,
#'   `mode` and `segmentation`.
resolveCalibrationParameters <- function(defaults = .config_defaults(),
                                         overrides = NULL,
                                         mode = c("single", "multi")) {
    mode <- match.arg(mode)
    defaults <- .calibration_defaults(defaults)
    calibration <- defaults$calibration
    if (is.null(calibration)) {
        calibration <- list()
    }
    if (!is.list(calibration)) {
        .mamamia_error(
            "MaMaMIACalibrationError",
            "`defaults$calibration` must be a mapping"
        )
    }
    if (!is.null(overrides) && !is.list(overrides)) {
        .mamamia_error(
            "MaMaMIACalibrationError", "`overrides` must be NULL or a list"
        )
    }
    if (is.null(overrides)) {
        overrides <- list()
    }

    length_bins_mb <- overrides$length_bins_mb %||%
        calibration$length_bins_mb %||% c(5, 10, 25, 50, 100, Inf)
    bins <- .calibration_bins(as.numeric(unlist(length_bins_mb)))
    n_bins <- length(bins$labels)

    fdr <- .calibration_number(
        overrides$fdr %||% calibration$fdr %||% 0.01,
        "fdr", min = .Machine$double.eps, max = 1 - .Machine$double.eps
    )

    controls <- overrides$controls %||% calibration$controls
    if (!is.null(controls)) {
        controls <- as.character(unlist(controls, use.names = FALSE))
    }

    single_sample_floors <- overrides$single_sample_floors %||%
        calibration$single_sample_floors
    if (is.null(single_sample_floors) || length(single_sample_floors) == 0L) {
        single_sample_floors <- c(2.0, 1.5, 1.2, 1.0, 0.5)
    }
    single_sample_floors <- .calibration_stretch_floors(
        single_sample_floors, n_bins
    )

    floors <- overrides$floors
    if (is.null(floors)) {
        floors <- if (identical(mode, "single")) {
            single_sample_floors
        } else {
            rep(0.25, n_bins)
        }
    } else {
        floors <- as.numeric(unlist(floors, use.names = FALSE))
        if (length(floors) != 1L && length(floors) != n_bins) {
            .mamamia_error(
                "MaMaMIACalibrationError",
                "`floors` must have one value per length bin (%d), got %d",
                n_bins, length(floors)
            )
        }
    }
    floors <- .calibration_floor(floors, bins$labels, mode)

    segmentation <- list(
        alpha = .calibration_number(
            overrides$alpha %||% 0.001, "alpha", min = 0, max = 1
        ),
        undo_SD = .calibration_number(
            overrides$undo_SD %||% 1, "undo_SD", min = 0
        ),
        min_width = as.integer(.calibration_number(
            overrides$min_width %||% 5L, "min_width", min = 1
        )),
        merge_gap = as.integer(.calibration_number(
            overrides$merge_gap %||% 100L, "merge_gap", min = 0
        )),
        merge_floor = .calibration_number(
            overrides$merge_floor %||% -0.5, "merge_floor",
            min = -Inf, max = Inf
        ),
        merge_min = .calibration_number(
            overrides$merge_min %||% 0.5, "merge_min",
            min = -Inf, max = Inf
        ),
        min_mb = .calibration_number(
            overrides$min_mb %||% 5, "min_mb", min = 0
        )
    )

    structure(
        list(
            bins = bins,
            length_bins_mb = bins$breaks,
            labels = bins$labels,
            n_bins = n_bins,
            floors = floors,
            single_sample_floors = stats::setNames(
                single_sample_floors, bins$labels
            ),
            controls = controls,
            fdr = fdr,
            mode = mode,
            segmentation = segmentation
        ),
        class = c("MaMaMIACalibrationParameters", "list")
    )
}
