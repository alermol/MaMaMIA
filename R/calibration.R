## Length-stratified empirical calibration of segment thresholds.

.calibration_bins <- function(length_bins_mb) {
    if (length(length_bins_mb) < 2L || any(diff(length_bins_mb) <= 0)) {
        .mamamia_error(
            "MaMaMIACalibrationError",
            "`length_bins_mb` must be an increasing numeric vector"
        )
    }
    labels <- paste0(
        format(length_bins_mb[-length(length_bins_mb)], trim = TRUE),
        "-",
        format(length_bins_mb[-1L], trim = TRUE)
    )
    labels[[length(labels)]] <- paste0(
        ">", format(length_bins_mb[[length(labels)]], trim = TRUE)
    )
    list(breaks = length_bins_mb, labels = labels)
}

.calibration_floor <- function(floors, labels, mode) {
    if (is.null(floors)) {
        floors <- if (identical(mode, "single")) {
            c(2.0, 1.5, 1.2, 1.0, 0.5)
        } else {
            rep(0.25, length(labels))
        }
    }
    if (!is.numeric(floors) || anyNA(floors) || any(floors < 0)) {
        .mamamia_error(
            "MaMaMIACalibrationError",
            "`floors` must be non-negative numbers"
        )
    }
    if (length(floors) == 1L) {
        floors <- rep(floors, length(labels))
    }
    if (length(floors) != length(labels)) {
        .mamamia_error(
            "MaMaMIACalibrationError",
            "`floors` must have one value per length bin (%d), got %d",
            length(labels), length(floors)
        )
    }
    stats::setNames(as.numeric(floors), labels)
}

.calibration_q99 <- function(values, fdr, mode, has_controls) {
    if (length(values) < 5L) {
        return(NA_real_)
    }
    if (identical(mode, "single") && !isTRUE(has_controls)) {
        ## A single sample can contain several real events, so the upper tail
        ## of all segment means is not a null. Split the segment means into a
        ## null and an event component and take the null quantile.
        fit <- .mask_em2(values)
        posterior_upper <- .mask_posterior(values, fit)
        null_values <- values[posterior_upper < 0.5]
        if (length(null_values) >= 5L) {
            return(as.numeric(stats::quantile(
                null_values, probs = 1 - fdr, na.rm = TRUE, names = FALSE
            )))
        }
        return(stats::qnorm(1 - fdr, fit$mu[[1L]], max(fit$sd[[1L]], 1e-6)))
    }
    ## Cohort calibration is robust when a few positives are present: trim
    ## the upper 5% before taking the empirical quantile.
    cutoff <- as.numeric(stats::quantile(
        values, 0.95, na.rm = TRUE, names = FALSE
    ))
    values <- values[values <= cutoff]
    if (length(values) == 0L) {
        return(NA_real_)
    }
    as.numeric(stats::quantile(
        values, probs = 1 - fdr, na.rm = TRUE, names = FALSE
    ))
}

#' Calibrate segment thresholds by length
#'
#' For each mode and segment-length bin, the threshold is the empirical
#' `1 - fdr` quantile of the segment means in the control set, raised to a
#' configurable floor. In a single-sample run the sample is its own control; in
#' a cohort every sample (or the explicitly supplied `controls`) contributes to
#' the null.
#'
#' @param calls A `Calls` data frame from [segmentScores()].
#' @param controls Optional sample names used as controls.
#' @param mode Optional mode filter (`"single"` or `"multi"`).
#' @param length_bins_mb Increasing length-bin breaks in Mb.
#' @param fdr Target false discovery rate.
#' @param floors Optional non-negative floor per bin.
#'
#' @return A list with `calibration` (`CalibrationTable`) and `calls` (the
#'   input plus `length_bin`, `threshold` and `pass`).
#' @export
calibrateSegments <- function(
    calls,
    controls = NULL,
    mode = NULL,
    length_bins_mb = c(5, 10, 25, 50, 100, Inf),
    fdr = 0.01,
    floors = NULL
) {
    if (!is.data.frame(calls)) {
        .mamamia_error(
            "MaMaMIACalibrationError", "`calls` must be a data frame"
        )
    }
    if (nrow(calls) == 0L) {
        return(list(
            calibration = data.frame(
                mode = character(0), length_bin = character(0),
                n_control = integer(0), q99 = numeric(0),
                threshold = numeric(0), fdr = numeric(0),
                version = character(0), stringsAsFactors = FALSE
            ),
            calls = calls
        ))
    }
    required <- c("sample", "mode", "width_bp", "score_mean")
    missing <- setdiff(required, names(calls))
    if (length(missing) > 0L) {
        .mamamia_error(
            "MaMaMIACalibrationError",
            "`calls` is missing column(s): %s", paste(missing, collapse = ", ")
        )
    }
    fdr <- as.numeric(fdr)[[1L]]
    if (is.na(fdr) || fdr <= 0 || fdr >= 1) {
        .mamamia_error(
            "MaMaMIACalibrationError", "`fdr` must be in (0, 1)"
        )
    }
    bins <- .calibration_bins(length_bins_mb)
    calls$length_bin <- cut(
        calls$width_bp / 1e6,
        breaks = bins$breaks,
        labels = bins$labels,
        include.lowest = TRUE,
        right = FALSE
    )
    modes <- unique(as.character(calls$mode))
    if (!is.null(mode)) {
        modes <- intersect(modes, as.character(mode))
    }
    rows <- list()
    for (current_mode in modes) {
        floor_values <- .calibration_floor(floors, bins$labels, current_mode)
        for (label in bins$labels) {
            pool <- calls[
                calls$mode == current_mode & calls$length_bin == label, ,
                drop = FALSE
            ]
            if (!is.null(controls)) {
                pool <- pool[pool$sample %in% controls, , drop = FALSE]
            }
            values <- pool$score_mean[is.finite(pool$score_mean)]
            q99 <- .calibration_q99(
                values, fdr = fdr, mode = current_mode,
                has_controls = !is.null(controls)
            )
            threshold <- max(floor_values[[label]], q99, na.rm = TRUE)
            if (!is.finite(threshold)) {
                threshold <- floor_values[[label]]
            }
            rows[[length(rows) + 1L]] <- data.frame(
                mode = current_mode,
                length_bin = label,
                n_control = length(values),
                q99 = q99,
                threshold = threshold,
                fdr = fdr,
                stringsAsFactors = FALSE
            )
        }
    }
    calibration <- do.call(rbind, rows)
    version <- hash_object(list(
        mode = calibration$mode,
        length_bin = calibration$length_bin,
        threshold = round(calibration$threshold, 8),
        fdr = fdr
    ))
    calibration$version <- version
    rownames(calibration) <- NULL

    key <- paste(calls$mode, calls$length_bin, sep = ":")
    lookup <- paste(calibration$mode, calibration$length_bin, sep = ":")
    position <- match(key, lookup)
    calls$threshold <- calibration$threshold[position]
    calls$pass <- is.finite(calls$score_mean) &
        is.finite(calls$threshold) & calls$score_mean > calls$threshold
    calls$threshold[is.na(calls$threshold)] <- Inf
    calls$pass[is.na(calls$pass)] <- FALSE
    list(calibration = calibration, calls = calls)
}
