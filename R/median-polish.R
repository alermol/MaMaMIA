## Robust two-way median polish used by the multi-sample score.
##
## The model is x[window, sample] = mu + a[sample] + b[window] + residual.
## Alternating row and column medians remove the shared window effect and the
## per-sample library/ploidy offset. The residual is the introgression score.
## Missing values are allowed; a row that is entirely missing stays missing.

#' Median polish of a matrix
#'
#' @param x A numeric matrix with one row per window and one column per sample.
#' @param max_iter Maximum number of row/column sweeps.
#' @param tol Convergence tolerance on the sum of absolute effects.
#' @param na.rm Ignore missing values when computing medians.
#'
#' @return A list with `residuals` (the same shape as `x`), `row` (window
#'   effects), `column` (sample effects), `mu`, `iterations` and `converged`.
#' @export
medianPolish <- function(x, max_iter = 30L, tol = 1e-6, na.rm = TRUE) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("`x` must be a numeric matrix", call. = FALSE)
    }
    max_iter <- as.integer(max_iter)[[1L]]
    if (is.na(max_iter) || max_iter < 1L) {
        stop("`max_iter` must be a positive integer", call. = FALSE)
    }
    residuals <- matrix(as.numeric(x), nrow = nrow(x), ncol = ncol(x),
        dimnames = dimnames(x))
    row_effects <- stats::setNames(numeric(nrow(residuals)), rownames(residuals))
    column_effects <- stats::setNames(
        numeric(ncol(residuals)), colnames(residuals)
    )
    mu <- 0
    converged <- FALSE
    iteration <- 0L
    for (iteration in seq_len(max_iter)) {
        old <- sum(abs(row_effects), na.rm = TRUE) +
            sum(abs(column_effects), na.rm = TRUE) + abs(mu)
        row_medians <- apply(residuals, 1L, stats::median, na.rm = na.rm)
        row_medians[!is.finite(row_medians)] <- 0
        residuals <- residuals - row_medians
        row_effects <- row_effects + row_medians
        column_medians <- apply(residuals, 2L, stats::median, na.rm = na.rm)
        column_medians[!is.finite(column_medians)] <- 0
        residuals <- sweep(residuals, 2L, column_medians, "-")
        column_effects <- column_effects + column_medians
        mu <- mu + stats::median(
            c(row_medians, column_medians), na.rm = TRUE
        )
        current <- sum(abs(row_effects), na.rm = TRUE) +
            sum(abs(column_effects), na.rm = TRUE) + abs(mu)
        if (is.finite(current) && abs(old - current) < tol) {
            converged <- TRUE
            break
        }
    }
    list(
        residuals = residuals,
        row = row_effects,
        column = column_effects,
        mu = mu,
        iterations = iteration,
        converged = converged
    )
}
