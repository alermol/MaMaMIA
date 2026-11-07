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
medianPolish <- function(x, max_iter = 30L, tol = 1e-6, na.rm = TRUE) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("`x` must be a numeric matrix", call. = FALSE)
    }
    max_iter <- as.integer(max_iter)[[1L]]
    if (is.na(max_iter) || max_iter < 1L) {
        stop("`max_iter` must be a positive integer", call. = FALSE)
    }
    result <- median_polish_cpp(as.matrix(x), max_iter, tol, na.rm)
    dimnames(result$residuals) <- dimnames(x)
    names(result$row) <- rownames(x)
    names(result$column) <- colnames(x)
    result
}
