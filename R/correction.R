## GC-bias correction of RCA coverage: fit a zero-inflated negative binomial
## model on the ideal windows, then scale every window by the model's expected
## response at the reference GC value of its subgenome.

#' @export
#' @importFrom rlang .data
correctReadCounts <- function(RCA, cores = 1L, verbose = TRUE, cache = NULL,
                              max_fit_windows = Inf) {
    stopifnot("Input is not RCA object" = inherits(RCA, "RCA"))

    cov_outlier <- 0.01
    gc_outlier <- 0.01

    max_fit_windows <- validate_max_fit_windows(max_fit_windows)
    if (identical(cache, FALSE)) {
        cache <- NULL
    }
    cores <- check_cores(cores)

    if (verbose) {
        message("Applying filter on data...")
    }
    RCA$data <- RCA$data |>
        dplyr::mutate(valid = (.data$cov >= 0) & (.data$gc > 0))

    n_valid <- sum(RCA$data$valid, na.rm = TRUE)
    if (n_valid < 2L) {
        stop("Fewer than 2 valid windows; cannot estimate GC/coverage quantiles",
            call. = FALSE
        )
    }

    gc_q <- stats::quantile(RCA$data$gc[RCA$data$valid], c(gc_outlier, 1 - gc_outlier))
    cov_q <- stats::quantile(RCA$data$cov[RCA$data$valid], 1 - cov_outlier)

    RCA$data <- dplyr::mutate(
        RCA$data,
        ideal = .data$valid &
            (.data$gc >= gc_q[1]) &
            (.data$gc <= gc_q[2]) &
            (.data$cov <= cov_q)
    )

    n_ideal <- sum(RCA$data$ideal, na.rm = TRUE)
    min_ideal <- 10L
    if (n_ideal < min_ideal) {
        stop(
            sprintf(
                "Fewer than %d ideal windows (%d); cannot fit GC-bias smooth",
                min_ideal,
                n_ideal
            ),
            call. = FALSE
        )
    }

    if (verbose) {
        message("Correcting for GC bias using ", cores, " cores...")
    }
    ideal_data <- dplyr::filter(RCA$data, .data$ideal)
    fit_call <- coverage_model_call(cores)
    sampling_seed <- if (is.finite(max_fit_windows) &&
        nrow(ideal_data) > max_fit_windows) {
        fit_sampling_seed(ideal_data, cores, fit_call, max_fit_windows)
    } else {
        NA_integer_
    }
    fit_sample <- select_fit_windows(ideal_data, max_fit_windows, sampling_seed)
    ideal_data <- fit_sample$data
    cache_key <- if (is.null(cache)) {
        NULL
    } else {
        fit_cache_key(
            ideal_data,
            cores,
            fit_call,
            max_fit_windows = max_fit_windows,
            selected_rows = ideal_data[c("chr_id", "iid")]
        )
    }
    fit <- read_cached_fit(cache, cache_key)
    if (is.null(fit)) {
        fit <- eval(fit_call)
        write_cached_fit(cache, cache_key, fit)
    } else if (verbose) {
        message("Reusing the cached model fit from ", cache)
    }

    gc_ref <- stats::median(RCA$data$gc[RCA$data$ideal], na.rm = TRUE)
    fitted_response <- predict_zinb_response(
        fit,
        gc = RCA$data$gc,
        subgenome = RCA$data$subgenome,
        gc_ref = gc_ref
    )
    RCA$data$cor.gc <- RCA$data$cov *
        (fitted_response$ref / (fitted_response$actual + 1e-8))

    RCA$data$ideal <- RCA$data$ideal &
        RCA$data$cor.gc < stats::quantile(RCA$data$cor.gc,
            probs = 1 - cov_outlier,
            na.rm = TRUE
        )

    n_ideal_post <- sum(RCA$data$ideal, na.rm = TRUE)
    if (n_ideal_post < 1L) {
        stop("No ideal windows remain after GC-corrected coverage filtering",
            call. = FALSE
        )
    }

    RCA$corrected <- TRUE
    RCA$fit <- fit
    RCA$fit_sampling <- fit_sample$sampling
    RCA$outliers <- list(
        gc_lower_bound = gc_q[1],
        gc_upper_bound = gc_q[2],
        cov_upper_bound = cov_q
    )

    return(RCA)
}

## Fit-window subsampling ----------------------------------------------------

#' Validate the fit-window cap
#' @noRd
validate_max_fit_windows <- function(max_fit_windows) {
    if (length(max_fit_windows) != 1L || is.na(max_fit_windows) ||
        !is.numeric(max_fit_windows) || max_fit_windows <= 0 ||
        (!is.infinite(max_fit_windows) &&
            max_fit_windows != floor(max_fit_windows))) {
        stop(
            "`max_fit_windows` must be a positive integer or Inf",
            call. = FALSE
        )
    }
    as.numeric(max_fit_windows)
}

#' Deterministic stratified subsample of the ideal windows
#'
#' Strata are the cross-classification of subgenome and GC decile. Every
#' non-empty stratum contributes at least ten windows (or all of its windows
#' when it is smaller). The remaining slots are allocated proportionally to the
#' remaining capacity. Sampling within a stratum uses a caller-supplied,
#' cache-key-derived seed, and the caller's RNG state is restored afterwards.
#'
#' @param ideal_data Ideal rows of `RCA$data`.
#' @param max_fit_windows Positive integer cap or `Inf`.
#' @param seed Integer seed; required when subsampling.
#'
#' @return A list with `data` (the rows used for fitting) and `sampling` (the
#'   `RCA$fit_sampling` metadata).
#' @noRd
select_fit_windows <- function(ideal_data, max_fit_windows,
                               seed = NA_integer_) {
    n_ideal <- nrow(ideal_data)
    max_fit_windows <- validate_max_fit_windows(max_fit_windows)

    gc_decile <- as.integer(dplyr::ntile(ideal_data$gc, 10L))
    stratum <- paste(
        as.character(ideal_data$subgenome), gc_decile,
        sep = "\r"
    )
    keys <- sort(unique(stratum))
    sizes <- as.integer(table(factor(stratum, levels = keys)))
    strata <- data.frame(
        subgenome = sub("\r.*$", "", keys),
        gc_decile = as.integer(sub("^.*\r", "", keys)),
        n_ideal = sizes,
        stringsAsFactors = FALSE
    )

    sampled <- is.finite(max_fit_windows) && n_ideal > max_fit_windows
    if (!sampled) {
        strata$n_fitted <- strata$n_ideal
        return(list(
            data = ideal_data,
            sampling = list(
                sampled = FALSE,
                max_fit_windows = max_fit_windows,
                n_ideal = n_ideal,
                n_fitted = n_ideal,
                seed = NA_integer_,
                strata = strata
            )
        ))
    }

    size <- as.integer(min(max_fit_windows, n_ideal))
    allocation <- allocate_fit_windows(sizes, size)
    strata$n_fitted <- allocation

    if (length(seed) != 1L || is.na(seed) || !is.numeric(seed)) {
        stop("A finite integer seed is required for fit subsampling",
            call. = FALSE
        )
    }
    selected <- with_preserved_seed(seed, {
        sort(unlist(lapply(seq_along(keys), function(j) {
            rows <- which(stratum == keys[[j]])
            if (allocation[[j]] >= length(rows)) {
                rows
            } else {
                sample(rows, allocation[[j]])
            }
        }), use.names = FALSE))
    })

    list(
        data = ideal_data[selected, , drop = FALSE],
        sampling = list(
            sampled = TRUE,
            max_fit_windows = max_fit_windows,
            n_ideal = n_ideal,
            n_fitted = length(selected),
            seed = as.integer(seed),
            strata = strata
        )
    )
}

#' Allocate the fit sample across strata
#'
#' Each stratum keeps at least `min(n, 10)` windows. The remaining slots are
#' allocated in proportion to the remaining capacity, with the largest
#' fractional remainders receiving the leftover slots. The result always sums to
#' `size` and never exceeds the stratum sizes.
#' @noRd
allocate_fit_windows <- function(sizes, size) {
    sizes <- as.integer(sizes)
    minimum <- pmin(sizes, 10L)
    required <- sum(minimum)
    if (size < required) {
        stop(
            sprintf(
                paste0(
                    "`max_fit_windows` (%s) is too small for %d ",
                    "subgenome x GC-decile strata; at least %s windows are ",
                    "required to fit at least 10 windows per stratum"
                ),
                format(size, scientific = FALSE),
                length(sizes),
                format(required, scientific = FALSE)
            ),
            call. = FALSE
        )
    }

    allocation <- minimum
    remaining <- size - sum(allocation)
    capacity <- sizes - allocation
    total_capacity <- sum(capacity)
    if (remaining <= 0L || total_capacity == 0L) {
        return(allocation)
    }

    quota <- remaining * capacity / total_capacity
    add <- floor(quota)
    allocation <- allocation + as.integer(add)
    remaining <- size - sum(allocation)
    if (remaining <= 0L) {
        return(allocation)
    }

    ## The fractional part is below one per stratum, so at most one extra
    ## window per open stratum is needed.
    fraction <- quota - floor(quota)
    open <- which(allocation < sizes)
    order_open <- open[order(-fraction[open], open)]
    for (j in order_open) {
        if (remaining <= 0L) {
            break
        }
        allocation[[j]] <- allocation[[j]] + 1L
        remaining <- remaining - 1L
    }
    allocation
}

#' Evaluate an expression under a temporary seed
#'
#' The caller's `.Random.seed` is restored on exit, so fit subsampling is
#' reproducible and does not perturb the surrounding analysis.
#' @noRd
with_preserved_seed <- function(seed, expr) {
    global <- globalenv()
    had_seed <- exists(".Random.seed", envir = global, inherits = FALSE)
    if (had_seed) {
        old_seed <- get(".Random.seed", envir = global, inherits = FALSE)
    }
    on.exit({
        if (had_seed) {
            assign(".Random.seed", old_seed, envir = global)
        } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
            rm(".Random.seed", envir = global)
        }
    }, add = TRUE)
    set.seed(seed)
    force(expr)
}

## The ZINB response is evaluated once for the distinct (GC, subgenome) pairs.
## The design matrices are assembled from the smooth stored in the fit instead
## of glmmTMB's `predict(debug = TRUE)`, which is two orders of magnitude
## slower; if that structure is unavailable, `predict()` is used as a fallback.

#' Fixed/random design columns of a stored mgcv smooth for new values.
#' @noRd
smooth_design <- function(sm, gc) {
    ev <- eigen(sm$S[[1L]], symmetric = TRUE)
    if (ev$vectors[1L, 1L] < 0) {
        ev$vectors <- -ev$vectors
    }
    p_rank <- min(sm$rank, ncol(sm$X))
    null_rank <- sm$df - sm$rank
    scaling <- 1 / sqrt(c(ev$values[seq_len(p_rank)], rep(1, null_rank)))
    basis <- mgcv::PredictMat(sm, data = data.frame(gc = gc)) %*%
        t(t(ev$vectors) * scaling)
    fixed <- if (p_rank < sm$df) {
        basis[, (p_rank + 1L):sm$df, drop = FALSE]
    } else {
        matrix(0, nrow(basis), 0L)
    }
    list(
        fixed = fixed,
        random = basis[, seq_len(p_rank), drop = FALSE],
        label = sm$label
    )
}

#' Fixed/random design matrices for the fitted coverage model.
#' @noRd
glmmtmb_design <- function(fit, newdata) {
    subgenomes <- levels(factor(fit$frame$subgenome))
    if (anyNA(factor(newdata$subgenome, levels = subgenomes))) {
        stop("unknown subgenome level", call. = FALSE)
    }
    dummies <- unname(stats::model.matrix(
        ~ factor(newdata$subgenome, levels = subgenomes)
    )[, -1L, drop = FALSE])

    build <- function(part, coef_names) {
        info <- fit$modelInfo$reTrms[[part]]$smooth_info
        if (length(info) != 1L) {
            stop("expected exactly one smooth in `", part, "`", call. = FALSE)
        }
        sd <- smooth_design(info[[1L]]$sm, newdata$gc)
        expected <- c(
            "(Intercept)",
            paste0("subgenome", subgenomes[-1L]),
            paste0(sd$label, seq_len(ncol(sd$fixed)))
        )
        if (!identical(expected, as.character(coef_names))) {
            stop("unexpected fixed-effect structure in `", part, "`", call. = FALSE)
        }
        list(X = cbind(1, dummies, sd$fixed), Z = sd$random)
    }

    list(
        cond = build("cond", names(glmmTMB::fixef(fit)$cond)),
        zi = build("zi", names(glmmTMB::fixef(fit)$zi))
    )
}

#' Design matrices from glmmTMB's prediction machinery.
#' @noRd
design_from_tmb <- function(fit, newdata) {
    tmb <- stats::predict(fit, newdata = newdata, debug = TRUE)$data.tmb
    if (is.null(dim(tmb$X))) {
        stop("Fitted model exposes no dense conditional design matrix", call. = FALSE)
    }
    n_aug <- nrow(tmb$X)
    n_new <- nrow(newdata)
    if (n_aug != nrow(fit$frame) + n_new) {
        stop("Unexpected augmented design size; cannot align predictions", call. = FALSE)
    }
    idx <- (n_aug - n_new + 1L):n_aug
    trim <- function(M) {
        as.matrix(M)[idx, , drop = FALSE]
    }
    list(
        cond = list(X = trim(tmb$X), Z = trim(tmb$Z)),
        zi = list(X = trim(tmb$Xzi), Z = trim(tmb$Zzi))
    )
}

#' Mean response of the fitted ZINB coverage model for new data.
#' @noRd
eta_zinb_response <- function(fit, newdata) {
    pars <- fit$fit$parfull
    get_par <- function(name) {
        out <- pars[names(pars) == name]
        if (length(out) == 0L) {
            stop("Fitted model has no `", name, "` parameter", call. = FALSE)
        }
        out
    }

    design <- tryCatch(
        glmmtmb_design(fit, newdata),
        error = function(e) NULL
    )
    if (is.null(design)) {
        design <- design_from_tmb(fit, newdata)
    }

    eta_cond <- as.numeric(design$cond$X %*% get_par("beta")) +
        as.numeric(design$cond$Z %*% get_par("b"))
    eta_zi <- as.numeric(design$zi$X %*% get_par("betazi")) +
        as.numeric(design$zi$Z %*% get_par("bzi"))

    (1 - stats::plogis(eta_zi)) * exp(eta_cond)
}

#' Expected ZINB coverage at observed GC and at reference GC.
#' @noRd
predict_zinb_response <- function(fit, gc, subgenome, gc_ref) {
    n <- length(gc)

    fallback <- function() {
        nd <- data.frame(gc = gc, subgenome = subgenome)
        list(
            actual = as.numeric(stats::predict(fit, newdata = nd, type = "response")),
            ref = as.numeric(stats::predict(fit,
                newdata = transform(nd, gc = gc_ref),
                type = "response"
            ))
        )
    }

    ok <- !is.na(gc) & !is.na(subgenome)
    if (!any(ok)) {
        na <- rep(NA_real_, n)
        return(list(actual = na, ref = na))
    }

    key <- function(g, s) paste(sprintf("%.17g", g), s, sep = "\r")
    rows <- data.frame(gc = gc[ok], subgenome = as.character(subgenome[ok]))
    row_keys <- key(rows$gc, rows$subgenome)
    unique_rows <- rows[!duplicated(row_keys), , drop = FALSE]
    unique_keys <- key(unique_rows$gc, unique_rows$subgenome)
    subgenomes <- unique(unique_rows$subgenome)
    ref_rows <- data.frame(gc = gc_ref, subgenome = subgenomes)

    pred <- tryCatch(
        eta_zinb_response(fit, rbind(unique_rows, ref_rows)),
        error = function(e) NULL
    )
    if (is.null(pred) || length(pred) != nrow(unique_rows) + nrow(ref_rows)) {
        return(fallback())
    }

    n_unique <- nrow(unique_rows)
    actual <- rep(NA_real_, n)
    actual[ok] <- pred[match(row_keys, unique_keys)]
    ref <- rep(NA_real_, n)
    ref[ok] <- pred[n_unique + match(rows$subgenome, subgenomes)]
    list(actual = actual, ref = ref)
}
