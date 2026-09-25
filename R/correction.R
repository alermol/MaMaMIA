#' @export
#' @importFrom rlang .data
RCA <- function(dcov, dgc, rcov, rgc, meta, presorted = TRUE) {
    dcov <- as_window_df(dcov, "cov", "dcov")
    dgc <- as_window_df(dgc, "gc", "dgc")
    rcov <- as_window_df(rcov, "cov", "rcov")
    rgc <- as_window_df(rgc, "gc", "rgc")

    if (!presorted) {
        message("Sorting...")
        dcov <- dplyr::arrange(dcov, dplyr::across(c("chr_id", "start", "stop")))
        rcov <- dplyr::arrange(rcov, dplyr::across(c("chr_id", "start", "stop")))
        dgc <- dplyr::arrange(dgc, dplyr::across(c("chr_id", "start", "stop")))
        rgc <- dplyr::arrange(rgc, dplyr::across(c("chr_id", "start", "stop")))
    }

    validate_windows_tiled(dcov, "dcov")
    validate_windows_tiled(dgc, "dgc")
    validate_windows_tiled(rcov, "rcov")
    validate_windows_tiled(rgc, "rgc")

    # Validate meta info
    meta <- as_meta_df(meta)

    don_chr_ids <- unique(dcov[["chr_id"]])
    rec_chr_ids <- unique(rcov[["chr_id"]])
    stopifnot(setequal(meta$chr_id, c(don_chr_ids, rec_chr_ids)), length(intersect(don_chr_ids, rec_chr_ids)) == 0)

    don_subgenomes <- unique(meta$subgenome[meta$chr_id %in% don_chr_ids])
    rec_subgenomes <- unique(meta$subgenome[meta$chr_id %in% rec_chr_ids])
    if (length(intersect(don_subgenomes, rec_subgenomes)) > 0L) {
        stop(
            "Donor and recipient subgenome labels must be disjoint; shared: ",
            paste(intersect(don_subgenomes, rec_subgenomes), collapse = ", "),
            call. = FALSE
        )
    }

    don_data <- dplyr::inner_join(dcov, dgc, by = c("chr_id", "start", "stop"))
    stopifnot(length(unique(c(
        nrow(dcov), nrow(dgc), nrow(don_data)
    ))) == 1)

    rec_data <- dplyr::inner_join(rcov, rgc, by = c("chr_id", "start", "stop"))
    stopifnot(length(unique(c(
        nrow(rcov), nrow(rgc), nrow(rec_data)
    ))) == 1)

    full_data <- dplyr::bind_rows(rec_data, don_data)
    validate_equal_n_windows(full_data$chr_id)

    full_data <- full_data |>
        dplyr::left_join(meta, by = "chr_id") |>
        dplyr::group_by(.data$chr_id) |>
        dplyr::mutate(iid = dplyr::row_number(), .after = "chr_id") |>
        dplyr::ungroup()

    structure(list(
        data = full_data,
        meta = list(
            meta = meta,
            don_chr_ids = don_chr_ids,
            rec_chr_ids = rec_chr_ids
        ),
        corrected = FALSE
    ), class = "RCA")
}

#' @export
#' @importFrom rlang .data
correctReadCounts <- function(RCA, cores = 1L, verbose = TRUE, cache = NULL) {
    stopifnot("Input is not RCA object" = inherits(RCA, "RCA"))

    cov_outlier <- 0.01
    gc_outlier <- 0.01

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
    cache_key <- if (is.null(cache)) NULL else fit_cache_key(ideal_data, cores, fit_call)
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
    RCA$outliers <- list(
        gc_lower_bound = gc_q[1],
        gc_upper_bound = gc_q[2],
        cov_upper_bound = cov_q
    )

    return(RCA)
}

#' @rdname reverseWindows
#' @export
reverseWindows.RCA <- function(x, chr_ids, ...) {
    stopifnot("Input is not RCA object" = inherits(x, "RCA"))

    chr_ids <- as.character(chr_ids)
    if (length(chr_ids) == 0L) {
        warning("No chromosomes requested; returning input unchanged",
            call. = FALSE
        )
        return(x)
    }
    unknown <- setdiff(chr_ids, unique(as.character(x$data$chr_id)))
    stopifnot(
        "Some chromosome IDs in `chr_ids` are not present in the object" =
            length(unknown) == 0L
    )

    value_cols <- intersect(
        names(x$data),
        c("cov", "gc", "cor.gc", "valid", "ideal")
    )
    x$data <- reverse_chr_value_cols(x$data, chr_ids, value_cols)

    return(x)
}


#' @exportS3Method base::plot
#' @importFrom rlang .data
plot.RCA <- function(x,
                     plot.type = c("orig_cov", "corr_cov"),
                     show_outliers = FALSE,
                     ...) {
    plot_type <- match.arg(plot.type)

    if (!x$corrected) {
        warning("Correction of GC-bias was not performed. See correctReadCounts().")
        warning("Fallback to original coverage with outliers")
        plot_type <- "orig_cov"
        show_outliers <- TRUE
    }

    if (show_outliers == TRUE) {
        plotting_df <- x$data
    } else {
        plotting_df <- dplyr::filter(x$data, .data$ideal == TRUE)
    }
    if (nrow(plotting_df) == 0L) {
        stop("No windows available to plot (empty data after filtering)",
            call. = FALSE
        )
    }
    plotting_df$color <- (match(plotting_df$chr_id, unique(plotting_df$chr_id)) - 1) %% 2
    plotting_df$x <- seq_len(nrow(plotting_df))

    breaks <- plotting_df |>
        dplyr::mutate(rn = .data$x) |>
        dplyr::group_by(.data$chr_name) |>
        dplyr::reframe(pos = min(.data$rn) + (diff(range(.data$rn)) / 2))

    plot <- ggplot2::ggplot(plotting_df, ggplot2::aes(x = .data$x, color = factor(.data$color))) +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.01, 0.01)),
            breaks = breaks$pos,
            labels = breaks$chr_name
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)), name = "Coverage") +
        ggplot2::scale_color_manual(values = c("grey60", "grey80")) +
        ggplot2::theme(
            aspect.ratio = 1 / 5,
            legend.position = "none",
            panel.grid = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank()
        )

    if (plot_type == "orig_cov") {
        plot + ggplot2::geom_point(
            size = 0.5,
            mapping = ggplot2::aes(y = .data$cov)
        )
    } else {
        plot + ggplot2::geom_point(
            size = 0.5,
            mapping = ggplot2::aes(y = .data$cor.gc)
        )
    }
}

#' @exportS3Method base::subset
subset.RCA <- function(x,
                       don_chromlist = NULL,
                       rec_chromlist = NULL,
                       don_subgenomes = NULL,
                       rec_subgenomes = NULL,
                       ...) {
    target_chrs <- resolve_chr_targets(
        x$meta$meta,
        don_chromlist,
        rec_chromlist,
        don_subgenomes,
        rec_subgenomes,
        don_chr_ids = x$meta$don_chr_ids,
        rec_chr_ids = x$meta$rec_chr_ids
    )

    x$data <- dplyr::filter(x$data, .data$chr_id %in% c(target_chrs[[1]], target_chrs[[2]]))
    x$meta$meta <- dplyr::filter(
        x$meta$meta,
        .data$chr_id %in% c(target_chrs[[1]], target_chrs[[2]])
    )
    x$meta$rec_chr_ids <- target_chrs[[2]]
    x$meta$don_chr_ids <- target_chrs[[1]]

    return(x)
}


#' @exportS3Method base::print
print.RCA <- function(x, ...) {
    n_chr <- length(unique(x$data$chr_id))
    cat(
        "RCA object: ",
        n_chr,
        " chromosomes (",
        length(x$meta$don_chr_ids),
        " donor / ",
        length(x$meta$rec_chr_ids),
        " recipient), ",
        nrow(x$data),
        " windows, GC-corrected: ",
        x$corrected,
        "\n",
        sep = ""
    )
    invisible(x)
}

#' @exportS3Method base::summary
summary.RCA <- function(object, ...) {
    subgenomes <- unique(object$meta$meta$subgenome)
    chromosomes_per_subgenome <- stats::setNames(vapply(subgenomes, \(d) {
        length(unique(object$meta$meta$chr_id[object$meta$meta$subgenome == d]))
    }, integer(1L)), subgenomes)

    out <- list(
        n_chromosomes = length(unique(object$data$chr_id)),
        n_subgenomes = length(subgenomes),
        chromosomes_per_subgenome = chromosomes_per_subgenome,
        corrected = object$corrected,
        outliers = if (isTRUE(object$corrected)) {
            object$outliers
        } else {
            NULL
        }
    )
    class(out) <- "summary.RCA"
    out
}

#' @exportS3Method base::print
print.summary.RCA <- function(x, ...) {
    cat(
        sprintf(
            "RCA Summary: %d chromosomes in %d subgenome(s): %s\n",
            x$n_chromosomes,
            x$n_subgenomes,
            paste(
                sprintf(
                    "%s: %d",
                    names(x$chromosomes_per_subgenome),
                    x$chromosomes_per_subgenome
                ),
                collapse = "; "
            )
        )
    )

    cat("GC-bias corrected:", if (isTRUE(x$corrected)) {
        "Yes"
    } else {
        "No"
    }, "\n")
    if (isTRUE(x$corrected) && !is.null(x$outliers)) {
        cat(
            sprintf(
                "Outlier thresholds | Coverage upper: %.2f, GC lower: %.3f, GC upper: %.3f\n",
                x$outliers$cov_upper_bound,
                x$outliers$gc_lower_bound,
                x$outliers$gc_upper_bound
            )
        )
    }
    invisible(x)
}
