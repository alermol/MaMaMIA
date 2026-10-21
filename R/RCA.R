## The Read Count Array (RCA): windowed donor and recipient coverage and GC
## values joined into one table, plus the methods that inspect or subset it.
## See correction.R for the GC-bias correction of an RCA object.

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

## Internal: resolve the donor and recipient chromosomes selected for
## `subset.RCA()` from chromosome IDs and subgenome labels.
#' @noRd
resolve_chr_targets <- function(meta,
                                don_chromlist,
                                rec_chromlist,
                                don_subgenomes,
                                rec_subgenomes,
                                don_chr_ids,
                                rec_chr_ids) {
    stopifnot(
        "Require at least one of don_chromlist, rec_chromlist, don_subgenomes, rec_subgenomes" = !all(vapply(
            list(
                don_chromlist,
                rec_chromlist,
                don_subgenomes,
                rec_subgenomes
            ),
            is.null,
            logical(1)
        ))
    )

    if (!is.null(don_chromlist) && length(don_chromlist) == 0L) {
        stop("`don_chromlist` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(rec_chromlist) && length(rec_chromlist) == 0L) {
        stop("`rec_chromlist` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(don_subgenomes) && length(don_subgenomes) == 0L) {
        stop("`don_subgenomes` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(rec_subgenomes) && length(rec_subgenomes) == 0L) {
        stop("`rec_subgenomes` must be NULL or non-empty", call. = FALSE)
    }

    if (!is.null(don_subgenomes)) {
        stopifnot(
            "Unknown donor subgenome label(s)" =
                all(don_subgenomes %in% meta$subgenome)
        )
    }
    if (!is.null(rec_subgenomes)) {
        stopifnot(
            "Unknown recipient subgenome label(s)" =
                all(rec_subgenomes %in% meta$subgenome)
        )
    }
    if (!is.null(don_chromlist)) {
        stopifnot(
            "Unknown donor chromosome ID(s)" =
                all(don_chromlist %in% don_chr_ids)
        )
    }
    if (!is.null(rec_chromlist)) {
        stopifnot(
            "Unknown recipient chromosome ID(s)" =
                all(rec_chromlist %in% rec_chr_ids)
        )
    }

    resolve_one <- function(chromlist, subgenomes, current_ids) {
        current_ids <- as.character(current_ids)
        if (is.null(chromlist) && is.null(subgenomes)) {
            return(current_ids)
        }
        ids <- character(0)
        if (!is.null(subgenomes)) {
            ids <- c(ids, intersect(as.character(meta$chr_id[meta$subgenome %in% subgenomes]), current_ids))
        }
        if (!is.null(chromlist)) {
            ids <- c(ids, as.character(chromlist))
        }
        unique(ids)
    }

    target_don <- resolve_one(don_chromlist, don_subgenomes, don_chr_ids)
    target_rec <- resolve_one(rec_chromlist, rec_subgenomes, rec_chr_ids)
    if (length(target_don) == 0L || length(target_rec) == 0L) {
        stop("Subset would leave donor or recipient with no chromosomes",
            call. = FALSE
        )
    }

    list(target_don, target_rec)
}

#' @exportS3Method base::print
print.RCA <- function(x, ...) {
    n_chr <- length(unique(x$data$chr_id))
    normalization <- if (is.null(x$normalization)) {
        ""
    } else {
        paste0(", normalized: ", x$normalization$mode)
    }
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
        normalization,
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
        normalization = object$normalization,
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
    if (!is.null(x$normalization)) {
        cat(
            sprintf(
                "Normalization: %s (by side: %s, pseudocount: %.3g)\n",
                x$normalization$mode,
                if (isTRUE(x$normalization$by_side)) "yes" else "no",
                x$normalization$pseudocount
            )
        )
    }
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
