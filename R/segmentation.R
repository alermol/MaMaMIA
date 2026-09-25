#' @export
#' @importFrom rlang .data
segments <- function(RCA,
                     target_pairs,
                     min_segment_median = 0L,
                     min_segment_mean = 0L,
                     threshold_type = c("median", "mean"),
                     alpha = 0.001,
                     min_width = 2,
                     undo_SD = 1,
                     seed = 1L,
                     verbose = TRUE) {
    stopifnot("Input is not RCA object" = inherits(RCA, "RCA"))
    stopifnot("Coverage must be corrected before segmentation" = RCA$corrected)

    thr_type <- match.arg(threshold_type)
    threshold <- switch(thr_type,
        "mean" = min_segment_mean,
        "median" = min_segment_median
    )
    colname <- switch(thr_type,
        "mean" = "seg.mean",
        "median" = "seg.median"
    )

    param <- validate_segmentation_params(alpha, min_width, undo_SD, seed)

    target_pairs_df <- validate_target_pairs(target_pairs, RCA$meta$don_chr_ids, RCA$meta$rec_chr_ids)

    by_chr_split <- split(RCA$data, RCA$data$chr_id)

    if (verbose) {
        message("Performing segmentation...")
    }

    prepared_data <- target_pairs_df |>
        dplyr::mutate(
            data = purrr::map2(.data$chr_id.don, .data$chr_id.rec, \(x, y) {
                don_data <- by_chr_split[[as.character(x)]]
                rec_data <- by_chr_split[[as.character(y)]]
                dplyr::inner_join(don_data,
                    rec_data,
                    by = "iid",
                    suffix = c(".don", ".rec")
                ) |>
                    dplyr::filter(.data$ideal.don &
                        .data$ideal.rec) |>
                    dplyr::mutate(diff = .data$cor.gc.don - .data$cor.gc.rec)
            }),
            segments = purrr::pmap(list(
                .data$data, .data$chr_id.don, .data$chr_id.rec
            ), segment_pair_data,
            meta = RCA$meta$meta,
            param = param,
            colname = colname,
            threshold = threshold,
            min_width = min_width
            ),
            .keep = "unused"
        )

    keep_don <- unique(target_pairs_df$chr_id.don)
    keep_rec <- unique(target_pairs_df$chr_id.rec)
    keep_chrs <- c(keep_don, keep_rec)

    structure(
        list(
            data = dplyr::filter(RCA$data, .data$chr_id %in% keep_chrs),
            out = dplyr::bind_rows(prepared_data$data),
            segments = dplyr::bind_rows(
                empty_segments(),
                purrr::compact(prepared_data$segments)
            ),
            param = param,
            meta = list(
                meta = dplyr::filter(RCA$meta$meta, .data$chr_id %in% keep_chrs),
                don_chr_ids = keep_don,
                rec_chr_ids = keep_rec,
                target_pairs = target_pairs_df,
                metrics = colname,
                metrics_threshold = threshold
            )
        ),
        class = "ISA"
    )
}

#' @rdname reverseWindows
#' @export
#' @importFrom rlang .data
reverseWindows.ISA <- function(x, chr_ids, ...) {
    stopifnot("Input is not ISA object" = inherits(x, "ISA"))

    chr_ids <- as.character(chr_ids)
    if (length(chr_ids) == 0L) {
        warning("No chromosomes requested; returning input unchanged",
            call. = FALSE
        )
        return(x)
    }
    known <- c(
        as.character(x$meta$don_chr_ids),
        as.character(x$meta$rec_chr_ids)
    )
    unknown <- setdiff(chr_ids, known)
    stopifnot(
        "Some chromosome IDs in `chr_ids` are not present in the object" =
            length(unknown) == 0L
    )

    don_targets <- intersect(chr_ids, as.character(x$meta$don_chr_ids))
    rec_targets <- intersect(chr_ids, as.character(x$meta$rec_chr_ids))

    don_cols <- intersect(
        names(x$out),
        paste0(c("cov", "gc", "cor.gc", "valid", "ideal"), ".don")
    )
    rec_cols <- intersect(
        names(x$out),
        paste0(c("cov", "gc", "cor.gc", "valid", "ideal"), ".rec")
    )

    # Reverse the donor-side value vectors for the requested donor chromosomes,
    # within each (donor, recipient) pair so row order is preserved.
    if (!is.null(x$out) && nrow(x$out) > 0L &&
        length(don_targets) > 0L && length(don_cols) > 0L) {
        x$out <- reverse_paired_side(
            x$out, chr_id_side = "chr_id.don",
            target_ids = don_targets, value_cols = don_cols
        )
    }

    # Reverse the recipient-side value vectors.
    if (!is.null(x$out) && nrow(x$out) > 0L &&
        length(rec_targets) > 0L && length(rec_cols) > 0L) {
        x$out <- reverse_paired_side(
            x$out, chr_id_side = "chr_id.rec",
            target_ids = rec_targets, value_cols = rec_cols
        )
    }

    # Recompute the derived difference so it stays consistent with the reversed
    # cor.gc vectors.
    if (!is.null(x$out) && nrow(x$out) > 0L) {
        x$out <- dplyr::mutate(
            x$out,
            diff = .data$cor.gc.don - .data$cor.gc.rec
        )
    }

    # Re-run CBS for every pair containing a reversed chromosome, since those
    # pairs' difference profiles have changed orientation. Unaffected pairs'
    # segments are kept as-is.
    affected_pairs <- dplyr::filter(
        x$meta$target_pairs,
        .data$chr_id.don %in% chr_ids | .data$chr_id.rec %in% chr_ids
    )
    if (!is.null(x$out) && nrow(x$out) > 0L && nrow(affected_pairs) > 0L) {
        affected_segments <- purrr::pmap(
            list(affected_pairs$chr_id.don, affected_pairs$chr_id.rec),
            \(don, rec) {
                pair_rows <- dplyr::filter(
                    x$out,
                    .data$chr_id.don == !!don & .data$chr_id.rec == !!rec
                )
                segment_pair_data(
                    pair_rows,
                    don = don,
                    rec = rec,
                    meta = x$meta$meta,
                    param = x$param,
                    colname = x$meta$metrics,
                    threshold = x$meta$metrics_threshold,
                    min_width = x$param$min.width
                )
            }
        )
        x$segments <- dplyr::bind_rows(
            dplyr::anti_join(
                x$segments,
                affected_pairs,
                by = c("chr_id.don", "chr_id.rec")
            ),
            purrr::compact(affected_segments)
        )
    }

    # Keep the underlying per-window table aligned as well.
    if (!is.null(x$data) && nrow(x$data) > 0L) {
        value_cols <- intersect(
            names(x$data),
            c("cov", "gc", "cor.gc", "valid", "ideal")
        )
        x$data <- reverse_chr_value_cols(x$data, chr_ids, value_cols)
    }

    return(x)
}


#' @exportS3Method base::subset
#' @importFrom rlang .data
subset.ISA <- function(x,
                       min_width = NULL,
                       min_segment_median = NULL,
                       min_segment_mean = NULL,
                       threshold_type = NULL,
                       target_pairs = NULL,
                       ...) {
    if (!is.null(threshold_type)) {
        stopifnot(
            "Threshold types must be 'mean' or 'median'" =
                !is.na(match(
                    threshold_type, c("mean", "median")
                ))
        )
        new_metrics <- ifelse(threshold_type == "mean", "seg.mean", "seg.median")
        if (!identical(new_metrics, x$meta$metrics)) {
            if (identical(threshold_type, "mean") &&
                is.null(min_segment_mean)) {
                stop(paste0(
                    c(
                        "When switching threshold_type to 'mean', ",
                        "supply min_segment_mean"
                    ),
                    collapse = ""
                ), call. = FALSE)
            }
            if (identical(threshold_type, "median") &&
                is.null(min_segment_median)) {
                stop(
                    "When switching threshold_type to 'median', supply min_segment_median",
                    call. = FALSE
                )
            }
        }
        x$meta$metrics <- new_metrics
    }

    stopifnot(
        "Only median or mean threshold could be used" = is.null(min_segment_mean) ||
            is.null(min_segment_median)
    )
    if (!is.null(min_segment_mean)) {
        stopifnot("Metrics in current ISA is not mean" = x$meta$metrics == "seg.mean")
        x$meta$metrics_threshold <- min_segment_mean
    }

    if (!is.null(min_segment_median)) {
        stopifnot("Metrics in current ISA is not median" = x$meta$metrics == "seg.median")
        x$meta$metrics_threshold <- min_segment_median
    }

    x$segments$putative_introgression <- x$segments[[x$meta$metrics]] >= x$meta$metrics_threshold

    if (!is.null(target_pairs)) {
        target_pairs_df <- validate_target_pairs(target_pairs, x$meta$don_chr_ids, x$meta$rec_chr_ids)
        x$meta$target_pairs <- target_pairs_df
        x$meta$don_chr_ids <- unique(target_pairs_df$chr_id.don)
        x$meta$rec_chr_ids <- unique(target_pairs_df$chr_id.rec)

        keep_chrs <- c(x$meta$don_chr_ids, x$meta$rec_chr_ids)
        x$meta$meta <- dplyr::filter(x$meta$meta, .data$chr_id %in% keep_chrs)
        x$data <- dplyr::filter(x$data, .data$chr_id %in% keep_chrs)

        x$out <- dplyr::semi_join(x$out,
            x$meta$target_pairs,
            by = c("chr_id.don", "chr_id.rec")
        )

        x$segments <- dplyr::semi_join(x$segments,
            x$meta$target_pairs,
            by = c("chr_id.don", "chr_id.rec")
        )
    }

    if (!is.null(min_width)) {
        stopifnot("min.width must be >=2" = min_width >= 2)
        x$segments <- dplyr::filter(
            x$segments,
            (.data$loc.end - .data$loc.start + 1) >= min_width
        )
        x$param$min.width <- min_width
    }

    return(x)
}


#' @exportS3Method base::plot
#' @importFrom rlang .data
plot.ISA <- function(x,
                     plot.type = c("mirror", "diff"),
                     all_segments = FALSE,
                     ...) {
    plot_type <- match.arg(plot.type)

    if (is.null(x$out) || nrow(x$out) == 0L) {
        stop("No pairwise windows available to plot", call. = FALSE)
    }

    # Add throughout number of windows
    x$out$iid.e2e <- seq_along(x$out$iid)

    # Convert windows coordinates into end-to-end scale
    x$segments <- dplyr::inner_join(
        x$segments,
        dplyr::select(x$out, dplyr::all_of(
            c("chr_id.don", "chr_id.rec", "iid", "iid.e2e")
        )),
        by = c(
            "chr_id.don" = "chr_id.don",
            "chr_id.rec" = "chr_id.rec",
            "loc.start" = "iid"
        )
    ) |>
        dplyr::inner_join(
            dplyr::select(x$out, dplyr::all_of(
                c("chr_id.don", "chr_id.rec", "iid", "iid.e2e")
            )),
            by = c(
                "chr_id.don" = "chr_id.don",
                "chr_id.rec" = "chr_id.rec",
                "loc.end" = "iid"
            ),
            suffix = c(".start", ".end")
        ) |>
        (\(x) {
            if (!all_segments) {
                dplyr::filter(x, .data$putative_introgression)
            } else {
                x
            }
        })()


    # Color of each chromosome
    x$out$color.don <- (match(x$out$chr_name.don, unique(x$out$chr_name.don)) - 1) %% 2
    x$out$color.rec <- (match(x$out$chr_name.rec, unique(x$out$chr_name.rec)) - 1) %% 2

    # Calculate break coordinates
    breaks <- x$out |>
        dplyr::group_by(.data$chr_name.don, .data$chr_name.rec) |>
        dplyr::reframe(pos = min(.data$iid.e2e) + (diff(range(.data$iid.e2e)) / 2))

    plot <- ggplot2::ggplot(data = x$out) +
        ggplot2::theme_bw() +
        ggplot2::scale_x_continuous(
            expand = ggplot2::expansion(mult = c(0.01, 0.01)),
            breaks = breaks$pos,
            labels = paste(breaks$chr_name.don, breaks$chr_name.rec, sep = "/")
        ) +
        ggplot2::theme(
            aspect.ratio = 1 / 5,
            axis.title.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            panel.grid = ggplot2::element_blank(),
            legend.position = "none"
        )

    if (plot_type == "mirror") {
        don_offset <- round(max(x$out$cor.gc.don, na.rm = TRUE) * 0.1)
        rec_offset <- round(max(x$out$cor.gc.rec, na.rm = TRUE) * 0.1)
        max_y <- max(c(x$out$cor.gc.don + don_offset,
                       x$out$cor.gc.rec + rec_offset))
        
        plot + ggplot2::geom_point(
            mapping = ggplot2::aes(
                x = .data$iid.e2e,
                y = .data$cor.gc.don + don_offset,
                color = factor(.data$color.don)
            ),
            size = 0.5
        ) +
            ggplot2::scale_color_manual(values = c("grey60", "grey80")) +
            ggnewscale::new_scale_color() +
            ggplot2::geom_point(
                mapping = ggplot2::aes(
                    x = .data$iid.e2e,
                    y = -.data$cor.gc.rec - rec_offset,
                    color = factor(.data$color.rec)
                ),
                size = 0.5
            ) +
            ggplot2::scale_color_manual(values = c("grey80", "grey60")) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.03, 0.03)),
                labels = \(x) gsub("-", "", x),
                name = "Coverage",
                limits = c(-max_y, max_y)
            )
    } else {
        max_y <- max(abs(x$out$diff))
        plot +
            ggplot2::geom_point(
                mapping = ggplot2::aes(
                    x = .data$iid.e2e,
                    y = .data$diff,
                    color = factor(.data$color.don)
                ),
                size = 0.5
            ) +
            ggplot2::geom_hline(
                yintercept = 0,
                color = "green",
                linetype = "dashed",
                linewidth = 1
            ) +
            ggplot2::geom_segment(
                data = x$segments,
                mapping = ggplot2::aes(
                    x = .data$iid.e2e.start,
                    y = .data[[x$meta$metrics]],
                    xend = .data$iid.e2e.end,
                    yend = .data[[x$meta$metrics]]
                ),
                color = "orangered",
                linewidth = 1
            ) +
            ggplot2::scale_color_manual(values = c("grey60", "grey80")) + 
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.03, 0.03)),
                labels = \(x) gsub("-", "", x),
                name = "Coverage",
                limits = c(-max_y, max_y)
            )
    }
}


#' @exportS3Method base::print
print.ISA <- function(x, ...) {
    n_putative <- if (!is.null(x$segments) &&
        "putative_introgression" %in% names(x$segments)) {
        sum(x$segments$putative_introgression, na.rm = TRUE)
    } else {
        0L
    }
    cat(
        "ISA object: ",
        nrow(x$meta$target_pairs),
        " pair(s), ",
        nrow(x$segments),
        " segment(s), ",
        n_putative,
        " putative introgression(s), metric: ",
        x$meta$metrics,
        "\n",
        sep = ""
    )
    invisible(x)
}

#' @exportS3Method base::summary
#' @importFrom rlang .data
summary.ISA <- function(object,
                        unit = c("bp", "Kb", "Mb", "Gb"),
                        digits = 2,
                        ...) {
    unit <- match.arg(unit)
    stopifnot("Argument digits must be non-negative" = digits >= 0)

    intro_data <- annotate_segment_coords(
        dplyr::filter(object$segments, .data$putative_introgression),
        object$out
    ) |>
        dplyr::select(dplyr::all_of(
            c(
                "chr_id.don",
                "chr_name.don",
                "subgenome.don",
                "start.don",
                "stop.don",
                "chr_id.rec",
                "chr_name.rec",
                "subgenome.rec",
                "start.rec",
                "stop.rec"
            )
        )) |>
        dplyr::mutate(
            length.don = .data$stop.don - .data$start.don,
            .after = "stop.don"
        ) |>
        dplyr::mutate(
            length.rec = .data$stop.rec - .data$start.rec,
            .after = "stop.rec"
        )

    out <- list(
        introgressions = intro_data,
        n_introgressions = nrow(intro_data),
        unit = unit,
        digits = digits
    )
    class(out) <- "summary.ISA"
    out
}

#' @exportS3Method base::print
print.summary.ISA <- function(x, ...) {
    scale_factor <- list(
        "bp" = 1,
        "Kb" = 1e3,
        "Mb" = 1e6,
        "Gb" = 1e9
    )

    if (x$n_introgressions == 0L) {
        cat("No putative introgressions found\n")
    } else {
        cat(
            "Putative introgressions found:",
            x$n_introgressions,
            "\n"
        )
        cat(
            "Coordinates and length reported in",
            x$unit,
            "scale",
            "\n\n"
        )
        for (i in seq_len(x$n_introgressions)) {
            intro_info <- as.list(x$introgressions[i, ])
            cat("Segment", i, "\n")
            coords <- vapply(
                c(
                    intro_info$start.don,
                    intro_info$stop.don,
                    intro_info$start.rec,
                    intro_info$stop.rec,
                    intro_info$length.don,
                    intro_info$length.rec
                ),
                \(v) round(v / scale_factor[[x$unit]], digits = x$digits),
                numeric(1)
            )

            data.frame(
                row.names = c("Donor", "Recipient"),
                "Chromosome" = c(
                    intro_info$chr_name.don,
                    intro_info$chr_name.rec
                ),
                "Subgenome" = c(
                    intro_info$subgenome.don,
                    intro_info$subgenome.rec
                ),
                "From" = coords[c(1, 3)],
                "To" = coords[c(2, 4)],
                "Length" = coords[c(5, 6)]
            ) |> print()
            cat("\n")
        }
    }
    invisible(x)
}


#' @export
#' @importFrom rlang .data
writeToBedpe <- function(ISA, file) {
    stopifnot("Input is not ISA object" = inherits(ISA, "ISA"))

    output_colnames <- c(
        "chrom1",
        "start1",
        "end1",
        "chrom2",
        "start2",
        "end2",
        "name",
        "score",
        "strand1",
        "strand2",
        "chrom_name1",
        "chrom_name2",
        "subgen_name1",
        "subgen_name2"
    )

    putative <- dplyr::filter(ISA$segments, .data$putative_introgression)

    if (nrow(putative) == 0L) {
        warning("No putative introgressions to write; writing header only",
            call. = FALSE
        )
        empty <- data.frame(matrix(
            ncol = length(output_colnames),
            nrow = 0L,
            dimnames = list(NULL, output_colnames)
        ))
        utils::write.table(
            empty,
            file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )
        message("Output was written in ", file)
        return(invisible(NULL))
    }

    annotate_segment_coords(putative, ISA$out) |>
        dplyr::select(dplyr::all_of(
            c(
                "chr_id.don",
                "start.don",
                "stop.don",
                "chr_id.rec",
                "start.rec",
                "stop.rec",
                ISA$meta$metrics,
                "chr_name.don",
                "chr_name.rec",
                "subgenome.don",
                "subgenome.rec"
            )
        )) |>
        dplyr::mutate(
            name = paste("PIS", seq_len(dplyr::n()), sep = ""),
            .after = "stop.rec"
        ) |>
        dplyr::mutate(
            strand1 = ".",
            strand2 = ".",
            .after = ISA$meta$metrics
        ) |>
        stats::setNames(output_colnames) |>
        utils::write.table(file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
        )

    message("Output was written in ", file)
    invisible(NULL)
}
