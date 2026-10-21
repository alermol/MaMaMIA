## The Introgression Segment Array (ISA): methods to inspect, subset and plot
## an existing segmentation result. See segmentation.R for `segments()`, which
## builds an ISA object.

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

    value_col <- isa_value_col(x)
    don_values <- x$out[[paste0(value_col, ".don")]]
    rec_values <- x$out[[paste0(value_col, ".rec")]]
    value_label <- if (identical(value_col, "dose")) "Dosage" else "Coverage"

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
        don_offset <- round(max(don_values, na.rm = TRUE) * 0.1)
        rec_offset <- round(max(rec_values, na.rm = TRUE) * 0.1)
        max_y <- max(c(don_values + don_offset,
                       rec_values + rec_offset))

        plot + ggplot2::geom_point(
            mapping = ggplot2::aes(
                x = .data$iid.e2e,
                y = don_values + don_offset,
                color = factor(.data$color.don)
            ),
            size = 0.5
        ) +
            ggplot2::scale_color_manual(values = c("grey60", "grey80")) +
            ggnewscale::new_scale_color() +
            ggplot2::geom_point(
                mapping = ggplot2::aes(
                    x = .data$iid.e2e,
                    y = -rec_values - rec_offset,
                    color = factor(.data$color.rec)
                ),
                size = 0.5
            ) +
            ggplot2::scale_color_manual(values = c("grey80", "grey60")) +
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.03, 0.03)),
                labels = \(x) gsub("-", "", x),
                name = value_label,
                limits = c(-max_y, max_y)
            )
    } else {
        if (is.null(x$out$diff)) {
            x$out$diff <- don_values - rec_values
        }
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
                name = paste(value_label, "difference"),
                limits = c(-max_y, max_y)
            )
    }
}

## Resolve the per-window value column of an ISA, falling back to `cor.gc`
## when the recorded column is absent (older ISA objects).
#' @noRd
isa_value_col <- function(x) {
    value_col <- x$meta$value_col
    if (is.null(value_col) || !is.character(value_col) ||
        length(value_col) != 1L || is.na(value_col) || !nzchar(value_col)) {
        value_col <- if (all(c("dose.don", "dose.rec") %in% names(x$out))) {
            "dose"
        } else {
            "cor.gc"
        }
    }
    if (!all(c(paste0(value_col, ".don"), paste0(value_col, ".rec")) %in% names(x$out))) {
        value_col <- "cor.gc"
    }
    value_col
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
        if (!is.null(x$meta$value_col)) {
            paste0(", value: ", x$meta$value_col)
        } else {
            ""
        },
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
