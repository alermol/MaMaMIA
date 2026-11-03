## Diagnostic plots for the downloadable result bundle. Every plot is built
## from the same renamed tables that are written to `qc/`, so the PNG files and
## the raw CSV data always describe the same rows.

bam_bundle_plot_theme <- function() {
    ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_line(
                colour = "#e6ece8", linewidth = 0.3
            ),
            plot.title = ggplot2::element_text(
                size = 12.5, face = "bold", colour = "#1f2a24"
            ),
            plot.subtitle = ggplot2::element_text(
                size = 9.5, colour = "#4f5f55", margin =
                    ggplot2::margin(b = 6)
            ),
            plot.caption = ggplot2::element_text(
                size = 8.5, colour = "#65736a", hjust = 0
            ),
            axis.title = ggplot2::element_text(
                size = 9.5, colour = "#33443a"
            ),
            axis.text = ggplot2::element_text(
                size = 8.5, colour = "#4f5f55"
            ),
            strip.text = ggplot2::element_text(
                size = 8.5, colour = "#33443a", face = "bold"
            ),
            legend.title = ggplot2::element_blank(),
            legend.position = "top",
            legend.text = ggplot2::element_text(
                size = 9, colour = "#33443a"
            ),
            plot.background = ggplot2::element_rect(
                fill = "white", colour = NA
            ),
            panel.background = ggplot2::element_rect(
                fill = "white", colour = NA
            ),
            plot.margin = ggplot2::margin(10, 12, 8, 10)
        )
}

bam_bundle_plot_groups <- function(x) {
    groups <- bam_bundle_group_order(x)
    factor(as.character(x), levels = groups)
}

bam_bundle_plot_coverage <- function(metrics) {
    if (!is.data.frame(metrics) || nrow(metrics) == 0L ||
        !all(c("chromosome_group", "role", "breadth") %in% names(metrics))) {
        return(NULL)
    }
    x <- metrics[is.finite(metrics$breadth), , drop = FALSE]
    if (nrow(x) == 0L) {
        return(NULL)
    }
    x$chromosome_group <- bam_bundle_plot_groups(x$chromosome_group)
    x$role <- factor(as.character(x$role), levels = c("recipient", "donor"))
    ggplot2::ggplot(
        x,
        ggplot2::aes(
            x = .data$chromosome_group, y = .data$breadth,
            fill = .data$role
        )
    ) +
        ggplot2::geom_boxplot(
            outlier.size = 0.25, outlier.alpha = 0.25, linewidth = 0.3
        ) +
        ggplot2::scale_fill_manual(
            values = c(recipient = "#2f6f5e", donor = "#c26b3f"),
            drop = FALSE
        ) +
        ggplot2::labs(
            x = "Chromosome group",
            y = "Window breadth",
            title = "Coverage breadth by chromosome group",
            subtitle = paste(
                "One box per role and group; each observation is a",
                "window."
            ),
            caption = paste(
                "Breadth is the fraction of the window covered by at least",
                "one high-MAPQ read. Raw values: qc/metrics-windows.csv.gz."
            )
        ) +
        bam_bundle_plot_theme()
}

bam_bundle_plot_scores <- function(scores) {
    if (!is.data.frame(scores) || nrow(scores) == 0L ||
        !all(c("chromosome_group", "introgression_score") %in%
            names(scores))) {
        return(NULL)
    }
    x <- scores[is.finite(scores$introgression_score), , drop = FALSE]
    if (nrow(x) == 0L) {
        return(NULL)
    }
    x$chromosome_group <- bam_bundle_plot_groups(x$chromosome_group)
    ggplot2::ggplot(
        x,
        ggplot2::aes(x = .data$chromosome_group, y = .data$introgression_score)
    ) +
        ggplot2::geom_hline(
            yintercept = 0, linewidth = 0.3, linetype = "dashed",
            colour = "#8a6d3b"
        ) +
        ggplot2::geom_boxplot(
            fill = "#2f6f5e", colour = "#1f4f42", alpha = 0.55,
            outlier.shape = NA, linewidth = 0.3, width = 0.68
        ) +
        ggplot2::labs(
            x = "Chromosome group",
            y = "Introgression score",
            title = "Introgression score by chromosome group",
            subtitle = paste(
                "Boxes show the median and interquartile range;",
                "whiskers reach 1.5 IQR."
            ),
            caption = paste(
                "Outliers are not drawn. Raw values:",
                "qc/scores-windows.csv.gz."
            )
        ) +
        bam_bundle_plot_theme()
}

bam_bundle_plot_score_profile <- function(scores) {
    if (!is.data.frame(scores) || nrow(scores) == 0L ||
        !all(c(
            "chromosome_group", "window_index", "introgression_score"
        ) %in% names(scores))) {
        return(NULL)
    }
    x <- scores[
        is.finite(scores$window_index) &
            is.finite(scores$introgression_score), ,
        drop = FALSE
    ]
    if (nrow(x) == 0L) {
        return(NULL)
    }
    x$chromosome_group <- bam_bundle_plot_groups(x$chromosome_group)
    n_groups <- length(unique(x$chromosome_group))
    n_col <- min(4L, max(1L, n_groups))
    ggplot2::ggplot(
        x,
        ggplot2::aes(x = .data$window_index, y = .data$introgression_score)
    ) +
        ggplot2::geom_hline(
            yintercept = 0, linewidth = 0.3, linetype = "dashed",
            colour = "#8a6d3b"
        ) +
        ggplot2::geom_line(
            linewidth = 0.3, colour = "#2f6f5e", alpha = 0.85
        ) +
        ggplot2::facet_wrap(
            ggplot2::vars(.data$chromosome_group),
            scales = "free_x", ncol = n_col
        ) +
        ggplot2::labs(
            x = "Window index",
            y = "Introgression score",
            title = "Introgression score along each chromosome group",
            subtitle = paste(
                "Every window score in the selected sample;",
                "each panel is one chromosome group."
            ),
            caption = "Raw values: qc/scores-windows.csv.gz."
        ) +
        bam_bundle_plot_theme() +
        ggplot2::theme(legend.position = "none")
}

bam_bundle_plot_residuals <- function(scores) {
    columns <- c(
        "donor_residual", "recipient_residual", "total_depth_residual"
    )
    if (!is.data.frame(scores) || nrow(scores) == 0L ||
        !all(columns %in% names(scores))) {
        return(NULL)
    }
    labels <- c(
        donor_residual = "Donor",
        recipient_residual = "Recipient",
        total_depth_residual = "Total depth"
    )
    parts <- lapply(columns, function(column) {
        values <- suppressWarnings(as.numeric(scores[[column]]))
        keep <- is.finite(values)
        if (!any(keep)) {
            return(NULL)
        }
        data.frame(
            series = unname(labels[[column]]),
            residual = values[keep],
            stringsAsFactors = FALSE
        )
    })
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0L) {
        return(NULL)
    }
    long <- do.call(rbind, parts)
    long$series <- factor(
        long$series, levels = unname(labels[columns])
    )
    limits <- stats::quantile(
        long$residual, c(0.005, 0.995), na.rm = TRUE
    )
    limits <- as.numeric(limits)
    if (!all(is.finite(limits)) || limits[[1L]] == limits[[2L]]) {
        limits <- range(long$residual, na.rm = TRUE)
    }
    if (!all(is.finite(limits)) || limits[[1L]] == limits[[2L]]) {
        limits <- limits + c(-1, 1)
    }
    ggplot2::ggplot(
        long,
        ggplot2::aes(
            x = .data$residual, colour = .data$series,
            fill = .data$series
        )
    ) +
        ggplot2::geom_vline(
            xintercept = 0, linewidth = 0.3, linetype = "dashed",
            colour = "#8a6d3b"
        ) +
        ggplot2::geom_density(
            linewidth = 0.45, alpha = 0.12, na.rm = TRUE
        ) +
        ggplot2::scale_colour_manual(
            values = c(
                Donor = "#c26b3f", Recipient = "#2f6f5e",
                `Total depth` = "#5b6b60"
            )
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                Donor = "#c26b3f", Recipient = "#2f6f5e",
                `Total depth` = "#5b6b60"
            )
        ) +
        ggplot2::coord_cartesian(xlim = limits) +
        ggplot2::labs(
            x = "Residual",
            y = "Density",
            title = "Depth residuals after GC and copy-number adjustment",
            subtitle = paste(
                "Densities are clipped to the 0.5th-99.5th percentile;",
                "the full range stays in the raw table."
            ),
            caption = "Raw values: qc/scores-windows.csv.gz."
        ) +
        bam_bundle_plot_theme()
}

bam_bundle_plot_masks <- function(masks_summary) {
    if (!is.data.frame(masks_summary) || nrow(masks_summary) == 0L ||
        !all(c("chromosome_group", "masked_fraction") %in%
            names(masks_summary))) {
        return(NULL)
    }
    x <- masks_summary[
        is.finite(masks_summary$masked_fraction), , drop = FALSE
    ]
    if (nrow(x) == 0L) {
        return(NULL)
    }
    x$chromosome_group <- bam_bundle_plot_groups(x$chromosome_group)
    masked <- if ("masked_windows" %in% names(x)) {
        sum(as.numeric(x$masked_windows), na.rm = TRUE)
    } else {
        sum(x$masked_fraction, na.rm = TRUE)
    }
    upper <- max(x$masked_fraction, na.rm = TRUE)
    if (!is.finite(upper) || upper <= 0) {
        upper <- 1
    }
    subtitle <- if (isTRUE(masked > 0)) {
        "Fraction of windows masked before scoring in each group."
    } else {
        "No windows were masked in this run; every bar is zero."
    }
    ggplot2::ggplot(
        x,
        ggplot2::aes(x = .data$chromosome_group, y = .data$masked_fraction)
    ) +
        ggplot2::geom_col(
            fill = "#d8c6a5", colour = "#8a6d3b", linewidth = 0.3,
            width = 0.7
        ) +
        ggplot2::scale_y_continuous(
            labels = function(value) paste0(round(value * 100), "%"),
            limits = c(0, upper)
        ) +
        ggplot2::labs(
            x = "Chromosome group",
            y = "Masked fraction",
            title = "Masked windows by chromosome group",
            subtitle = subtitle,
            caption = paste(
                "Raw values: qc/masks-summary.csv and",
                "qc/masks-windows.csv.gz."
            )
        ) +
        bam_bundle_plot_theme()
}

bam_bundle_plot_calls <- function(calls_summary) {
    if (!is.data.frame(calls_summary) || nrow(calls_summary) == 0L ||
        !all(c(
            "call_class", "calls", "passing_calls", "confident_calls"
        ) %in% names(calls_summary))) {
        return(NULL)
    }
    x <- calls_summary
    classes <- unique(as.character(x$call_class))
    classes <- classes[!is.na(classes) & nzchar(classes)]
    if (length(classes) == 0L) {
        return(NULL)
    }
    series <- c("All calls", "Passing calibration", "Confident")
    long <- do.call(rbind, lapply(seq_along(classes), function(i) {
        data.frame(
            call_class = classes[[i]],
            series = series,
            value = c(
                x$calls[[i]], x$passing_calls[[i]], x$confident_calls[[i]]
            ),
            stringsAsFactors = FALSE
        )
    }))
    long$call_class <- factor(long$call_class, levels = classes)
    long$series <- factor(long$series, levels = series)
    ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$call_class, y = .data$value, fill = .data$series)
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.82),
            width = 0.74
        ) +
        ggplot2::scale_fill_manual(
            values = c(
                `All calls` = "#b9c6bf",
                `Passing calibration` = "#2f6f5e",
                Confident = "#8a6d3b"
            ),
            drop = FALSE
        ) +
        ggplot2::labs(
            x = "Call class",
            y = "Calls",
            title = "Calls by class and review tier",
            subtitle = paste(
                "Every call, calls passing calibration, and confident",
                "calls."
            ),
            caption = "Raw values: calls-summary.csv."
        ) +
        bam_bundle_plot_theme()
}

bam_bundle_plot_cohort <- function(cohort) {
    if (!is.data.frame(cohort) || nrow(cohort) < 2L ||
        !all(c("sample", "n_all_calls") %in% names(cohort))) {
        return(NULL)
    }
    x <- cohort
    confident <- if ("n_class_confident" %in% names(x)) {
        as.numeric(x$n_class_confident)
    } else {
        rep(NA_real_, nrow(x))
    }
    long <- rbind(
        data.frame(
            sample = as.character(x$sample),
            series = "All calls",
            value = as.numeric(x$n_all_calls),
            stringsAsFactors = FALSE
        ),
        data.frame(
            sample = as.character(x$sample),
            series = "Confident",
            value = confident,
            stringsAsFactors = FALSE
        )
    )
    long <- long[is.finite(long$value), , drop = FALSE]
    if (nrow(long) == 0L) {
        return(NULL)
    }
    samples <- unique(as.character(x$sample))
    long$sample <- factor(long$sample, levels = samples)
    ggplot2::ggplot(
        long,
        ggplot2::aes(x = .data$sample, y = .data$value, fill = .data$series)
    ) +
        ggplot2::geom_col(
            position = ggplot2::position_dodge(width = 0.8),
            width = 0.7
        ) +
        ggplot2::scale_fill_manual(
            values = c(`All calls` = "#b9c6bf", Confident = "#2f6f5e"),
            drop = FALSE
        ) +
        ggplot2::labs(
            x = NULL,
            y = "Calls",
            title = "Calls per sample in this run",
            subtitle = paste(
                "All samples in the selected BAM, not only the exported",
                "sample."
            ),
            caption = "Raw values: qc/cohort-summary.csv."
        ) +
        bam_bundle_plot_theme() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(
                angle = 30, hjust = 1, vjust = 1
            )
        )
}

bam_bundle_save_plot <- function(plot, path, width = 9, height = 5.5) {
    if (is.null(plot)) {
        return(FALSE)
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    attempts <- list(
        function() {
            ggplot2::ggsave(
                filename = path, plot = plot,
                width = width, height = height, units = "in",
                dpi = 120, bg = "white", device = "png", type = "cairo"
            )
        },
        function() {
            ggplot2::ggsave(
                filename = path, plot = plot,
                width = width, height = height, units = "in",
                dpi = 120, bg = "white", device = "png"
            )
        }
    )
    for (attempt in attempts) {
        ok <- tryCatch(
            {
                suppressMessages(attempt())
                TRUE
            },
            error = function(e) FALSE
        )
        if (ok && file.exists(path) &&
            isTRUE(file.info(path)$size > 0)) {
            return(TRUE)
        }
        unlink(path, force = TRUE)
    }
    FALSE
}

bam_bundle_plot_files <- function(bundle, qc, calls_summary) {
    plot_dir <- file.path(bundle, "qc", "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    specs <- list(
        list(
            file = "coverage-by-group.png",
            title = "Coverage breadth by chromosome group",
            caption = paste(
                "Window breadth for donor and recipient contigs.",
                "Raw values: qc/metrics-windows.csv.gz."
            ),
            plot = bam_bundle_plot_coverage(qc$metrics),
            width = 9, height = 5.5
        ),
        list(
            file = "score-by-group.png",
            title = "Introgression score by chromosome group",
            caption = paste(
                "Score distribution for every group window.",
                "Raw values: qc/scores-windows.csv.gz."
            ),
            plot = bam_bundle_plot_scores(qc$scores),
            width = 9, height = 5.5
        ),
        list(
            file = "score-profile.png",
            title = "Introgression score along each chromosome group",
            caption = paste(
                "Score track across window positions.",
                "Raw values: qc/scores-windows.csv.gz."
            ),
            plot = bam_bundle_plot_score_profile(qc$scores),
            width = 10, height = 6.8
        ),
        list(
            file = "depth-residuals.png",
            title = "Depth residuals",
            caption = paste(
                "Donor, recipient and total-depth residuals after",
                "adjustment. Raw values: qc/scores-windows.csv.gz."
            ),
            plot = bam_bundle_plot_residuals(qc$scores),
            width = 9, height = 5.5
        ),
        list(
            file = "masked-fraction-by-group.png",
            title = "Masked windows by chromosome group",
            caption = paste(
                "Masked fraction per group.",
                "Raw values: qc/masks-summary.csv."
            ),
            plot = bam_bundle_plot_masks(qc$masks_summary),
            width = 9, height = 5.0
        ),
        list(
            file = "calls-by-class.png",
            title = "Calls by class and review tier",
            caption = paste(
                "All calls, passing calls and confident calls by class.",
                "Raw values: calls-summary.csv."
            ),
            plot = bam_bundle_plot_calls(calls_summary),
            width = 9, height = 5.5
        ),
        list(
            file = "cohort-calls.png",
            title = "Calls per sample in this run",
            caption = paste(
                "Calls for every sample in the BAM.",
                "Raw values: qc/cohort-summary.csv."
            ),
            plot = bam_bundle_plot_cohort(qc$cohort),
            width = 10, height = 5.2
        )
    )
    rows <- lapply(specs, function(spec) {
        path <- file.path(plot_dir, spec$file)
        if (!bam_bundle_save_plot(
            spec$plot, path, width = spec$width, height = spec$height
        )) {
            return(NULL)
        }
        data.frame(
            path = file.path("qc", "plots", spec$file),
            title = spec$title,
            caption = spec$caption,
            stringsAsFactors = FALSE
        )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) {
        return(data.frame(
            path = character(), title = character(), caption = character(),
            stringsAsFactors = FALSE
        ))
    }
    do.call(rbind, rows)
}

bam_bundle_png_data_uri <- function(path) {
    if (!file.exists(path)) {
        return(NA_character_)
    }
    size <- file.info(path)$size
    if (!is.finite(size) || size <= 0) {
        return(NA_character_)
    }
    encoded <- jsonlite::base64_enc(readBin(path, "raw", n = size))
    paste0(
        "data:image/png;base64,",
        gsub("[[:space:]]", "", encoded, useBytes = TRUE)
    )
}

bam_bundle_plot_html <- function(plots, bundle, file) {
    if (!is.data.frame(plots) || nrow(plots) == 0L ||
        !all(c("path", "title", "caption") %in% names(plots))) {
        return(character())
    }
    wanted <- file.path("qc", "plots", file)
    index <- match(wanted, as.character(plots$path))
    if (is.na(index)) {
        return(character())
    }
    uri <- bam_bundle_png_data_uri(file.path(bundle, wanted))
    if (is.na(uri)) {
        return(character())
    }
    c(
        "<figure>",
        paste0(
            "<img src=\"", uri, "\" alt=\"",
            bam_bundle_html_escape(plots$title[[index]]), "\">"
        ),
        paste0(
            "<figcaption><strong>",
            bam_bundle_html_escape(plots$title[[index]]),
            "</strong>",
            if (nzchar(plots$caption[[index]])) {
                paste0(
                    "<br>", bam_bundle_html_escape(plots$caption[[index]])
                )
            } else {
                ""
            },
            "</figcaption>"
        ),
        "</figure>"
    )
}
