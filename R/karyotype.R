## Interactive vertical karyotype for the introgression explorer.
##
## The widget shows the recipient karyotype only. A recipient chromosome is
## kept in place when it has no coverage, drawn as a dashed "deleted" bar. A
## full-chromosome donor addition adds that donor chromosome (one shape per
## homolog copy) to the same group column. The payload is plain lists so
## htmlwidgets can serialise it without a second data-frame layer.

karyotype_palette <- function() {
    list(
        negative = "#2166ac",
        zero = "#f7f7f7",
        positive = "#b2182b",
        substitution = "#d95f02",
        addition = "#1f78b4",
        other = "#7570b3",
        deleted = "#d9d9d9",
        added = "#a6cee3",
        grid = "#e5e5e5",
        text = "#333333"
    )
}

#' Resolve the karyotype sample
#' @noRd
karyotype_sample <- function(result, sample = NULL) {
    if (!inherits(result, "BamAnalysisResult")) {
        .mamamia_error(
            "MaMaMIAKaryotypeError",
            "`result` must be a BamAnalysisResult"
        )
    }
    samples <- unique(as.character(result$samples$sample))
    if (length(samples) == 0L) {
        samples <- unique(as.character(result$scores$sample))
    }
    if (length(samples) == 0L) {
        .mamamia_error(
            "MaMaMIAKaryotypeError",
            "The result contains no samples"
        )
    }
    if (is.null(sample)) {
        return(samples[[1L]])
    }
    if (length(sample) != 1L || is.na(sample) || !nzchar(sample) ||
        !sample %in% samples) {
        .mamamia_error(
            "MaMaMIAKaryotypeError",
            "Unknown sample '%s'",
            paste(sample, collapse = ", ")
        )
    }
    as.character(sample)
}

#' Sort chromosome groups naturally
#' @noRd
karyotype_groups <- function(groups) {
    groups <- unique(as.character(groups))
    groups <- groups[!is.na(groups) & nzchar(groups)]
    numeric <- suppressWarnings(as.numeric(groups))
    groups[order(numeric, groups, na.last = TRUE)]
}

#' Recipient chromosomes with no coverage
#'
#' A chromosome is deleted when its mean breadth is below `deletion_breadth`
#' and below `deletion_relative` times the sample's median recipient breadth.
#' The relative rule keeps a globally shallow sample from losing every
#' chromosome.
#' @noRd
karyotype_deleted <- function(result, sample, deletion_breadth = 0.001,
                              deletion_relative = 0.1) {
    metrics <- result$metrics
    if (is.null(metrics) || nrow(metrics) == 0L) {
        return(character())
    }
    keep <- as.character(metrics$sample) == sample &
        as.character(metrics$role) == "recipient"
    metrics <- metrics[keep, , drop = FALSE]
    if (nrow(metrics) == 0L) {
        return(character())
    }
    breadth <- split(metrics$breadth, as.character(metrics$chr_id))
    mean_breadth <- vapply(breadth, function(x) {
        x <- x[is.finite(x)]
        if (length(x) == 0L) NA_real_ else mean(x)
    }, numeric(1))
    mean_breadth <- mean_breadth[is.finite(mean_breadth)]
    if (length(mean_breadth) == 0L) {
        return(character())
    }
    median_breadth <- stats::median(mean_breadth)
    if (!is.finite(median_breadth) || median_breadth <= 0) {
        return(character())
    }
    names(mean_breadth)[
        mean_breadth < deletion_breadth &
            mean_breadth < deletion_relative * median_breadth
    ]
}

#' Full-chromosome calls
#'
#' A call is full-chromosome when it passes calibration, starts at the first
#' window, ends at the last window and covers at least 95% of the group.
#' @noRd
karyotype_full_calls <- function(result, sample, n_per = NULL) {
    calls <- result$calls
    if (is.null(calls) || nrow(calls) == 0L) {
        return(calls)
    }
    calls <- calls[as.character(calls$sample) == sample, , drop = FALSE]
    if (nrow(calls) == 0L) {
        return(calls)
    }
    if ("pass" %in% names(calls)) {
        calls <- calls[isTRUE_vec(calls$pass), , drop = FALSE]
    }
    if (nrow(calls) == 0L) {
        return(calls)
    }
    if (is.null(n_per) || length(n_per) == 0L || is.na(n_per[[1L]])) {
        n_per <- as.integer(result$windows$n_per_chromosome)[[1L]]
    } else {
        n_per <- as.integer(n_per[[1L]])
    }
    scores <- result$scores
    scores <- scores[as.character(scores$sample) == sample, , drop = FALSE]
    group_length <- tapply(
        scores$length, as.character(scores$group),
        function(x) {
            x <- x[is.finite(x)]
            if (length(x) == 0L) NA_real_ else max(x)
        }
    )
    length_bp <- group_length[as.character(calls$group)]
    fraction <- calls$width_bp / unname(length_bp)
    ## A full-chromosome event starts at the first window and covers at least
    ## 95% of the group. The end is allowed to fall a few windows short: real
    ## terminal segments (for example the 7S addition in the study data) end
    ## at 99% of the profile rather than exactly at the last window.
    full <- calls$start_iid == 0L & is.finite(fraction) & fraction >= 0.95
    calls$group_length <- unname(length_bp)
    calls$fraction <- fraction
    calls$full <- full
    calls[full, , drop = FALSE]
}

#' Build the karyotype payload
#'
#' @param result A `BamAnalysisResult`.
#' @param sample Sample to display; defaults to the first.
#' @param mode `"heatmap"` for per-window scores or `"calls"` for segments.
#' @param calls Calls to overlay; defaults to every call in the result. Pass a
#'   filtered table to respect the explorer controls.
#' @param selected_groups Groups currently selected in the explorer. Unselected
#'   groups are dimmed rather than removed.
#' @param deletion_breadth Absolute mean-breadth cutoff for a deletion.
#' @param deletion_relative Relative-to-median cutoff for a deletion.
#'
#' @return A list suitable for `mamamia_karyotype()`.
#' @noRd
karyotype_payload <- function(result, sample = NULL,
                              mode = c("heatmap", "calls"),
                              calls = NULL, selected_groups = NULL,
                              deletion_breadth = 0.001,
                              deletion_relative = 0.1) {
    mode <- match.arg(mode)
    sample <- karyotype_sample(result, sample)
    assembly <- result$assembly
    recipient <- assembly$contigs$recipient
    donor <- assembly$contigs$donor
    if (is.null(recipient) || nrow(recipient) == 0L) {
        .mamamia_error(
            "MaMaMIAKaryotypeError",
            "The assembly contains no recipient contigs"
        )
    }
    groups <- karyotype_groups(recipient$group)
    recipient_subgenomes <- unique(as.character(recipient$subgenome))
    donor_subgenomes <- if (is.null(donor)) {
        character()
    } else {
        unique(as.character(donor$subgenome))
    }
    n_per <- as.integer(result$windows$n_per_chromosome)[[1L]]
    scores <- result$scores
    scores <- scores[as.character(scores$sample) == sample, , drop = FALSE]
    deleted <- karyotype_deleted(
        result, sample,
        deletion_breadth = deletion_breadth,
        deletion_relative = deletion_relative
    )
    full <- karyotype_full_calls(result, sample, n_per)
    if (is.null(full) || nrow(full) == 0L) {
        full <- data.frame(
            group = character(),
            donor_mean = numeric(),
            stringsAsFactors = FALSE
        )
    }
    overlay <- if (is.null(calls)) {
        result$calls
    } else {
        calls
    }
    if (is.null(overlay)) {
        overlay <- full[0, , drop = FALSE]
    } else {
        overlay <- overlay[as.character(overlay$sample) == sample, , drop = FALSE]
    }
    selected_groups <- if (is.null(selected_groups) ||
        length(selected_groups) == 0L) {
        groups
    } else {
        intersect(as.character(selected_groups), groups)
    }

    ## A cross-group assignment is displayed on the recipient group, not on
    ## the donor group that produced the call. The donor group stays in the
    ## call payload so the tooltip can explain the translocation.
    call_group <- as.character(overlay$group)
    if (all(c("recipient_group", "recipient_cross_group") %in%
        names(overlay))) {
        cross <- as.logical(overlay$recipient_cross_group) %in% TRUE
        recipient_group <- as.character(overlay$recipient_group)
        move <- cross & !is.na(recipient_group) & nzchar(recipient_group)
        call_group[move] <- recipient_group[move]
    }

    group_payload <- lapply(groups, function(group) {
        rec <- recipient[as.character(recipient$group) == group, , drop = FALSE]
        rec <- rec[order(match(
            as.character(rec$subgenome), recipient_subgenomes
        )), , drop = FALSE]
        rec_payload <- lapply(seq_len(nrow(rec)), function(i) {
            list(
                id = as.character(rec$id[[i]]),
                name = as.character(rec$name[[i]]),
                subgenome = as.character(rec$subgenome[[i]]),
                length = as.numeric(rec$length[[i]]),
                deleted = as.character(rec$id[[i]]) %in% deleted
            )
        })

        add <- full[
            as.character(full$group) == group & full$donor_mean > 0.5,
            , drop = FALSE
        ]
        add_payload <- list()
        if (nrow(add) > 0L && !is.null(donor) && nrow(donor) > 0L) {
            don <- donor[as.character(donor$group) == group, , drop = FALSE]
            don <- don[order(match(
                as.character(don$subgenome), donor_subgenomes
            )), , drop = FALSE]
            for (i in seq_len(nrow(don))) {
                copies <- suppressWarnings(as.integer(don$copy_number[[i]]))
                if (length(copies) != 1L || is.na(copies) || copies < 1L) {
                    copies <- 1L
                }
                for (copy in seq_len(copies)) {
                    add_payload[[length(add_payload) + 1L]] <- list(
                        id = as.character(don$id[[i]]),
                        name = as.character(don$name[[i]]),
                        subgenome = as.character(don$subgenome[[i]]),
                        length = as.numeric(don$length[[i]]),
                        copy = copy,
                        copies = copies
                    )
                }
            }
        }

        group_scores <- scores[
            as.character(scores$group) == group, , drop = FALSE
        ]
        group_scores <- group_scores[order(group_scores$iid), , drop = FALSE]
        score_values <- if (nrow(group_scores) > 0L) {
            as.numeric(group_scores$score)
        } else {
            numeric()
        }
        group_calls <- overlay[call_group == group, , drop = FALSE]
        call_payload <- lapply(seq_len(nrow(group_calls)), function(i) {
            start <- if (n_per > 0L) {
                as.numeric(group_calls$start_iid[[i]]) / n_per
            } else {
                0
            }
            end <- if (n_per > 0L) {
                (as.numeric(group_calls$end_iid[[i]]) + 1) / n_per
            } else {
                1
            }
            list(
                start = max(0, min(1, start)),
                end = max(0, min(1, end)),
                start_bp = as.numeric(group_calls$start_bp[[i]]),
                end_bp = as.numeric(group_calls$end_bp[[i]]),
                width_bp = as.numeric(group_calls$width_bp[[i]]),
                class = as.character(group_calls$class[[i]]),
                confidence = as.character(group_calls$confidence[[i]]),
                score = as.numeric(group_calls$score_mean[[i]]),
                pass = if ("pass" %in% names(group_calls)) {
                    isTRUE(group_calls$pass[[i]])
                } else {
                    TRUE
                },
                recipient_subgenome = if ("recipient_subgenome" %in%
                    names(group_calls)) {
                    as.character(group_calls$recipient_subgenome[[i]])
                } else {
                    NA_character_
                },
                recipient_chr_id = if ("recipient_chr_id" %in%
                    names(group_calls)) {
                    as.character(group_calls$recipient_chr_id[[i]])
                } else {
                    NA_character_
                },
                recipient_chr_name = if ("recipient_chr_name" %in%
                    names(group_calls)) {
                    as.character(group_calls$recipient_chr_name[[i]])
                } else {
                    NA_character_
                },
                recipient_start_bp = if ("recipient_start_bp" %in%
                    names(group_calls)) {
                    as.numeric(group_calls$recipient_start_bp[[i]])
                } else {
                    NA_real_
                },
                recipient_end_bp = if ("recipient_end_bp" %in%
                    names(group_calls)) {
                    as.numeric(group_calls$recipient_end_bp[[i]])
                } else {
                    NA_real_
                },
                recipient_assignment = if ("recipient_assignment" %in%
                    names(group_calls)) {
                    as.character(group_calls$recipient_assignment[[i]])
                } else {
                    NA_character_
                },
                recipient_weight = if ("recipient_weight" %in%
                    names(group_calls)) {
                    as.numeric(group_calls$recipient_weight[[i]])
                } else {
                    NA_real_
                },
                donor_subgenome = if ("donor_subgenome" %in%
                    names(group_calls)) {
                    as.character(group_calls$donor_subgenome[[i]])
                } else {
                    NA_character_
                },
                donor_chr_id = if ("donor_chr_id" %in% names(group_calls)) {
                    as.character(group_calls$donor_chr_id[[i]])
                } else {
                    NA_character_
                },
                donor_chr_name = if ("donor_chr_name" %in%
                    names(group_calls)) {
                    as.character(group_calls$donor_chr_name[[i]])
                } else {
                    NA_character_
                },
                donor_group = as.character(group_calls$group[[i]]),
                recipient_group = if ("recipient_group" %in%
                    names(group_calls)) {
                    as.character(group_calls$recipient_group[[i]])
                } else {
                    NA_character_
                },
                cross_group = if ("recipient_cross_group" %in%
                    names(group_calls)) {
                    isTRUE(as.logical(
                        group_calls$recipient_cross_group[[i]]
                    ))
                } else {
                    FALSE
                },
                recipient_alternative = if ("recipient_alternative" %in%
                    names(group_calls)) {
                    as.character(group_calls$recipient_alternative[[i]])
                } else {
                    NA_character_
                },
                recipient_mean_depletion = if ("recipient_mean_depletion" %in%
                    names(group_calls)) {
                    as.numeric(group_calls$recipient_mean_depletion[[i]])
                } else {
                    NA_real_
                }
            )
        })

        list(
            group = group,
            label = paste("Group", group),
            length = if (nrow(rec) > 0L) {
                max(rec$length, na.rm = TRUE)
            } else {
                NA_real_
            },
            selected = group %in% selected_groups,
            recipient = rec_payload,
            added = add_payload,
            scores = score_values,
            calls = call_payload,
            deleted = any(vapply(
                rec_payload, function(x) isTRUE(x$deleted), logical(1)
            ))
        )
    })

    limit <- suppressWarnings(stats::quantile(
        abs(scores$score), 0.99, na.rm = TRUE, names = FALSE
    ))
    if (length(limit) != 1L || !is.finite(limit) || limit <= 0) {
        limit <- 0.5
    }
    list(
        mode = mode,
        sample = sample,
        score_limit = max(0.5, as.numeric(limit)),
        palette = karyotype_palette(),
        groups = group_payload
    )
}

#' Create the karyotype htmlwidget
#' @noRd
mamamia_karyotype <- function(payload, width = NULL, height = NULL) {
    css <- htmltools::htmlDependency(
        name = "mamamiaKaryotype-css",
        version = as.character(utils::packageVersion("MaMaMIA")),
        src = c(file = system.file("htmlwidgets", package = "MaMaMIA")),
        stylesheet = "mamamiaKaryotype.css"
    )
    htmlwidgets::createWidget(
        name = "mamamiaKaryotype",
        x = payload,
        width = width,
        height = height,
        package = "MaMaMIA",
        dependencies = css,
        sizingPolicy = htmlwidgets::sizingPolicy(
            viewer.fill = TRUE,
            browser.fill = TRUE,
            padding = 0
        )
    )
}

#' Shiny output for the karyotype
#' @noRd
mamamia_karyotype_output <- function(outputId, width = "100%",
                                     height = "480px") {
    htmlwidgets::shinyWidgetOutput(
        outputId, "mamamiaKaryotype", width, height, package = "MaMaMIA"
    )
}

#' Shiny renderer for the karyotype
#' @noRd
render_mamamia_karyotype <- function(expr, env = parent.frame(),
                                     quoted = FALSE) {
    if (!quoted) {
        expr <- substitute(expr)
    }
    htmlwidgets::shinyRenderWidget(
        expr, mamamia_karyotype_output, env, quoted = TRUE
    )
}
