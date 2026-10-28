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

    ## Calls are attached to the chromosome they belong to, not to a group
    ## column. This lets the layout use one column per subgenome and scale
    ## every chromosome bar by its configured length.
    recipient_ids <- if ("recipient_chr_id" %in% names(overlay)) {
        as.character(overlay$recipient_chr_id)
    } else {
        rep(NA_character_, nrow(overlay))
    }
    recipient_subs <- if ("recipient_subgenome" %in% names(overlay)) {
        as.character(overlay$recipient_subgenome)
    } else {
        rep(NA_character_, nrow(overlay))
    }
    recipient_start <- if ("recipient_start_bp" %in% names(overlay)) {
        as.numeric(overlay$recipient_start_bp)
    } else {
        rep(NA_real_, nrow(overlay))
    }
    recipient_end <- if ("recipient_end_bp" %in% names(overlay)) {
        as.numeric(overlay$recipient_end_bp)
    } else {
        rep(NA_real_, nrow(overlay))
    }
    recipient_assignment <- if ("recipient_assignment" %in% names(overlay)) {
        as.character(overlay$recipient_assignment)
    } else {
        rep(NA_character_, nrow(overlay))
    }
    donor_ids <- if ("donor_chr_id" %in% names(overlay)) {
        as.character(overlay$donor_chr_id)
    } else {
        rep(NA_character_, nrow(overlay))
    }
    donor_subs <- if ("donor_subgenome" %in% names(overlay)) {
        as.character(overlay$donor_subgenome)
    } else {
        rep(NA_character_, nrow(overlay))
    }
    donor_start <- if ("donor_start_bp" %in% names(overlay)) {
        as.numeric(overlay$donor_start_bp)
    } else {
        rep(NA_real_, nrow(overlay))
    }
    donor_end <- if ("donor_end_bp" %in% names(overlay)) {
        as.numeric(overlay$donor_end_bp)
    } else {
        rep(NA_real_, nrow(overlay))
    }
    call_groups <- if ("group" %in% names(overlay)) {
        as.character(overlay$group)
    } else {
        rep(NA_character_, nrow(overlay))
    }
    call_value <- function(name, i, default = NA) {
        if (!name %in% names(overlay)) {
            return(default)
        }
        overlay[[name]][[i]]
    }
    score_values <- lapply(
        split(scores, as.character(scores$group)),
        function(x) {
            x <- x[order(x$iid), , drop = FALSE]
            as.numeric(x$score)
        }
    )
    make_call <- function(i, start_bp, end_bp, length_bp, n_per) {
        start <- NA_real_
        end <- NA_real_
        if (is.finite(start_bp) && is.finite(length_bp) && length_bp > 0) {
            start <- start_bp / length_bp
        }
        if (is.finite(end_bp) && is.finite(length_bp) && length_bp > 0) {
            end <- end_bp / length_bp
        }
        if (!is.finite(start) || !is.finite(end) || end <= start) {
            start <- as.numeric(call_value("start_iid", i, 0)) / n_per
            end <- (as.numeric(call_value("end_iid", i, n_per - 1L)) + 1) /
                n_per
            if (is.finite(length_bp) && length_bp > 0) {
                start_bp <- start * length_bp
                end_bp <- end * length_bp
            }
        }
        start <- max(0, min(1, start))
        end <- max(0, min(1, end))
        list(
            start = start,
            end = end,
            start_bp = as.numeric(start_bp),
            end_bp = as.numeric(end_bp),
            width_bp = as.numeric(call_value("width_bp", i, NA_real_)),
            class = as.character(call_value("class", i, NA_character_)),
            confidence = as.character(
                call_value("confidence", i, NA_character_)
            ),
            score = as.numeric(call_value("score_mean", i, NA_real_)),
            pass = isTRUE(as.logical(call_value("pass", i, TRUE))),
            recipient_subgenome = as.character(
                call_value("recipient_subgenome", i, NA_character_)
            ),
            recipient_chr_id = as.character(
                call_value("recipient_chr_id", i, NA_character_)
            ),
            recipient_chr_name = as.character(
                call_value("recipient_chr_name", i, NA_character_)
            ),
            recipient_start_bp = as.numeric(
                call_value("recipient_start_bp", i, NA_real_)
            ),
            recipient_end_bp = as.numeric(
                call_value("recipient_end_bp", i, NA_real_)
            ),
            recipient_assignment = as.character(
                call_value("recipient_assignment", i, NA_character_)
            ),
            recipient_weight = as.numeric(
                call_value("recipient_weight", i, NA_real_)
            ),
            recipient_group = as.character(
                call_value("recipient_group", i, NA_character_)
            ),
            donor_subgenome = as.character(
                call_value("donor_subgenome", i, NA_character_)
            ),
            donor_chr_id = as.character(
                call_value("donor_chr_id", i, NA_character_)
            ),
            donor_chr_name = as.character(
                call_value("donor_chr_name", i, NA_character_)
            ),
            donor_group = as.character(
                call_value("group", i, NA_character_)
            ),
            cross_group = isTRUE(as.logical(
                call_value("recipient_cross_group", i, FALSE)
            )),
            recipient_alternative = as.character(
                call_value("recipient_alternative", i, NA_character_)
            ),
            recipient_mean_depletion = as.numeric(
                call_value("recipient_mean_depletion", i, NA_real_)
            )
        )
    }
    recipient_calls <- function(chr_id, group, subgenome, length_bp) {
        if (nrow(overlay) == 0L) {
            return(list())
        }
        own <- !is.na(recipient_ids) & recipient_ids == chr_id
        fallback <- is.na(recipient_ids) & !is.na(recipient_subs) &
            recipient_subs == subgenome & call_groups == group
        keep <- is.na(recipient_assignment) |
            recipient_assignment != "no_recipient_loss"
        rows <- which((own | fallback) & keep)
        lapply(rows, function(i) {
            make_call(i, recipient_start[[i]], recipient_end[[i]],
                length_bp, n_per
            )
        })
    }
    donor_calls <- function(chr_id, group, subgenome, length_bp) {
        if (nrow(overlay) == 0L) {
            return(list())
        }
        own <- !is.na(donor_ids) & donor_ids == chr_id
        fallback <- is.na(donor_ids) & !is.na(donor_subs) &
            donor_subs == subgenome & call_groups == group
        rows <- which(own | fallback)
        lapply(rows, function(i) {
            make_call(i, donor_start[[i]], donor_end[[i]],
                length_bp, n_per
            )
        })
    }

    added_groups <- if (nrow(full) > 0L) {
        keep <- !is.na(full$donor_mean) & full$donor_mean > 0.5
        unique(as.character(full$group[keep]))
    } else {
        character()
    }
    if (is.null(donor) || nrow(donor) == 0L) {
        donor <- data.frame(
            id = character(), name = character(), subgenome = character(),
            group = character(), length = numeric(),
            copy_number = numeric(), stringsAsFactors = FALSE
        )
    }
    added_donor <- donor[
        as.character(donor$group) %in% added_groups, , drop = FALSE
    ]
    added_donor <- added_donor[order(
        suppressWarnings(as.numeric(added_donor$group)),
        as.character(added_donor$name)
    ), , drop = FALSE]
    subgenome_levels <- unique(c(
        as.character(recipient$subgenome),
        as.character(added_donor$subgenome)
    ))
    subgenome_payload <- lapply(subgenome_levels, function(level) {
        rec <- recipient[
            as.character(recipient$subgenome) == level, , drop = FALSE
        ]
        rec <- rec[order(
            suppressWarnings(as.numeric(rec$group)),
            as.character(rec$name)
        ), , drop = FALSE]
        add <- added_donor[
            as.character(added_donor$subgenome) == level, , drop = FALSE
        ]
        chromosomes <- list()
        for (i in seq_len(nrow(rec))) {
            chr_id <- as.character(rec$id[[i]])
            chr_name <- as.character(rec$name[[i]])
            chr_group <- as.character(rec$group[[i]])
            chr_length <- as.numeric(rec$length[[i]])
            chromosomes[[length(chromosomes) + 1L]] <- list(
                id = chr_id,
                name = chr_name,
                label = chr_name,
                subgenome = level,
                group = chr_group,
                length = chr_length,
                deleted = chr_id %in% deleted,
                added = FALSE,
                copy = 1L,
                copies = 1L,
                selected = chr_group %in% selected_groups,
                scores = if (is.null(score_values[[chr_group]])) {
                    numeric()
                } else {
                    score_values[[chr_group]]
                },
                calls = recipient_calls(
                    chr_id, chr_group, level, chr_length
                )
            )
        }
        for (i in seq_len(nrow(add))) {
            chr_id <- as.character(add$id[[i]])
            chr_name <- as.character(add$name[[i]])
            chr_group <- as.character(add$group[[i]])
            chr_length <- as.numeric(add$length[[i]])
            copies <- suppressWarnings(as.integer(add$copy_number[[i]]))
            if (length(copies) != 1L || is.na(copies) || copies < 1L) {
                copies <- 1L
            }
            for (copy in seq_len(copies)) {
                label <- if (copies > 1L) {
                    paste(chr_name, copy)
                } else {
                    chr_name
                }
                chromosomes[[length(chromosomes) + 1L]] <- list(
                    id = chr_id,
                    name = chr_name,
                    label = label,
                    subgenome = level,
                    group = chr_group,
                    length = chr_length,
                    deleted = FALSE,
                    added = TRUE,
                    copy = copy,
                    copies = copies,
                    selected = chr_group %in% selected_groups,
                    scores = if (is.null(score_values[[chr_group]])) {
                        numeric()
                    } else {
                        score_values[[chr_group]]
                    },
                    calls = donor_calls(
                        chr_id, chr_group, level, chr_length
                    )
                )
            }
        }
        list(
            subgenome = level,
            label = paste("Subgenome", level),
            selected = any(vapply(
                chromosomes,
                function(x) isTRUE(x$selected),
                logical(1)
            )),
            chromosomes = chromosomes
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
        subgenomes = subgenome_payload
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
