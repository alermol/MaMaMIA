## Interactive vertical karyotype for the introgression explorer.
##
## The widget shows the recipient karyotype only. A recipient chromosome is
## kept in place when it has no coverage, drawn as a dashed "deleted" bar. A
## full-chromosome donor addition adds that donor chromosome (one shape per
## homolog copy) as another diploid pair. The payload is plain lists so
## htmlwidgets can serialise it without a second data-frame layer.

karyotype_palette <- function() {
    list(
        negative = "#2166ac",
        zero = "#f1efe7",
        positive = "#b2182b",
        substitution = "#d95f02",
        addition = "#1f78b4",
        deletion = "#e31a1c",
        complex = "#e6ab02",
        other = "#7570b3",
        deleted = "#d7d3c6",
        added = "#9cc3d5",
        grid = "#e8e4da",
        text = "#3a3d35"
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
    full <- karyotype_full_calls(result, sample, n_per)
    if (is.null(full) || nrow(full) == 0L) {
        full <- data.frame(
            group = character(),
            donor_mean = numeric(),
            stringsAsFactors = FALSE
        )
    }
    deleted <- karyotype_deleted(
        result, sample,
        deletion_breadth = deletion_breadth,
        deletion_relative = deletion_relative
    )
    ## A full balanced substitution replaces the recipient chromosome. Show
    ## that as the loss of the recipient copy plus an added donor chromosome,
    ## so the plot does not caption the donor signal as the recipient.
    full_substitution <- if ("class" %in% names(full)) {
        unique(as.character(full$group[
            as.character(full$class) == "substitution"
        ]))
    } else {
        character()
    }
    if (length(full_substitution) > 0L) {
        substituted_ids <- as.character(recipient$id[
            as.character(recipient$group) %in% full_substitution
        ])
        deleted <- union(deleted, substituted_ids)
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
    ## column. This lets the layout place chromosomes in config order and
    ## scale every chromosome bar by its configured length.
    recipient_ids <- if ("recipient_chr_id" %in% names(overlay)) {
        as.character(overlay$recipient_chr_id)
    } else {
        rep(NA_character_, nrow(overlay))
    }
    recipient_subs <- if ("recipient_component" %in% names(overlay)) {
        as.character(overlay$recipient_component)
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
    donor_subs <- if ("donor_component" %in% names(overlay)) {
        as.character(overlay$donor_component)
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
            class_group = as.character(
                call_value("class_group", i, NA_character_)
            ),
            confidence = as.character(
                call_value("confidence", i, NA_character_)
            ),
            score = as.numeric(call_value("score_mean", i, NA_real_)),
            donor_mean = as.numeric(call_value("donor_mean", i, NA_real_)),
            recipient_mean = as.numeric(
                call_value("recipient_mean", i, NA_real_)
            ),
            pass = isTRUE(as.logical(call_value("pass", i, TRUE))),
            recipient_component = as.character(
                call_value("recipient_component", i, NA_character_)
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
            recipient_run_start_bp = as.numeric(
                call_value("recipient_run_start_bp", i, NA_real_)
            ),
            recipient_run_end_bp = as.numeric(
                call_value("recipient_run_end_bp", i, NA_real_)
            ),
            recipient_deletion_width_bp = as.numeric(
                call_value("recipient_deletion_width_bp", i, NA_real_)
            ),
            recipient_assignment = as.character(
                call_value("recipient_assignment", i, NA_character_)
            ),
            recipient_weight = as.numeric(
                call_value("recipient_weight", i, NA_real_)
            ),
            recipient_margin = as.numeric(
                call_value("recipient_margin", i, NA_real_)
            ),
            recipient_group = as.character(
                call_value("recipient_group", i, NA_character_)
            ),
            donor_component = as.character(
                call_value("donor_component", i, NA_character_)
            ),
            donor_chr_id = as.character(
                call_value("donor_chr_id", i, NA_character_)
            ),
            donor_chr_name = as.character(
                call_value("donor_chr_name", i, NA_character_)
            ),
            donor_start_bp = as.numeric(
                call_value("donor_start_bp", i, NA_real_)
            ),
            donor_end_bp = as.numeric(
                call_value("donor_end_bp", i, NA_real_)
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
            ),
            class_evidence = as.character(
                call_value("class_evidence", i, NA_character_)
            ),
            class_confidence = as.character(
                call_value("class_confidence", i, NA_character_)
            ),
            class_status = as.character(
                call_value("class_status", i, NA_character_)
            ),
            copy_state = as.character(
                call_value("copy_state", i, NA_character_)
            ),
            total_evidence = as.character(
                call_value("total_evidence", i, NA_character_)
            ),
            component_status = as.character(
                call_value("component_status", i, NA_character_)
            ),
            mappability_status = as.character(
                call_value("mappability_status", i, NA_character_)
            ),
            mapq_fraction = as.numeric(
                call_value("mapq_fraction", i, NA_real_)
            ),
            total_mean = as.numeric(
                call_value("total_mean", i, NA_real_)
            ),
            donor_copies = as.numeric(
                call_value("donor_copies", i, NA_real_)
            ),
            recipient_depletion_se = as.numeric(
                call_value("recipient_depletion_se", i, NA_real_)
            ),
            recipient_width_overlap = as.numeric(
                call_value("recipient_width_overlap", i, NA_real_)
            ),
            recipient_fraction_ratio = as.numeric(
                call_value("recipient_fraction_ratio", i, NA_real_)
            ),
            recipient_copy_number = as.numeric(
                call_value("recipient_copy_number", i, NA_real_)
            ),
            affected_fraction = as.numeric(
                call_value("affected_fraction", i, NA_real_)
            ),
            affected_copies = as.numeric(
                call_value("affected_copies", i, NA_real_)
            ),
            dosage_confidence = as.character(
                call_value("dosage_confidence", i, NA_character_)
            )
        )
    }
    recipient_calls <- function(chr_id, group, component, length_bp) {
        if (nrow(overlay) == 0L || chr_id %in% deleted) {
            return(list())
        }
        own <- !is.na(recipient_ids) & recipient_ids == chr_id
        fallback <- is.na(recipient_ids) & !is.na(recipient_subs) &
            recipient_subs == component & call_groups == group
        keep <- is.na(recipient_assignment) |
            recipient_assignment != "no_recipient_loss"
        rows <- which((own | fallback) & keep)
        lapply(rows, function(i) {
            make_call(i, recipient_start[[i]], recipient_end[[i]],
                length_bp, n_per
            )
        })
    }
    donor_calls <- function(chr_id, group, component, length_bp) {
        if (nrow(overlay) == 0L) {
            return(list())
        }
        own <- !is.na(donor_ids) & donor_ids == chr_id
        fallback <- is.na(donor_ids) & !is.na(donor_subs) &
            donor_subs == component & call_groups == group
        rows <- which(own | fallback)
        lapply(rows, function(i) {
            make_call(i, donor_start[[i]], donor_end[[i]],
                length_bp, n_per
            )
        })
    }

    added_groups <- if (nrow(full) > 0L) {
        keep <- !is.na(full$donor_mean) & full$donor_mean > 0.5
        if ("class" %in% names(full)) {
            keep <- keep & as.character(full$class) %in%
                c("addition", "substitution")
        }
        unique(as.character(full$group[keep]))
    } else {
        character()
    }
    if (is.null(donor) || nrow(donor) == 0L) {
        donor <- data.frame(
            id = character(), name = character(), component = character(),
            group = character(), length = numeric(),
            copy_number = numeric(), stringsAsFactors = FALSE
        )
    }
    ## Preserve the registry order for donor chromosomes as well as for the
    ## recipient chromosomes above.
    added_donor <- donor[
        as.character(donor$group) %in% added_groups, , drop = FALSE
    ]
    chromosome_payload <- list()
    for (i in seq_len(nrow(recipient))) {
        chr_id <- as.character(recipient$id[[i]])
        chr_name <- as.character(recipient$name[[i]])
        chr_component <- as.character(recipient$component[[i]])
        chr_group <- as.character(recipient$group[[i]])
        chr_length <- as.numeric(recipient$length[[i]])
        copies <- if ("copy_number" %in% names(recipient)) {
            suppressWarnings(as.integer(recipient$copy_number[[i]]))
        } else {
            1L
        }
        if (length(copies) != 1L || is.na(copies) || copies < 1L) {
            copies <- 1L
        }
        chr_calls <- recipient_calls(
            chr_id, chr_group, chr_component, chr_length
        )
        for (copy in seq_len(copies)) {
            chromosome_payload[[length(chromosome_payload) + 1L]] <- list(
                id = chr_id,
                name = chr_name,
                label = chr_name,
                component = chr_component,
                group = chr_group,
                length = chr_length,
                deleted = chr_id %in% deleted,
                added = FALSE,
                copy = copy,
                copies = copies,
                selected = chr_group %in% selected_groups,
                scores = if (is.null(score_values[[chr_group]])) {
                    numeric()
                } else {
                    score_values[[chr_group]]
                },
                calls = chr_calls
            )
        }
    }
    for (i in seq_len(nrow(added_donor))) {
        chr_id <- as.character(added_donor$id[[i]])
        chr_name <- as.character(added_donor$name[[i]])
        chr_component <- as.character(added_donor$component[[i]])
        chr_group <- as.character(added_donor$group[[i]])
        chr_length <- as.numeric(added_donor$length[[i]])
        copies <- suppressWarnings(as.integer(added_donor$copy_number[[i]]))
        if (length(copies) != 1L || is.na(copies) || copies < 1L) {
            copies <- 1L
        }
        chr_calls <- donor_calls(
            chr_id, chr_group, chr_component, chr_length
        )
        for (copy in seq_len(copies)) {
            chromosome_payload[[length(chromosome_payload) + 1L]] <- list(
                id = chr_id,
                name = chr_name,
                label = chr_name,
                component = chr_component,
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
                calls = chr_calls
            )
        }
    }

    ## Keep each added donor chromosome next to its recipient group instead of
    ## appending it after every recipient chromosome. Donor groups without a
    ## recipient contig stay at the end of the layout.
    if (length(chromosome_payload) > 0L) {
        payload_group <- vapply(
            chromosome_payload,
            function(x) as.character(x$group),
            character(1)
        )
        payload_added <- vapply(
            chromosome_payload,
            function(x) isTRUE(x$added),
            logical(1)
        )
        group_order <- unique(as.character(recipient$group))
        group_order <- c(
            group_order,
            setdiff(unique(payload_group), group_order)
        )
        chromosome_payload <- unlist(lapply(group_order, function(group) {
            entries <- chromosome_payload[payload_group == group]
            added <- payload_added[payload_group == group]
            c(entries[!added], entries[added])
        }), recursive = FALSE)
    }

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
        chromosomes = chromosome_payload
    )
}

#' Create the karyotype htmlwidget
#' @noRd
mamamia_karyotype <- function(payload, width = NULL, height = NULL) {
    css <- htmltools::htmlDependency(
        name = "mamamiaKaryotypeCard-css",
        version = "0.1.0",
        src = c(file = system.file("htmlwidgets", package = "MaMaMIA")),
        stylesheet = "mamamiaKaryotypeCard.css"
    )
    ## htmlwidgets names its default binding dependency after the widget and
    ## versions it with the package version. Declaring the same name here
    ## keeps a single script URL; raise this asset version for a widget-only
    ## change so the browser cannot keep serving the previous file.
    js <- htmltools::htmlDependency(
        name = "mamamiaKaryotypeCard-binding",
        version = "0.1.0",
        src = "htmlwidgets",
        package = "MaMaMIA",
        script = "mamamiaKaryotypeCard.js",
        all_files = FALSE
    )
    htmlwidgets::createWidget(
        name = "mamamiaKaryotypeCard",
        x = payload,
        width = width,
        height = height,
        package = "MaMaMIA",
        dependencies = list(css, js),
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
        outputId, "mamamiaKaryotypeCard", width, height,
        package = "MaMaMIA"
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
