## Interactive vertical karyotype for the introgression explorer.
##
## The widget shows the recipient karyotype only. A recipient chromosome is
## kept in place when it has no coverage, drawn as a dashed "deleted" bar. A
## full-chromosome donor addition adds that donor chromosome as a single bar.
## Every chromosome is drawn once, without diploid copy pairs. The payload is
## plain lists so htmlwidgets can serialise it without a second data-frame layer.

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
#' window and covers at least 95% of the donor chromosome. A substitution is
#' judged on the recipient loss instead when that measurement is available, so
#' a terminus that the donor read stops short of does not hide the replacement.
#' @noRd
karyotype_full_calls <- function(result, sample, n_per = NULL,
                                 calls = NULL) {
    if (is.null(calls)) {
        calls <- result$calls
    }
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
    score_length <- tapply(
        scores$length, as.character(scores$group),
        function(x) {
            x <- x[is.finite(x)]
            if (length(x) == 0L) NA_real_ else max(x)
        }
    )
    ## `width_bp` comes from the donor windows (segmentScores() uses the donor
    ## group length), so the full-chromosome fraction must use the same
    ## coordinate system. Dividing by the recipient group maximum makes a full
    ## donor chromosome look partial whenever the donor group is shorter.
    donor_length <- if (!is.null(result$windows) &&
        is.data.frame(result$windows$table) &&
        nrow(result$windows$table) > 0L &&
        all(c("role", "group", "stop") %in%
            names(result$windows$table))) {
        donor_windows <- result$windows$table[
            as.character(result$windows$table$role) == "donor", ,
            drop = FALSE
        ]
        tapply(
            as.numeric(donor_windows$stop),
            as.character(donor_windows$group),
            function(x) {
                x <- x[is.finite(x)]
                if (length(x) == 0L) NA_real_ else max(x)
            }
        )
    } else {
        numeric()
    }
    length_bp <- donor_length[as.character(calls$group)]
    missing_length <- !is.finite(length_bp)
    if (any(missing_length)) {
        length_bp[missing_length] <- score_length[
            as.character(calls$group)[missing_length]
        ]
    }
    donor_width_fraction <- calls$width_bp / unname(length_bp)

    ## A substitution removes a recipient chromosome. When assignment measured
    ## that loss, use its fraction of the assigned recipient contig: the donor
    ## segment can stop a few percent short at a low-mappability terminus while
    ## the recipient chromosome is still fully replaced. Additions and older
    ## calls without recipient coordinates keep the donor fraction.
    recipient_length_bp <- rep(NA_real_, nrow(calls))
    recipient_contigs <- result$assembly$contigs$recipient
    if (is.data.frame(recipient_contigs) &&
        nrow(recipient_contigs) > 0L &&
        all(c("id", "length") %in% names(recipient_contigs))) {
        if ("recipient_chr_id" %in% names(calls)) {
            position <- match(
                as.character(calls$recipient_chr_id),
                as.character(recipient_contigs$id)
            )
            recipient_length_bp <- as.numeric(
                recipient_contigs$length[position]
            )
        }
        if (anyNA(recipient_length_bp) &&
            "recipient_chr_name" %in% names(calls) &&
            "group" %in% names(recipient_contigs) &&
            "name" %in% names(recipient_contigs)) {
            missing <- which(is.na(recipient_length_bp))
            for (i in missing) {
                name <- as.character(calls$recipient_chr_name[[i]])
                group <- as.character(calls$group[[i]])
                if (length(name) != 1L || is.na(name) || !nzchar(name)) {
                    next
                }
                hit <- which(
                    as.character(recipient_contigs$name) == name &
                        as.character(recipient_contigs$group) == group
                )
                if (length(hit) == 1L) {
                    recipient_length_bp[[i]] <-
                        as.numeric(recipient_contigs$length[[hit]])
                }
            }
        }
    }
    missing_recipient <- !is.finite(recipient_length_bp)
    if (any(missing_recipient)) {
        recipient_length_bp[missing_recipient] <- score_length[
            as.character(calls$group)[missing_recipient]
        ]
    }
    recipient_fraction <- rep(NA_real_, nrow(calls))
    if ("recipient_deletion_width_bp" %in% names(calls)) {
        recipient_fraction <- as.numeric(calls$recipient_deletion_width_bp) /
            unname(recipient_length_bp)
    }
    class <- if ("class" %in% names(calls)) {
        as.character(calls$class)
    } else {
        rep(NA_character_, nrow(calls))
    }
    use_recipient <- class == "substitution" & is.finite(recipient_fraction)
    fraction <- donor_width_fraction
    fraction[use_recipient] <- recipient_fraction[use_recipient]

    ## A full-chromosome event starts at the first window and covers at least
    ## 95% of the chromosome. The end is allowed to fall a few windows short:
    ## real terminal segments (for example the 7S addition in the study data)
    ## end at 99% of the profile rather than exactly at the last window.
    full <- calls$start_iid == 0L & is.finite(fraction) & fraction >= 0.95
    calls$group_length <- unname(length_bp)
    calls$fraction <- fraction
    calls$full <- full
    calls[full, , drop = FALSE]
}

#' Donor chromosomes covered by passed call segments
#'
#' A chromosome-scale addition can be split into several adjacent calls when
#' the score profile dips between them. No single call then starts at the first
#' window and reaches the coverage threshold, so the full-chromosome test in
#' [karyotype_full_calls()] misses the event and the donor chromosome is never
#' added to the plot. Measure the union of the passed donor intervals instead.
#'
#' @param calls Calls to inspect, already filtered to one sample.
#' @param donor Donor contig table with `id` and `length`.
#' @param minimum Fraction of the donor chromosome that must be covered.
#' @param start_max Largest allowed start position, as a fraction of the
#'   chromosome. `NULL` allows a segment anywhere on the chromosome.
#' @param classes Call classes to include. Defaults to the two classes that
#'   add donor chromosome material.
#' @noRd
karyotype_covered_donor_ids <- function(calls, donor, minimum = 0.95,
                                        start_max = 0.01,
                                        classes = c(
                                            "addition", "substitution"
                                        )) {
    if (is.null(calls) || !is.data.frame(calls) || nrow(calls) == 0L ||
        is.null(donor) || !is.data.frame(donor) || nrow(donor) == 0L ||
        !all(c("donor_chr_id", "donor_start_bp", "donor_end_bp") %in%
            names(calls)) ||
        !all(c("id", "length") %in% names(donor))) {
        return(character())
    }
    keep <- rep(TRUE, nrow(calls))
    if ("pass" %in% names(calls)) {
        keep <- keep & isTRUE_vec(calls$pass)
    }
    if ("class" %in% names(calls)) {
        keep <- keep & as.character(calls$class) %in% classes
    }
    if ("donor_mean" %in% names(calls)) {
        keep <- keep & is.finite(as.numeric(calls$donor_mean)) &
            as.numeric(calls$donor_mean) > 0.5
    }
    ids <- as.character(calls$donor_chr_id)
    starts <- suppressWarnings(as.numeric(calls$donor_start_bp))
    ends <- suppressWarnings(as.numeric(calls$donor_end_bp))
    keep <- keep & !is.na(ids) & nzchar(ids) &
        is.finite(starts) & is.finite(ends) & ends > starts
    if (!any(keep)) {
        return(character())
    }
    position <- match(ids[keep], as.character(donor$id))
    length_bp <- suppressWarnings(as.numeric(donor$length[position]))
    usable <- is.finite(length_bp) & length_bp > 0
    if (!any(usable)) {
        return(character())
    }
    ids <- ids[keep][usable]
    starts <- pmax(0, starts[keep][usable]) / length_bp[usable]
    ends <- pmin(length_bp[usable], ends[keep][usable]) / length_bp[usable]
    covered <- character()
    for (id in unique(ids)) {
        at <- which(ids == id)
        from <- starts[at]
        to <- ends[at]
        order_at <- order(from, to)
        from <- from[order_at]
        to <- to[order_at]
        span_start <- from[[1L]]
        span_end <- to[[1L]]
        total <- 0
        if (length(from) > 1L) {
            for (i in 2:length(from)) {
                if (from[[i]] <= span_end) {
                    span_end <- max(span_end, to[[i]])
                } else {
                    total <- total + span_end - span_start
                    span_start <- from[[i]]
                    span_end <- to[[i]]
                }
            }
        }
        total <- total + span_end - span_start
        starts_near_zero <- is.null(start_max) || min(from) <= start_max
        if (total >= minimum && starts_near_zero) {
            covered <- c(covered, id)
        }
    }
    covered
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
    overlay <- if (is.null(calls)) {
        result$calls
    } else {
        calls
    }
    if (!is.null(overlay)) {
        overlay <- overlay[
            as.character(overlay$sample) == sample, , drop = FALSE
        ]
    }
    full <- karyotype_full_calls(
        result, sample, n_per, calls = overlay
    )
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
    ## A full balanced substitution replaces one recipient chromosome. Show
    ## that as the loss of that recipient copy plus an added donor chromosome,
    ## so the plot does not caption the donor signal as the recipient.
    ##
    ## A group can contain several homoeologs (7A, 7B and 7D), and a call is
    ## assigned to only one of them. Marking every contig in the group would
    ## dash 7A and 7B when only 7D was replaced.
    full_substitution <- if ("class" %in% names(full)) {
        full[as.character(full$class) == "substitution", , drop = FALSE]
    } else {
        full[0, , drop = FALSE]
    }
    if (nrow(full_substitution) > 0L) {
        substituted_ids <- character()
        for (i in seq_len(nrow(full_substitution))) {
            group <- if ("recipient_group" %in%
                names(full_substitution)) {
                as.character(full_substitution$recipient_group[[i]])
            } else {
                NA_character_
            }
            if (length(group) != 1L || is.na(group) || !nzchar(group)) {
                group <- as.character(full_substitution$group[[i]])
            }
            candidates <- character()
            if ("recipient_chr_id" %in% names(full_substitution)) {
                id <- as.character(
                    full_substitution$recipient_chr_id[[i]]
                )
                if (length(id) == 1L && !is.na(id) && nzchar(id) &&
                    id %in% as.character(recipient$id)) {
                    candidates <- id
                }
            }
            if (length(candidates) == 0L &&
                "recipient_chr_name" %in% names(full_substitution)) {
                name <- as.character(
                    full_substitution$recipient_chr_name[[i]]
                )
                if (length(name) == 1L && !is.na(name) && nzchar(name)) {
                    candidates <- as.character(recipient$id[
                        as.character(recipient$name) == name &
                            as.character(recipient$group) == group
                    ])
                }
            }
            if (length(candidates) == 0L &&
                "recipient_component" %in% names(full_substitution) &&
                "component" %in% names(recipient)) {
                component <- as.character(
                    full_substitution$recipient_component[[i]]
                )
                if (length(component) == 1L && !is.na(component) &&
                    nzchar(component)) {
                    candidates <- as.character(recipient$id[
                        as.character(recipient$component) == component &
                            as.character(recipient$group) == group
                    ])
                }
            }
            if (length(candidates) == 0L) {
                group_ids <- as.character(recipient$id[
                    as.character(recipient$group) == group
                ])
                ## A single-contig group is unambiguous; a multi-contig group
                ## without an assigned recipient must not lose all homoeologs.
                if (length(group_ids) == 1L) {
                    candidates <- group_ids
                }
            }
            substituted_ids <- union(substituted_ids, candidates)
        }
        deleted <- union(deleted, substituted_ids)
    }
    if (is.null(overlay)) {
        overlay <- full[0, , drop = FALSE]
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
    ## A deletion or a recipient-only event has no donor segment. The donor
    ## contig fields still describe the homologous donor interval used for
    ## placement, but the tooltip must not present that as donor evidence.
    call_donor_supported <- function(i) {
        evidence <- as.character(
            call_value("class_evidence", i, NA_character_)
        )
        if (length(evidence) == 1L && !is.na(evidence) &&
            nzchar(evidence)) {
            return(grepl("donor", evidence, fixed = TRUE))
        }
        as.character(call_value("class", i, "")) %in%
            c("substitution", "addition", "complex")
    }
    ## The reciprocal case: an addition has donor support but no recipient
    ## loss. The recipient contig fields are the homologous host interval used
    ## for placement, not an observed recipient segment.
    call_recipient_supported <- function(i) {
        evidence <- as.character(
            call_value("class_evidence", i, NA_character_)
        )
        if (length(evidence) == 1L && !is.na(evidence) &&
            nzchar(evidence)) {
            return(grepl("recipient", evidence, fixed = TRUE))
        }
        as.character(call_value("class", i, "")) %in%
            c("substitution", "deletion", "complex")
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
            donor_supported = call_donor_supported(i),
            recipient_supported = call_recipient_supported(i),
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
            recipient_depletion_se = as.numeric(
                call_value("recipient_depletion_se", i, NA_real_)
            ),
            recipient_width_overlap = as.numeric(
                call_value("recipient_width_overlap", i, NA_real_)
            ),
            recipient_fraction_ratio = as.numeric(
                call_value("recipient_fraction_ratio", i, NA_real_)
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
        ## A partial addition that has its own donor bar is drawn there.
        ## Smaller additions stay on the recipient bar.
        rows <- which((own | fallback) & keep)
        if (length(rows) > 0L && length(partial_donor_ids) > 0L) {
            pass <- if ("pass" %in% names(overlay)) {
                isTRUE_vec(overlay$pass[rows])
            } else {
                rep(TRUE, length(rows))
            }
            moved <- as.character(overlay$class[rows]) == "addition" &
                pass & !is.na(donor_ids[rows]) &
                donor_ids[rows] %in% partial_donor_ids
            rows <- rows[!moved]
        }
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
        supported <- vapply(
            seq_len(nrow(overlay)), call_donor_supported, logical(1)
        )
        rows <- which((own | fallback) & supported)
        lapply(rows, function(i) {
            make_call(i, donor_start[[i]], donor_end[[i]],
                length_bp, n_per
            )
        })
    }

    if (is.null(donor) || nrow(donor) == 0L) {
        donor <- data.frame(
            id = character(), name = character(), component = character(),
            group = character(), length = numeric(),
            copy_number = numeric(), stringsAsFactors = FALSE
        )
    }
    added_calls <- if (nrow(full) > 0L) {
        keep <- !is.na(full$donor_mean) & full$donor_mean > 0.5
        if ("class" %in% names(full)) {
            keep <- keep & as.character(full$class) %in%
                c("addition", "substitution")
        }
        full[keep, , drop = FALSE]
    } else {
        full[0, , drop = FALSE]
    }
    ## Add only the donor chromosome named by the call. A donor group can hold
    ## several components (for example J02, St02 and V02); adding every contig
    ## in the group would invent the components the caller did not observe.
    full_donor_ids <- character()
    for (i in seq_len(nrow(added_calls))) {
        group <- as.character(added_calls$group[[i]])
        candidates <- character()
        if ("donor_chr_id" %in% names(added_calls)) {
            id <- as.character(added_calls$donor_chr_id[[i]])
            if (length(id) == 1L && !is.na(id) && nzchar(id) &&
                id %in% as.character(donor$id)) {
                candidates <- id
            }
        }
        if (length(candidates) == 0L &&
            "donor_chr_name" %in% names(added_calls)) {
            name <- as.character(added_calls$donor_chr_name[[i]])
            if (length(name) == 1L && !is.na(name) && nzchar(name)) {
                candidates <- as.character(donor$id[
                    as.character(donor$name) == name &
                        as.character(donor$group) == group
                ])
            }
        }
        if (length(candidates) == 0L &&
            "donor_component" %in% names(added_calls)) {
            component <- as.character(added_calls$donor_component[[i]])
            if (length(component) == 1L && !is.na(component) &&
                nzchar(component)) {
                candidates <- as.character(donor$id[
                    as.character(donor$component) == component &
                        as.character(donor$group) == group
                ])
            }
        }
        if (length(candidates) == 0L) {
            candidates <- as.character(donor$id[
                as.character(donor$group) == group
            ])
        }
        full_donor_ids <- union(full_donor_ids, candidates)
    }
    ## A chromosome-scale addition that was split into adjacent calls never
    ## satisfies the single-call full-chromosome test. Add the donor
    ## chromosome when the union of its passed segments covers it.
    full_donor_ids <- union(
        full_donor_ids,
        karyotype_covered_donor_ids(overlay, donor)
    )
    ## A partial addition gets its own donor chromosome once its passed
    ## segments cover at least 20% of that donor, whether or not the call also
    ## has a recipient placement. Only the added part is filled; the rest of
    ## the donor bar is drawn as a dashed donor deletion.
    partial_donor_ids <- setdiff(
        karyotype_covered_donor_ids(
            overlay, donor,
            minimum = 0.2, start_max = NULL,
            classes = "addition"
        ),
        full_donor_ids
    )
    added_donor_ids <- union(full_donor_ids, partial_donor_ids)
    ## Preserve the registry order for donor chromosomes as well as for the
    ## recipient chromosomes above.
    added_donor <- donor[
        as.character(donor$id) %in% added_donor_ids, , drop = FALSE
    ]
    chromosome_payload <- list()
    for (i in seq_len(nrow(recipient))) {
        chr_id <- as.character(recipient$id[[i]])
        chr_name <- as.character(recipient$name[[i]])
        chr_component <- as.character(recipient$component[[i]])
        chr_group <- as.character(recipient$group[[i]])
        chr_length <- as.numeric(recipient$length[[i]])
        chr_calls <- recipient_calls(
            chr_id, chr_group, chr_component, chr_length
        )
        chromosome_payload[[length(chromosome_payload) + 1L]] <- list(
            id = chr_id,
            name = chr_name,
            label = chr_name,
            component = chr_component,
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
            calls = chr_calls
        )
    }
    for (i in seq_len(nrow(added_donor))) {
        chr_id <- as.character(added_donor$id[[i]])
        chr_name <- as.character(added_donor$name[[i]])
        chr_component <- as.character(added_donor$component[[i]])
        chr_group <- as.character(added_donor$group[[i]])
        chr_length <- as.numeric(added_donor$length[[i]])
        chr_calls <- donor_calls(
            chr_id, chr_group, chr_component, chr_length
        )
        chromosome_payload[[length(chromosome_payload) + 1L]] <- list(
            id = chr_id,
            name = chr_name,
            label = chr_name,
            component = chr_component,
            group = chr_group,
            length = chr_length,
            deleted = FALSE,
            added = TRUE,
            added_partial = chr_id %in% partial_donor_ids,
            copy = 1L,
            copies = 1L,
            selected = chr_group %in% selected_groups,
            scores = if (is.null(score_values[[chr_group]])) {
                numeric()
            } else {
                score_values[[chr_group]]
            },
            calls = chr_calls
        )
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

#' Asset dependencies for the karyotype widget
#' @noRd
mamamia_karyotype_dependencies <- function() {
    src <- c(file = mamamia_path("www", "htmlwidgets"))
    list(
        htmltools::htmlDependency(
            name = "mamamiaKaryotypeCard-css",
            version = "0.1.0",
            src = src,
            stylesheet = "mamamiaKaryotypeCard.css"
        ),
        ## htmlwidgets names its default binding dependency after the widget
        ## and versions it with the application version. Declaring the same
        ## name here keeps a single script URL; raise this asset version for a
        ## widget-only change so the browser cannot keep serving the previous
        ## file.
        htmltools::htmlDependency(
            name = "mamamiaKaryotypeCard-binding",
            version = "0.1.1",
            src = src,
            script = "mamamiaKaryotypeCard.js",
            all_files = FALSE
        )
    )
}

#' Create the karyotype htmlwidget
#' @noRd
mamamia_karyotype <- function(payload, width = NULL, height = NULL) {
    htmlwidgets::createWidget(
        name = "mamamiaKaryotypeCard",
        x = payload,
        width = width,
        height = height,
        dependencies = mamamia_karyotype_dependencies(),
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
    tag <- htmlwidgets::shinyWidgetOutput(
        outputId, "mamamiaKaryotypeCard", width, height
    )
    ## `shinyWidgetOutput()` already carries the htmlwidgets runtime
    ## dependency. Append the widget's own CSS/binding so the binding script
    ## loads after `HTMLWidgets` exists; replacing the dependencies here makes
    ## the browser evaluate `HTMLWidgets.widget()` before htmlwidgets.js.
    htmltools::attachDependencies(
        tag, mamamia_karyotype_dependencies(), append = TRUE
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
