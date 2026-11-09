## Backend-neutral display scene for the karyotype payload.
##
## The payload carries data. This module resolves the display policy once so
## the browser widget and the server-side static SVG cannot drift apart:
## legend entries, chromosome base styles, call colours/opacity/stroke/dash,
## labels, heatmap bins and the omitted-segment note all come from one place.
## Adapters only turn these primitives into SVG/HTML nodes.

#' Merge a payload palette with the application defaults
#' @noRd
karyotype_scene_palette <- function(payload = NULL) {
    palette <- karyotype_palette()
    supplied <- if (is.list(payload)) payload$palette else NULL
    if (is.list(supplied)) {
        for (name in names(supplied)) {
            value <- supplied[[name]]
            if (is.character(value) && length(value) == 1L &&
                !is.na(value) && nzchar(value)) {
                palette[[name]] <- value
            }
        }
    }
    palette
}

#' RGB components of a hex colour
#' @noRd
karyotype_scene_rgb <- function(hex) {
    hex <- as.character(hex)
    if (length(hex) != 1L || is.na(hex) || !nzchar(hex)) {
        return(c(238, 238, 238))
    }
    value <- sub("^#", "", hex)
    if (nchar(value) == 3L) {
        value <- paste0(
            substr(value, 1L, 1L), substr(value, 1L, 1L),
            substr(value, 2L, 2L), substr(value, 2L, 2L),
            substr(value, 3L, 3L), substr(value, 3L, 3L)
        )
    }
    if (nchar(value) != 6L) {
        return(c(238, 238, 238))
    }
    component <- function(from) {
        parsed <- suppressWarnings(strtoi(substr(value, from, from + 1L),
            base = 16L
        ))
        if (is.na(parsed)) 238L else parsed
    }
    c(component(1L), component(3L), component(5L))
}

#' Mix two colours like the browser renderer
#' @noRd
karyotype_scene_mix <- function(from, to, t) {
    a <- karyotype_scene_rgb(from)
    b <- karyotype_scene_rgb(to)
    t <- max(0, min(1, as.numeric(t)[[1L]]))
    mixed <- round(a + (b - a) * t)
    sprintf("rgb(%d,%d,%d)", mixed[[1L]], mixed[[2L]], mixed[[3L]])
}

#' Heatmap colour for a score
#' @noRd
karyotype_scene_score_colour <- function(value, limit, palette,
                                         missing = "#eeeeee") {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) {
        return(missing)
    }
    limit <- suppressWarnings(as.numeric(limit))
    if (length(limit) != 1L || !is.finite(limit) || limit <= 0) {
        limit <- 0.5
    }
    t <- max(-1, min(1, value / limit))
    if (t < 0) {
        return(karyotype_scene_mix(palette$zero, palette$negative, -t))
    }
    karyotype_scene_mix(palette$zero, palette$positive, t)
}

#' Call legend entries resolved from the payload policy
#' @noRd
karyotype_scene_legend_entries <- function(palette, hide_deletions) {
    entries <- list(
        list(key = "substitution", label = "Substitution"),
        list(key = "addition", label = "Addition"),
        list(key = "deletion", label = "Deletion"),
        list(key = "complex", label = "Complex"),
        list(key = "other", label = "Other"),
        list(key = "deleted", label = "Deleted recipient"),
        list(key = "added", label = "Added donor"),
        list(key = "unplaced", label = "Unplaced donor")
    )
    if (isTRUE(hide_deletions)) {
        entries <- Filter(
            function(entry) {
                !entry$key %in% c("deletion", "deleted")
            },
            entries
        )
    }
    lapply(entries, function(entry) {
        deleted <- identical(entry$key, "deleted")
        list(
            key = entry$key,
            label = entry$label,
            fill = if (is.null(palette[[entry$key]])) {
                palette$other
            } else {
                palette[[entry$key]]
            },
            stroke = if (deleted) palette$negative else palette$text,
            stroke_width = if (deleted) 1.2 else 0.6,
            dash = if (deleted) "3 2" else NULL
        )
    })
}

#' Heatmap legend model
#' @noRd
karyotype_scene_heatmap_legend <- function(palette, limit) {
    list(
        negative = palette$negative,
        zero = palette$zero,
        positive = palette$positive,
        score_limit = limit,
        labels = list(
            min = sprintf("-%.1f", limit),
            zero = "0",
            max = sprintf("+%.1f", limit)
        )
    )
}

#' Merge overlapping call intervals for heatmap coverage
#' @noRd
karyotype_scene_intervals <- function(calls) {
    if (length(calls) == 0L) {
        return(list())
    }
    starts <- vapply(calls, function(call) {
        value <- suppressWarnings(as.numeric(call$start))
        if (length(value) != 1L || !is.finite(value)) 0 else value
    }, numeric(1))
    ends <- vapply(calls, function(call) {
        value <- suppressWarnings(as.numeric(call$end))
        if (length(value) != 1L || !is.finite(value)) 1 else value
    }, numeric(1))
    starts <- pmax(0, pmin(1, starts))
    ends <- pmax(0, pmin(1, ends))
    keep <- ends > starts
    if (!any(keep)) {
        return(list())
    }
    order_at <- order(starts[keep], ends[keep])
    list(start = starts[keep][order_at], end = ends[keep][order_at])
}

#' Fixed heatmap bin policy
#'
#' The browser and the static SVG must bin the same score profile the same
#' way. Both consume the bins already resolved in the payload scene; this
#' count is the one policy that decides their resolution, independent of the
#' rendered bar height.
#' @noRd
karyotype_scene_heatmap_bin_count <- function() {
    120L
}

#' Heatmap bins for one chromosome
#'
#' A bin is kept only when a non-cross-group call covers its centre. Bins
#' carry colour so the browser and the static SVG use the same scale.
#' @noRd
karyotype_scene_bins <- function(scores, calls, bins, palette, limit) {
    scores <- suppressWarnings(as.numeric(scores))
    scores <- scores[seq_len(length(scores))]
    if (length(scores) == 0L || length(calls) == 0L) {
        return(list())
    }
    overlaps <- Filter(
        function(call) !isTRUE(call$cross_group),
        calls
    )
    intervals <- karyotype_scene_intervals(overlaps)
    if (length(intervals$start) == 0L) {
        return(list())
    }
    n_bins <- suppressWarnings(as.integer(bins))
    if (length(n_bins) != 1L || is.na(n_bins) || n_bins < 1L) {
        n_bins <- karyotype_scene_heatmap_bin_count()
    }
    n_bins <- min(length(scores), max(40L, n_bins))
    bin_height_fraction <- 1 / n_bins
    starts <- intervals$start
    ends <- intervals$end
    pointer <- 1L
    active_end <- -Inf
    out <- list()
    for (i in seq_len(n_bins)) {
        centre <- (i - 0.5) / n_bins
        while (pointer <= length(starts) && starts[[pointer]] <= centre) {
            active_end <- max(active_end, ends[[pointer]])
            pointer <- pointer + 1L
        }
        if (centre > active_end) {
            next
        }
        from <- floor((i - 1L) * length(scores) / n_bins) + 1L
        to <- max(from, floor(i * length(scores) / n_bins))
        values <- scores[from:to]
        values <- values[is.finite(values)]
        value <- if (length(values) == 0L) NA_real_ else mean(values)
        out[[length(out) + 1L]] <- list(
            index = i - 1L,
            start = (i - 1L) * bin_height_fraction,
            end = i * bin_height_fraction,
            fill = karyotype_scene_score_colour(value, limit, palette)
        )
    }
    out
}

#' Resolved style for one call segment
#' @noRd
karyotype_scene_call <- function(call, index, palette) {
    confidence <- call$confidence
    if (is.null(confidence) || length(confidence) != 1L ||
        is.na(confidence) || !nzchar(as.character(confidence))) {
        confidence <- call$class_confidence
    }
    confidence <- as.character(confidence)
    opacity <- if (identical(confidence, "high")) {
        0.95
    } else if (identical(confidence, "moderate")) {
        0.72
    } else {
        0.45
    }
    if (identical(call$pass, FALSE)) {
        opacity <- opacity * 0.45
    }
    assignment <- as.character(call$recipient_assignment)
    if (identical(assignment, "ambiguous")) {
        opacity <- opacity * 0.55
    }
    if (identical(assignment, "no_recipient_loss")) {
        opacity <- opacity * 0.35
    }
    cross <- isTRUE(call$cross_group)
    class <- as.character(call$class)
    start <- suppressWarnings(as.numeric(call$start))
    end <- suppressWarnings(as.numeric(call$end))
    if (length(start) != 1L || !is.finite(start)) start <- 0
    if (length(end) != 1L || !is.finite(end)) end <- 1
    start <- max(0, min(1, start))
    end <- max(0, min(1, end))
    status <- as.character(call$class_status)
    if (length(status) != 1L || is.na(status)) {
        status <- as.character(call$confidence)
    }
    if (length(status) != 1L || is.na(status)) {
        status <- ""
    }
    list(
        index = as.integer(index),
        start = start,
        end = end,
        class = class,
        class_status = status,
        score = suppressWarnings(as.numeric(call$score)),
        start_bp = suppressWarnings(as.numeric(call$start_bp)),
        end_bp = suppressWarnings(as.numeric(call$end_bp)),
        cross_group = cross,
        fill = if (is.null(palette[[class]])) palette$other else palette[[class]],
        opacity = opacity,
        stroke = if (cross) "#333333" else "#ffffff",
        stroke_width = if (cross) 1.3 else 0.5,
        dash = if (cross) "3 2" else NULL
    )
}

#' First non-empty label from a chromosome
#' @noRd
karyotype_scene_label <- function(chromosome) {
    for (name in c("label", "name", "component", "id")) {
        value <- chromosome[[name]]
        if (is.character(value) && length(value) >= 1L &&
            !is.na(value[[1L]]) && nzchar(value[[1L]])) {
            return(as.character(value[[1L]]))
        }
    }
    ""
}

#' Resolve one chromosome's geometry and styles
#' @noRd
karyotype_scene_chromosome <- function(chromosome, mode, heatmap_bins,
                                       palette, limit) {
    added <- isTRUE(chromosome$added)
    added_partial <- added && isTRUE(chromosome$added_partial)
    unplaced <- isTRUE(chromosome$unplaced)
    deleted <- !added && !unplaced && isTRUE(chromosome$deleted)
    calls <- chromosome$calls
    if (!is.list(calls)) {
        calls <- list()
    }
    length_bp <- suppressWarnings(as.numeric(chromosome$length))
    if (length(length_bp) != 1L || !is.finite(length_bp)) length_bp <- 0
    display_length <- suppressWarnings(as.numeric(chromosome$display_length))
    if (length(display_length) != 1L || !is.finite(display_length) ||
        display_length <= 0) {
        display_length <- if (length_bp > 0) length_bp else 0
    }
    ## An added donor chromosome is a chromosome state. In calls mode its base
    ## uses the added-donor colour; in heatmap mode it uses the zero-coverage
    ## colour so regions without score bins read as missing rather than as a
    ## state block. Calls drawn on it keep their class colour. Using the sole
    ## call's class here painted the whole bar, including bases the call never
    ## covered, in the call class colour.
    added_fill <- as.character(palette$added %||% "#9cc3d5")[[1L]]
    added_text <- as.character(palette$added_text %||% added_fill)[[1L]]
    unplaced_fill <- as.character(palette$unplaced %||% "#bdbdbd")[[1L]]
    unplaced_text <- as.character(
        palette$unplaced_text %||% unplaced_fill
    )[[1L]]
    zero_fill <- as.character(palette$zero %||% "#f1efe7")[[1L]]
    heatmap <- identical(mode, "heatmap")
    base_fill <- if (deleted) {
        palette$deleted
    } else if (unplaced) {
        if (heatmap) zero_fill else unplaced_fill
    } else if (added) {
        if (heatmap) {
            zero_fill
        } else if (added_partial) {
            "none"
        } else {
            added_fill
        }
    } else {
        "#f4f4f4"
    }
    base_stroke <- if (deleted) {
        palette$negative
    } else if (unplaced) {
        unplaced_text
    } else if (added) {
        if (added_partial) palette$negative else added_text
    } else {
        "#666666"
    }
    base_stroke_width <- if (deleted) {
        1.4
    } else if (unplaced) {
        1.4
    } else if (added) {
        if (added_partial) 1.4 else 1.8
    } else {
        0.7
    }
    base_dash <- if (deleted || added_partial || unplaced) {
        if (unplaced) "3 2" else "4 3"
    } else {
        NULL
    }
    copies <- suppressWarnings(as.integer(chromosome$copies))
    if (length(copies) != 1L || is.na(copies) || copies < 1L) copies <- 1L
    copy <- suppressWarnings(as.integer(chromosome$copy))
    if (length(copy) != 1L || is.na(copy) || copy < 1L) copy <- 1L
    caption <- if (deleted && copy == 1L) {
        list(text = "Deleted", fill = palette$negative)
    } else if (unplaced && copy == 1L) {
        list(text = "Unplaced", fill = unplaced_text)
    } else if (added && copy == 1L) {
        list(text = "Added", fill = added_text)
    } else {
        NULL
    }
    list(
        id = as.character(chromosome$id %||% "")[[1L]],
        name = as.character(chromosome$name %||% "")[[1L]],
        label = karyotype_scene_label(chromosome),
        component = as.character(chromosome$component %||% "")[[1L]],
        group = as.character(chromosome$group %||% "")[[1L]],
        added = added,
        added_partial = added_partial,
        unplaced = unplaced,
        deleted = deleted,
        selected = !identical(chromosome$selected, FALSE),
        copy = copy,
        copies = copies,
        length = length_bp,
        display_length = display_length,
        base = list(
            fill = base_fill,
            stroke = base_stroke,
            stroke_width = base_stroke_width,
            dash = base_dash
        ),
        caption = caption,
        show_label = copy == 1L || copies == 1L,
        calls = lapply(seq_along(calls), function(i) {
            karyotype_scene_call(calls[[i]], i, palette)
        }),
        heatmap_bins = if (identical(mode, "heatmap") && !deleted) {
            karyotype_scene_bins(
                chromosome$scores, calls, heatmap_bins, palette, limit
            )
        } else {
            list()
        }
    )
}

#' Build the backend-neutral karyotype scene
#'
#' @param payload A karyotype payload from [karyotype_payload()].
#' @param heatmap_bins Number of score bins. Defaults to the shared
#'   [karyotype_scene_heatmap_bin_count()] policy so every caller that renders
#'   the heatmap uses the same resolution.
#'
#' @return A list with `mode`, `palette`, `legend`, `omitted_note` and one
#'   resolved entry per chromosome.
#' @noRd
karyotype_scene <- function(payload,
                            heatmap_bins = karyotype_scene_heatmap_bin_count()) {
    if (!is.list(payload)) {
        payload <- list()
    }
    palette <- karyotype_scene_palette(payload)
    mode <- as.character(payload$mode %||% "calls")[[1L]]
    if (!mode %in% c("heatmap", "calls")) {
        mode <- "calls"
    }
    hide_deletions <- isTRUE(payload$hide_deletions)
    chromosomes <- payload$chromosomes
    if (!is.list(chromosomes)) {
        chromosomes <- list()
    }
    limit <- suppressWarnings(as.numeric(payload$score_limit))
    if (length(limit) != 1L || !is.finite(limit) || limit <= 0) {
        limit <- 0.5
    }
    hidden_segments <- suppressWarnings(as.integer(
        payload$hidden_deletions$segments
    ))
    if (length(hidden_segments) != 1L || is.na(hidden_segments)) {
        hidden_segments <- 0L
    }
    omitted_note <- if (hide_deletions && hidden_segments > 0L) {
        list(
            segments = hidden_segments,
            text = paste0(
                hidden_segments, " segment",
                if (hidden_segments == 1L) "" else "s", " omitted"
            )
        )
    } else {
        NULL
    }
    list(
        mode = mode,
        sample = as.character(payload$sample %||% "")[[1L]],
        score_limit = limit,
        hide_deletions = hide_deletions,
        palette = palette,
        omitted_note = omitted_note,
        legend = list(
            kind = mode,
            entries = if (identical(mode, "calls")) {
                karyotype_scene_legend_entries(palette, hide_deletions)
            } else {
                list()
            },
            heatmap = karyotype_scene_heatmap_legend(palette, limit)
        ),
        heatmap = list(
            missing = "#eeeeee",
            negative = palette$negative,
            zero = palette$zero,
            positive = palette$positive,
            score_limit = limit
        ),
        chromosomes = lapply(chromosomes, karyotype_scene_chromosome,
            mode = mode,
            heatmap_bins = heatmap_bins,
            palette = palette,
            limit = limit
        )
    )
}
