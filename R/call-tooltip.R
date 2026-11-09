## Resolved tooltip projection for one introgression call.
##
## The call schema declares the fields, their display labels and their
## descriptions. This module turns one call record into the display model the
## browser tooltip needs: header text, status, subtitle, sections, rows,
## metrics and chips. Formatting lives here so the JavaScript adapter renders
## whatever it is given and never inspects canonical call fields itself.

#' Display text for a scalar call value
#' @noRd
call_tooltip_text <- function(value) {
    if (is.null(value) || length(value) != 1L || is.na(value)) {
        return("")
    }
    as.character(value)
}

#' Whether a resolved tooltip value should be shown
#' @noRd
call_tooltip_has <- function(value) {
    text <- call_tooltip_text(value)
    nzchar(text) && !identical(text, "NA")
}

#' Schema label for one field
#' @noRd
call_tooltip_label <- function(name, fallback = name) {
    labels <- call_schema_labels()
    value <- labels[[name]]
    if (is.null(value) || length(value) != 1L || is.na(value) ||
        !nzchar(value)) {
        fallback
    } else {
        as.character(value)
    }
}

#' Short schema label for a metric such as `score_mean`
#' @noRd
call_tooltip_metric_label <- function(name) {
    sub(" mean$", "", call_tooltip_label(name))
}

#' Format a signed finite number like the browser tooltip did
#' @noRd
call_tooltip_signed <- function(value, digits = 2L) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) {
        return("NA")
    }
    prefix <- if (value > 0) "+" else ""
    paste0(prefix, formatC(value, format = "f", digits = digits))
}

#' Format a base-pair count
#' @noRd
call_tooltip_bp <- function(value) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value) || value < 0) {
        return("NA")
    }
    if (value >= 1e6) {
        return(paste0(formatC(value / 1e6, format = "f", digits = 1L), " Mb"))
    }
    if (value >= 1e3) {
        return(paste0(formatC(value / 1e3, format = "f", digits = 1L), " kb"))
    }
    paste0(formatC(round(value), format = "f", digits = 0L), " bp")
}

#' Format a base-pair span
#' @noRd
call_tooltip_span <- function(start, end) {
    start <- suppressWarnings(as.numeric(start))
    end <- suppressWarnings(as.numeric(end))
    start_ok <- length(start) == 1L && is.finite(start)
    end_ok <- length(end) == 1L && is.finite(end)
    if (start_ok && end_ok) {
        return(paste0(call_tooltip_bp(start), " \u2013 ", call_tooltip_bp(end)))
    }
    if (start_ok) {
        return(paste0("from ", call_tooltip_bp(start)))
    }
    if (end_ok) {
        return(paste0("to ", call_tooltip_bp(end)))
    }
    ""
}

#' Format a fraction as a percentage
#' @noRd
call_tooltip_percent <- function(value) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) {
        return("NA")
    }
    paste0(formatC(round(100 * value), format = "f", digits = 0L), "%")
}

#' Sign class used by the tooltip metric styling
#' @noRd
call_tooltip_sign_class <- function(value) {
    value <- suppressWarnings(as.numeric(value))
    if (length(value) != 1L || !is.finite(value)) {
        return("muted")
    }
    if (value > 0) {
        return("pos")
    }
    if (value < 0) {
        return("neg")
    }
    "muted"
}

#' Resolved status text and CSS class for a call
#' @noRd
call_tooltip_status <- function(call) {
    status <- call_tooltip_text(call$class_status)
    if (!nzchar(status)) {
        confidence <- call_tooltip_text(call$class_confidence)
        status <- if (identical(confidence, "high")) {
            "confident"
        } else if (identical(confidence, "moderate")) {
            "ambiguous"
        } else if (identical(confidence, "low")) {
            "needs_review"
        } else {
            "call"
        }
    }
    key <- tolower(gsub("_", " ", status))
    css <- if (identical(key, "confident")) {
        "mamamia-tip-status-confident"
    } else if (identical(key, "ambiguous")) {
        "mamamia-tip-status-ambiguous"
    } else if (key %in% c("needs review", "needs_review")) {
        "mamamia-tip-status-review"
    } else {
        "mamamia-tip-status-neutral"
    }
    list(
        text = gsub("_", " ", status),
        class = css
    )
}

#' Donor support flag resolved from the call record
#' @noRd
call_tooltip_donor_supported <- function(call) {
    value <- call$donor_supported
    if (is.logical(value) && length(value) == 1L && !is.na(value)) {
        return(isTRUE(value))
    }
    evidence <- call_tooltip_text(call$class_evidence)
    if (nzchar(evidence)) {
        return(grepl("donor", evidence, fixed = TRUE))
    }
    call_tooltip_text(call$class) %in%
        c("substitution", "addition", "complex")
}

#' Recipient support flag resolved from the call record
#' @noRd
call_tooltip_recipient_supported <- function(call) {
    value <- call$recipient_supported
    if (is.logical(value) && length(value) == 1L && !is.na(value)) {
        return(isTRUE(value))
    }
    evidence <- call_tooltip_text(call$class_evidence)
    if (nzchar(evidence)) {
        return(grepl("recipient", evidence, fixed = TRUE))
    }
    call_tooltip_text(call$class) %in%
        c("substitution", "deletion", "complex")
}

#' Display text for the placement outcome
#' @noRd
call_tooltip_assignment <- function(value) {
    key <- call_tooltip_text(value)
    if (identical(key, "assigned")) {
        return("assigned")
    }
    if (identical(key, "cross_group")) {
        return("cross-group candidate")
    }
    if (identical(key, "ambiguous")) {
        return("ambiguous")
    }
    if (identical(key, "no_recipient_loss")) {
        return("no host")
    }
    if (nzchar(key)) key else "group only"
}

#' One tooltip row
#' @noRd
call_tooltip_row <- function(label, value) {
    if (!call_tooltip_has(value)) {
        return(NULL)
    }
    list(
        label = as.character(label),
        value = as.character(value)
    )
}

#' Remove empty tooltip rows
#' @noRd
call_tooltip_rows <- function(rows) {
    Filter(Negate(is.null), rows)
}

#' Resolve the tooltip model for one call
#'
#' @param call A karyotype call record from the payload.
#' @param sample Optional sample name for the header.
#' @param chromosome_group Optional group of the chromosome bar carrying the
#'   call. This is the donor group for a donor bar and the recipient group for
#'   a recipient bar.
#' @param palette Optional palette, defaults to the application palette.
#'
#' @return A list with `accent`, `title`, `status`, `subtitle` and `sections`.
#'   Every string is already formatted for display.
#' @noRd
call_tooltip <- function(call, sample = NULL, chromosome_group = NULL,
                         palette = karyotype_palette()) {
    if (!is.list(call)) {
        call <- list()
    }
    class <- call_tooltip_text(call$class)
    if (!nzchar(class)) {
        class <- "call"
    }
    accent <- palette[[class]]
    if (is.null(accent) || length(accent) != 1L || is.na(accent) ||
        !nzchar(accent)) {
        accent <- "#6c757d"
    }
    donor_supported <- call_tooltip_donor_supported(call)
    recipient_supported <- call_tooltip_recipient_supported(call)

    subtitle <- character()
    sample <- call_tooltip_text(sample)
    if (nzchar(sample)) {
        subtitle <- c(subtitle, sample)
    }
    chromosome_group <- call_tooltip_text(chromosome_group)
    if (nzchar(chromosome_group)) {
        subtitle <- c(subtitle, paste0("group ", chromosome_group))
    }
    span <- if (donor_supported) {
        call_tooltip_span(call$donor_start_bp, call$donor_end_bp)
    } else {
        ""
    }
    span_label <- "donor "
    if (!nzchar(span)) {
        span <- call_tooltip_span(call$start_bp, call$end_bp)
        span_label <- ""
    }
    if (nzchar(span)) {
        subtitle <- c(subtitle, paste0(span_label, span))
    }

    placement <- list()
    if (recipient_supported) {
        recipient <- call_tooltip_text(call$recipient_chr_name)
        if (!nzchar(recipient)) {
            recipient <- if (identical(
                call_tooltip_text(call$recipient_assignment),
                "no_recipient_loss"
            )) {
                "no host"
            } else {
                "unassigned"
            }
        }
        recipient_group <- call_tooltip_text(call$recipient_group)
        if (nzchar(recipient_group)) {
            recipient <- paste0(recipient, " \u00b7 group ", recipient_group)
        }
        recipient_span <- call_tooltip_span(
            call$recipient_start_bp, call$recipient_end_bp
        )
        if (nzchar(recipient_span)) {
            recipient <- paste0(recipient, " \u00b7 ", recipient_span)
        }
        placement <- c(placement, list(call_tooltip_row(
            call_tooltip_label("recipient_chr_name", "Recipient"), recipient
        )))
    }
    if (donor_supported) {
        donor_name <- call_tooltip_text(call$donor_chr_name)
        donor_group <- call_tooltip_text(call$donor_group)
        donor <- if (nzchar(donor_name)) {
            donor_name
        } else if (nzchar(donor_group)) {
            paste0("group ", donor_group)
        } else {
            "unknown"
        }
        if (nzchar(donor_name) && nzchar(donor_group)) {
            donor <- paste0(donor, " \u00b7 group ", donor_group)
        }
        donor_span <- call_tooltip_span(call$donor_start_bp, call$donor_end_bp)
        if (nzchar(donor_span)) {
            donor <- paste0(donor, " \u00b7 ", donor_span)
        }
        placement <- c(placement, list(call_tooltip_row(
            call_tooltip_label("donor_chr_name", "Donor"), donor
        )))
    }
    if (recipient_supported && donor_supported) {
        relation <- if (identical(
            call_tooltip_text(call$recipient_assignment),
            "no_recipient_loss"
        )) {
            "not placed"
        } else if (isTRUE(call$cross_group)) {
            "cross-group"
        } else {
            "same group"
        }
        placement <- c(placement, list(call_tooltip_row("Relation", relation)))
    }
    if (recipient_supported) {
        assignment <- call_tooltip_assignment(call$recipient_assignment)
        margin <- suppressWarnings(as.numeric(call$recipient_margin))
        if (length(margin) == 1L && is.finite(margin)) {
            assignment <- paste0(
                assignment, " \u00b7 margin ", call_tooltip_signed(margin)
            )
        }
        placement <- c(placement, list(call_tooltip_row(
            "Assignment", assignment
        )))
    }
    class_group <- call_tooltip_text(call$class_group)
    if (nzchar(class_group) && !identical(class_group, class)) {
        placement <- c(placement, list(call_tooltip_row(
            call_tooltip_label("class_group", "Group class"),
            paste0(class_group, " \u2192 ", class)
        )))
    }
    if (recipient_supported &&
        identical(call_tooltip_text(call$recipient_assignment), "ambiguous") &&
        call_tooltip_has(call$recipient_alternative)) {
        placement <- c(placement, list(call_tooltip_row(
            "Alternative",
            call_tooltip_text(call$recipient_alternative)
        )))
    }
    placement <- call_tooltip_rows(placement)

    metric_fields <- c(
        "score_mean", "donor_mean", "recipient_mean", "total_mean"
    )
    metric_values <- list(
        score_mean = call$score,
        donor_mean = call$donor_mean,
        recipient_mean = call$recipient_mean,
        total_mean = call$total_mean
    )
    metrics <- lapply(metric_fields, function(name) {
        value <- metric_values[[name]]
        list(
            label = call_tooltip_metric_label(name),
            value = call_tooltip_signed(value),
            class = call_tooltip_sign_class(value)
        )
    })

    chips <- list()
    chip <- function(text, kind) {
        if (!call_tooltip_has(text)) {
            return(NULL)
        }
        list(text = as.character(text), kind = as.character(kind))
    }
    confidence <- call_tooltip_text(call$class_confidence)
    if (nzchar(confidence)) {
        chips <- c(chips, list(chip(
            paste("class confidence", confidence), "confidence"
        )))
    }
    evidence <- call_tooltip_text(call$class_evidence)
    if (nzchar(evidence)) {
        chips <- c(chips, list(chip(paste("evidence", evidence), "evidence")))
    }
    copy_state <- call_tooltip_text(call$copy_state)
    if (nzchar(copy_state)) {
        chips <- c(chips, list(chip(paste("copy", copy_state), "copy")))
    }
    total_evidence <- call_tooltip_text(call$total_evidence)
    if (nzchar(total_evidence) && !identical(total_evidence, "neutral")) {
        chips <- c(chips, list(chip(paste("total", total_evidence), "copy")))
    }
    if (identical(call_tooltip_text(call$component_status), "heterogeneous")) {
        chips <- c(chips, list(chip("components heterogeneous", "warn")))
    }
    if (identical(call_tooltip_text(call$mappability_status), "low")) {
        chips <- c(chips, list(chip("low mappability", "warn")))
    }
    if (identical(call$pass, FALSE)) {
        chips <- c(chips, list(chip("not passing", "warn")))
    }
    chips <- Filter(Negate(is.null), chips)

    geometry <- list()
    width <- suppressWarnings(as.numeric(call$width_bp))
    if (length(width) == 1L && is.finite(width)) {
        geometry <- c(geometry, list(call_tooltip_row(
            if (donor_supported) "Donor width" else "Call width",
            call_tooltip_bp(width)
        )))
    }
    recipient_width <- suppressWarnings(
        as.numeric(call$recipient_deletion_width_bp)
    )
    if (recipient_supported && length(recipient_width) == 1L &&
        is.finite(recipient_width)) {
        geometry <- c(geometry, list(call_tooltip_row(
            "Recipient run", call_tooltip_bp(recipient_width)
        )))
    }
    overlap <- suppressWarnings(as.numeric(call$recipient_width_overlap))
    if (donor_supported && recipient_supported &&
        length(overlap) == 1L && is.finite(overlap)) {
        geometry <- c(geometry, list(call_tooltip_row(
            "Overlap", call_tooltip_percent(overlap)
        )))
    }
    ratio <- suppressWarnings(as.numeric(call$recipient_fraction_ratio))
    if (donor_supported && recipient_supported &&
        length(ratio) == 1L && is.finite(ratio)) {
        geometry <- c(geometry, list(call_tooltip_row(
            call_tooltip_label("recipient_fraction_ratio", "Width ratio"),
            formatC(ratio, format = "f", digits = 2L)
        )))
    }
    geometry <- call_tooltip_rows(geometry)

    placement_section <- list(
        kind = "rows",
        label = "Placement",
        rows = placement
    )
    if (length(placement) == 0L) {
        placement_section$empty <-
            "No donor or recipient segment assigned."
    }
    evidence_section <- list(
        kind = "evidence",
        label = "Evidence",
        metrics = metrics,
        chips = chips
    )
    sections <- list(placement_section, evidence_section)
    if (length(geometry) > 0L) {
        sections <- c(sections, list(list(
            kind = "rows",
            label = "Geometry",
            rows = geometry
        )))
    }

    list(
        accent = as.character(accent),
        title = class,
        status = call_tooltip_status(call),
        subtitle = paste(subtitle, collapse = " \u00b7 "),
        sections = sections
    )
}
