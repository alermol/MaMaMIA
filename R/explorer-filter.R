## Pure explorer filtering.
##
## The Shiny explorer adapts `input` into the plain `filters` list accepted
## here. Keeping the selection rules free of session state means the table,
## karyotype payloads, BEDPE selection and bundle payload all share exactly
## one implementation.

#' Group a call is displayed under
#'
#' A cross-group call is shown under the recipient group it was assigned to;
#' every other call keeps its own group. The recipient columns are optional,
#' so older call tables without them display under `group`.
#'
#' @param calls A call data frame, or `NULL`.
#' @return A character vector as long as `calls`.
#' @noRd
explorer_display_group <- function(calls) {
    if (is.null(calls)) {
        return(character())
    }
    groups <- as.character(calls$group)
    if (all(c("recipient_group", "recipient_cross_group") %in%
        names(calls))) {
        cross <- as.logical(calls$recipient_cross_group) %in% TRUE
        recipient <- as.character(calls$recipient_group)
        move <- cross & !is.na(recipient) & nzchar(recipient)
        groups[move] <- recipient[move]
    }
    groups
}

#' Calls with no recipient host
#'
#' A call that ends with `recipient_assignment = "no_recipient_loss"` has no
#' recipient interval, so the karyotype cannot place it. It remains a visible
#' call for review and export unless the explorer hides unplaced calls.
#'
#' @param calls A call data frame, or `NULL`.
#' @return A logical vector as long as `calls`.
#' @noRd
explorer_unplaced <- function(calls) {
    if (is.null(calls) || nrow(calls) == 0L) {
        return(logical())
    }
    if (!"recipient_assignment" %in% names(calls)) {
        return(rep(FALSE, nrow(calls)))
    }
    assignment <- as.character(calls$recipient_assignment)
    !is.na(assignment) & assignment == "no_recipient_loss"
}

#' Visible explorer calls
#'
#' Applies the explorer filters in one place. Every supplied filter follows
#' the same precedence as the dashboard: pass-only first, then group
#' membership, call class, class status, copy state, confidence, minimum
#' width, and finally hiding unplaced calls. `NULL` and zero-row calls pass
#' through unchanged.
#'
#' @param calls A call data frame, or `NULL`.
#' @param filters A plain list with any of:
#'   `pass_only` (logical), `groups` (character), `class`, `class_status`,
#'   `copy_state`, `confidence` (character), `min_width_bp` (numeric,
#'   already in base pairs) and `hide_unplaced` (logical). Missing or empty
#'   entries do not filter.
#' @return A list with `calls` (the visible rows), `display_group` (the
#'   display group for each visible row) and `unplaced` (how many rows that
#'   passed the other filters have no recipient host, including rows hidden
#'   by `hide_unplaced`).
#' @noRd
explorer_calls <- function(calls, filters = list()) {
    if (is.null(calls) || nrow(calls) == 0L) {
        return(list(
            calls = calls,
            display_group = explorer_display_group(calls),
            unplaced = 0L
        ))
    }
    if (is.null(filters)) {
        filters <- list()
    }
    keep <- rep(TRUE, nrow(calls))
    if (isTRUE(filters[["pass_only"]]) && "pass" %in% names(calls)) {
        keep <- keep & isTRUE_vec(calls$pass)
    }
    groups <- filters[["groups"]]
    if (length(groups) > 0L) {
        keep <- keep & (
            as.character(calls$group) %in% groups |
                explorer_display_group(calls) %in% groups
        )
    }
    if (length(filters[["class"]]) > 0L) {
        keep <- keep & calls$class %in% filters[["class"]]
    }
    if (length(filters[["class_status"]]) > 0L) {
        keep <- keep & calls$class_status %in% filters[["class_status"]]
    }
    if (length(filters[["copy_state"]]) > 0L &&
        "copy_state" %in% names(calls)) {
        keep <- keep & calls$copy_state %in% filters[["copy_state"]]
    }
    if (length(filters[["confidence"]]) > 0L) {
        keep <- keep & calls$confidence %in% filters[["confidence"]]
    }
    width <- suppressWarnings(as.numeric(filters[["min_width_bp"]]))
    if (length(width) == 1L && is.finite(width) && width > 0) {
        keep <- keep & calls$width_bp >= width
    }
    unplaced <- explorer_unplaced(calls)
    unplaced_count <- sum(keep & unplaced, na.rm = TRUE)
    if (isTRUE(filters[["hide_unplaced"]])) {
        keep <- keep & !unplaced
    }
    list(
        calls = calls[keep, , drop = FALSE],
        display_group = explorer_display_group(calls)[keep],
        unplaced = as.integer(unplaced_count)
    )
}
