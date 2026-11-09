## Shared call-status rules.
##
## Component consistency and mappability are group-level diagnostics that both
## `classifyCalls()` and `assignCallPositions()` apply to a call. Declaring the
## rule once keeps the two stages from drifting apart.

#' Component status for a call
#'
#' A large donor or recipient component spread, or more than one supported
#' donor component, marks a group-level call as heterogeneous rather than
#' failing it. Non-finite inputs count as zero.
#'
#' @param donor_component_spread,recipient_component_spread Log-depth spread
#'   across donor or recipient components.
#' @param n_donor_supported Number of donor components that support the call.
#' @param spread_cutoff Spread above which the call is heterogeneous.
#' @noRd
call_component_status <- function(donor_component_spread = NA_real_,
                                  recipient_component_spread = NA_real_,
                                  n_donor_supported = NA_real_,
                                  spread_cutoff = 2.0) {
    clean <- function(x) {
        values <- as.numeric(x)
        values[!is.finite(values)] <- 0
        values
    }
    donor <- clean(donor_component_spread)
    recipient <- clean(recipient_component_spread)
    supported <- clean(n_donor_supported)
    ifelse(
        donor > spread_cutoff |
            recipient > spread_cutoff |
            supported > 1,
        "heterogeneous", "consistent"
    )
}

#' Component status from a calls table or a single call
#' @noRd
call_component_status_from <- function(calls, spread_cutoff = 2.0) {
    n <- if (is.data.frame(calls)) nrow(calls) else 1L
    column <- function(name) {
        if (is.data.frame(calls) && name %in% names(calls)) {
            calls[[name]]
        } else if (!is.data.frame(calls) && name %in% names(calls)) {
            calls[[name]]
        } else {
            rep(0, n)
        }
    }
    call_component_status(
        column("donor_component_spread"),
        column("recipient_component_spread"),
        column("n_donor_supported"),
        spread_cutoff = spread_cutoff
    )
}

#' Mappability status for a call
#'
#' `low` when the high-MAPQ read fraction is below `cutoff`, `ok` when it is
#' known and at least the cutoff, and `NA` when it is unavailable.
#' @noRd
call_mappability_status <- function(mapq_fraction, cutoff = 0.5) {
    values <- as.numeric(mapq_fraction)
    ifelse(
        !is.finite(values),
        NA_character_,
        ifelse(values < cutoff, "low", "ok")
    )
}

#' Mappability status from a calls table
#' @noRd
call_mappability_status_from <- function(calls, cutoff = 0.5) {
    values <- if (is.data.frame(calls)) {
        if ("mapq_fraction" %in% names(calls)) {
            calls$mapq_fraction
        } else {
            rep(NA_real_, nrow(calls))
        }
    } else {
        calls$mapq_fraction %||% NA_real_
    }
    call_mappability_status(values, cutoff = cutoff)
}

#' Whether a call's class confidence should be downgraded
#'
#' A heterogeneous component pattern or low mappability turns a high class
#' confidence into moderate. The returned logical vector is `FALSE` when the
#' confidence is not high.
#' @noRd
call_status_downgrade <- function(component_status, mappability_status,
                                  class_confidence) {
    high <- !is.na(class_confidence) & class_confidence == "high"
    heterogeneous <- !is.na(component_status) &
        component_status == "heterogeneous"
    low_mappability <- !is.na(mappability_status) &
        mappability_status == "low"
    high & (heterogeneous | low_mappability)
}
