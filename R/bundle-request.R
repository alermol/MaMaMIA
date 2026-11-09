## One description of a result-bundle download.
##
## The dashboard, the bundle cache and the archive writer all depend on the
## same inputs. Building this request once keeps the cache key and the archive
## from drifting apart.

#' Analysis identity carried by a result
#'
#' New results carry a pipeline-generated `provenance$analysis_id`. Older
#' results fall back to the provenance fields that identify one run.
#' @noRd
bundle_analysis_id <- function(result) {
    provenance <- if (is.list(result)) result$provenance else NULL
    if (is.list(provenance)) {
        id <- provenance$analysis_id
        if (!is.null(id) && length(id) == 1L && !is.na(id) &&
            nzchar(id)) {
            return(as.character(id))
        }
        fields <- c(
            "started_at", "duration_secs", "assembly", "n_samples",
            "calibration_version", "mask_version", "config_hash"
        )
        fallback <- vapply(fields, function(name) {
            value <- provenance[[name]]
            if (is.null(value) || length(value) == 0L) {
                ""
            } else {
                as.character(value[[1L]])
            }
        }, character(1))
        if (any(nzchar(fallback))) {
            return(paste(fallback, collapse = "|"))
        }
    }
    assembly <- if (is.list(result$assembly)) result$assembly$key else NULL
    paste(
        as.character(assembly %||% ""),
        nrow(result$calls %||% data.frame()),
        paste(unique(as.character(result$samples$sample)), collapse = ","),
        sep = "|"
    )
}

#' Build a result-bundle request
#'
#' @param result A `BamAnalysisResult`.
#' @param sample Sample to export.
#' @param mode `"heatmap"` or `"calls"`.
#' @param calls Rendered calls carrying `.call_id`.
#' @param groups Selected chromosome groups.
#' @param hide_deletions, hide_unplaced Current view switches.
#' @param current_payload Optional current karyotype payload, built once by
#'   the dashboard and reused on a cache miss.
#'
#' @return A request list with key material and the current payload.
#' @noRd
bundle_request <- function(result, sample, mode, calls = NULL,
                           groups = NULL, hide_deletions = FALSE,
                           hide_unplaced = FALSE,
                           current_payload = NULL) {
    if (!inherits(result, "BamAnalysisResult")) {
        stop("`result` must be a BamAnalysisResult", call. = FALSE)
    }
    if (!is.character(sample) || length(sample) != 1L || is.na(sample) ||
        !nzchar(sample)) {
        stop("`sample` must be a single sample name", call. = FALSE)
    }
    mode <- match.arg(mode, c("heatmap", "calls"))
    if (!is.null(calls) && nrow(calls) > 0L &&
        !".call_id" %in% names(calls)) {
        calls$.call_id <- seq_len(nrow(calls))
    }
    list(
        result = result,
        analysis_id = bundle_analysis_id(result),
        bundle_version = bam_bundle_version(),
        sample = as.character(sample),
        mode = mode,
        calls = calls,
        groups = as.character(groups %||% character()),
        hide_deletions = isTRUE(hide_deletions),
        hide_unplaced = isTRUE(hide_unplaced),
        current_payload = current_payload
    )
}

#' Key that identifies one generated result bundle
#'
#' The key is a structured list, not a delimiter-joined string. Rendering
#' inputs cannot collide when a value happens to contain a separator.
#' @noRd
bundle_request_key <- function(request) {
    if (!is.list(request)) {
        stop("`request` must come from bundle_request()", call. = FALSE)
    }
    calls <- request$calls
    call_ids <- if (is.null(calls) || nrow(calls) == 0L) {
        integer()
    } else if (".call_id" %in% names(calls)) {
        sort(unique(as.integer(calls$.call_id)))
    } else {
        integer()
    }
    list(
        analysis_id = as.character(request$analysis_id %||% ""),
        bundle_version = as.character(
            request$bundle_version %||% ""
        ),
        sample = as.character(request$sample %||% ""),
        mode = as.character(request$mode %||% ""),
        call_ids = call_ids,
        groups = sort(unique(as.character(
            request$groups %||% character()
        ))),
        hide_deletions = isTRUE(request$hide_deletions),
        hide_unplaced = isTRUE(request$hide_unplaced)
    )
}
