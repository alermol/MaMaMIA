## Per-call component evidence.
##
## Group-level calls answer "where is the introgression?"; this table answers
## "which donor component or components support that group?" It is diagnostic:
## it does not change the group-level call or the calibrated pass flag.

.component_evidence_empty <- function() {
    data.frame(
        sample = character(),
        call_number = integer(),
        call_id = character(),
        chromosome_group = character(),
        call_start_window = integer(),
        call_end_window = integer(),
        call_class = character(),
        call_pass = logical(),
        donor_component = character(),
        component_rank = integer(),
        component_selected = logical(),
        n_component_contigs = integer(),
        n_windows = integer(),
        mean_component_lambda = numeric(),
        mean_recipient_lambda = numeric(),
        mean_component_score = numeric(),
        median_component_score = numeric(),
        max_component_score = numeric(),
        supported_fraction = numeric(),
        component_supported = logical(),
        score_type = character(),
        metric = character(),
        stringsAsFactors = FALSE
    )
}

.component_evidence_lambda <- function(metrics, metric, eps) {
    if (identical(metric, "breadth")) {
        return(as.numeric(metrics$lambda_breadth))
    }
    if (identical(metric, "count")) {
        return(as.numeric(metrics$lambda_count))
    }
    breadth <- pmax(as.numeric(metrics$lambda_breadth), 0)
    count <- pmax(as.numeric(metrics$lambda_count), 0)
    sqrt((breadth + eps) * (count + eps))
}

#' Build per-call component evidence
#'
#' @param result A `BamAnalysisResult`.
#' @param sample Sample name. Defaults to the first sample.
#' @param calls Optional call rows; defaults to `result$calls` for the sample.
#' @param support_threshold Component-score threshold used for the supported
#'   fraction. The default is a log-scale threshold of 0.5.
#'
#' @return A data frame with one row per call and donor component. This is a
#'   diagnostic table; it does not alter the primary group-level calls.
#' @noRd
componentEvidence <- function(result, sample = NULL, calls = NULL,
                              support_threshold = 0.5) {
    if (!inherits(result, "BamAnalysisResult")) {
        .mamamia_error(
            "MaMaMIAComponentEvidenceError",
            "`result` must be a BamAnalysisResult"
        )
    }
    if (!is.data.frame(result$metrics) || !is.data.frame(result$scores) ||
        !is.data.frame(result$calls)) {
        .mamamia_error(
            "MaMaMIAComponentEvidenceError",
            "`result` must contain metrics, scores and calls"
        )
    }
    metrics <- result$metrics
    if (is.null(sample)) {
        samples <- unique(as.character(metrics$sample))
        sample <- if (length(samples) == 0L) NULL else samples[[1L]]
    }
    sample <- as.character(sample)[[1L]]
    if (is.na(sample) || !nzchar(sample)) {
        return(.component_evidence_empty())
    }
    metrics <- metrics[as.character(metrics$sample) == sample, , drop = FALSE]
    scores <- result$scores[
        as.character(result$scores$sample) == sample, , drop = FALSE
    ]
    if (is.null(calls)) {
        calls <- result$calls[
            as.character(result$calls$sample) == sample, , drop = FALSE
        ]
    } else {
        calls <- calls[
            as.character(calls$sample) == sample, , drop = FALSE
        ]
    }
    if (nrow(metrics) == 0L || nrow(calls) == 0L) {
        return(.component_evidence_empty())
    }
    metric <- if ("metric" %in% names(scores)) {
        unique(as.character(scores$metric))
    } else {
        "breadth"
    }
    if (length(metric) != 1L || is.na(metric) || !nzchar(metric)) {
        metric <- "breadth"
    }
    score_type <- if ("score_type" %in% names(scores)) {
        unique(as.character(scores$score_type))
    } else {
        "paired_log_ratio"
    }
    if (length(score_type) != 1L) {
        score_type <- "paired_log_ratio"
    }
    eps <- .scoring_eps
    metrics$component_lambda <- .component_evidence_lambda(
        metrics, metric, eps
    )
    donor <- metrics[metrics$role == "donor", , drop = FALSE]
    recipient <- metrics[metrics$role == "recipient", , drop = FALSE]
    if (nrow(donor) == 0L || nrow(recipient) == 0L) {
        return(.component_evidence_empty())
    }
    recipient_frame <- stats::aggregate(
        component_lambda ~ group + iid,
        data = recipient,
        FUN = mean, na.rm = TRUE
    )
    names(recipient_frame) <- c(
        "chromosome_group", "iid", "recipient_lambda"
    )
    recipient_frame$iid <- as.integer(recipient_frame$iid)
    score_frame <- scores[, c("group", "iid", "score", "donor_mean"),
                          drop = FALSE]
    names(score_frame) <- c(
        "chromosome_group", "iid", "group_score", "group_donor_lambda"
    )
    score_frame$iid <- as.integer(score_frame$iid)

    contigs <- result$assembly$contigs$donor
    contig_count <- function(group, component) {
        if (!is.data.frame(contigs) ||
            !all(c("group", "component") %in% names(contigs))) {
            return(NA_integer_)
        }
        sum(
            as.character(contigs$group) == as.character(group) &
                as.character(contigs$component) == as.character(component)
        )
    }

    rows <- list()
    for (i in seq_len(nrow(calls))) {
        call <- calls[i, , drop = FALSE]
        group <- as.character(call$group[[1L]])
        start <- as.integer(call$start_iid[[1L]])
        end <- as.integer(call$end_iid[[1L]])
        components <- sort(unique(as.character(
            donor$component[as.character(donor$group) == group]
        )))
        if (length(components) == 0L || is.na(start) || is.na(end)) {
            next
        }
        rec <- recipient_frame[
            recipient_frame$chromosome_group == group &
                recipient_frame$iid >= start &
                recipient_frame$iid <= end, ,
            drop = FALSE
        ]
        if (nrow(rec) == 0L) next
        component_rows <- list()
        for (component in components) {
            d <- donor[
                as.character(donor$group) == group &
                    as.character(donor$component) == component &
                    donor$iid >= start & donor$iid <= end, ,
                drop = FALSE
            ]
            if (nrow(d) == 0L) next
            d_lambda <- tapply(d$component_lambda, d$iid, mean, na.rm = TRUE)
            d_frame <- data.frame(
                iid = as.integer(names(d_lambda)),
                component_lambda = as.numeric(d_lambda),
                stringsAsFactors = FALSE
            )
            merged <- merge(
                d_frame, rec, by = "iid", all = FALSE, sort = TRUE
            )
            merged <- merge(
                merged, score_frame, by = "iid", all = FALSE, sort = TRUE
            )
            if (nrow(merged) == 0L) next
            ## Transfer the group-level GC/centring correction to each
            ## component by expressing it relative to the selected donor
            ## component at the same window.
            component_score <- merged$group_score +
                log(merged$component_lambda + eps) -
                log(merged$group_donor_lambda + eps)
            component_rows[[length(component_rows) + 1L]] <- data.frame(
                sample = sample,
                call_number = i,
                call_id = paste(
                    sample, group, start, end, sep = "_"
                ),
                chromosome_group = group,
                call_start_window = start,
                call_end_window = end,
                call_class = as.character(call$class[[1L]]),
                call_pass = isTRUE(call$pass[[1L]]),
                donor_component = component,
                n_component_contigs = contig_count(group, component),
                n_windows = nrow(merged),
                mean_component_lambda = mean(
                    merged$component_lambda, na.rm = TRUE
                ),
                mean_recipient_lambda = mean(
                    merged$recipient_lambda, na.rm = TRUE
                ),
                mean_component_score = mean(component_score, na.rm = TRUE),
                median_component_score = stats::median(
                    component_score, na.rm = TRUE
                ),
                max_component_score = max(component_score, na.rm = TRUE),
                supported_fraction = mean(
                    component_score > support_threshold, na.rm = TRUE
                ),
                score_type = score_type,
                metric = metric,
                stringsAsFactors = FALSE
            )
        }
        if (length(component_rows) > 0L) {
            component_rows <- do.call(rbind, component_rows)
            component_rows$component_rank <- rank(
                -component_rows$mean_component_score,
                ties.method = "min", na.last = "keep"
            )
            selected <- as.character(call$donor_component[[1L]])
            component_rows$component_selected <- if (is.na(selected)) {
                rep(FALSE, nrow(component_rows))
            } else {
                component_rows$donor_component == selected
            }
            component_rows$component_supported <-
                component_rows$supported_fraction >= 0.5
            rows[[length(rows) + 1L]] <- component_rows
        }
    }
    if (length(rows) == 0L) {
        return(.component_evidence_empty())
    }
    out <- do.call(rbind, rows)
    out <- out[order(
        out$sample, out$chromosome_group, out$call_start_window,
        out$component_rank
    ), , drop = FALSE]
    rownames(out) <- NULL
    out[, names(.component_evidence_empty()), drop = FALSE]
}
