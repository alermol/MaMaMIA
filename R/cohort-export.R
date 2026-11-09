## Cohort export table, column picker and writer.
## The launcher keeps the Shiny wiring; the table and export logic lives
## here so it can be tested without a session.

#' Brief cohort-table columns for the Export tab
#' @noRd
bam_cohort_brief_columns <- function() {
    c(
        "sample", "n_windows", "reads_donor", "reads_recipient",
        "reads_total", "n_calls", "n_high_confidence", "masked_fraction"
    )
}

#' One row per sample with reads, call counts and mask burden
#'
#' Every column returned here is available to the Export table. The brief
#' columns are shown by default; the remaining diagnostics can be added from
#' the column picker.
#' @noRd
bam_cohort_summary <- function(result) {
    empty <- data.frame(
        sample = character(),
        mode = character(),
        n_windows = integer(),
        n_contigs = integer(),
        n_groups = integer(),
        n_recipient_contigs = integer(),
        n_donor_contigs = integer(),
        reads_donor = numeric(),
        reads_recipient = numeric(),
        reads_total = numeric(),
        reads_donor_per_window = numeric(),
        reads_recipient_per_window = numeric(),
        reads_total_per_window = numeric(),
        donor_read_share = numeric(),
        recipient_read_share = numeric(),
        n_calls = integer(),
        n_high_confidence = integer(),
        n_all_calls = integer(),
        n_class_confident = integer(),
        n_class_ambiguous = integer(),
        n_class_needs_review = integer(),
        n_substitution = integer(),
        n_addition = integer(),
        n_deletion = integer(),
        n_complex = integer(),
        n_other = integer(),
        n_groups_with_calls = integer(),
        calls_per_1000_windows = numeric(),
        high_confidence_per_1000_windows = numeric(),
        pass_rate = numeric(),
        class_confident_rate = numeric(),
        high_confidence_share = numeric(),
        masked_fraction = numeric(),
        masked_windows = numeric(),
        unmasked_fraction = numeric(),
        total_width_mb = numeric(),
        mean_width_mb = numeric(),
        max_width_mb = numeric(),
        mean_call_score = numeric(),
        max_call_score = numeric(),
        stringsAsFactors = FALSE
    )
    if (is.null(result) || is.null(result$samples) ||
        !is.data.frame(result$samples) || nrow(result$samples) == 0L) {
        return(empty)
    }
    samples <- bam_sample_choices(result)
    if (length(samples) == 0L) {
        return(empty)
    }

    meta <- result$samples
    position <- match(samples, as.character(meta$sample))
    column <- function(name, type = "numeric", default = NA) {
        value <- if (name %in% names(meta)) {
            meta[[name]][position]
        } else {
            rep(default, length(samples))
        }
        switch(
            type,
            integer = as.integer(value),
            numeric = as.numeric(value),
            character = as.character(value),
            value
        )
    }
    n_windows <- column("n_windows", "integer")
    n_calls <- column("n_calls", "integer")
    n_high_confidence <- column("n_high_confidence", "integer")
    masked_fraction <- column("masked_fraction", "numeric")
    mode <- column("mode", "character", NA_character_)

    metrics <- result$metrics
    reads <- function(roles) {
        out <- matrix(
            0, nrow = length(samples), ncol = length(roles),
            dimnames = list(samples, roles)
        )
        if (!is.data.frame(metrics) ||
            !all(c("sample", "role", "count") %in% names(metrics))) {
            return(out)
        }
        sample_value <- as.character(metrics$sample)
        role_value <- as.character(metrics$role)
        keep <- sample_value %in% samples & role_value %in% roles
        if (!any(keep)) {
            return(out)
        }
        ## One grouped pass over the window table. The previous per-sample
        ## vapply() rescanned every sample's windows once per role.
        totals <- tapply(
            as.numeric(metrics$count[keep]),
            list(
                factor(sample_value[keep], levels = samples),
                factor(role_value[keep], levels = roles)
            ),
            sum, na.rm = TRUE
        )
        totals[is.na(totals)] <- 0
        out[, colnames(totals)] <- totals
        out
    }
    read_totals <- reads(c("donor", "recipient"))
    reads_donor <- read_totals[, "donor"]
    reads_recipient <- read_totals[, "recipient"]
    reads_total <- reads_donor + reads_recipient

    windows <- result$windows$table
    if (is.data.frame(windows) && nrow(windows) > 0L &&
        all(c("chr_id", "role", "group") %in% names(windows))) {
        n_contigs <- length(unique(as.character(windows$chr_id)))
        n_groups <- length(unique(as.character(windows$group)))
        recipient <- as.character(windows$role) == "recipient"
        donor <- as.character(windows$role) == "donor"
        n_recipient_contigs <- length(unique(
            as.character(windows$chr_id[recipient])
        ))
        n_donor_contigs <- length(unique(
            as.character(windows$chr_id[donor])
        ))
    } else {
        n_contigs <- NA_integer_
        n_groups <- NA_integer_
        n_recipient_contigs <- NA_integer_
        n_donor_contigs <- NA_integer_
    }

    call_stats <- data.frame(
        n_all_calls = integer(length(samples)),
        n_class_confident = integer(length(samples)),
        n_class_ambiguous = integer(length(samples)),
        n_class_needs_review = integer(length(samples)),
        n_substitution = integer(length(samples)),
        n_addition = integer(length(samples)),
        n_deletion = integer(length(samples)),
        n_complex = integer(length(samples)),
        n_other = integer(length(samples)),
        n_groups_with_calls = integer(length(samples)),
        total_width_mb = rep(NA_real_, length(samples)),
        mean_width_mb = rep(NA_real_, length(samples)),
        max_width_mb = rep(NA_real_, length(samples)),
        mean_call_score = rep(NA_real_, length(samples)),
        max_call_score = rep(NA_real_, length(samples)),
        stringsAsFactors = FALSE
    )
    calls <- result$calls
    if (is.data.frame(calls) && nrow(calls) > 0L &&
        "sample" %in% names(calls)) {
        for (i in seq_along(samples)) {
            one <- calls[as.character(calls$sample) == samples[[i]], ,
                drop = FALSE
            ]
            if (nrow(one) == 0L) {
                next
            }
            call_stats$n_all_calls[[i]] <- nrow(one)
            status <- if ("class_status" %in% names(one)) {
                as.character(one$class_status)
            } else {
                rep(NA_character_, nrow(one))
            }
            call_stats$n_class_confident[[i]] <-
                sum(status == "confident", na.rm = TRUE)
            call_stats$n_class_ambiguous[[i]] <-
                sum(status == "ambiguous", na.rm = TRUE)
            call_stats$n_class_needs_review[[i]] <-
                sum(status == "needs_review", na.rm = TRUE)
            passing <- if ("pass" %in% names(one)) {
                as.logical(one$pass)
            } else {
                rep(TRUE, nrow(one))
            }
            passing[is.na(passing)] <- FALSE
            selected <- one[passing, , drop = FALSE]
            if (nrow(selected) == 0L) {
                next
            }
            classes <- if ("class" %in% names(selected)) {
                as.character(selected$class)
            } else {
                rep(NA_character_, nrow(selected))
            }
            call_stats$n_substitution[[i]] <-
                sum(classes == "substitution", na.rm = TRUE)
            call_stats$n_addition[[i]] <-
                sum(classes == "addition", na.rm = TRUE)
            call_stats$n_deletion[[i]] <-
                sum(classes == "deletion", na.rm = TRUE)
            call_stats$n_complex[[i]] <-
                sum(classes == "complex", na.rm = TRUE)
            call_stats$n_other[[i]] <-
                sum(classes == "other", na.rm = TRUE)
            if ("group" %in% names(selected)) {
                groups <- as.character(selected$group)
                call_stats$n_groups_with_calls[[i]] <-
                    length(unique(groups[!is.na(groups)]))
            }
            if ("width_bp" %in% names(selected)) {
                widths <- as.numeric(selected$width_bp)
                widths <- widths[is.finite(widths)]
                if (length(widths) > 0L) {
                    call_stats$total_width_mb[[i]] <- sum(widths) / 1e6
                    call_stats$mean_width_mb[[i]] <- mean(widths) / 1e6
                    call_stats$max_width_mb[[i]] <- max(widths) / 1e6
                }
            }
            if ("score_mean" %in% names(selected)) {
                scores <- as.numeric(selected$score_mean)
                scores <- scores[is.finite(scores)]
                if (length(scores) > 0L) {
                    call_stats$mean_call_score[[i]] <- mean(scores)
                    call_stats$max_call_score[[i]] <- max(scores)
                }
            }
        }
    }

    safe_div <- function(numerator, denominator) {
        out <- rep(NA_real_, length(numerator))
        ok <- is.finite(numerator) & is.finite(denominator) &
            denominator > 0
        out[ok] <- numerator[ok] / denominator[ok]
        out
    }

    data.frame(
        sample = samples,
        mode = mode,
        n_windows = n_windows,
        n_contigs = rep(n_contigs, length(samples)),
        n_groups = rep(n_groups, length(samples)),
        n_recipient_contigs = rep(n_recipient_contigs, length(samples)),
        n_donor_contigs = rep(n_donor_contigs, length(samples)),
        reads_donor = reads_donor,
        reads_recipient = reads_recipient,
        reads_total = reads_total,
        reads_donor_per_window = safe_div(reads_donor, n_windows),
        reads_recipient_per_window = safe_div(reads_recipient, n_windows),
        reads_total_per_window = safe_div(reads_total, n_windows),
        donor_read_share = safe_div(reads_donor, reads_total),
        recipient_read_share = safe_div(reads_recipient, reads_total),
        n_calls = n_calls,
        n_high_confidence = n_high_confidence,
        n_all_calls = call_stats$n_all_calls,
        n_class_confident = call_stats$n_class_confident,
        n_class_ambiguous = call_stats$n_class_ambiguous,
        n_class_needs_review = call_stats$n_class_needs_review,
        n_substitution = call_stats$n_substitution,
        n_addition = call_stats$n_addition,
        n_deletion = call_stats$n_deletion,
        n_complex = call_stats$n_complex,
        n_other = call_stats$n_other,
        n_groups_with_calls = call_stats$n_groups_with_calls,
        calls_per_1000_windows = safe_div(n_calls, n_windows) * 1000,
        high_confidence_per_1000_windows =
            safe_div(n_high_confidence, n_windows) * 1000,
        pass_rate = safe_div(n_calls, call_stats$n_all_calls),
        class_confident_rate =
            safe_div(call_stats$n_class_confident, call_stats$n_all_calls),
        high_confidence_share = safe_div(n_high_confidence, n_calls),
        masked_fraction = masked_fraction,
        masked_windows = masked_fraction * n_windows,
        unmasked_fraction = 1 - masked_fraction,
        total_width_mb = call_stats$total_width_mb,
        mean_width_mb = call_stats$mean_width_mb,
        max_width_mb = call_stats$max_width_mb,
        mean_call_score = call_stats$mean_call_score,
        max_call_score = call_stats$max_call_score,
        stringsAsFactors = FALSE
    )
}

#' Human-readable labels for cohort-table columns
#' @noRd
bam_cohort_column_labels <- function(columns = NULL) {
    labels <- c(
        sample = "Sample",
        mode = "Mode",
        n_windows = "Windows",
        n_contigs = "Contigs",
        n_groups = "Groups",
        n_recipient_contigs = "Recipient contigs",
        n_donor_contigs = "Donor contigs",
        reads_donor = "Donor reads",
        reads_recipient = "Recipient reads",
        reads_total = "Total reads",
        reads_donor_per_window = "Donor reads/window",
        reads_recipient_per_window = "Recipient reads/window",
        reads_total_per_window = "Reads/window",
        donor_read_share = "Donor share (%)",
        recipient_read_share = "Recipient share (%)",
        n_calls = "Passing calls",
        n_high_confidence = "High-conf. calls",
        n_all_calls = "All calls",
        n_class_confident = "Confident calls",
        n_class_ambiguous = "Ambiguous calls",
        n_class_needs_review = "Needs-review calls",
        n_substitution = "Substitutions",
        n_addition = "Additions",
        n_deletion = "Deletions",
        n_complex = "Complex calls",
        n_other = "Other calls",
        n_groups_with_calls = "Groups with calls",
        calls_per_1000_windows = "Calls/1,000 windows",
        high_confidence_per_1000_windows = "High-conf./1,000 windows",
        pass_rate = "Pass rate (%)",
        class_confident_rate = "Confident rate (%)",
        high_confidence_share = "High-conf. share (%)",
        masked_fraction = "Masked %",
        masked_windows = "Masked windows",
        unmasked_fraction = "Unmasked %",
        total_width_mb = "Total width (Mb)",
        mean_width_mb = "Mean width (Mb)",
        max_width_mb = "Largest call (Mb)",
        mean_call_score = "Mean call score",
        max_call_score = "Max call score"
    )
    if (is.null(columns)) {
        return(labels)
    }
    out <- unname(labels[columns])
    out[is.na(out)] <- columns[is.na(out)]
    out
}

#' Column metadata for the Export cohort table
#' @noRd
bam_cohort_column_metadata <- function(columns = NULL) {
    labels <- bam_cohort_column_labels()
    descriptions <- c(
        sample = "Sample name from the BAM @RG SM tag.",
        mode = "Analysis mode: single for one sample, multi for a cohort.",
        n_windows = "Analysis windows tiled across all analysed contigs.",
        n_contigs = "Analysed contigs in the donor + recipient assembly.",
        n_groups = "Chromosome groups in the assembly.",
        n_recipient_contigs = "Recipient contigs in the assembly.",
        n_donor_contigs = "Donor contigs in the assembly.",
        reads_donor = "Reads assigned to donor chromosomes in analysis windows.",
        reads_recipient = "Reads assigned to recipient chromosomes in analysis windows.",
        reads_total = "Donor plus recipient reads in analysis windows.",
        reads_donor_per_window = "Mean donor reads per analysis window.",
        reads_recipient_per_window = "Mean recipient reads per analysis window.",
        reads_total_per_window = "Mean donor + recipient reads per analysis window.",
        donor_read_share = "Donor reads as a share of all analysed reads.",
        recipient_read_share = "Recipient reads as a share of all analysed reads.",
        n_calls = "Calls that pass calibration filters.",
        n_high_confidence = "Passing calls with high confidence.",
        n_all_calls = "All calls, including those that fail calibration.",
        n_class_confident = "Calls whose class status is confident.",
        n_class_ambiguous = "Calls whose class status is ambiguous.",
        n_class_needs_review = "Calls whose class status needs review.",
        n_substitution = "Passing calls classified as substitutions.",
        n_addition = "Passing calls classified as additions.",
        n_deletion = "Passing calls classified as deletions.",
        n_complex = "Passing calls classified as complex.",
        n_other = "Passing calls classified as other.",
        n_groups_with_calls = "Chromosome groups with at least one passing call.",
        calls_per_1000_windows = "Passing calls per 1,000 analysis windows.",
        high_confidence_per_1000_windows =
            "High-confidence passing calls per 1,000 analysis windows.",
        pass_rate = "Passing calls as a share of all calls.",
        class_confident_rate = "Confident class status as a share of all calls.",
        high_confidence_share = "High-confidence calls as a share of passing calls.",
        masked_fraction = "Share of analysis windows masked by low mappability.",
        masked_windows = "Number of analysis windows masked by low mappability.",
        unmasked_fraction = "Share of analysis windows not masked.",
        total_width_mb = "Total width of passing calls in megabases.",
        mean_width_mb = "Mean width of passing calls in megabases.",
        max_width_mb = "Width of the largest passing call in megabases.",
        mean_call_score = "Mean score of passing calls.",
        max_call_score = "Maximum score among passing calls."
    )
    formats <- c(
        sample = "text", mode = "text", n_windows = "integer",
        n_contigs = "integer", n_groups = "integer",
        n_recipient_contigs = "integer", n_donor_contigs = "integer",
        reads_donor = "integer", reads_recipient = "integer",
        reads_total = "integer", reads_donor_per_window = "decimal1",
        reads_recipient_per_window = "decimal1",
        reads_total_per_window = "decimal1", donor_read_share = "percent",
        recipient_read_share = "percent", n_calls = "integer",
        n_high_confidence = "integer", n_all_calls = "integer",
        n_class_confident = "integer", n_class_ambiguous = "integer",
        n_class_needs_review = "integer", n_substitution = "integer",
        n_addition = "integer", n_deletion = "integer",
        n_complex = "integer", n_other = "integer",
        n_groups_with_calls = "integer",
        calls_per_1000_windows = "decimal1",
        high_confidence_per_1000_windows = "decimal1",
        pass_rate = "percent", class_confident_rate = "percent",
        high_confidence_share = "percent", masked_fraction = "percent",
        masked_windows = "integer", unmasked_fraction = "percent",
        total_width_mb = "decimal2", mean_width_mb = "decimal2",
        max_width_mb = "decimal2", mean_call_score = "signif",
        max_call_score = "signif"
    )
    if (is.null(columns)) {
        columns <- names(labels)
    }
    index <- match(columns, names(labels))
    label <- unname(labels[index])
    label[is.na(label)] <- columns[is.na(label)]
    description <- unname(descriptions[index])
    description[is.na(description)] <- "No description available."
    format <- unname(formats[index])
    format[is.na(format)] <- "text"
    data.frame(
        column = columns,
        label = label,
        description = description,
        format = format,
        stringsAsFactors = FALSE
    )
}

#' Selectable extra columns for the Export cohort table
#' @noRd
bam_cohort_extra_column_choices <- function() {
    metadata <- bam_cohort_column_metadata()
    extra <- metadata[
        !metadata$column %in% bam_cohort_brief_columns(),
        ,
        drop = FALSE
    ]
    stats::setNames(extra$column, extra$label)
}

#' Additional-columns picker for the Export cohort table
#' @noRd
bam_cohort_columns_ui <- function() {
    shiny::div(
        class = "mamamia-call-columns",
        shinyWidgets::pickerInput(
            "cohort_columns",
            "Additional columns",
            choices = bam_cohort_extra_column_choices(),
            selected = character(0),
            multiple = TRUE,
            width = "fit",
            options = shinyWidgets::pickerOptions(
                actionsBox = TRUE,
                liveSearch = TRUE,
                title = "Add columns",
                noneSelectedText = "Add columns",
                selectedTextFormat = "count > 0",
                countSelectedText = "{0} added"
            )
        )
    )
}

#' Header for the Export cohort table
#' @noRd
bam_cohort_table_header <- function(columns) {
    bam_table_header(bam_cohort_column_metadata(columns))
}

#' Sortable Export cohort table
#' @noRd
bam_cohort_table <- function(summary, selected_columns = character(0)) {
    defaults <- bam_cohort_brief_columns()
    selected_columns <- as.character(selected_columns)
    selected_columns <- selected_columns[!is.na(selected_columns)]
    extras <- names(summary)[names(summary) %in% selected_columns]
    columns <- unique(c(defaults, extras))
    columns <- columns[columns %in% names(summary)]
    display <- summary[, columns, drop = FALSE]
    table <- DT::datatable(
        display,
        container = bam_cohort_table_header(columns),
        rownames = FALSE,
        selection = "none",
        options = list(
            pageLength = 10,
            scrollX = TRUE,
            orderMulti = TRUE,
            columnDefs = list(list(
                orderSequence = c("asc", "desc", ""),
                targets = "_all"
            )),
            initComplete = DT::JS(
                "function(settings) {
                    var api = this.api();
                    api.on('order.dt', function() {
                        var order = api.order();
                        var unsorted = order.some(function(item) {
                            return item[1] === '';
                        });
                        if (unsorted) {
                            api.order([]).draw();
                        }
                    });
                }"
            )
        )
    )
    if (nrow(display) == 0L) {
        return(table)
    }
    metadata <- bam_cohort_column_metadata(columns)
    for (kind in unique(metadata$format)) {
        cols <- which(metadata$format == kind)
        if (length(cols) == 0L) {
            next
        }
        table <- switch(
            kind,
            integer = DT::formatRound(table, cols, digits = 0),
            decimal1 = DT::formatRound(table, cols, digits = 1),
            decimal2 = DT::formatRound(table, cols, digits = 2),
            percent = DT::formatPercentage(table, cols, digits = 1),
            signif = DT::formatSignif(table, cols, digits = 3),
            table
        )
    }
    table
}

#' Export-tab cohort card
#'
#' Builds the cohort/sample summary card for the launcher. The Shiny
#' registration stays in `bam_launcher_server()`; this function only owns the
#' markup, so the card can be checked without a session.
#'
#' @param result A `BamAnalysisResult` or `NULL` before a run.
#' @return A `shiny.tag` describing the cohort export card.
#' @noRd
bam_cohort_export_ui <- function(result) {
    if (is.null(result)) {
        return(shiny::p(
            class = "text-muted small",
            "Run the analysis to enable the exports."
        ))
    }
    samples <- bam_sample_choices(result)
    multiple_samples <- length(samples) > 1L
    shiny::div(
        class = paste(
            "mamamia-export-section mamamia-export-card",
            "mamamia-export-cohort"
        ),
        shiny::div(
            class = "mamamia-export-card-head",
            shiny::h5(
                if (multiple_samples) "Cohort" else "Sample summary"
            ),
            bam_cohort_columns_ui()
        ),
        shiny::p(
            class = "mamamia-call-note",
            paste(
                "Hover a column heading for its description;",
                "use Additional columns to add diagnostics.",
                "Shift-click headings to sort by several",
                "columns; click again to clear the sort."
            )
        ),
        DT::DTOutput("cohort_table"),
        shiny::div(
            class = "mamamia-export-card-foot",
            shiny::p(
                class = "mamamia-export-note",
                if (multiple_samples) {
                    paste(
                        "The results bundle above contains one folder",
                        "per sample plus cohort summaries at the top",
                        "level. The cohort BEDPE adds a `sample` column",
                        "for every sample in the BAM."
                    )
                } else {
                    paste(
                        "The results bundle above contains a folder",
                        "for this sample plus run-wide summaries at",
                        "the top level."
                    )
                }
            ),
            if (multiple_samples) {
                explorer_download_button(
                    "download_cohort",
                    "Download all samples (BEDPE)",
                    class = "btn-sm",
                    busy = "Preparing cohort BEDPE..."
                )
            }
        )
    )
}

#' Write the cohort BEDPE download
#'
#' @param result A `BamAnalysisResult` or `NULL`.
#' @param file Output path supplied by the download handler.
#' @return The output path, invisibly.
#' @noRd
bam_cohort_write_bedpe <- function(result, file) {
    if (is.null(result)) {
        stop("Run the analysis before downloading", call. = FALSE)
    }
    writeCallsToBedpe(result, file, pass_only = TRUE)
}
