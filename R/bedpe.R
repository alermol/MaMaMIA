## BEDPE export, writing the segments flagged as putative introgressions with
## both chromosome coordinates.

.bedpe_columns <- c(
    "chrom1",
    "start1",
    "end1",
    "chrom2",
    "start2",
    "end2",
    "name",
    "score",
    "strand1",
    "strand2",
    "chrom_name1",
    "chrom_name2",
    "subgen_name1",
    "subgen_name2"
)

#' Write putative introgressions to a BEDPE file
#'
#' Exports putative introgressed segments from an `ISA` object to a BEDPE
#' (Browser Extensible Data Paired-End) file. Each row is one donor-recipient
#' segment pair.
#'
#' Only segments with `putative_introgression = TRUE` are written. Coordinates
#' are taken from the donor and recipient windows that span each segment
#' (`start`/`stop` of the first and last ideal windows). If no putative
#' introgressions are present, a warning is issued and a header-only file is
#' written.
#'
#' Output columns, in order: `chrom1`, `start1`, `end1`, `chrom2`, `start2`,
#' `end2`, `name`, `score`, `strand1`, `strand2`, `chrom_name1`, `chrom_name2`,
#' `subgen_name1` and `subgen_name2`. With `sample` supplied, a `sample` column
#' is appended as the last column. The sample column is the frozen multi-sample
#' BEDPE format (D11, resolved question 5).
#'
#' @param ISA An object of class `ISA`, as returned by [segments()].
#' @param file Path to the output file. Written as a tab-separated table with a
#'   header row.
#' @param sample Optional sample name. When supplied, it is written to an extra
#'   `sample` column so that single-sample and cohort exports can be merged
#'   without losing the sample identity. `NULL` preserves the legacy column
#'   layout byte for byte.
#'
#' @return Called for its side effect of writing `file`. Invisibly returns
#'   `NULL`. Prints a message with the output path when finished.
#' @seealso [writeCohortToBedpe()]
#' @examples
#' data(triticum)
#' rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
#' \donttest{
#' ## Note: correctReadCounts() can take a while on the full dataset
#' rca <- correctReadCounts(rca, cores = 4L, verbose = FALSE)
#' }
#'
#' ## For demonstration purposes, assume that GC bias correction was performed,
#' ## by copying the 'cov' column to 'cor.gc'.
#' rca$data$valid <- TRUE
#' gc_q <- quantile(rca$data$gc[rca$data$valid], c(0.01, 0.99))
#' cov_q <- quantile(rca$data$cov[rca$data$valid], 0.99)
#' rca$data$ideal <- rca$data$valid &
#'                     (rca$data$gc >= gc_q[1]) &
#'                     (rca$data$gc <= gc_q[2]) &
#'                     (rca$data$cov <= cov_q)
#' rca$corrected <- TRUE
#' rca$outliers <- list(
#'         gc_lower_bound = gc_q[1],
#'         gc_upper_bound = gc_q[2],
#'         cov_upper_bound = cov_q
#'     )
#' rca$data$cor.gc <- rca$data$cov
#'
#' isa <- segments(
#'   rca,
#'   target_pairs = c("OY997261.1" = "CP155614.1"),
#'   verbose = FALSE
#' )
#' \donttest{
#' writeToBedpe(isa, file = 'introgressions.bedpe')
#' }
#' @export
writeToBedpe <- function(ISA, file, sample = NULL) {
    stopifnot("Input is not ISA object" = inherits(ISA, "ISA"))
    sample <- validate_bedpe_sample(sample)

    table <- bedpe_table(ISA, sample = sample)
    write_bedpe_table(table, file)
}

#' Write a combined BEDPE for a cohort
#'
#' Combines putative introgressions from one or more samples of a
#' `BamAnalysisResult` (or a named list of `ISA` objects) into one BEDPE file.
#' The `sample` column is always written as the last column, so the result
#' follows the frozen multi-sample format (D11, resolved question 5).
#'
#' By default every putative introgression of every exported sample is written.
#' Pass `selected` to export specific segment rows instead, as the dashboard
#' does: it is a named list mapping a sample name to 1-based row ids in that
#' sample's `ISA$segments` (the `.row_id` values of the explorer table). Samples
#' absent from the list contribute no rows. Row ids are validated against each
#' sample's segment table; a selected row is exported even when it is not
#' flagged as a putative introgression.
#'
#' When no rows remain, the header is still written and a warning is issued.
#'
#' @param result A `BamAnalysisResult` from [analyzeBamInput()], or a named list
#'   of `ISA` objects.
#' @param file Path to the output file.
#' @param samples Optional character vector of sample names to export. Defaults
#'   to all samples in `result`, in their original order.
#' @param selected Optional named list of segment row ids per sample. `NULL`
#'   (default) exports every putative introgression. Samples missing from the
#'   list export no rows, `integer(0)` also exports none for that sample, and a
#'   `NULL` entry exports every putative introgression for that sample.
#'
#' @return Called for its side effect of writing `file`. Invisibly returns
#'   `NULL`. Prints a message with the output path when finished.
#' @seealso [writeToBedpe()]
#' @examples
#' \dontrun{
#' result <- analyzeBamInput(
#'     c("sample1.bam", "sample2.bam"),
#'     config,
#'     normalization = "median"
#' )
#' writeCohortToBedpe(result, "cohort.bedpe")
#' writeCohortToBedpe(
#'     result, "selected.bedpe",
#'     selected = list(sample1 = c(1L, 4L), sample2 = 2L)
#' )
#' }
#' @export
writeCohortToBedpe <- function(result, file, samples = NULL, selected = NULL) {
    isas <- cohort_isa_list(result)
    sample_names <- names(isas)

    if (is.null(samples)) {
        samples <- sample_names
    } else {
        if (!is.character(samples) || length(samples) == 0L || anyNA(samples)) {
            stop("`samples` must be NULL or a non-empty character vector", call. = FALSE)
        }
        unknown <- setdiff(samples, sample_names)
        if (length(unknown) > 0L) {
            stop(
                "Unknown sample(s): ",
                paste(unknown, collapse = ", "),
                call. = FALSE
            )
        }
        samples <- unique(samples)
    }

    selected <- validate_bedpe_selection(selected, samples, isas)

    tables <- lapply(samples, function(sample) {
        rows <- if (is.null(selected)) {
            NULL
        } else if (sample %in% names(selected)) {
            selected[[sample]]
        } else {
            integer(0L)
        }
        bedpe_table(isas[[sample]], sample = sample, selected_rows = rows)
    })
    table <- dplyr::bind_rows(tables)

    write_bedpe_table(table, file)
}

## Internal BEDPE helpers -----------------------------------------------------

#' @noRd
validate_bedpe_sample <- function(sample) {
    if (is.null(sample)) {
        return(NULL)
    }
    if (!is.character(sample) || length(sample) != 1L ||
        is.na(sample) || !nzchar(sample)) {
        stop(
            "`sample` must be NULL or a single non-empty string",
            call. = FALSE
        )
    }
    sample
}

#' @noRd
empty_bedpe_table <- function(sample = NULL) {
    out <- data.frame(
        chrom1 = character(),
        start1 = integer(),
        end1 = integer(),
        chrom2 = character(),
        start2 = integer(),
        end2 = integer(),
        name = character(),
        score = numeric(),
        strand1 = character(),
        strand2 = character(),
        chrom_name1 = character(),
        chrom_name2 = character(),
        subgen_name1 = character(),
        subgen_name2 = character(),
        stringsAsFactors = FALSE
    )
    if (!is.null(sample)) {
        out$sample <- character()
    }
    out
}

## Build the legacy-compatible table for one ISA. `selected_rows` is either
## NULL (use the ISA's putative flags) or an already validated integer vector.
#' @noRd
bedpe_table <- function(ISA, sample = NULL, selected_rows = NULL) {
    if (is.null(selected_rows)) {
        rows <- dplyr::filter(ISA$segments, .data$putative_introgression)
    } else {
        rows <- ISA$segments[as.integer(selected_rows), , drop = FALSE]
        if (nrow(rows) > 0L) {
            rows$putative_introgression <- TRUE
        }
    }

    if (nrow(rows) == 0L) {
        return(empty_bedpe_table(sample))
    }

    out <- annotate_segment_coords(rows, ISA$out) |>
        dplyr::select(dplyr::all_of(
            c(
                "chr_id.don",
                "start.don",
                "stop.don",
                "chr_id.rec",
                "start.rec",
                "stop.rec",
                ISA$meta$metrics,
                "chr_name.don",
                "chr_name.rec",
                "subgenome.don",
                "subgenome.rec"
            )
        )) |>
        dplyr::mutate(
            name = paste("PIS", seq_len(dplyr::n()), sep = ""),
            .after = "stop.rec"
        ) |>
        dplyr::mutate(
            strand1 = ".",
            strand2 = ".",
            .after = ISA$meta$metrics
        ) |>
        stats::setNames(.bedpe_columns)

    if (!is.null(sample)) {
        out$sample <- sample
    }
    out
}

#' @noRd
write_bedpe_table <- function(table, file, warn_empty = TRUE) {
    if (nrow(table) == 0L && isTRUE(warn_empty)) {
        warning("No putative introgressions to write; writing header only",
            call. = FALSE
        )
    }
    utils::write.table(
        table,
        file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
    message("Output was written in ", file)
    invisible(NULL)
}

#' @noRd
cohort_isa_list <- function(result) {
    isas <- if (inherits(result, "BamAnalysisResult")) {
        result$isa
    } else if (is.list(result) && !inherits(result, "ISA")) {
        result
    } else {
        stop(
            "`result` must be a BamAnalysisResult or a named list of ISA objects",
            call. = FALSE
        )
    }

    valid_names <- is.list(isas) && length(isas) > 0L &&
        !is.null(names(isas)) && all(nzchar(names(isas)))
    if (!valid_names) {
        stop(
            "`result` must contain a non-empty named list of ISA objects",
            call. = FALSE
        )
    }
    bad <- !vapply(isas, inherits, logical(1L), what = "ISA")
    if (any(bad)) {
        stop(
            "`result` element(s) are not ISA objects: ",
            paste(names(isas)[bad], collapse = ", "),
            call. = FALSE
        )
    }
    isas
}

#' @noRd
validate_bedpe_selection <- function(selected, samples, isas) {
    if (is.null(selected)) {
        return(NULL)
    }
    if (!is.list(selected)) {
        stop(
            "`selected` must be NULL or a named list of segment row ids",
            call. = FALSE
        )
    }
    if (length(selected) > 0L &&
        (is.null(names(selected)) || any(!nzchar(names(selected))))) {
        stop(
            "`selected` must be NULL or a named list of segment row ids",
            call. = FALSE
        )
    }

    unknown <- setdiff(names(selected), samples)
    if (length(unknown) > 0L) {
        stop(
            "`selected` names sample(s) that are not exported: ",
            paste(unknown, collapse = ", "),
            call. = FALSE
        )
    }

    for (sample in names(selected)) {
        rows <- selected[[sample]]
        if (is.null(rows)) {
            next
        }
        if (!is.numeric(rows) || anyNA(rows) || any(!is.finite(rows)) ||
            any(rows != trunc(rows))) {
            stop(
                sprintf(
                    "Sample '%s': selected rows must be whole numbers",
                    sample
                ),
                call. = FALSE
            )
        }
        n <- nrow(isas[[sample]]$segments)
        if (any(rows < 1 | rows > n)) {
            stop(
                sprintf(
                    "Sample '%s': selected rows must be between 1 and %d",
                    sample, n
                ),
                call. = FALSE
            )
        }
        selected[[sample]] <- unique(as.integer(rows))
    }
    selected
}
