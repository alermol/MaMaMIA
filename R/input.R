## Reading the on-disk input layout behind runMaMaMIA()
##
## analysis_dir/<sample>/<reference>.w<N>.dcov.bed
## analysis_dir/<sample>/<reference>.w<N>.rcov.bed
## assemblies_dir/<reference>/donor.w<N>.gc.bed
## assemblies_dir/<reference>/recipient.w<N>.gc.bed
## assemblies_dir/<reference>/chromosomes.tsv
##
## Coverage and GC tables are tab-separated with four columns (chromosome ID,
## start, end, value); GC tables may carry the bedtools nuc header, which starts
## with "#". chromosomes.tsv holds three columns: chromosome ID, short display
## name, and subgenome label. <N> labels the windowing scheme: the number of
## windows each chromosome was tiled into.

#' @noRd
read_window_table <- function(path, value_col, what = value_col) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !file.exists(path)) {
        stop(sprintf("Missing %s file: %s", what, path), call. = FALSE)
    }
    if (isTRUE(file.size(path) == 0L)) {
        stop(sprintf("%s file %s contains no data rows", what, path), call. = FALSE)
    }

    table <- tryCatch(
        utils::read.delim(
            path,
            header = FALSE,
            comment.char = "#",
            colClasses = "character",
            stringsAsFactors = FALSE,
            check.names = FALSE
        ),
        error = function(e) {
            stop(
                sprintf("Could not read %s file %s: %s", what, path, conditionMessage(e)),
                call. = FALSE
            )
        }
    )

    if (nrow(table) == 0L) {
        stop(sprintf("%s file %s contains no data rows", what, path), call. = FALSE)
    }
    if (ncol(table) < 4L) {
        stop(
            sprintf(
                "%s file %s must have at least 4 columns (chr_id, start, stop, %s), found %d",
                what, path, value_col, ncol(table)
            ),
            call. = FALSE
        )
    }

    table <- table[, 1:4, drop = FALSE]
    names(table) <- c("chr_id", "start", "stop", value_col)

    table$chr_id <- trimws(table$chr_id)
    bad_chr <- is.na(table$chr_id) | !nzchar(table$chr_id)
    if (any(bad_chr)) {
        stop(
            sprintf(
                "%s file %s has %d empty chromosome ID(s) (first at row %d)",
                what, path, sum(bad_chr), which(bad_chr)[[1L]]
            ),
            call. = FALSE
        )
    }

    for (col in c("start", "stop", value_col)) {
        raw <- trimws(table[[col]])
        values <- suppressWarnings(as.numeric(raw))
        bad <- is.na(values)
        if (any(bad)) {
            stop(
                sprintf(
                    "Column `%s` of %s file %s has %d non-numeric or missing value(s) (first at row %d)",
                    col, what, path, sum(bad), which(bad)[[1L]]
                ),
                call. = FALSE
            )
        }
        table[[col]] <- values
    }

    table
}

#' @noRd
read_chromosome_meta <- function(path) {
    if (!is.character(path) || length(path) != 1L || is.na(path) || !file.exists(path)) {
        stop(sprintf("Missing chromosome metadata file: %s", path), call. = FALSE)
    }
    if (isTRUE(file.size(path) == 0L)) {
        stop(sprintf("Chromosome metadata %s contains no rows", path), call. = FALSE)
    }

    table <- tryCatch(
        utils::read.delim(
            path,
            header = FALSE,
            comment.char = "#",
            colClasses = "character",
            stringsAsFactors = FALSE,
            check.names = FALSE
        ),
        error = function(e) {
            stop(
                sprintf("Could not read chromosome metadata %s: %s", path, conditionMessage(e)),
                call. = FALSE
            )
        }
    )

    if (ncol(table) != 3L) {
        stop(
            sprintf(
                "Chromosome metadata %s must have exactly 3 columns (chr_id, chr_name, subgenome), found %d",
                path, ncol(table)
            ),
            call. = FALSE
        )
    }

    table[] <- lapply(table, trimws)
    names(table) <- c("chr_id", "chr_name", "subgenome")

    if (nrow(table) == 0L) {
        stop(sprintf("Chromosome metadata %s contains no rows", path), call. = FALSE)
    }
    for (col in names(table)) {
        if (anyNA(table[[col]]) || any(!nzchar(table[[col]]))) {
            stop(
                sprintf("Chromosome metadata %s must not contain empty `%s` values", path, col),
                call. = FALSE
            )
        }
    }
    if (anyDuplicated(table$chr_id) > 0L) {
        stop(
            sprintf(
                "Chromosome metadata %s has duplicated chromosome IDs (e.g. %s)",
                path, table$chr_id[anyDuplicated(table$chr_id)]
            ),
            call. = FALSE
        )
    }

    table
}

#' @noRd
check_single_string <- function(value, arg_name) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
        stop(sprintf("`%s` must be a single non-empty string", arg_name), call. = FALSE)
    }
    invisible(value)
}

#' @noRd
check_directory <- function(path, arg_name) {
    check_single_string(path, arg_name)
    if (!dir.exists(path)) {
        stop(sprintf("`%s` does not exist: %s", arg_name, path), call. = FALSE)
    }
    invisible(path)
}

#' List samples found in an analysis directory
#'
#' Scans the analysis directory of the [readMaMaMIAInput()] layout and returns
#' the sample names, that is the subdirectories that hold at least one
#' `*.dcov.bed` coverage table.
#'
#' @param analysis_dir Path to the analysis directory.
#'
#' @return A sorted character vector of sample names.
#' @seealso [readMaMaMIAInput()]
#' @export
listMaMaMIASamples <- function(analysis_dir) {
    check_directory(analysis_dir, "analysis_dir")
    entries <- list.dirs(analysis_dir, recursive = FALSE, full.names = FALSE)
    has_coverage <- vapply(
        file.path(analysis_dir, entries),
        function(dir) length(list.files(dir, pattern = "\\.dcov\\.bed$")) > 0L,
        logical(1L)
    )
    sort(entries[has_coverage])
}

#' List references found in an assemblies directory
#'
#' Scans the assemblies directory of the [readMaMaMIAInput()] layout and
#' returns the reference names, that is the subdirectories that hold a
#' `chromosomes.tsv` metadata table.
#'
#' @param assemblies_dir Path to the assemblies directory.
#'
#' @return A sorted character vector of reference names.
#' @seealso [readMaMaMIAInput()]
#' @export
listMaMaMIAReferences <- function(assemblies_dir) {
    check_directory(assemblies_dir, "assemblies_dir")
    entries <- list.dirs(assemblies_dir, recursive = FALSE, full.names = FALSE)
    has_meta <- file.exists(file.path(assemblies_dir, entries, "chromosomes.tsv"))
    sort(entries[has_meta])
}

#' @noRd
coverage_labels <- function(dir, reference, kind) {
    files <- list.files(dir, pattern = paste0("\\.", kind, "\\.bed$"))
    prefix <- paste0(reference, ".w")
    files <- files[startsWith(files, prefix)]
    rest <- substring(files, nchar(prefix) + 1L)
    rest <- sub(paste0("\\.", kind, "\\.bed$"), "", rest)
    ok <- grepl("^[0-9]+$", rest)
    as.integer(rest[ok])
}

#' @noRd
list_window_labels <- function(analysis_dir, assemblies_dir, sample, reference) {
    check_directory(analysis_dir, "analysis_dir")
    check_directory(assemblies_dir, "assemblies_dir")
    check_single_string(sample, "sample")
    check_single_string(reference, "reference")

    sample_dir <- file.path(analysis_dir, sample)
    reference_dir <- file.path(assemblies_dir, reference)
    if (!dir.exists(sample_dir)) {
        stop(sprintf("Sample directory does not exist: %s", sample_dir), call. = FALSE)
    }
    if (!dir.exists(reference_dir)) {
        stop(sprintf("Reference directory does not exist: %s", reference_dir), call. = FALSE)
    }

    labels <- sort(intersect(
        coverage_labels(sample_dir, reference, "dcov"),
        coverage_labels(sample_dir, reference, "rcov")
    ))
    have_assemblies <- vapply(labels, function(label) {
        all(file.exists(
            file.path(reference_dir, sprintf("donor.w%d.gc.bed", label)),
            file.path(reference_dir, sprintf("recipient.w%d.gc.bed", label))
        ))
    }, logical(1L))
    labels[have_assemblies]
}

#' Read windowed input data for one sample and reference
#'
#' Reads the on-disk layout used by [runMaMaMIA()]. An analysis directory holds
#' one subdirectory per sample, named `<sample>`, with donor and recipient
#' coverage tables `<reference>.w<N>.dcov.bed` and `<reference>.w<N>.rcov.bed`.
#' An assemblies directory holds one subdirectory per donor reference, named
#' `<reference>`, with the donor and recipient GC-content tables
#' `donor.w<N>.gc.bed` and `recipient.w<N>.gc.bed`, plus the chromosome
#' metadata `chromosomes.tsv`. All tables are tab-separated four-column
#' (chromosome ID, start, end, value) tables, except `chromosomes.tsv`, which
#' holds chromosome ID, short display name, and subgenome label. GC tables may
#' carry the `bedtools nuc` header, which starts with `#`.
#'
#' The `w<N>` label denotes the windowing scheme: the number of windows each
#' chromosome was tiled into. When several schemes are available for a sample
#' and reference, choose one with the `window` argument.
#'
#' @param analysis_dir Path to the analysis directory.
#' @param assemblies_dir Path to the assemblies directory.
#' @param sample Sample name (a subdirectory of `analysis_dir`).
#' @param reference Donor reference name (a subdirectory of `assemblies_dir`).
#' @param window Windowing scheme label, an integer such as `1000` for
#'   `w1000` files. Defaults to the only scheme available.
#'
#' @return A list with the donor and recipient coverage and GC tables
#'   (`dcov`, `dgc`, `rcov`, `rgc`) ready for [RCA()], the chromosome
#'   `meta` table restricted to the chromosomes present in the coverage
#'   tables, the donor and recipient chromosome IDs, and the resolved
#'   `sample`, `reference`, and `window`.
#' @seealso [listMaMaMIASamples()], [listMaMaMIAReferences()], [RCA()]
#' @export
readMaMaMIAInput <- function(analysis_dir, assemblies_dir, sample, reference, window = NULL) {
    labels <- list_window_labels(analysis_dir, assemblies_dir, sample, reference)
    if (length(labels) == 0L) {
        stop(
            sprintf(
                paste0(
                    "No windowing scheme found for sample '%s' and reference '%s'.\n",
                    "Expected %s/<sample>/%s.w<N>.dcov.bed and .rcov.bed, and %s/%s/donor.w<N>.gc.bed and recipient.w<N>.gc.bed"
                ),
                sample, reference, analysis_dir, reference, assemblies_dir, reference
            ),
            call. = FALSE
        )
    }

    if (is.null(window)) {
        if (length(labels) > 1L) {
            stop(
                sprintf(
                    "Multiple windowing schemes available (%s); choose one with `window`",
                    paste(sprintf("w%d", labels), collapse = ", ")
                ),
                call. = FALSE
            )
        }
        window <- labels[[1L]]
    }
    window <- as.character(window)[[1L]]
    if (!(window %in% as.character(labels))) {
        stop(
            sprintf(
                "Window '%s' is not available for sample '%s' and reference '%s'; available: %s",
                window, sample, reference,
                paste(sprintf("w%d", labels), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    sample_dir <- file.path(analysis_dir, sample)
    reference_dir <- file.path(assemblies_dir, reference)

    dcov <- read_window_table(
        file.path(sample_dir, sprintf("%s.w%s.dcov.bed", reference, window)),
        "cov", "donor coverage"
    )
    rcov <- read_window_table(
        file.path(sample_dir, sprintf("%s.w%s.rcov.bed", reference, window)),
        "cov", "recipient coverage"
    )
    dgc <- read_window_table(
        file.path(reference_dir, sprintf("donor.w%s.gc.bed", window)),
        "gc", "donor GC"
    )
    rgc <- read_window_table(
        file.path(reference_dir, sprintf("recipient.w%s.gc.bed", window)),
        "gc", "recipient GC"
    )
    meta <- read_chromosome_meta(file.path(reference_dir, "chromosomes.tsv"))

    don_chr_ids <- unique(dcov$chr_id)
    rec_chr_ids <- unique(rcov$chr_id)
    shared <- intersect(don_chr_ids, rec_chr_ids)
    if (length(shared) > 0L) {
        stop(
            sprintf(
                "Chromosome IDs %s appear in both donor and recipient coverage",
                paste(utils::head(shared, 5L), collapse = ", ")
            ),
            call. = FALSE
        )
    }

    present <- c(don_chr_ids, rec_chr_ids)
    missing_meta <- setdiff(present, meta$chr_id)
    if (length(missing_meta) > 0L) {
        stop(
            sprintf(
                "chromosomes.tsv in %s is missing: %s",
                reference_dir,
                paste(utils::head(missing_meta, 5L), collapse = ", ")
            ),
            call. = FALSE
        )
    }
    meta <- meta[meta$chr_id %in% present, , drop = FALSE]

    list(
        dcov = dcov,
        dgc = dgc,
        rcov = rcov,
        rgc = rgc,
        meta = meta,
        don_chr_ids = don_chr_ids,
        rec_chr_ids = rec_chr_ids,
        sample = sample,
        reference = reference,
        window = window
    )
}

#' @noRd
chr_number <- function(chr_name) {
    has_number <- grepl("[0-9]", chr_name)
    number <- sub("^[^0-9]*([0-9]+).*$", "\\1", chr_name)
    number[!has_number] <- NA_character_
    number
}

#' Pair donor and recipient chromosomes by chromosome number
#'
#' Builds a named `target_pairs` vector for [segments()] by matching donor and
#' recipient chromosomes with the same number in their short display names
#' (`chr_name`, for example `"1At"` with `"1A"`). Donor chromosomes whose
#' number matches no recipient chromosome, or more than one, are left out.
#' This is the default pairing used by [runMaMaMIA()]; it keeps only one
#' recipient chromosome per donor chromosome.
#'
#' @param meta Chromosome metadata with columns `chr_id`, `chr_name`, and
#'   `subgenome`, as read by [readMaMaMIAInput()].
#' @param don_chr_ids Donor chromosome IDs present in the data.
#' @param rec_chr_ids Recipient chromosome IDs present in the data.
#' @param don_subgenomes Optional donor subgenome labels to restrict the
#'   donor side to.
#' @param rec_subgenomes Optional recipient subgenome labels to restrict the
#'   recipient side to.
#'
#' @return A named character vector: names are donor chromosome IDs, values
#'   are recipient chromosome IDs.
#' @seealso [segments()], [readMaMaMIAInput()]
#' @export
defaultTargetPairs <- function(meta, don_chr_ids, rec_chr_ids,
                               don_subgenomes = NULL, rec_subgenomes = NULL) {
    meta <- as_meta_df(meta)
    don_chr_ids <- as.character(don_chr_ids)
    rec_chr_ids <- as.character(rec_chr_ids)

    don_keep <- meta$chr_id %in% don_chr_ids
    if (!is.null(don_subgenomes)) {
        don_keep <- don_keep & meta$subgenome %in% as.character(don_subgenomes)
    }
    rec_keep <- meta$chr_id %in% rec_chr_ids
    if (!is.null(rec_subgenomes)) {
        rec_keep <- rec_keep & meta$subgenome %in% as.character(rec_subgenomes)
    }

    don <- meta[don_keep, , drop = FALSE]
    rec <- meta[rec_keep, , drop = FALSE]
    if (nrow(don) == 0L || nrow(rec) == 0L) {
        stop("No donor or recipient chromosomes match the requested subgenomes", call. = FALSE)
    }

    don$number <- chr_number(don$chr_name)
    rec$number <- chr_number(rec$chr_name)

    pairs <- character()
    for (i in seq_len(nrow(don))) {
        if (is.na(don$number[[i]])) {
            next
        }
        matched <- rec$chr_id[!is.na(rec$number) & rec$number == don$number[[i]]]
        if (length(matched) == 1L) {
            pairs[[don$chr_id[[i]]]] <- matched
        }
    }
    if (length(pairs) == 0L) {
        stop(
            paste0(
                "No donor->recipient pairs could be matched by chromosome number; ",
                "check the chr_name values or restrict the search with ",
                "`don_subgenomes` and `rec_subgenomes`"
            ),
            call. = FALSE
        )
    }

    pairs
}
