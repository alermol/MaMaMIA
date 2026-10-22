## Rhtslib BAM header reading and sample resolution (D2, D8).
##
## The header is parsed exactly once per BAM, before any index or count work.
## Samples are read from the @RG SM tags, read groups from the @RG ID tags. A
## BAM may contain one or several samples; every read group must be complete
## and every read group ID must be unique. A BAM without any @RG line keeps
## the legacy single-sample fallback: the file basename without .bam.

#' Construct a MaMaMIA error condition
#'
#' Internal helper shared by the BAM-first modules. Every error carries the
#' base `MaMaMIAError` class plus a specific subclass, and may carry structured
#' data (for example candidate assembly matches) in `$data`.
#'
#' @param class Specific condition class.
#' @param fmt `sprintf()` format string.
#' @param ... Values for `fmt`.
#' @param data Optional structured payload.
#'
#' @return Never returns; raises a condition.
#' @noRd
.mamamia_error <- function(class, fmt, ..., data = NULL) {
    message <- sprintf(fmt, ...)
    condition <- structure(
        class = c(class, "MaMaMIAError", "error", "condition"),
        list(message = message, call = NULL, data = data)
    )
    stop(condition)
}

.bam_header_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIABamHeaderError", fmt, ...)
}

#' Convert a path to an absolute path without resolving a final symlink
#'
#' `normalizePath()` resolves symlinks, which is undesirable for staged BAMs:
#' htslib must write the CSI next to the symlink in session temp. This helper
#' only resolves the containing directory.
#'
#' @param path A file path.
#'
#' @return An absolute path, or the input when the directory does not exist.
#' @noRd
.absolute_path <- function(path) {
    directory <- dirname(path)
    resolved <- tryCatch(
        normalizePath(directory, winslash = "/", mustWork = TRUE),
        error = function(e) NULL
    )
    if (is.null(resolved)) {
        return(path)
    }
    file.path(resolved, basename(path))
}

.validate_bam_path <- function(path, class = "MaMaMIABamHeaderError") {
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(path)) {
        .mamamia_error(class, "`path` must be a single non-empty file path")
    }
    path <- path.expand(path)
    if (!file.exists(path)) {
        .mamamia_error(class, "BAM file does not exist: %s", path)
    }
    if (dir.exists(path)) {
        .mamamia_error(class, "BAM path is a directory, not a file: %s", path)
    }
    size <- file.info(path)$size
    if (length(size) != 1L || is.na(size) || size <= 0) {
        .mamamia_error(class, "BAM file is empty: %s", path)
    }
    .absolute_path(path)
}

.header_lines <- function(header_text) {
    if (length(header_text) != 1L || is.na(header_text) ||
        !nzchar(header_text)) {
        return(character(0))
    }
    lines <- strsplit(header_text, "\n", fixed = TRUE)[[1L]]
    lines <- sub("\r$", "", lines)
    lines[nzchar(lines)]
}

.header_records <- function(lines, tag) {
    prefix <- paste0("@", tag)
    lines[lines == prefix | startsWith(lines, paste0(prefix, "\t"))]
}

.header_fields <- function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
    if (length(fields) <= 1L) {
        return(stats::setNames(character(0), character(0)))
    }
    fields <- fields[-1L]
    keys <- sub(":.*$", "", fields)
    values <- sub("^[^:]*:", "", fields)
    stats::setNames(values, keys)
}

.header_field <- function(fields, key) {
    if (!key %in% names(fields)) {
        return(NA_character_)
    }
    value <- unname(fields[[key]])
    if (length(value) != 1L || is.na(value)) {
        return(NA_character_)
    }
    value
}

.parse_read_groups <- function(lines) {
    records <- .header_records(lines, "RG")
    if (length(records) == 0L) {
        return(data.frame(
            id = character(0),
            sample = character(0),
            stringsAsFactors = FALSE
        ))
    }
    rows <- lapply(records, function(line) {
        fields <- .header_fields(line)
        data.frame(
            id = .header_field(fields, "ID"),
            sample = .header_field(fields, "SM"),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

.parse_programs <- function(lines) {
    records <- .header_records(lines, "PG")
    if (length(records) == 0L) {
        return(data.frame(
            id = character(0),
            name = character(0),
            version = character(0),
            command = character(0),
            stringsAsFactors = FALSE
        ))
    }
    rows <- lapply(records, function(line) {
        fields <- .header_fields(line)
        data.frame(
            id = .header_field(fields, "ID"),
            name = .header_field(fields, "PN"),
            version = .header_field(fields, "VN"),
            command = .header_field(fields, "CL"),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

.sort_order <- function(lines, path) {
    records <- .header_records(lines, "HD")
    if (length(records) == 0L) {
        .bam_header_error(
            "BAM '%s' has no @HD header line; the file must be coordinate-sorted before analysis",
            path
        )
    }
    fields <- .header_fields(records[[1L]])
    sort_order <- .header_field(fields, "SO")
    if (!is.na(sort_order)) {
        sort_order <- trimws(sort_order)
    }
    if (is.na(sort_order) || !nzchar(sort_order)) {
        .bam_header_error(
            "BAM '%s' has no @HD SO field; the file must be coordinate-sorted before analysis",
            path
        )
    }
    if (!identical(tolower(sort_order), "coordinate")) {
        .bam_header_error(
            "BAM '%s' is not coordinate-sorted (@HD SO:%s); coordinate-sort it with 'samtools sort' first",
            path, sort_order
        )
    }
    sort_order
}

.resolve_samples <- function(path, read_groups) {
    if (nrow(read_groups) == 0L) {
        sample <- sub("\\.bam$", "", basename(path), ignore.case = TRUE)
        if (!nzchar(sample)) {
            .mamamia_error(
                "MaMaMIASampleError",
                paste(
                    "Could not resolve a sample name for BAM '%s':",
                    "no @RG SM value and the file name yields an empty basename"
                ),
                path
            )
        }
        return(list(
            samples = sample,
            source = "basename",
            read_groups = data.frame(
                id = "<basename>",
                sample = sample,
                stringsAsFactors = FALSE
            )
        ))
    }

    ids <- trimws(read_groups$id)
    samples <- trimws(read_groups$sample)
    missing_id <- is.na(ids) | !nzchar(ids)
    if (any(missing_id)) {
        .mamamia_error(
            "MaMaMIASampleError",
            paste(
                "BAM '%s' has an @RG line without a read group ID (@RG ID);",
                "every read group needs both ID and SM tags"
            ),
            path
        )
    }
    missing_sample <- is.na(samples) | !nzchar(samples)
    if (any(missing_sample)) {
        .mamamia_error(
            "MaMaMIASampleError",
            paste(
                "BAM '%s' has an @RG line without a sample ID (@RG SM);",
                "every read group needs both ID and SM tags"
            ),
            path
        )
    }
    duplicated_ids <- unique(ids[duplicated(ids)])
    if (length(duplicated_ids) > 0L) {
        .mamamia_error(
            "MaMaMIASampleError",
            "BAM '%s' has duplicate @RG ID value%s: %s",
            path,
            if (length(duplicated_ids) == 1L) "" else "s",
            paste(sprintf("'%s'", duplicated_ids), collapse = ", ")
        )
    }

    read_groups$id <- ids
    read_groups$sample <- samples
    list(
        samples = unique(samples),
        source = "rg",
        read_groups = read_groups
    )
}

#' Read BAM header metadata
#'
#' Parses the header of a coordinate-sorted BAM with Rhtslib, without reading
#' any alignments or requiring an index. Samples come from the non-empty
#' `@RG SM` tags and read groups from `@RG ID`. A BAM may contain one or several
#' samples, and several read groups may share one sample. Every `@RG` line must
#' contain both `ID` and `SM`, and the IDs must be unique; otherwise a BAM
#' without any `@RG` line falls back to `basename(path)` without the `.bam`
#' suffix.
#'
#' @param path Path to a BAM file.
#'
#' @return A list with `path`, `format`, `targets` (a data frame with `name`
#'   and `length`), `samples`, `read_groups` (`id`, `sample`), `programs`
#'   (`id`, `name`, `version`, `command`), `header_text`, `sample` (the only
#'   sample, or `NA_character_` when the BAM contains several samples),
#'   `sample_source` (`"rg"` or `"basename"`), `multi_sample`, `sort_order`
#'   and `n_targets`.
#'
#' @examples
#' demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
#' header <- bamHeaderInfo(file.path(demo, "demo.bam"))
#' header$samples
#' header$targets
#'
#' @export
bamHeaderInfo <- function(path) {
    path <- .validate_bam_path(path)

    info <- tryCatch(
        bam_header_info_cpp(path),
        error = function(e) {
            message <- conditionMessage(e)
            if (grepl("^Unsupported input format", message)) {
                .mamamia_error(
                    "MaMaMIAUnsupportedFormatError", "%s", message
                )
            }
            .bam_header_error("%s", message)
        }
    )

    header_text <- info$header_text
    lines <- .header_lines(header_text)
    sort_order <- .sort_order(lines, path)

    targets <- info$targets
    if (!is.data.frame(targets) || nrow(targets) == 0L) {
        .bam_header_error(
            "BAM '%s' contains no target sequences in its header", path
        )
    }
    if (anyDuplicated(targets$name)) {
        duplicate <- unique(targets$name[duplicated(targets$name)])[[1L]]
        .bam_header_error(
            "BAM '%s' has duplicate target sequence name '%s' in its header",
            path, duplicate
        )
    }

    read_groups <- .parse_read_groups(lines)
    programs <- .parse_programs(lines)
    resolved <- .resolve_samples(path, read_groups)
    samples <- resolved$samples

    list(
        path = path,
        format = "BAM",
        targets = targets,
        samples = samples,
        read_groups = resolved$read_groups,
        programs = programs,
        header_text = header_text,
        sample = if (length(samples) == 1L) samples[[1L]] else NA_character_,
        sample_source = resolved$source,
        multi_sample = length(samples) > 1L,
        sort_order = sort_order,
        n_targets = as.integer(nrow(targets))
    )
}
