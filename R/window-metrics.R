## WindowMetrics: breadth and counts in the exact-N analysis windows.
##
## `countWindowMetrics()` is the R wrapper around bam_metrics_rg_cpp(). It
## returns one row per sample and window with both the alignment count and the
## number of covered bases, so breadth and count depth come from a single pass
## over the BAM. `buildWindowMetrics()` joins those measurements with the
## window metadata and the reference GC and derives the two depth scores.

.window_metrics_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAWindowMetricsError", fmt, ...)
}

.window_metrics_integer <- function(x, what, minimum, maximum) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
        !is.finite(x) || x != floor(x) || x < minimum || x > maximum) {
        .window_metrics_error(
            "`%s` must be a single integer between %s and %s",
            what,
            format(minimum, scientific = FALSE),
            format(maximum, scientific = FALSE)
        )
    }
    as.integer(x)
}

.validate_metric_windows <- function(windows) {
    if (!is.list(windows) || is.null(windows$table) ||
        !is.data.frame(windows$table) || is.null(windows$n_per_chromosome)) {
        .window_metrics_error(
            "`windows` must be a list returned by makeAnalysisWindows()"
        )
    }
    table <- windows$table
    required <- c("chr_id", "start", "stop", "role", "subgenome", "group")
    missing <- setdiff(required, names(table))
    if (length(missing) > 0L) {
        .window_metrics_error(
            "`windows$table` is missing column(s): %s",
            paste(missing, collapse = ", ")
        )
    }
    n_per <- .window_metrics_integer(
        windows$n_per_chromosome, "n_per_chromosome", 1L, .Machine$integer.max
    )
    chr_ids <- unique(as.character(table$chr_id))
    if (length(chr_ids) == 0L) {
        .window_metrics_error("`windows$table` contains no contigs")
    }
    if (!identical(
        as.character(table$chr_id), rep(chr_ids, each = n_per)
    )) {
        .window_metrics_error(
            "`windows$table` must list each contig as one contiguous block"
        )
    }
    lengths <- vapply(chr_ids, function(chr_id) {
        rows <- which(table$chr_id == chr_id)
        if (length(rows) != n_per) {
            .window_metrics_error(
                "Contig '%s' has %d windows, expected %d",
                chr_id, length(rows), n_per
            )
        }
        start <- as.numeric(table$start[rows])
        stop <- as.numeric(table$stop[rows])
        if (start[[1L]] != 0 || any(start[-1L] != stop[-n_per]) ||
            any(stop <= start)) {
            .window_metrics_error(
                "Windows for contig '%s' are not a contiguous half-open tiling",
                chr_id
            )
        }
        stop[[n_per]]
    }, numeric(1L))
    list(table = table, n_per = n_per, chr_ids = chr_ids, lengths = lengths)
}

#' Count breadth and read counts per sample in a BAM
#'
#' Counts, for every exact-N analysis window, the number of alignments with at
#' least one aligned base in the window and the number of reference bases
#' covered by at least one alignment. Both values come from one htslib pass.
#' A one-sample BAM counts every read and ignores the `RG` tag; a multi-sample
#' BAM separates samples through `@RG` (`ID` -> `SM`) exactly like
#' [preflightBams()].
#'
#' @param bam_path Path to a coordinate-sorted, indexed BAM.
#' @param windows A window list from [makeAnalysisWindows()].
#' @param samples Sample names in the order of the returned rows.
#' @param read_groups A data frame with `id` and `sample` columns from
#'   [bamHeaderInfo()].
#' @param mapq Minimum mapping quality (0-255).
#' @param exclude_flags SAM flag mask; a read is excluded when any bit is set.
#' @param read_length Known read length, or `0` to infer it from the first
#'   passing alignment.
#' @param count_deletions Treat `D` CIGAR blocks as covered bases.
#'
#' @return A list with `metrics` (a long data frame with `sample`, `chr_id`,
#'   `iid`, `count` and `covered`), `counts` and `covered` matrices (samples x
#'   windows), `read_length` and `n_per_chromosome`.
#' @export
countWindowMetrics <- function(bam_path, windows, samples, read_groups,
                               mapq = 20L, exclude_flags = 1796L,
                               read_length = 0L, count_deletions = TRUE) {
    bam_path <- .validate_bam_path(bam_path)
    mapq <- .window_metrics_integer(mapq, "mapq", 0L, 255L)
    exclude_flags <- .window_metrics_integer(
        exclude_flags, "exclude_flags", 0L, 4095L
    )
    read_length <- .window_metrics_integer(
        read_length, "read_length", 0L, .Machine$integer.max
    )
    if (length(count_deletions) != 1L || is.na(count_deletions) ||
        !is.logical(count_deletions)) {
        .window_metrics_error("`count_deletions` must be true or false")
    }
    validated <- .validate_metric_windows(windows)
    table <- validated$table
    n_per <- validated$n_per
    chr_ids <- validated$chr_ids
    lengths <- validated$lengths

    if (!is.character(samples) || length(samples) == 0L || anyNA(samples) ||
        any(!nzchar(samples))) {
        .window_metrics_error("`samples` must be a non-empty character vector")
    }
    samples <- as.character(samples)
    if (anyDuplicated(samples)) {
        .window_metrics_error("`samples` must be unique")
    }
    if (!is.data.frame(read_groups) ||
        !all(c("id", "sample") %in% names(read_groups))) {
        .window_metrics_error("`read_groups` must have id and sample columns")
    }
    read_groups$id <- as.character(read_groups$id)
    read_groups$sample <- as.character(read_groups$sample)
    read_groups <- read_groups[
        !is.na(read_groups$id) & nzchar(read_groups$id) &
            !is.na(read_groups$sample) & nzchar(read_groups$sample),
        ,
        drop = FALSE
    ]
    if (length(samples) == 1L) {
        read_group_ids <- character()
        sample_index <- integer()
    } else {
        missing_samples <- setdiff(samples, read_groups$sample)
        if (length(missing_samples) > 0L) {
            .window_metrics_error(
                "No @RG line maps sample%s %s to a read group",
                if (length(missing_samples) == 1L) "" else "s",
                paste(sprintf("'%s'", missing_samples), collapse = ", ")
            )
        }
        read_groups <- read_groups[
            read_groups$sample %in% samples, , drop = FALSE
        ]
        read_group_ids <- read_groups$id
        sample_index <- match(read_groups$sample, samples)
    }

    result <- tryCatch(
        bam_metrics_rg_cpp(
            bam_path,
            contig_names = chr_ids,
            contig_lengths = lengths,
            n_per_chromosome = n_per,
            read_group_ids = read_group_ids,
            sample_index = sample_index,
            sample_names = samples,
            mapq = mapq,
            exclude_flags = exclude_flags,
            read_length = read_length,
            count_deletions = count_deletions
        ),
        error = function(e) {
            .window_metrics_error(
                "Could not measure reads per sample in BAM '%s': %s",
                bam_path, conditionMessage(e)
            )
        }
    )
    counts <- result$counts
    covered <- result$covered
    if (!is.matrix(counts) || !is.matrix(covered) ||
        !identical(dim(counts), dim(covered)) ||
        nrow(counts) != length(samples) ||
        ncol(counts) != nrow(table)) {
        .window_metrics_error(
            "Per-sample counting returned matrices of the wrong shape"
        )
    }
    dimnames(counts) <- list(samples, NULL)
    dimnames(covered) <- list(samples, NULL)
    storage.mode(counts) <- "integer"
    storage.mode(covered) <- "integer"

    n_windows <- nrow(table)
    ## The label columns are sample-major (one sample's whole window grid,
    ## then the next). The matrices are windows x samples once transposed, so
    ## flattening them column-major gives the same ordering. Flattening the
    ## untransposed matrix would interleave samples per window and attach
    ## every value to the wrong sample.
    metrics <- data.frame(
        sample = rep(samples, each = n_windows),
        chr_id = rep(as.character(table$chr_id), times = length(samples)),
        iid = rep(rep(seq.int(0L, n_per - 1L), times = length(chr_ids)),
            times = length(samples)),
        count = as.integer(t(counts)),
        covered = as.integer(t(covered)),
        stringsAsFactors = FALSE
    )

    list(
        metrics = metrics,
        counts = counts,
        covered = covered,
        read_length = as.numeric(result$read_length),
        n_per_chromosome = n_per
    )
}

#' Build the WindowMetrics table
#'
#' Joins per-sample count/covered measurements with the window metadata and the
#' reference GC, then derives breadth, `lambda_breadth = -log(1 - breadth)` and
#' `lambda_count = count * read_length / length`.
#'
#' @param measured A result from [countWindowMetrics()], or a data frame with
#'   `sample`, `chr_id`, `iid`, `count`, `covered` and an optional
#'   `read_length` column.
#' @param windows A window list from [makeAnalysisWindows()].
#' @param gc A GC table from `extractWindowGc()` with `chr_id`, `start`,
#'   `stop` and `gc`.
#' @param read_length Fallback read length when the measurements do not carry
#'   one.
#'
#' @return A `WindowMetrics` data frame with one row per sample and window.
#' @export
buildWindowMetrics <- function(measured, windows, gc, read_length = NULL) {
    validated <- .validate_metric_windows(windows)
    table <- validated$table
    n_per <- validated$n_per
    chr_ids <- validated$chr_ids

    measured_read_length <- NULL
    if (is.list(measured) && !is.data.frame(measured) &&
        is.data.frame(measured$metrics)) {
        measured_read_length <- measured$read_length
        measured <- measured$metrics
        if (!"read_length" %in% names(measured) &&
            !is.null(measured_read_length) &&
            is.finite(measured_read_length) && measured_read_length > 0) {
            measured$read_length <- as.numeric(measured_read_length)
        }
    }
    if (!is.data.frame(measured)) {
        .window_metrics_error(
            "`measured` must be a countWindowMetrics() result or a data frame"
        )
    }
    required <- c("sample", "chr_id", "iid", "count", "covered")
    missing <- setdiff(required, names(measured))
    if (length(missing) > 0L) {
        .window_metrics_error(
            "`measured` is missing column(s): %s", paste(missing, collapse = ", ")
        )
    }
    if (nrow(measured) == 0L) {
        .window_metrics_error("`measured` has no rows")
    }
    if (!"read_length" %in% names(measured)) {
        if (is.null(read_length) || length(read_length) != 1L ||
            is.na(read_length) || !is.finite(read_length) || read_length <= 0) {
            .window_metrics_error(
                "A positive `read_length` is required to derive lambda_count"
            )
        }
        measured$read_length <- as.numeric(read_length)
    }

    if (!is.data.frame(gc) || !all(c("chr_id", "start", "stop", "gc") %in%
        names(gc))) {
        .window_metrics_error(
            "`gc` must be a data frame with chr_id, start, stop and gc"
        )
    }
    gc_key <- paste(gc$chr_id, gc$start, gc$stop, sep = ":")
    window_key <- paste(table$chr_id, table$start, table$stop, sep = ":")
    gc_position <- match(window_key, gc_key)
    if (anyNA(gc_position)) {
        .window_metrics_error(
            "GC is missing for %d window(s); extract GC with the same window table",
            sum(is.na(gc_position))
        )
    }
    table$gc <- as.numeric(gc$gc[gc_position])

    measured$sample <- as.character(measured$sample)
    measured$chr_id <- as.character(measured$chr_id)
    measured$iid <- as.integer(measured$iid)
    measured$count <- as.numeric(measured$count)
    measured$covered <- as.numeric(measured$covered)
    measured$read_length <- as.numeric(measured$read_length)
    if (anyNA(measured$sample) || any(!nzchar(measured$sample))) {
        .window_metrics_error("`measured$sample` contains an empty value")
    }
    if (anyNA(measured$iid) || any(measured$iid < 0L) ||
        any(measured$iid >= n_per)) {
        .window_metrics_error("`measured$iid` is outside the window range")
    }
    if (anyNA(measured$count) || any(measured$count < 0) ||
        anyNA(measured$covered) || any(measured$covered < 0)) {
        .window_metrics_error("Counts and covered bases must be non-negative")
    }

    position <- match(
        paste(measured$chr_id, measured$iid, sep = ":"),
        paste(table$chr_id, rep(seq.int(0L, n_per - 1L), times = length(chr_ids)),
            sep = ":")
    )
    if (anyNA(position)) {
        .window_metrics_error(
            "`measured` contains %d window(s) that are not in `windows`",
            sum(is.na(position))
        )
    }
    meta <- table[position, , drop = FALSE]
    out <- data.frame(
        sample = measured$sample,
        chr_id = meta$chr_id,
        chr_name = meta$chr_name,
        group = meta$group,
        role = meta$role,
        subgenome = meta$subgenome,
        copy_number = as.numeric(meta$copy_number),
        iid = measured$iid,
        start = as.numeric(meta$start),
        stop = as.numeric(meta$stop),
        length = as.numeric(meta$stop) - as.numeric(meta$start),
        count = as.integer(measured$count),
        covered = as.integer(measured$covered),
        gc = as.numeric(meta$gc),
        read_length = measured$read_length,
        stringsAsFactors = FALSE
    )
    out$breadth <- out$covered / out$length
    out$breadth[!is.finite(out$breadth)] <- 0
    out$breadth <- pmin(pmax(out$breadth, 0), 1)
    epsilon <- 1e-12
    out$lambda_breadth <- -log(pmax(1 - out$breadth, epsilon))
    out$lambda_count <- out$count * out$read_length / out$length
    out$lambda_count[!is.finite(out$lambda_count)] <- 0
    out <- out[order(
        out$sample,
        match(out$chr_id, chr_ids),
        out$iid
    ), , drop = FALSE]
    rownames(out) <- NULL
    out[, c(
        "sample", "chr_id", "chr_name", "group", "role", "subgenome",
        "copy_number", "iid", "start", "stop", "length", "count", "covered", "breadth",
        "gc", "lambda_breadth", "lambda_count", "read_length"
    )]
}
