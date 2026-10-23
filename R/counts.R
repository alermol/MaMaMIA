## Read counting through htslib.
##
## countReads() counts arbitrary GRanges with bam_count_ranges_cpp();
## .count_reads_by_sample() counts the pipeline's exact-N windows with
## bam_count_rg_cpp(). Both live in src/htslib-utils.cpp and share the same
## filtering rules.
##
## `exclude_flags` follows samtools view -F: a read is excluded when *any* bit
## in the mask is set. The htslib loops apply the mask directly, so a combined
## mask costs one pass rather than an inclusion-exclusion expansion.
##
## Reads are assigned to the window containing their 5' end, not to every
## window they overlap: a forward-strand read uses its leftmost coordinate and
## a reverse-strand read uses its rightmost coordinate (the alignment end), so
## a reverse read can belong to the window after the one containing its SAM
## POS. This is documented as the deliberate difference from
## bedtools coverage -counts.

.count_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIACountError", fmt, ...)
}

.count_integer <- function(x, what, minimum, maximum) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
        !is.finite(x) || x != floor(x) || x < minimum || x > maximum) {
        .count_error(
            "`%s` must be a single integer between %s and %s",
            what,
            format(minimum, scientific = FALSE),
            format(maximum, scientific = FALSE)
        )
    }
    as.integer(x)
}

.count_logical <- function(x, what) {
    if (length(x) != 1L || is.na(x) || !is.logical(x)) {
        .count_error("`%s` must be true or false", what)
    }
    x
}

.validate_count_gr <- function(gr) {
    if (!inherits(gr, "GRanges")) {
        .count_error("`gr` must be a GenomicRanges::GRanges object")
    }
    if (length(gr) == 0L) {
        return(invisible(gr))
    }

    seqnames <- as.character(GenomicRanges::seqnames(gr))
    if (anyNA(seqnames) || any(!nzchar(seqnames))) {
        .count_error("`gr` contains an empty or missing seqname")
    }
    start <- IRanges::start(gr)
    end <- IRanges::end(gr)
    if (anyNA(start) || anyNA(end)) {
        .count_error("`gr` contains missing window coordinates")
    }
    if (any(IRanges::width(gr) <= 0L)) {
        bad <- which(IRanges::width(gr) <= 0L)[[1L]]
        .count_error(
            "`gr` contains an empty window: %s:%s-%s",
            seqnames[[bad]], start[[bad]], end[[bad]]
        )
    }

    info <- GenomicRanges::seqinfo(gr)
    seqlengths <- stats::setNames(
        as.numeric(info@seqlengths),
        as.character(GenomicRanges::seqnames(info))
    )
    used <- unique(seqnames)
    missing <- used[!used %in% names(seqlengths)]
    if (length(missing) > 0L) {
        .count_error(
            "`gr` has no seqinfo entry for: %s", paste(missing, collapse = ", ")
        )
    }
    lengths <- as.numeric(seqlengths[used])
    if (anyNA(lengths) || any(!is.finite(lengths)) || any(lengths <= 0)) {
        bad <- used[which(is.na(lengths) | !is.finite(lengths) | lengths <= 0)[[1L]]]
        .count_error(
            "`gr` must define a positive seqlength for '%s'", bad
        )
    }
    length_for_window <- lengths[match(seqnames, used)]
    if (any(start < 1L)) {
        bad <- which(start < 1L)[[1L]]
        .count_error(
            "`gr` window starts before position 1: %s:%s-%s",
            seqnames[[bad]], start[[bad]], end[[bad]]
        )
    }
    if (any(end > length_for_window)) {
        bad <- which(end > length_for_window)[[1L]]
        .count_error(
            "`gr` window %s:%s-%s exceeds the seqlength %s",
            seqnames[[bad]], start[[bad]], end[[bad]],
            format(length_for_window[[bad]], scientific = FALSE)
        )
    }
    invisible(gr)
}

## Convert a GRanges to the arrays consumed by bam_count_ranges_cpp().
.count_reads_cpp <- function(bam_path, gr, mapq, exclude_flags) {
    seqnames <- as.character(GenomicRanges::seqnames(gr))
    contigs <- unique(seqnames)
    info <- GenomicRanges::seqinfo(gr)
    seqlengths <- stats::setNames(
        as.numeric(info@seqlengths),
        as.character(GenomicRanges::seqnames(info))
    )
    tryCatch(
        bam_count_ranges_cpp(
            bam_path,
            contig_names = contigs,
            contig_lengths = as.numeric(seqlengths[contigs]),
            range_contig = match(seqnames, contigs),
            range_start = as.numeric(IRanges::start(gr)) - 1,
            range_end = as.numeric(IRanges::end(gr)),
            mapq = mapq,
            exclude_flags = exclude_flags
        ),
        error = function(e) {
            message <- conditionMessage(e)
            hint <- if (grepl("index|idx", message, ignore.case = TRUE)) {
                "; check the CSI index with ensureBamIndex() or stage the BAM in session temp"
            } else {
                ""
            }
            .count_error(
                "Could not count reads in BAM '%s': %s%s",
                bam_path, message, hint
            )
        }
    )
}

#' Count reads in genomic ranges
#'
#' Counts reads whose 5' end falls inside each range using the package's
#' htslib iterator. A valid CSI index is ensured first. This is the public
#' helper for arbitrary `GRanges`; [analyzeBamInput()] uses the same htslib
#' rules through `.count_reads_by_sample()`.
#'
#' `exclude_flags` follows `samtools view -F`: a read is excluded when **any**
#' bit in the mask is set. The default
#' `2308 = 0x904 = unmapped + secondary + supplementary` matches the legacy
#' `samtools view -F 2308` filter.
#'
#' Counting is 5'-end based: a forward-strand read belongs to the range
#' containing its leftmost coordinate and a reverse-strand read to the range
#' containing its rightmost coordinate (the alignment end). A reverse read can
#' therefore belong to the range after the one containing its SAM `POS`. This
#' differs from `bedtools coverage -counts`, which counts a read in every range
#' it overlaps; the difference is confined to reads whose 5' end lies outside
#' the range while the alignment span overlaps it.
#'
#' For very large range sets, the option `MaMaMIA.count_chunk_size` splits the
#' GRanges into chunks before counting. It is disabled (`0`) by default and is
#' only useful to bound the size of the returned count vector.
#'
#' @param bam_path Path to a coordinate-sorted BAM.
#' @param gr A `GenomicRanges::GRanges` object with positive seqlengths for
#'   every seqname.
#' @param mapq Minimum mapping quality (0-255).
#' @param exclude_flags SAM flag mask; a read is excluded when any bit is set
#'   (0-4095).
#' @param verbose Print a start message. Accepted for compatibility with the
#'   previous backend.
#'
#' @return An integer vector aligned to `gr`.
#'
#' @examples
#' demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
#' bam <- file.path(tempdir(), "demo.bam")
#' file.copy(file.path(demo, "demo.bam"), bam, overwrite = TRUE)
#' gr <- GenomicRanges::GRanges(
#'     c("d1", "r1"), IRanges::IRanges(1L, 20000L),
#'     seqlengths = c(d1 = 20000L, r1 = 20000L)
#' )
#' countReads(bam, gr)
#'
#' @export
countReads <- function(bam_path, gr, mapq = 20L, exclude_flags = 2308L,
                       verbose = FALSE) {
    bam_path <- .validate_bam_path(bam_path)
    mapq <- .count_integer(mapq, "mapq", 0L, 255L)
    exclude_flags <- .count_integer(exclude_flags, "exclude_flags", 0L, 4095L)
    verbose <- .count_logical(verbose, "verbose")
    .validate_count_gr(gr)
    if (length(gr) == 0L) {
        return(integer(0))
    }

    ensureBamIndex(bam_path)
    if (verbose) {
        message(sprintf(
            "Counting reads in %d range(s) from %s",
            length(gr), bam_path
        ))
    }

    chunk_size <- getOption("MaMaMIA.count_chunk_size", 0L)
    if (!is.numeric(chunk_size) || length(chunk_size) != 1L ||
        is.na(chunk_size) || chunk_size < 1 ||
        chunk_size > .Machine$integer.max) {
        chunk_size <- 0L
    } else {
        chunk_size <- as.integer(chunk_size)
    }
    if (chunk_size > 0L && length(gr) > chunk_size) {
        chunks <- split(
            seq_along(gr),
            ceiling(seq_along(gr) / chunk_size)
        )
        counts <- integer(length(gr))
        for (idx in chunks) {
            counts[idx] <- .count_reads_cpp(
                bam_path, gr[idx], mapq, exclude_flags
            )
        }
        return(counts)
    }

    .count_reads_cpp(bam_path, gr, mapq, exclude_flags)
}

#' Count reads per sample in a BAM
#'
#' Counts the 5' end of every alignment in the windows of `windows`. A
#' one-sample BAM counts every read and ignores the `RG` tag; a multi-sample BAM
#' separates the counts by the read group (`RG`) tag of each read. For several
#' samples every read group must be declared in the BAM header, every read group
#' must map to one of `samples`, and reads without an `RG` tag are rejected.
#'
#' This is the counting backend for the BAM-first pipeline, whether the file
#' contains one sample or several.
#'
#' @param bam_path Path to a coordinate-sorted, indexed BAM.
#' @param windows A window list from [makeAnalysisWindows()].
#' @param samples Sample names, in the order of the returned matrix rows.
#' @param read_groups A data frame with `id` and `sample` columns, as returned
#'   by [bamHeaderInfo()].
#' @param mapq Minimum mapping quality (0-255).
#' @param exclude_flags SAM flag mask; a read is excluded when any bit is set.
#'
#' @return An integer matrix with one row per sample and one column per window.
#' @noRd
.count_reads_by_sample <- function(bam_path, windows, samples, read_groups,
                                   mapq = 20L, exclude_flags = 2308L) {
    bam_path <- .validate_bam_path(bam_path)
    mapq <- .count_integer(mapq, "mapq", 0L, 255L)
    exclude_flags <- .count_integer(exclude_flags, "exclude_flags", 0L, 4095L)

    if (!is.list(windows) || is.null(windows$table) ||
        !is.data.frame(windows$table) || is.null(windows$n_per_chromosome)) {
        .count_error("`windows` must be a list returned by makeAnalysisWindows()")
    }
    table <- windows$table
    required <- c("chr_id", "start", "stop")
    if (!all(required %in% names(table)) || nrow(table) == 0L) {
        .count_error("`windows$table` must contain chr_id, start and stop")
    }
    n_per <- .count_integer(
        windows$n_per_chromosome, "n_per_chromosome", 1L, .Machine$integer.max
    )
    chr_ids <- unique(as.character(table$chr_id))
    if (length(chr_ids) == 0L) {
        .count_error("`windows$table` contains no contigs")
    }
    if (!identical(
        as.character(table$chr_id), rep(chr_ids, each = n_per)
    )) {
        .count_error(
            "`windows$table` must list each contig as one contiguous block"
        )
    }
    lengths <- vapply(chr_ids, function(chr_id) {
        rows <- which(table$chr_id == chr_id)
        if (length(rows) != n_per) {
            .count_error(
                "Contig '%s' has %d windows, expected %d",
                chr_id, length(rows), n_per
            )
        }
        start <- as.numeric(table$start[rows])
        stop <- as.numeric(table$stop[rows])
        if (start[[1L]] != 0 || any(start[-1L] != stop[-n_per]) ||
            any(stop <= start)) {
            .count_error(
                "Windows for contig '%s' are not a contiguous half-open tiling",
                chr_id
            )
        }
        stop[[n_per]]
    }, numeric(1L))

    if (!is.character(samples) || length(samples) == 0L || anyNA(samples) ||
        any(!nzchar(samples))) {
        .count_error("`samples` must be a non-empty character vector")
    }
    samples <- as.character(samples)
    if (anyDuplicated(samples)) {
        .count_error("`samples` must be unique")
    }
    if (!is.data.frame(read_groups) ||
        !all(c("id", "sample") %in% names(read_groups))) {
        .count_error("`read_groups` must have id and sample columns")
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
        ## The C++ counter ignores RG tags when there is only one sample, so
        ## a BAM without @RG (or with unrelated read groups) still works.
        read_group_ids <- character()
        sample_index <- integer()
    } else {
        missing_samples <- setdiff(samples, read_groups$sample)
        if (length(missing_samples) > 0L) {
            .count_error(
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

    counts <- tryCatch(
        bam_count_rg_cpp(
            bam_path,
            contig_names = chr_ids,
            contig_lengths = lengths,
            n_per_chromosome = n_per,
            read_group_ids = read_group_ids,
            sample_index = sample_index,
            sample_names = samples,
            mapq = mapq,
            exclude_flags = exclude_flags
        ),
        error = function(e) {
            .count_error(
                "Could not count reads per sample in BAM '%s': %s",
                bam_path, conditionMessage(e)
            )
        }
    )
    if (!is.matrix(counts) || nrow(counts) != length(samples) ||
        ncol(counts) != nrow(table)) {
        .count_error(
            "Per-sample counting returned a matrix of the wrong shape"
        )
    }
    dimnames(counts) <- list(samples, NULL)
    storage.mode(counts) <- "integer"
    counts
}
