## Read counting through htslib.
##
## countReads() counts arbitrary GRanges with bam_count_ranges_cpp(). The
## pipeline measures exact-N windows with countWindowMetrics(), which uses
## bam_metrics_rg_cpp() in the same file and shares the filtering rules.
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
#' filters through [countWindowMetrics()].
#'
#' `exclude_flags` follows `samtools view -F`: a read is excluded when **any**
#' bit in the mask is set. The default
#' `2308 = 0x904 = unmapped + secondary + supplementary` drops the alignments
#' that cannot contribute coverage.
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
