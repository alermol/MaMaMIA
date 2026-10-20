## Read counting through bamsignals with samtools-compatible filtering (T3.1).
##
## bamsignals::bamCount(filteredFlag = mask) excludes a read only when *all*
## bits in mask are set, whereas samtools view -F mask excludes a read when
## *any* bit is set. A combined mask therefore cannot be passed directly.
##
## For exclude_flags with bits b1..bk, inclusion-exclusion over the non-empty
## subsets gives the samtools-compatible count:
##
##   pass = sum_{S non-empty} (-1)^(|S| + 1) * C(mask(S))
##
## where C(mask) is bamsignals' count with filteredFlag = mask. For the default
## 2308 = 4 + 256 + 2048 this is seven bamsignals passes. The exponent is
## bounded because every extra bit doubles the work.
##
## Reads are assigned to the window containing their alignment start (5' end),
## not to every window they overlap. This matches bamsignals and is documented
## as a deliberate difference from bedtools coverage -counts.

.MaMaMIACountMaxBits <- 6L

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

.flag_bits <- function(mask) {
    bits <- as.integer(2^(0:11))
    bits[bitwAnd(as.integer(mask), bits) != 0L]
}

.flag_subsets <- function(bits) {
    out <- vector("list", 0L)
    for (k in seq_along(bits)) {
        for (index in utils::combn(seq_along(bits), k, simplify = FALSE)) {
            out[[length(out) + 1L]] <- bits[index]
        }
    }
    out
}

.bam_count_call <- function(bam_path, gr, mapq, filtered_flag, verbose) {
    tryCatch(
        as.numeric(bamsignals::bamCount(
            bam_path, gr,
            mapqual = mapq,
            filteredFlag = filtered_flag,
            ss = FALSE,
            verbose = verbose
        )),
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

.count_reads_bamsignals <- function(bam_path, gr, mapq, exclude_flags,
                                    verbose) {
    if (exclude_flags == 0L) {
        counts <- .bam_count_call(bam_path, gr, mapq, -1L, verbose)
    } else {
        bits <- .flag_bits(exclude_flags)
        if (length(bits) == 0L) {
            counts <- .bam_count_call(bam_path, gr, mapq, -1L, verbose)
        } else {
            if (length(bits) > .MaMaMIACountMaxBits) {
                .count_error(
                    "`exclude_flags` = %d uses %d flag bits; the bamsignals backend supports at most %d bits (each additional bit doubles the number of passes)",
                    exclude_flags, length(bits), .MaMaMIACountMaxBits
                )
            }

            counts <- numeric(length(gr))
            for (subset in .flag_subsets(bits)) {
                mask <- as.integer(sum(subset))
                sign <- if (length(subset) %% 2L == 1L) 1 else -1
                counts <- counts + sign *
                    .bam_count_call(bam_path, gr, mapq, mask, verbose)
            }
        }
    }

    if (any(counts < -1e-6)) {
        .count_error(
            "Internal counting error for BAM '%s': negative count", bam_path
        )
    }
    counts[counts < 0] <- 0
    if (any(abs(counts - round(counts)) > 1e-6)) {
        .count_error(
            "Internal counting error for BAM '%s': non-integer count", bam_path
        )
    }
    if (any(counts > .Machine$integer.max)) {
        .count_error(
            "Counts for BAM '%s' exceed R's integer range", bam_path
        )
    }
    as.integer(round(counts))
}

#' Count reads in genomic windows
#'
#' Counts reads whose alignment start falls inside each window using
#' [bamsignals::bamCount()]. A valid CSI index is ensured first.
#'
#' `exclude_flags` follows `samtools view -F`: a read is excluded when **any**
#' bit in the mask is set. bamsignals excludes a read only when **all** bits are
#' set, so a combined mask is decomposed with inclusion-exclusion. The default
#' `2308 = 0x904 = unmapped + secondary + supplementary` matches the legacy
#' `samtools view -F 2308` filter.
#'
#' Counting is start-based: a read belongs to the window containing its
#' alignment start (5' end), not to every window it overlaps. This differs from
#' `bedtools coverage -counts`, which counts overlaps; the difference is
#' confined to reads that start in the previous window and cross the boundary.
#'
#' For very large window sets, the option `MaMaMIA.count_chunk_size` splits the
#' GRanges into chunks before counting. It is disabled (`0`) by default; each
#' chunk repeats the bamsignals passes, so enable it only to bound memory.
#'
#' @param bam_path Path to a coordinate-sorted BAM.
#' @param gr A `GenomicRanges::GRanges` object with positive seqlengths for
#'   every seqname.
#' @param mapq Minimum mapping quality (0-255).
#' @param exclude_flags SAM flag mask; a read is excluded when any bit is set
#'   (0-4095, at most six bits).
#' @param verbose Print bamsignals progress messages.
#'
#' @return An integer vector aligned to `gr`.
#' @export
countReads <- function(bam_path, gr, mapq = 20L, exclude_flags = 2308L,
                       verbose = FALSE) {
    bam_path <- .validate_bam_path(bam_path)
    mapq <- .count_integer(mapq, "mapq", 0L, 255L)
    exclude_flags <- .count_integer(exclude_flags, "exclude_flags", 0L, 4095L)
    flag_bits <- .flag_bits(exclude_flags)
    if (length(flag_bits) > .MaMaMIACountMaxBits) {
        .count_error(
            "`exclude_flags` = %d uses %d flag bits; the bamsignals backend supports at most %d bits (each additional bit doubles the number of passes)",
            exclude_flags, length(flag_bits), .MaMaMIACountMaxBits
        )
    }
    verbose <- .count_logical(verbose, "verbose")
    .validate_count_gr(gr)
    if (length(gr) == 0L) {
        return(integer(0))
    }

    ensureBamIndex(bam_path)

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
            counts[idx] <- .count_reads_bamsignals(
                bam_path, gr[idx], mapq, exclude_flags, verbose
            )
        }
        return(counts)
    }

    .count_reads_bamsignals(bam_path, gr, mapq, exclude_flags, verbose)
}
