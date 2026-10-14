## reverseWindows(): flip the effective orientation of chromosomes that were
## assembled in opposite orientations between the donor and the recipient
## assembly. Methods live next to the generic because they share the value
## vectors and the pair bookkeeping they touch.

#' Reverse the effective orientation of selected chromosomes
#'
#' @description
#' When comparing two assemblies, homologous chromosomes can be assembled in
#' opposite orientations. As a result, the per-window coverage and GC vectors of
#' such chromosomes run in opposite directions along the genome.
#'
#' `reverseWindows()` flips the effective orientation of the requested
#' chromosomes so that windows can be compared position by position. It keeps
#' the window indexes (`iid`, `start`, `stop`, and the chromosome identifier)
#' unchanged, and instead reverses the value vectors (coverage, GC-corrected
#' coverage when present, and the per-window validity flags) within each
#' selected chromosome. The (modified) object is returned.
#'
#' For an `ISA` object, segmentation is **recomputed** for every
#' donor-recipient pair that contains a reversed chromosome, because these
#' pairs' difference profiles have changed orientation. Segments of unaffected
#' pairs are preserved unchanged.
#'
#' @param x An object of class `RCA` or `ISA`.
#' @param chr_ids Character vector of chromosome IDs in the object's meta table
#'   (`chr_id` column) whose orientation should be reversed. The short display
#'   name (`chr_name`), used mainly for plotting labels, must not be used here.
#' @param ... Further arguments passed to methods.
#'
#' @return The modified object, of the same class as `x`.
#'
#' @examples
#' \dontrun{
#' rev_rca <- reverseWindows(rca, chr_ids = c("OY997261.1", "CP155614.1"))
#' rev_isa <- reverseWindows(isa, chr_ids = "OY997261.1")
#' }
#' @export
reverseWindows <- function(x, chr_ids, ...) {
    UseMethod("reverseWindows")
}

#' @rdname reverseWindows
#' @export
reverseWindows.RCA <- function(x, chr_ids, ...) {
    stopifnot("Input is not RCA object" = inherits(x, "RCA"))

    chr_ids <- as.character(chr_ids)
    if (length(chr_ids) == 0L) {
        warning("No chromosomes requested; returning input unchanged",
            call. = FALSE
        )
        return(x)
    }
    unknown <- setdiff(chr_ids, unique(as.character(x$data$chr_id)))
    stopifnot(
        "Some chromosome IDs in `chr_ids` are not present in the object" =
            length(unknown) == 0L
    )

    value_cols <- intersect(
        names(x$data),
        c("cov", "gc", "cor.gc", "valid", "ideal")
    )
    x$data <- reverse_chr_value_cols(x$data, chr_ids, value_cols)

    return(x)
}

#' @rdname reverseWindows
#' @export
#' @importFrom rlang .data
reverseWindows.ISA <- function(x, chr_ids, ...) {
    stopifnot("Input is not ISA object" = inherits(x, "ISA"))

    chr_ids <- as.character(chr_ids)
    if (length(chr_ids) == 0L) {
        warning("No chromosomes requested; returning input unchanged",
            call. = FALSE
        )
        return(x)
    }
    known <- c(
        as.character(x$meta$don_chr_ids),
        as.character(x$meta$rec_chr_ids)
    )
    unknown <- setdiff(chr_ids, known)
    stopifnot(
        "Some chromosome IDs in `chr_ids` are not present in the object" =
            length(unknown) == 0L
    )

    don_targets <- intersect(chr_ids, as.character(x$meta$don_chr_ids))
    rec_targets <- intersect(chr_ids, as.character(x$meta$rec_chr_ids))

    don_cols <- intersect(
        names(x$out),
        paste0(c("cov", "gc", "cor.gc", "valid", "ideal"), ".don")
    )
    rec_cols <- intersect(
        names(x$out),
        paste0(c("cov", "gc", "cor.gc", "valid", "ideal"), ".rec")
    )

    # Reverse the donor-side value vectors for the requested donor chromosomes,
    # within each (donor, recipient) pair so row order is preserved.
    if (!is.null(x$out) && nrow(x$out) > 0L &&
        length(don_targets) > 0L && length(don_cols) > 0L) {
        x$out <- reverse_paired_side(
            x$out, chr_id_side = "chr_id.don",
            target_ids = don_targets, value_cols = don_cols
        )
    }

    # Reverse the recipient-side value vectors.
    if (!is.null(x$out) && nrow(x$out) > 0L &&
        length(rec_targets) > 0L && length(rec_cols) > 0L) {
        x$out <- reverse_paired_side(
            x$out, chr_id_side = "chr_id.rec",
            target_ids = rec_targets, value_cols = rec_cols
        )
    }

    # Recompute the derived difference so it stays consistent with the reversed
    # cor.gc vectors.
    if (!is.null(x$out) && nrow(x$out) > 0L) {
        x$out <- dplyr::mutate(
            x$out,
            diff = .data$cor.gc.don - .data$cor.gc.rec
        )
    }

    # Re-run CBS for every pair containing a reversed chromosome, since those
    # pairs' difference profiles have changed orientation. Unaffected pairs'
    # segments are kept as-is.
    affected_pairs <- dplyr::filter(
        x$meta$target_pairs,
        .data$chr_id.don %in% chr_ids | .data$chr_id.rec %in% chr_ids
    )
    if (!is.null(x$out) && nrow(x$out) > 0L && nrow(affected_pairs) > 0L) {
        affected_segments <- purrr::pmap(
            list(affected_pairs$chr_id.don, affected_pairs$chr_id.rec),
            \(don, rec) {
                pair_rows <- dplyr::filter(
                    x$out,
                    .data$chr_id.don == !!don & .data$chr_id.rec == !!rec
                )
                segment_pair_data(
                    pair_rows,
                    don = don,
                    rec = rec,
                    meta = x$meta$meta,
                    param = x$param,
                    colname = x$meta$metrics,
                    threshold = x$meta$metrics_threshold,
                    min_width = x$param$min.width
                )
            }
        )
        x$segments <- dplyr::bind_rows(
            dplyr::anti_join(
                x$segments,
                affected_pairs,
                by = c("chr_id.don", "chr_id.rec")
            ),
            purrr::compact(affected_segments)
        )
    }

    # Keep the underlying per-window table aligned as well.
    if (!is.null(x$data) && nrow(x$data) > 0L) {
        value_cols <- intersect(
            names(x$data),
            c("cov", "gc", "cor.gc", "valid", "ideal")
        )
        x$data <- reverse_chr_value_cols(x$data, chr_ids, value_cols)
    }

    return(x)
}

## Internal reversal helpers ---------------------------------------------------

#' @noRd
reverse_chr_value_cols <- function(df, chr_ids, cols) {
    if (length(cols) == 0L) {
        return(df)
    }
    for (id in chr_ids) {
        idx <- which(df$chr_id == id)
        if (length(idx) < 2L) {
            next
        }
        df[idx, cols] <- df[rev(idx), cols]
    }
    df
}

#' @noRd
reverse_paired_side <- function(df, chr_id_side, target_ids, value_cols) {
    if (length(value_cols) == 0L) {
        return(df)
    }
    target_rows <- which(as.character(df[[chr_id_side]]) %in% target_ids)
    if (length(target_rows) == 0L) {
        return(df)
    }

    pair <- paste(
        as.character(df$chr_id.don[target_rows]),
        as.character(df$chr_id.rec[target_rows])
    )
    for (g in unique(pair)) {
        gidx <- target_rows[pair == g]
        if (length(gidx) > 1L) {
            df[gidx, value_cols] <- df[rev(gidx), value_cols]
        }
    }
    df
}
