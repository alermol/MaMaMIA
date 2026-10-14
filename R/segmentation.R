## Pairwise segmentation: join each donor chromosome with its recipient
## counterpart on the corrected coverage difference, cut the difference profile
## into segments with circular binary segmentation (CBS), and flag the segments
## whose summary metric passes the call threshold.

#' @export
#' @importFrom rlang .data
segments <- function(RCA,
                     target_pairs,
                     min_segment_median = 0L,
                     min_segment_mean = 0L,
                     threshold_type = c("median", "mean"),
                     alpha = 0.001,
                     min_width = 2,
                     undo_SD = 1,
                     seed = 1L,
                     verbose = TRUE) {
    stopifnot("Input is not RCA object" = inherits(RCA, "RCA"))
    stopifnot("Coverage must be corrected before segmentation" = RCA$corrected)

    thr_type <- match.arg(threshold_type)
    threshold <- switch(thr_type,
        "mean" = min_segment_mean,
        "median" = min_segment_median
    )
    colname <- switch(thr_type,
        "mean" = "seg.mean",
        "median" = "seg.median"
    )

    param <- validate_segmentation_params(alpha, min_width, undo_SD, seed)

    target_pairs_df <- validate_target_pairs(target_pairs, RCA$meta$don_chr_ids, RCA$meta$rec_chr_ids)

    by_chr_split <- split(RCA$data, RCA$data$chr_id)

    if (verbose) {
        message("Performing segmentation...")
    }

    ## Windows of one donor-recipient pair, matched by window index (`iid`),
    ## restricted to the windows that are ideal on both sides.
    pair_data <- function(don_chr, rec_chr) {
        dplyr::inner_join(
            by_chr_split[[as.character(don_chr)]],
            by_chr_split[[as.character(rec_chr)]],
            by = "iid",
            suffix = c(".don", ".rec")
        ) |>
            dplyr::filter(.data$ideal.don & .data$ideal.rec) |>
            dplyr::mutate(diff = .data$cor.gc.don - .data$cor.gc.rec)
    }

    prepared_data <- target_pairs_df |>
        dplyr::mutate(
            data = purrr::map2(.data$chr_id.don, .data$chr_id.rec, pair_data),
            segments = purrr::pmap(
                list(.data$data, .data$chr_id.don, .data$chr_id.rec),
                segment_pair_data,
                meta = RCA$meta$meta,
                param = param,
                colname = colname,
                threshold = threshold,
                min_width = min_width
            ),
            .keep = "unused"
        )

    keep_don <- unique(target_pairs_df$chr_id.don)
    keep_rec <- unique(target_pairs_df$chr_id.rec)
    keep_chrs <- c(keep_don, keep_rec)

    structure(
        list(
            data = dplyr::filter(RCA$data, .data$chr_id %in% keep_chrs),
            out = dplyr::bind_rows(prepared_data$data),
            segments = dplyr::bind_rows(
                empty_segments(),
                purrr::compact(prepared_data$segments)
            ),
            param = param,
            meta = list(
                meta = dplyr::filter(RCA$meta$meta, .data$chr_id %in% keep_chrs),
                don_chr_ids = keep_don,
                rec_chr_ids = keep_rec,
                target_pairs = target_pairs_df,
                metrics = colname,
                metrics_threshold = threshold
            )
        ),
        class = "ISA"
    )
}

## Internal segmentation helpers ----------------------------------------------

#' @noRd
validate_segmentation_params <- function(alpha, min_width, undo_SD, seed = 1L) {
    if (alpha <= 0 || alpha > 1) {
        stop("1 >= alpha > 0 is not satisfied", call. = FALSE)
    }

    if (min_width < 2) {
        stop("min_width must be >=2", call. = FALSE)
    }

    if (undo_SD <= 0 || undo_SD > 10) {
        stop("10 >= undo_SD > 0 is not satisfied", call. = FALSE)
    }

    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
        stop("seed must be a single number", call. = FALSE)
    }

    list(
        alpha = alpha,
        undo.SD = undo_SD,
        min.width = min_width,
        seed = as.integer(seed)
    )
}

#' @noRd
validate_target_pairs <- function(target_pairs, don_chr_ids, rec_chr_ids) {
    if (length(target_pairs) < 1L || is.null(names(target_pairs)) ||
        any(!nzchar(names(target_pairs)))) {
        stop(
            "`target_pairs` must be a non-empty named vector of donor -> recipient chr_id",
            call. = FALSE
        )
    }

    target_don <- as.character(names(target_pairs))
    target_rec <- as.character(unname(target_pairs))
    stopifnot(
        "Some donor chromosome IDs in target_pairs are not present" =
            all(target_don %in% don_chr_ids)
    )
    stopifnot(
        "Some recipient chromosome IDs in target_pairs are not present" =
            all(target_rec %in% rec_chr_ids)
    )

    pairs <- data.frame(
        chr_id.don = target_don,
        chr_id.rec = target_rec,
        stringsAsFactors = FALSE
    )
    pairs <- pairs[!duplicated(pairs), , drop = FALSE]
    rownames(pairs) <- seq_len(nrow(pairs))
    pairs
}

#' Run code under a fixed seed, then restore the caller's RNG state.
#' @noRd
with_local_seed <- function(seed, code) {
    has_state <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    previous <- if (has_state) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(
        {
            if (has_state) {
                assign(".Random.seed", previous, envir = globalenv())
            } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
                rm(".Random.seed", envir = globalenv())
            }
        },
        add = TRUE
    )

    set.seed(as.integer(seed))
    force(code)
}

#' @noRd
getSegments <- function(counts, chrom, maploc, alpha, undo.SD, seed = 1L) {
    with_local_seed(seed, {
        CNA.object <- DNAcopy::CNA(
            genomdat = counts,
            chrom = chrom,
            maploc = maploc,
            data.type = "logratio"
        ) |> DNAcopy::smooth.CNA()
        DNAcopy::segment(
            CNA.object,
            verbose = 0,
            undo.splits = "sdundo",
            alpha = alpha,
            min.width = 2,
            undo.SD = undo.SD
        )
    })
}

#' @noRd
segment_pair_data <- function(pair_data, don, rec, meta, param, colname, threshold, min_width) {
    if (nrow(pair_data) < 2L) {
        warning(
            sprintf(
                "Skipping donor %s / recipient %s: need at least 2 ideal windows for CBS, got %d",
                don,
                rec,
                nrow(pair_data)
            ),
            call. = FALSE
        )
        return(NULL)
    }
    getSegments(
        counts = pair_data$diff,
        chrom = rep_len("1", nrow(pair_data)),
        maploc = pair_data$iid,
        alpha = param$alpha,
        undo.SD = param$undo.SD,
        seed = if (is.null(param$seed)) 1L else param$seed
    ) |>
        DNAcopy::segments.summary() |>
        dplyr::select(!dplyr::all_of(c("ID", "chrom"))) |>
        dplyr::mutate(
            putative_introgression = .data[[colname]] >= threshold,
            chr_id.don = don,
            chr_id.rec = rec,
            .before = 1L
        ) |>
        dplyr::inner_join(meta, by = c("chr_id.don" = "chr_id")) |>
        dplyr::inner_join(
            meta,
            by = c("chr_id.rec" = "chr_id"),
            suffix = c(".don", ".rec")
        ) |>
        dplyr::filter((.data$loc.end - .data$loc.start + 1) >= min_width)
}

#' @noRd
empty_segments <- function() {
    data.frame(
        chr_id.don = character(),
        chr_id.rec = character(),
        loc.start = integer(),
        loc.end = integer(),
        num.mark = integer(),
        seg.mean = numeric(),
        seg.sd = numeric(),
        seg.median = numeric(),
        seg.mad = numeric(),
        putative_introgression = logical(),
        chr_name.don = character(),
        subgenome.don = character(),
        chr_name.rec = character(),
        subgenome.rec = character(),
        stringsAsFactors = FALSE
    )
}

#' @noRd
annotate_segment_coords <- function(segments, out) {
    merge_data <- dplyr::select(out, dplyr::all_of(
        c(
            "iid",
            "chr_id.don",
            "start.don",
            "stop.don",
            "chr_id.rec",
            "start.rec",
            "stop.rec"
        )
    ))
    segments |>
        dplyr::inner_join(
            y = dplyr::select(merge_data, !dplyr::all_of(c(
                "stop.don", "stop.rec"
            ))),
            by = c("chr_id.don", "chr_id.rec", "loc.start" = "iid")
        ) |>
        dplyr::inner_join(
            y = dplyr::select(merge_data, !dplyr::all_of(c(
                "start.don", "start.rec"
            ))),
            by = c("chr_id.don", "chr_id.rec", "loc.end" = "iid")
        )
}
