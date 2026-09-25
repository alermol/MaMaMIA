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
as_window_df <- function(x, value_col, arg_name) {
    if (inherits(x, "matrix")) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    if (!is.data.frame(x)) {
        stop(sprintf("`%s` must be a data.frame, tibble, or matrix", arg_name),
            call. = FALSE
        )
    }
    if (ncol(x) != 4L) {
        stop(
            sprintf(
                "`%s` must have exactly 4 columns (chr_id, start, stop, %s)",
                arg_name,
                value_col
            ),
            call. = FALSE
        )
    }

    x <- as.data.frame(x, stringsAsFactors = FALSE)
    names(x) <- c("chr_id", "start", "stop", value_col)

    if (is.factor(x$chr_id)) {
        x$chr_id <- as.character(x$chr_id)
    }
    if (!is.atomic(x$chr_id) || is.list(x$chr_id)) {
        stop(
            sprintf(
                "`%s$chr_id` must be atomic (character/integer/factor)",
                arg_name
            ),
            call. = FALSE
        )
    }
    x$chr_id <- as.character(x$chr_id)

    for (col in c("start", "stop", value_col)) {
        if (!is.numeric(x[[col]])) {
            char_col <- as.character(x[[col]])
            non_na <- !is.na(x[[col]])
            looks_numeric <- !non_na | grepl(
                "^\\s*[-+]?(([0-9]+\\.?[0-9]*)|(\\.[0-9]+))([eE][-+]?[0-9]+)?\\s*$",
                char_col
            )
            if (!all(looks_numeric)) {
                stop(sprintf("`%s$%s` must be numeric", arg_name, col),
                    call. = FALSE
                )
            }
            x[[col]] <- as.numeric(char_col)
        }
    }

    if (anyNA(x$start) || anyNA(x$stop)) {
        stop(sprintf("`%s` start/stop must not contain NA", arg_name),
            call. = FALSE
        )
    }
    if (anyNA(x[[value_col]])) {
        stop(sprintf("`%s$%s` must not contain NA", arg_name, value_col),
            call. = FALSE
        )
    }
    bad <- which(!(x$start < x$stop))
    if (length(bad) > 0L) {
        stop(
            sprintf(
                "`%s` has intervals with start >= stop (first at row %d)",
                arg_name,
                bad[[1L]]
            ),
            call. = FALSE
        )
    }

    return(x)
}

#' @noRd
validate_windows_tiled <- function(x, arg_name) {
    for (chr in unique(x$chr_id)) {
        idx <- which(x$chr_id == chr)
        if (length(idx) < 2L) {
            next
        }
        win_start <- x$start[idx]
        win_stop <- x$stop[idx]
        next_start <- win_start[-1L]
        prev_stop <- win_stop[-length(win_stop)]
        if (any(next_start < prev_stop)) {
            stop(
                sprintf(
                    "`%s` has overlapping windows on chromosome %s",
                    arg_name,
                    chr
                ),
                call. = FALSE
            )
        }
        if (any(next_start > prev_stop)) {
            stop(
                sprintf(
                    "`%s` has gaps between windows on chromosome %s (expected contiguous half-open intervals)",
                    arg_name,
                    chr
                ),
                call. = FALSE
            )
        }
    }
    invisible(x)
}


#' @noRd
as_meta_df <- function(x, arg_name = "meta") {
    if (inherits(x, "matrix")) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    if (!is.data.frame(x)) {
        stop(sprintf("`%s` must be a data.frame, tibble, or matrix", arg_name),
            call. = FALSE
        )
    }
    if (ncol(x) != 3L) {
        stop(
            sprintf(
                "`%s` must have exactly 3 columns (chr_id, chr_name, subgenome)",
                arg_name
            ),
            call. = FALSE
        )
    }

    x <- as.data.frame(x, stringsAsFactors = FALSE)
    names(x) <- c("chr_id", "chr_name", "subgenome")

    for (col in c("chr_id", "chr_name", "subgenome")) {
        if (is.factor(x[[col]])) {
            x[[col]] <- as.character(x[[col]])
        }
        if (!is.atomic(x[[col]]) || is.list(x[[col]])) {
            stop(
                sprintf(
                    "`%s$%s` must be atomic (character/integer/factor)",
                    arg_name,
                    col
                ),
                call. = FALSE
            )
        }
        x[[col]] <- as.character(x[[col]])
    }

    if (anyNA(x)) {
        stop(sprintf("`%s` must not contain NA", arg_name), call. = FALSE)
    }
    if (anyDuplicated(x$chr_id)) {
        stop(sprintf("`%s$chr_id` must be unique", arg_name), call. = FALSE)
    }

    return(x)
}

#' @noRd
resolve_chr_targets <- function(meta,
                                don_chromlist,
                                rec_chromlist,
                                don_subgenomes,
                                rec_subgenomes,
                                don_chr_ids,
                                rec_chr_ids) {
    stopifnot(
        "Require at least one of don_chromlist, rec_chromlist, don_subgenomes, rec_subgenomes" = !all(vapply(
            list(
                don_chromlist,
                rec_chromlist,
                don_subgenomes,
                rec_subgenomes
            ),
            is.null,
            logical(1)
        ))
    )

    if (!is.null(don_chromlist) && length(don_chromlist) == 0L) {
        stop("`don_chromlist` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(rec_chromlist) && length(rec_chromlist) == 0L) {
        stop("`rec_chromlist` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(don_subgenomes) && length(don_subgenomes) == 0L) {
        stop("`don_subgenomes` must be NULL or non-empty", call. = FALSE)
    }
    if (!is.null(rec_subgenomes) && length(rec_subgenomes) == 0L) {
        stop("`rec_subgenomes` must be NULL or non-empty", call. = FALSE)
    }

    if (!is.null(don_subgenomes)) {
        stopifnot(
            "Unknown donor subgenome label(s)" =
                all(don_subgenomes %in% meta$subgenome)
        )
    }
    if (!is.null(rec_subgenomes)) {
        stopifnot(
            "Unknown recipient subgenome label(s)" =
                all(rec_subgenomes %in% meta$subgenome)
        )
    }
    if (!is.null(don_chromlist)) {
        stopifnot(
            "Unknown donor chromosome ID(s)" =
                all(don_chromlist %in% don_chr_ids)
        )
    }
    if (!is.null(rec_chromlist)) {
        stopifnot(
            "Unknown recipient chromosome ID(s)" =
                all(rec_chromlist %in% rec_chr_ids)
        )
    }

    resolve_one <- function(chromlist, subgenomes, current_ids) {
        current_ids <- as.character(current_ids)
        if (is.null(chromlist) && is.null(subgenomes)) {
            return(current_ids)
        }
        ids <- character(0)
        if (!is.null(subgenomes)) {
            ids <- c(ids, intersect(as.character(meta$chr_id[meta$subgenome %in% subgenomes]), current_ids))
        }
        if (!is.null(chromlist)) {
            ids <- c(ids, as.character(chromlist))
        }
        unique(ids)
    }

    target_don <- resolve_one(don_chromlist, don_subgenomes, don_chr_ids)
    target_rec <- resolve_one(rec_chromlist, rec_subgenomes, rec_chr_ids)
    if (length(target_don) == 0L || length(target_rec) == 0L) {
        stop("Subset would leave donor or recipient with no chromosomes",
            call. = FALSE
        )
    }

    list(target_don, target_rec)
}


#' @noRd
validate_equal_n_windows <- function(chr_id) {
    counts <- table(chr_id, useNA = "no")
    if (length(counts) == 0L) {
        stop("Input must contain at least one chromosome", call. = FALSE)
    }
    n_unique <- length(unique(as.integer(counts)))
    if (n_unique != 1L) {
        detail <- paste(sprintf("%s=%d", names(counts), as.integer(counts)), collapse = ", ")
        stop("All chromosomes must have the same number of windows; found: ",
            detail,
            call. = FALSE
        )
    }
    invisible(as.integer(counts[[1L]]))
}

#' @noRd
check_cores <- function(x) {
    x <- as.integer(x)[1L]
    max_cores <- parallel::detectCores()
    if (is.na(max_cores) || max_cores < 1L) {
        max_cores <- 1L
    }
    if (is.na(x) || x < 1L || x > max_cores) {
        warning("Cores must be between 1 and ",
            max_cores,
            ". Fallback to 1 core.",
            call. = FALSE
        )
        return(1L)
    }
    return(x)
}

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

#' @noRd
hash_object <- function(x) {
    path <- tempfile()
    on.exit(unlink(path), add = TRUE)
    saveRDS(x, path, version = 2)
    unname(tools::md5sum(path))
}

#' Coverage-model call shared by fitting and cache-keying.
#' @noRd
coverage_model_call <- function(cores) {
    substitute(
        glmmTMB::glmmTMB(
            cov ~ s(gc, k = 10) + subgenome,
            ziformula = ~ s(gc, k = 10) + subgenome,
            family = glmmTMB::nbinom2(),
            data = ideal_data,
            REML = TRUE,
            control = glmmTMB::glmmTMBControl(parallel = list(n = Cores))
        ),
        list(Cores = cores)
    )
}

#' Cache key for a fitted coverage model.
#' @noRd
fit_cache_key <- function(ideal_data, cores, model_call) {
    spec <- as.list(model_call)
    spec$data <- NULL
    hash_object(list(
        model = paste(deparse(as.call(spec)), collapse = " "),
        chr_id = as.character(ideal_data$chr_id),
        iid = as.integer(ideal_data$iid),
        cov = as.numeric(ideal_data$cov),
        gc = as.numeric(ideal_data$gc),
        subgenome = as.character(ideal_data$subgenome),
        cores = as.integer(cores),
        package_version = as.character(utils::packageVersion("MaMaMIA")),
        glmmTMB_version = as.character(utils::packageVersion("glmmTMB")),
        R_version = as.character(getRversion())
    ))
}

#' Cached-fit entry layout version.
#' @noRd
CACHE_LAYOUT_VERSION <- 1L

#' @noRd
read_cached_fit <- function(cache, key) {
    if (is.null(cache) || !file.exists(cache)) {
        return(NULL)
    }
    cached <- tryCatch(readRDS(cache), error = function(e) NULL)
    if (!is.list(cached) ||
        !identical(cached$format, CACHE_LAYOUT_VERSION) ||
        !identical(cached$key, key)) {
        return(NULL)
    }
    cached$fit
}

#' @noRd
write_cached_fit <- function(cache, key, fit) {
    if (is.null(cache)) {
        return(invisible(NULL))
    }
    dir <- dirname(cache)
    if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
    tmp <- paste0(cache, ".tmp")
    written <- tryCatch(
        {
            saveRDS(
                list(format = CACHE_LAYOUT_VERSION, key = key, fit = fit),
                tmp,
                version = 2
            )
            TRUE
        },
        error = function(e) FALSE
    )
    if (written) {
        file.rename(tmp, cache)
    }
    invisible(NULL)
}

#' Fixed/random design columns of a stored mgcv smooth for new values.
#' @noRd
smooth_design <- function(sm, gc) {
    ev <- eigen(sm$S[[1L]], symmetric = TRUE)
    if (ev$vectors[1L, 1L] < 0) {
        ev$vectors <- -ev$vectors
    }
    p_rank <- min(sm$rank, ncol(sm$X))
    null_rank <- sm$df - sm$rank
    scaling <- 1 / sqrt(c(ev$values[seq_len(p_rank)], rep(1, null_rank)))
    basis <- mgcv::PredictMat(sm, data = data.frame(gc = gc)) %*%
        t(t(ev$vectors) * scaling)
    fixed <- if (p_rank < sm$df) {
        basis[, (p_rank + 1L):sm$df, drop = FALSE]
    } else {
        matrix(0, nrow(basis), 0L)
    }
    list(
        fixed = fixed,
        random = basis[, seq_len(p_rank), drop = FALSE],
        label = sm$label
    )
}

#' Fixed/random design matrices for the fitted coverage model.
#' @noRd
glmmtmb_design <- function(fit, newdata) {
    subgenomes <- levels(factor(fit$frame$subgenome))
    if (anyNA(factor(newdata$subgenome, levels = subgenomes))) {
        stop("unknown subgenome level", call. = FALSE)
    }
    dummies <- unname(stats::model.matrix(
        ~ factor(newdata$subgenome, levels = subgenomes)
    )[, -1L, drop = FALSE])

    build <- function(part, coef_names) {
        info <- fit$modelInfo$reTrms[[part]]$smooth_info
        if (length(info) != 1L) {
            stop("expected exactly one smooth in `", part, "`", call. = FALSE)
        }
        sd <- smooth_design(info[[1L]]$sm, newdata$gc)
        expected <- c(
            "(Intercept)",
            paste0("subgenome", subgenomes[-1L]),
            paste0(sd$label, seq_len(ncol(sd$fixed)))
        )
        if (!identical(expected, as.character(coef_names))) {
            stop("unexpected fixed-effect structure in `", part, "`", call. = FALSE)
        }
        list(X = cbind(1, dummies, sd$fixed), Z = sd$random)
    }

    list(
        cond = build("cond", names(glmmTMB::fixef(fit)$cond)),
        zi = build("zi", names(glmmTMB::fixef(fit)$zi))
    )
}

#' Design matrices from glmmTMB's prediction machinery.
#' @noRd
design_from_tmb <- function(fit, newdata) {
    tmb <- stats::predict(fit, newdata = newdata, debug = TRUE)$data.tmb
    if (is.null(dim(tmb$X))) {
        stop("Fitted model exposes no dense conditional design matrix", call. = FALSE)
    }
    n_aug <- nrow(tmb$X)
    n_new <- nrow(newdata)
    if (n_aug != nrow(fit$frame) + n_new) {
        stop("Unexpected augmented design size; cannot align predictions", call. = FALSE)
    }
    idx <- (n_aug - n_new + 1L):n_aug
    trim <- function(M) {
        as.matrix(M)[idx, , drop = FALSE]
    }
    list(
        cond = list(X = trim(tmb$X), Z = trim(tmb$Z)),
        zi = list(X = trim(tmb$Xzi), Z = trim(tmb$Zzi))
    )
}

#' Mean response of the fitted ZINB coverage model for new data.
#' @noRd
eta_zinb_response <- function(fit, newdata) {
    pars <- fit$fit$parfull
    get_par <- function(name) {
        out <- pars[names(pars) == name]
        if (length(out) == 0L) {
            stop("Fitted model has no `", name, "` parameter", call. = FALSE)
        }
        out
    }

    design <- tryCatch(
        glmmtmb_design(fit, newdata),
        error = function(e) NULL
    )
    if (is.null(design)) {
        design <- design_from_tmb(fit, newdata)
    }

    eta_cond <- as.numeric(design$cond$X %*% get_par("beta")) +
        as.numeric(design$cond$Z %*% get_par("b"))
    eta_zi <- as.numeric(design$zi$X %*% get_par("betazi")) +
        as.numeric(design$zi$Z %*% get_par("bzi"))

    (1 - stats::plogis(eta_zi)) * exp(eta_cond)
}

#' Expected ZINB coverage at observed GC and at reference GC.
#' @noRd
predict_zinb_response <- function(fit, gc, subgenome, gc_ref) {
    n <- length(gc)

    fallback <- function() {
        nd <- data.frame(gc = gc, subgenome = subgenome)
        list(
            actual = as.numeric(stats::predict(fit, newdata = nd, type = "response")),
            ref = as.numeric(stats::predict(fit,
                newdata = transform(nd, gc = gc_ref),
                type = "response"
            ))
        )
    }

    ok <- !is.na(gc) & !is.na(subgenome)
    if (!any(ok)) {
        na <- rep(NA_real_, n)
        return(list(actual = na, ref = na))
    }

    key <- function(g, s) paste(sprintf("%.17g", g), s, sep = "\r")
    rows <- data.frame(gc = gc[ok], subgenome = as.character(subgenome[ok]))
    row_keys <- key(rows$gc, rows$subgenome)
    unique_rows <- rows[!duplicated(row_keys), , drop = FALSE]
    unique_keys <- key(unique_rows$gc, unique_rows$subgenome)
    subgenomes <- unique(unique_rows$subgenome)
    ref_rows <- data.frame(gc = gc_ref, subgenome = subgenomes)

    pred <- tryCatch(
        eta_zinb_response(fit, rbind(unique_rows, ref_rows)),
        error = function(e) NULL
    )
    if (is.null(pred) || length(pred) != nrow(unique_rows) + nrow(ref_rows)) {
        return(fallback())
    }

    n_unique <- nrow(unique_rows)
    actual <- rep(NA_real_, n)
    actual[ok] <- pred[match(row_keys, unique_keys)]
    ref <- rep(NA_real_, n)
    ref[ok] <- pred[n_unique + match(rows$subgenome, subgenomes)]
    list(actual = actual, ref = ref)
}
