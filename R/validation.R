## Validation of the window tables, chromosome metadata and numeric arguments
## shared by the constructors (`RCA()`, `readMaMaMIAInput()`).

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
