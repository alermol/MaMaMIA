## The GC-bias fit is by far the most expensive step, so it can be memoised to
## an `.rds` file. The cache key covers the fitted windows, the model
## specification, `cores` and the R, package and glmmTMB versions; any mismatch
## silently falls back to refitting.

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
