## Per-session cache for generated result bundles.
##
## Bundle generation is expensive and can take minutes for a cohort. A session
## keeps one archive for the current request so repeated download-button
## presses reuse it. A hit serves the archive from its first generation,
## including that build's generated timestamp. The cache is deliberately
## small: one entry, cleared when the Shiny session ends.

#' New empty bundle cache
#' @noRd
bundle_cache_state <- function() {
    cache <- new.env(parent = emptyenv())
    cache$key <- NULL
    cache$path <- NULL
    cache
}

#' Cached archive for this key, or `NULL`
#' @noRd
bundle_cache_get <- function(cache, key) {
    if (!is.environment(cache) || is.null(key) || is.null(cache$key) ||
        !identical(cache$key, key) || is.null(cache$path) ||
        !file.exists(cache$path)) {
        return(NULL)
    }
    cache$path
}

#' Store one archive and remove the previous one
#' @noRd
bundle_cache_store <- function(cache, key, path) {
    if (!is.environment(cache)) {
        stop("`cache` must come from bundle_cache_state()", call. = FALSE)
    }
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !file.exists(path) || is.na(file.info(path)$size) ||
        file.info(path)$size <= 0) {
        stop("Cannot cache an empty result bundle", call. = FALSE)
    }
    previous <- cache$path
    cache$key <- key
    cache$path <- path
    if (!is.null(previous) && !identical(previous, path) &&
        file.exists(previous)) {
        unlink(previous)
    }
    invisible(path)
}

#' Delete the cached archive
#' @noRd
bundle_cache_clear <- function(cache) {
    if (is.environment(cache) && !is.null(cache$path) &&
        file.exists(cache$path)) {
        unlink(cache$path)
    }
    if (is.environment(cache)) {
        cache$key <- NULL
        cache$path <- NULL
    }
    invisible(NULL)
}

#' Serve a bundle from cache or build it once
#'
#' @param cache Cache environment from `bundle_cache_state()`.
#' @param key Structured key from `bundle_request_key()`.
#' @param file Download destination.
#' @param build Function taking one path and writing the archive there.
#'
#' @return `file`, invisibly.
#' @noRd
bundle_cache_serve <- function(cache, key, file, build) {
    if (!is.function(build)) {
        stop("`build` must be a function", call. = FALSE)
    }
    path <- bundle_cache_get(cache, key)
    if (is.null(path)) {
        temporary <- tempfile("mamamia-bundle-cache-", fileext = ".zip")
        cached <- FALSE
        on.exit({
            if (!cached && file.exists(temporary)) {
                unlink(temporary)
            }
        }, add = TRUE)
        build(temporary)
        if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
            file.info(temporary)$size <= 0) {
            stop("Bundle generation produced no archive", call. = FALSE)
        }
        bundle_cache_store(cache, key, temporary)
        cached <- TRUE
        path <- temporary
    }
    if (!file.copy(path, file, overwrite = TRUE)) {
        stop("Could not copy the cached result bundle", call. = FALSE)
    }
    invisible(file)
}
