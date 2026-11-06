## Source-based runtime for the MaMaMIA web service.
##
## `mamamia_load()` sources the application modules into the global
## environment and loads the compiled native library. The Shiny entrypoint
## (app.R) and the test runner call it once per process.

#' Root of the deployed application
#'
#' The root is fixed by `MAMAMIA_ROOT` in the container and otherwise by the
#' current working directory. Asset lookups are relative to this root.
#' @noRd
mamamia_root <- function() {
    root <- getOption("MaMaMIA.root", NULL)
    if (is.null(root) || length(root) != 1L || is.na(root) || !nzchar(root)) {
        root <- Sys.getenv("MAMAMIA_ROOT", unset = getwd())
    }
    normalizePath(root, winslash = "/", mustWork = FALSE)
}

#' Path inside the application tree
#' @noRd
mamamia_path <- function(...) {
    file.path(mamamia_root(), ...)
}

#' Application version
#' @noRd
mamamia_version <- function() {
    path <- mamamia_path("VERSION")
    if (!file.exists(path)) {
        return("0.1.0")
    }
    value <- trimws(readLines(path, n = 1L, warn = FALSE))
    if (length(value) != 1L || is.na(value) || !nzchar(value)) {
        return("0.1.0")
    }
    value
}

#' Deployment mode for upload limits
#'
#' `MAMAMIA_DEPLOYMENT_MODE=shinyproxy` enables the upload caps from the
#' registry and `MAMAMIA_MAX_REQUEST_SIZE`. Any other value, including the
#' default, is a direct container run without a size cap.
#' @noRd
mamamia_deployment_mode <- function() {
    mode <- tolower(trimws(Sys.getenv(
        "MAMAMIA_DEPLOYMENT_MODE", unset = "direct"
    )))
    if (identical(mode, "shinyproxy")) "shinyproxy" else "direct"
}

#' Shiny request-body limit for the current deployment
#' @noRd
mamamia_max_request_size <- function() {
    if (!identical(mamamia_deployment_mode(), "shinyproxy")) {
        return(Inf)
    }
    value <- suppressWarnings(as.numeric(Sys.getenv(
        "MAMAMIA_MAX_REQUEST_SIZE",
        unset = as.character(5 * 1024^3)
    )))
    if (length(value) != 1L || !is.finite(value) || value <= 0) {
        value <- 5 * 1024^3
    }
    value
}

#' Path to the compiled native library
#' @noRd
mamamia_native_lib <- function(root = mamamia_root()) {
    name <- paste0("MaMaMIA", .Platform$dynlib.ext)
    candidates <- c(
        file.path(root, "lib", name),
        file.path(root, "src", name)
    )
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates) == 0L) {
        return(file.path(root, "lib", name))
    }
    normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE)
}

#' Load the application code and native library
#'
#' Sources every module under `R/` except this loader, exposes the rlang
#' `.data` pronoun used by ggplot2 aesthetics, and loads the compiled htslib
#' and gcrank bridge.
#'
#' @param root Application root.
#' @param native Whether to require and load the compiled native library.
#' @return The normalized application root, invisibly.
#' @noRd
mamamia_load <- function(root = getwd(), native = TRUE) {
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    options(MaMaMIA.root = root)

    ## Rcpp's C++ API resolves its callables from the loaded Rcpp namespace.
    loadNamespace("Rcpp")
    assign(".data", rlang::.data, envir = globalenv())

    files <- list.files(
        file.path(root, "R"), pattern = "[.]R$", full.names = TRUE
    )
    files <- files[basename(files) != "load.R"]
    first <- file.path(root, "R", "utils.R")
    if (file.exists(first)) {
        files <- c(first, setdiff(files, first))
    }
    for (file in files) {
        sys.source(file, envir = globalenv(), keep.source = FALSE)
    }

    if (isTRUE(native)) {
        library <- mamamia_native_lib(root)
        if (!file.exists(library)) {
            stop(
                "The native library is missing. Run ",
                "`Rscript tools/build-native.R` first.",
                call. = FALSE
            )
        }
        loaded <- vapply(
            getLoadedDLLs(),
            function(info) normalizePath(info[["path"]], mustWork = FALSE),
            character(1)
        )
        if (!library %in% loaded) {
            dyn.load(library)
        }
    }

    invisible(root)
}
