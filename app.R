## Shiny entrypoint for the MaMaMIA web service.
##
## ShinyProxy and direct Docker runs start this file through deployment/start.R.
## For local development, run `Rscript -e 'shiny::runApp(".")'` from the
## repository root after `Rscript tools/build-native.R`.

local({
    root <- Sys.getenv("MAMAMIA_ROOT", unset = "")
    if (!nzchar(root)) {
        args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", args, value = TRUE)
        if (length(file_arg) > 0L) {
            root <- dirname(normalizePath(
                sub("^--file=", "", file_arg[[1L]]),
                winslash = "/", mustWork = TRUE
            ))
        }
    }
    if (!nzchar(root)) {
        root <- getwd()
    }
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)

    source(file.path(root, "R", "load.R"), local = FALSE)
    mamamia_load(root)

    max_request <- suppressWarnings(as.numeric(Sys.getenv(
        "MAMAMIA_MAX_REQUEST_SIZE",
        unset = as.character(5 * 1024^3)
    )))
    if (length(max_request) != 1L || !is.finite(max_request) ||
        max_request <= 0) {
        max_request <- 5 * 1024^3
    }
    options(shiny.maxRequestSize = max_request)

    mamamia_app()
})
