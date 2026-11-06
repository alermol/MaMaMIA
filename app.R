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

    options(shiny.maxRequestSize = mamamia_max_request_size())

    mamamia_app()
})
