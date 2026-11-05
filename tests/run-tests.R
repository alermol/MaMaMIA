#!/usr/bin/env Rscript

## Run the internal test suite without installing an R package.
##
## Usage:
##   Rscript tests/run-tests.R [filter]
##
## The native library is rebuilt when it is missing. MAMAMIA_ROOT selects the
## application root; by default it is derived from this script's location.

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0L) {
    dirname(normalizePath(
        sub("^--file=", "", script_arg[[1L]]),
        winslash = "/", mustWork = TRUE
    ))
} else {
    getwd()
}
root <- Sys.getenv(
    "MAMAMIA_ROOT",
    unset = normalizePath(
        file.path(script_dir, ".."), winslash = "/", mustWork = TRUE
    )
)
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
if (!nzchar(Sys.getenv("MAMAMIA_ROOT", unset = ""))) {
    Sys.setenv(MAMAMIA_ROOT = root)
}

source(file.path(root, "R", "load.R"), local = FALSE)
if (!file.exists(mamamia_native_lib(root))) {
    message("Native library missing; building it first")
    status <- system2(
        file.path(R.home("bin"), "Rscript"),
        c(file.path(root, "tools", "build-native.R"), root)
    )
    if (!identical(status, 0L)) {
        stop("The native build failed", call. = FALSE)
    }
}
mamamia_load(root)

testthat::test_dir(
    file.path(root, "tests", "testthat"),
    filter = if (length(args) > 0L) args[[1L]] else NULL,
    reporter = "summary",
    stop_on_failure = TRUE
)
