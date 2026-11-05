#!/usr/bin/env Rscript

## Compile the htslib/gcrank bridge used by the Shiny service.
##
## Usage:
##   Rscript tools/build-native.R [app-root]
##
## The Rhtslib flags are substituted into src/Makevars, then R CMD SHLIB
## compiles every source into lib/MaMaMIA.so. The image build and the test
## runner call this script; no R package installation is involved.

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

root <- if (length(args) >= 1L && nzchar(args[[1L]])) {
    args[[1L]]
} else {
    normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
}
root <- normalizePath(root, winslash = "/", mustWork = TRUE)

src <- file.path(root, "src")
template <- file.path(src, "Makevars.in")
if (!file.exists(template)) {
    stop("Missing ", template, call. = FALSE)
}
if (!requireNamespace("Rhtslib", quietly = TRUE)) {
    stop("Rhtslib is required to compile the native library", call. = FALSE)
}

cppflags <- paste(
    capture.output(Rhtslib::pkgconfig("PKG_CPPFLAGS")), collapse = ""
)
libs <- paste(
    capture.output(Rhtslib::pkgconfig("PKG_LIBS")), collapse = ""
)
include_dirs <- c(
    system.file("include", package = "Rcpp"),
    system.file("include", package = "Rhtslib")
)
include_dirs <- include_dirs[nzchar(include_dirs)]
if (length(include_dirs) > 0L) {
    cppflags <- paste(
        c(sprintf('-I"%s"', include_dirs), cppflags),
        collapse = " "
    )
}
makevars <- readLines(template, warn = FALSE)
makevars <- gsub("@RHTLIBS_CPPFLAGS@", cppflags, makevars, fixed = TRUE)
makevars <- gsub("@RHTLIBS@", libs, makevars, fixed = TRUE)
writeLines(makevars, file.path(src, "Makevars"))

lib_dir <- file.path(root, "lib")
dir.create(lib_dir, showWarnings = FALSE, recursive = TRUE)
output <- file.path(lib_dir, paste0("MaMaMIA", .Platform$dynlib.ext))
unlink(output)

sources <- c(
    "RcppExports.cpp",
    "htslib-utils.cpp",
    "gcrank-wrapper.cpp",
    list.files(
        file.path(src, "gcrank", "src"),
        pattern = "[.]c$", full.names = FALSE
    )
)
sources <- sources[file.exists(file.path(src, sources))]
if (length(sources) == 0L) {
    stop("No native sources found under ", src, call. = FALSE)
}

old <- setwd(src)
on.exit(setwd(old), add = TRUE)
objects <- sub("[.][^.]+$", ".o", sources)
unlink(file.path(src, objects)[file.exists(file.path(src, objects))])
status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "SHLIB", "-o", output, sources)
)
if (!identical(status, 0L)) {
    stop("R CMD SHLIB failed with status ", status, call. = FALSE)
}
if (!file.exists(output) || file.info(output)$size <= 0) {
    stop("The native build did not produce ", output, call. = FALSE)
}

message("Built ", output)
