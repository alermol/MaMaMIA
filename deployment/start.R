#!/usr/bin/env Rscript

## ShinyProxy and Docker entrypoint. Run through uvr so the pinned R and
## .uvr/library are active:
##
##   uvr run /app/deployment/start.R
##
## The registry is mounted read-only and passed by MAMAMIA_CONFIG (see
## deployment/shinyproxy.yml). The app binds all interfaces because ShinyProxy
## connects to the container from the host network.

root <- Sys.getenv("MAMAMIA_ROOT", unset = "/app")
if (!dir.exists(root)) {
    stop("MAMAMIA_ROOT does not exist: ", root, call. = FALSE)
}
root <- normalizePath(root, winslash = "/", mustWork = TRUE)

config <- Sys.getenv("MAMAMIA_CONFIG", unset = "")
if (!nzchar(config)) {
    stop(
        "MAMAMIA_CONFIG is not set; mount the assembly registry and point ",
        "MAMAMIA_CONFIG at it",
        call. = FALSE
    )
}
if (!file.exists(config)) {
    stop("MAMAMIA_CONFIG does not exist: ", config, call. = FALSE)
}

version_file <- file.path(root, "VERSION")
version <- if (file.exists(version_file)) {
    trimws(readLines(version_file, n = 1L, warn = FALSE))
} else {
    "unknown"
}
message("MaMaMIA ", version)
message("Application root: ", root)
message("Assembly registry: ", config)
message("R: ", R.version.string)
message("Library: ", paste(.libPaths(), collapse = ", "))

## The uvr-managed R bundles cairo and fontconfig but not the system font
## configuration or fonts. Without them the first renderPlot() aborts the whole
## process inside grid's string metrics, so warn while startup can still report
## it. /etc/fonts/fonts.conf is the fontconfig default when no override is set.
fontconfig <- Sys.getenv("FONTCONFIG_FILE", unset = "")
if (!nzchar(fontconfig)) {
    fontconfig <- file.path(
        Sys.getenv("FONTCONFIG_PATH", unset = "/etc/fonts"),
        "fonts.conf"
    )
}
if (capabilities("cairo") && !file.exists(fontconfig)) {
    warning(
        "fontconfig has no configuration (", fontconfig, "); install ",
        "fontconfig and fonts-dejavu-core or the first plot will crash",
        call. = FALSE
    )
}

## Uploads are bounded by MAMAMIA_MAX_REQUEST_SIZE (read by app.R) and by the
## ShinyProxy multipart limits. Raise them together for larger cohorts.
shiny::runApp(
    appDir = root,
    host = "0.0.0.0",
    port = 3838L,
    launch.browser = FALSE
)
