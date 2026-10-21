#!/usr/bin/env Rscript

## ShinyProxy entrypoint. Run through uvr so the pinned R and .uvr/library are
## active:
##
##   uvr run /app/deployment/start.R
##
## The registry is mounted read-only and passed by MAMAMIA_CONFIG (see
## deployment/shinyproxy.yml). Shiny listens on all interfaces because
## ShinyProxy connects to the container from the host network.

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

message("MaMaMIA ", utils::packageVersion("MaMaMIA"))
message("Assembly registry: ", config)
message("R: ", R.version.string)
message("Library: ", paste(.libPaths(), collapse = ", "))

## 5 GiB matches the default server.max_total_bytes; raise both together when
## the registry allows larger cohorts. Uploads are also bounded by ShinyProxy
## and the reverse proxy (deployment/nginx.conf).
MaMaMIA::runMaMaMIA(
    host = "0.0.0.0",
    port = 3838L,
    launch.browser = FALSE,
    config = config,
    max_request_size = 5 * 1024^3
)
