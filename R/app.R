## The BAM-first application object. The Shiny entrypoint is the repository
## root app.R; this module only assembles the UI and server.

#' Build the app object without starting a server
#'
#' The `mamamia_input` attribute is always `"bam"`.
#' @noRd
mamamia_app <- function(config = NULL) {
    app <- bam_launcher_app(config = config)
    attr(app, "mamamia_input") <- "bam"
    app
}
