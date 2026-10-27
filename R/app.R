## The BAM-first launcher: pick a BAM and a species pair, run the pipeline,
## then hand over to the explorer module defined in dashboard.R.

#' Run MaMaMIA as a local Shiny dashboard
#'
#' `runMaMaMIA()` is the main entry point of the package and opens the
#' BAM-first Shiny application. The app reads BAM headers, matches them to the
#' configured species registry, stages the file in session temporary storage,
#' and runs [analyzeBamInput()] before handing the `BamAnalysisResult` to the
#' introgression explorer of [exploreIntrogressions()]. The dashboard accepts
#' one BAM per run; a multi-sample BAM is read through its `@RG SM` tags.
#'
#' The app reads its assembly registry from `config` or, when that is `NULL`,
#' from the `MAMAMIA_CONFIG` environment variable (the ShinyProxy deployment
#' path). A missing or invalid registry is rendered as a configuration error
#' card rather than a crash.
#'
#' The uploaded BAM is staged into a private per-session directory (symlinked
#' for large files, copied for small files and when symlinks are unsupported)
#' and indexed there. The source upload directory is never written to; the
#' session directory, the staged BAM and its CSI index are removed when the
#' session ends. GC index handles are process-wide and are not cleared at
#' session end.
#'
#' The app also offers a bundled synthetic example. It loads a pre-computed
#' result instead of reading a BAM, but replays the normal pipeline stages, so
#' the user sees the same progress messages as a real run without mounting an
#' example BAM in the container.
#'
#' The analysis runs in the R session, so the session blocks while it runs; the
#' progress bar reports the pipeline stages (headers, indexes, windows, GC,
#' breadth/counts, scoring, masking, segmentation, calibration). Errors are
#' reported with the failing stage and a suggested fix.
#'
#' The explorer shows the score profile, the call table with the BEDPE export,
#' and, for a cohort, a per-sample selector, a cohort summary and a combined
#' BEDPE download carrying the `sample` column.
#'
#' Dashboard packages are optional MaMaMIA dependencies. Install them with
#' `install.packages(c("shiny", "bslib", "shinyWidgets", "DT", "htmlwidgets"))`
#' before calling
#' this function.
#'
#' @param config Optional assembly registry: a path to a YAML configuration or
#'   a validated configuration list from [readMaMaMIAConfig()]. The app falls
#'   back to the `MAMAMIA_CONFIG` environment variable when this is `NULL`.
#' @param max_request_size Maximum size in bytes of a single Shiny upload,
#'   passed to `options(shiny.maxRequestSize = )`. Defaults to 5 GiB.
#' @param host Host passed to [shiny::runApp()].
#' @param port Port passed to [shiny::runApp()]. Use `NULL` to let Shiny choose.
#' @param launch.browser Logical passed to [shiny::runApp()].
#'
#' @return The return value of [shiny::runApp()]. Called for its side effect of
#'   launching the dashboard.
#' @seealso [exploreIntrogressions()], [analyzeBamInput()], [readMaMaMIAConfig()]
#' @examples
#' \dontrun{
#' runMaMaMIA(config = "assemblies.yml")
#' }
#' @export
runMaMaMIA <- function(config = NULL,
                       max_request_size = 5 * 1024^3,
                       host = "127.0.0.1", port = NULL,
                       launch.browser = interactive()) {
    require_dashboard_packages()
    if (!is.numeric(max_request_size) || length(max_request_size) != 1L ||
        is.na(max_request_size) || max_request_size <= 0) {
        stop(
            "`max_request_size` must be a positive number of bytes",
            call. = FALSE
        )
    }
    options(shiny.maxRequestSize = max_request_size)

    app <- mamamia_app(config = config)
    args <- list(appDir = app, host = host, launch.browser = launch.browser)
    if (!is.null(port)) {
        args$port <- port
    }
    do.call(shiny::runApp, args)
}

#' Build the app object without starting a server
#'
#' Kept internal so tests can inspect the BAM-first app object without
#' launching Shiny. The `mamamia_input` attribute is always `"bam"`.
#'
#' @noRd
mamamia_app <- function(config = NULL) {
    app <- bam_launcher_app(config = config)
    attr(app, "mamamia_input") <- "bam"
    app
}
