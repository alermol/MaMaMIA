## The launcher app: pick input, run RCA -> GC correction -> segmentation, then
## hand over to the explorer module defined in dashboard.R.

#' @noRd
default_cores <- function() {
    cores <- parallel::detectCores()
    if (is.na(cores) || cores < 1L) {
        return(1L)
    }
    as.integer(max(1L, min(4L, cores)))
}

#' @noRd
default_fit_cache <- function() {
    dir <- tryCatch(tools::R_user_dir("MaMaMIA", "cache"), error = function(e) NULL)
    if (is.null(dir)) {
        return(NULL)
    }
    file.path(dir, "zinb-fit.rds")
}

#' @noRd
subgenome_labels <- function(meta, chr_ids) {
    sort(unique(meta$subgenome[meta$chr_id %in% chr_ids]))
}

#' @noRd
default_subgenome <- function(meta, chr_ids) {
    subgenomes <- meta$subgenome[meta$chr_id %in% chr_ids]
    if (length(subgenomes) == 0L) {
        return(NULL)
    }
    counts <- table(subgenomes)
    names(counts)[[which.max(counts)]]
}

#' @noRd
example_input_data <- function(envir = parent.frame()) {
    utils::data("triticum", package = "MaMaMIA", envir = envir)
    dcov <- get("don_cov", envir = envir)
    rcov <- get("rec_cov", envir = envir)
    meta <- get("meta", envir = envir)
    names(meta) <- c("chr_id", "chr_name", "subgenome")
    list(
        dcov = dcov,
        dgc = get("don_gc", envir = envir),
        rcov = rcov,
        rgc = get("rec_gc", envir = envir),
        meta = meta,
        don_chr_ids = unique(dcov[[1L]]),
        rec_chr_ids = unique(rcov[[1L]]),
        sample = "triticum (packaged example)",
        reference = "Triticum timopheevii",
        window = "1000"
    )
}

#' @noRd
analyze_input_data <- function(input, cores = 1L, don_subgenomes = NULL,
                               rec_subgenomes = NULL, cache = NULL, verbose = TRUE) {
    stopifnot(is.list(input))
    missing <- setdiff(c("dcov", "dgc", "rcov", "rgc", "meta"), names(input))
    if (length(missing) > 0L) {
        stop("Input data is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }

    rca <- RCA(input$dcov, input$dgc, input$rcov, input$rgc, input$meta)
    corrected <- correctReadCounts(rca, cores = cores, verbose = verbose, cache = cache)
    if (!is.null(don_subgenomes) || !is.null(rec_subgenomes)) {
        corrected <- subset(corrected,
            don_subgenomes = don_subgenomes,
            rec_subgenomes = rec_subgenomes
        )
    }
    pairs <- defaultTargetPairs(
        corrected$meta$meta,
        corrected$meta$don_chr_ids,
        corrected$meta$rec_chr_ids
    )
    if (verbose) {
        message(sprintf("Segmenting %d donor-recipient pair(s)...", length(pairs)))
    }
    segments(corrected, target_pairs = pairs, verbose = verbose)
}

#' @noRd
keep_selection <- function(selected, choices) {
    if (length(choices) == 0L) {
        return(character())
    }
    if (!is.null(selected) && length(selected) == 1L && !is.na(selected) && selected %in% choices) {
        return(selected)
    }
    choices[[1L]]
}

#' @noRd
scan_value <- function(expr) {
    tryCatch(
        list(ok = TRUE, value = expr),
        error = function(e) list(ok = FALSE, message = conditionMessage(e))
    )
}

#' @noRd
scan_windows <- function(analysis_dir, assemblies_dir, sample, reference) {
    if (length(analysis_dir) != 1L || is.na(analysis_dir) || !nzchar(analysis_dir) ||
        length(assemblies_dir) != 1L || is.na(assemblies_dir) || !nzchar(assemblies_dir) ||
        length(sample) != 1L || is.na(sample) || !nzchar(sample) ||
        length(reference) != 1L || is.na(reference) || !nzchar(reference)) {
        return(integer())
    }
    tryCatch(
        list_window_labels(analysis_dir, assemblies_dir, sample, reference),
        error = function(e) integer()
    )
}

#' @noRd
launcher_app <- function(analysis_dir = NULL, assemblies_dir = NULL, cache = NULL) {
    analysis_dir <- if (is.null(analysis_dir)) "" else as.character(analysis_dir)[[1L]]
    assemblies_dir <- if (is.null(assemblies_dir)) "" else as.character(assemblies_dir)[[1L]]

    initial_samples <- tryCatch(
        listMaMaMIASamples(analysis_dir),
        error = function(e) character()
    )
    initial_references <- tryCatch(
        listMaMaMIAReferences(assemblies_dir),
        error = function(e) character()
    )
    initial_windows <- scan_windows(
        analysis_dir, assemblies_dir,
        initial_samples[1L], initial_references[1L]
    )
    max_cores <- parallel::detectCores()
    if (is.na(max_cores) || max_cores < 1L) {
        max_cores <- 1L
    }

    ui <- bslib::page_sidebar(
        title = "MaMaMIA",
        sidebar = bslib::sidebar(
            width = 400,
            shiny::radioButtons(
                "source", "Input source",
                choices = c(
                    "Directory layout" = "files",
                    "Packaged example (triticum)" = "example"
                ),
                selected = if (length(initial_samples) > 0L &&
                    length(initial_references) > 0L) "files" else "example"
            ),
            shiny::conditionalPanel(
                condition = "input.source == 'files'",
                shiny::textInput("analysis_dir", "Analysis directory", value = analysis_dir),
                shiny::textInput("assemblies_dir", "Assemblies directory", value = assemblies_dir),
                shiny::actionButton("scan", "Scan directories", class = "btn-sm"),
                shiny::selectInput("sample", "Sample", choices = initial_samples),
                shiny::selectInput("reference", "Reference", choices = initial_references),
                shiny::selectInput("window", "Windowing scheme", choices = initial_windows)
            ),
            shiny::numericInput(
                "cores", "Cores for the GC-bias fit",
                value = default_cores(), min = 1L, max = max_cores
            ),
            shiny::selectInput("don_subgenome", "Donor subgenome", choices = NULL),
            shiny::selectInput("rec_subgenome", "Recipient subgenome", choices = NULL),
            shiny::actionButton("run", "Run analysis", class = "btn-primary"),
            shiny::hr(),
            shiny::h6("Explorer filters"),
            explorer_controls_ui("explorer")
        ),
        bslib::card(
            bslib::card_header("Analysis"),
            shiny::textOutput("input_summary")
        ),
        explorer_results_ui("explorer")
    )

    server <- function(input, output, session) {
        input_data <- shiny::reactiveVal(NULL)
        isa <- shiny::reactiveVal(NULL)

        update_subgenome_choices <- function(data) {
            if (is.null(data)) {
                shiny::updateSelectInput(session, "don_subgenome", choices = character())
                shiny::updateSelectInput(session, "rec_subgenome", choices = character())
                return()
            }
            shiny::updateSelectInput(
                session, "don_subgenome",
                choices = subgenome_labels(data$meta, data$don_chr_ids),
                selected = default_subgenome(data$meta, data$don_chr_ids)
            )
            shiny::updateSelectInput(
                session, "rec_subgenome",
                choices = subgenome_labels(data$meta, data$rec_chr_ids),
                selected = default_subgenome(data$meta, data$rec_chr_ids)
            )
        }

        update_window_choices <- function() {
            choices <- scan_windows(
                input$analysis_dir, input$assemblies_dir,
                input$sample, input$reference
            )
            shiny::updateSelectInput(
                session, "window",
                choices = choices,
                selected = keep_selection(input$window, choices)
            )
        }

        shiny::observeEvent(input$scan, {
            samples <- scan_value(listMaMaMIASamples(input$analysis_dir))
            references <- scan_value(listMaMaMIAReferences(input$assemblies_dir))
            for (result in list(samples, references)) {
                if (!result$ok) {
                    shiny::showNotification(result$message, type = "error", duration = 10)
                }
            }
            if (samples$ok) {
                shiny::updateSelectInput(
                    session, "sample",
                    choices = samples$value,
                    selected = keep_selection(input$sample, samples$value)
                )
            }
            if (references$ok) {
                shiny::updateSelectInput(
                    session, "reference",
                    choices = references$value,
                    selected = keep_selection(input$reference, references$value)
                )
            }
            update_window_choices()
        })

        shiny::observeEvent(c(input$sample, input$reference), {
            if (identical(input$source, "files")) {
                update_window_choices()
            }
        }, ignoreInit = TRUE)

        load_input <- function() {
            if (identical(input$source, "example")) {
                return(example_input_data())
            }
            if (is.null(input$window) || !nzchar(input$window)) {
                return(NULL)
            }
            readMaMaMIAInput(
                input$analysis_dir, input$assemblies_dir,
                sample = input$sample,
                reference = input$reference,
                window = input$window
            )
        }

        shiny::observeEvent(list(input$source, input$sample, input$reference, input$window), {
            data <- tryCatch(load_input(), error = function(e) {
                shiny::showNotification(
                    paste("Could not load input:", conditionMessage(e)),
                    type = "error", duration = 10
                )
                NULL
            })
            input_data(data)
            update_subgenome_choices(data)
        }, ignoreNULL = FALSE)

        pairs_preview <- shiny::reactive({
            data <- input_data()
            if (is.null(data) || is.null(input$don_subgenome) || is.null(input$rec_subgenome)) {
                return(character())
            }
            tryCatch(
                defaultTargetPairs(
                    data$meta, data$don_chr_ids, data$rec_chr_ids,
                    don_subgenomes = input$don_subgenome,
                    rec_subgenomes = input$rec_subgenome
                ),
                error = function(e) character()
            )
        })

        output$input_summary <- shiny::renderText({
            data <- input_data()
            if (is.null(data)) {
                return(paste(
                    "Pick an input source in the sidebar (a directory layout or the",
                    "packaged example), then press Run analysis."
                ))
            }
            pairs <- pairs_preview()
            sprintf(
                "%s | reference: %s | window label: w%s | %d windows | donor: %d chromosome(s) in %s | recipient: %d chromosome(s) in %s | %d donor-recipient pair(s) will be segmented (%s -> %s)",
                data$sample,
                data$reference,
                data$window,
                nrow(data$dcov) + nrow(data$rcov),
                length(data$don_chr_ids),
                paste(subgenome_labels(data$meta, data$don_chr_ids), collapse = ", "),
                length(data$rec_chr_ids),
                paste(subgenome_labels(data$meta, data$rec_chr_ids), collapse = ", "),
                length(pairs),
                if (is.null(input$don_subgenome)) "?" else input$don_subgenome,
                if (is.null(input$rec_subgenome)) "?" else input$rec_subgenome
            )
        })

        shiny::observeEvent(input$run, {
            data <- input_data()
            if (is.null(data)) {
                shiny::showNotification("Load input data first", type = "warning", duration = 8)
                return()
            }
            result <- tryCatch(
                shiny::withProgress(message = "MaMaMIA analysis", value = 0, {
                    shiny::incProgress(0.05, detail = "Building the read count array")
                    ## The GC-bias fit dominates the runtime and blocks the session, so
                    ## the progress detail is announced before the call starts.
                    shiny::incProgress(0.15, detail = "Fitting the GC-bias model (can take minutes)")
                    isa_new <- analyze_input_data(
                        data,
                        cores = input$cores,
                        don_subgenomes = input$don_subgenome,
                        rec_subgenomes = input$rec_subgenome,
                        cache = cache,
                        verbose = FALSE
                    )
                    shiny::incProgress(0.8, detail = "Segmenting pairwise differences")
                    isa_new
                }),
                error = function(e) {
                    shiny::showNotification(
                        paste("Analysis failed:", conditionMessage(e)),
                        type = "error", duration = 15
                    )
                    NULL
                }
            )
            if (!is.null(result)) {
                isa(result)
                shiny::showNotification(
                    sprintf(
                        "%d segment(s), %d flagged as putative",
                        nrow(result$segments),
                        sum(result$segments$putative_introgression)
                    ),
                    type = "message", duration = 6
                )
            }
        })

        explorer_server("explorer", isa)
    }

    shiny::shinyApp(ui, server)
}

#' Run MaMaMIA as a local Shiny dashboard
#'
#' `runMaMaMIA()` is the main entry point of the package. It opens a locally
#' running Shiny application that reads windowed coverage and GC-content
#' tables, runs the MaMaMIA pipeline, and then hands the result to the
#' interactive introgression explorer of [exploreIntrogressions()].
#'
#' The input is either the packaged `triticum` example or a directory layout
#' read by [readMaMaMIAInput()]: an analysis directory with per-sample coverage
#' tables and an assemblies directory with per-reference GC-content tables and
#' chromosome metadata. Pick the sample, reference, donor and recipient
#' subgenomes in the app, then press *Run analysis*. Donor-recipient chromosome
#' pairs are matched by chromosome number with [defaultTargetPairs()].
#'
#' The app runs the pipeline in the R session, so the GC-bias fit keeps the
#' session busy while it runs; a fit cache (see the `cache` argument of
#' [correctReadCounts()]) makes re-runs of the same windows much faster.
#'
#' The explorer currently shows the segment table with the BEDPE export; its
#' BioCircos circos view is temporarily disabled together with the BioCircos
#' dependency.
#'
#' Dashboard packages are optional MaMaMIA dependencies. Install them with
#' `install.packages(c("shiny", "bslib", "shinyWidgets", "DT"))`
#' before calling this function.
#'
#' @param analysis_dir Optional path to an analysis directory, prefilling the
#'   corresponding input in the app.
#' @param assemblies_dir Optional path to an assemblies directory, prefilling
#'   the corresponding input in the app.
#' @param cache Optional path to an `.rds` file used to memoise the GC-bias
#'   model fit, passed to [correctReadCounts()]. Defaults to a file in the
#'   user cache directory.
#' @param host Host passed to [shiny::runApp()].
#' @param port Port passed to [shiny::runApp()]. Use `NULL` to let Shiny choose.
#' @param launch.browser Logical passed to [shiny::runApp()].
#'
#' @return The return value of [shiny::runApp()]. Called for its side effect of
#'   launching the dashboard.
#' @seealso [exploreIntrogressions()], [readMaMaMIAInput()]
#' @examples
#' \dontrun{
#' runMaMaMIA(
#'   analysis_dir = "mamamia_analysis",
#'   assemblies_dir = "assemblies"
#' )
#' }
#' @export
runMaMaMIA <- function(analysis_dir = NULL, assemblies_dir = NULL, cache = NULL,
                       host = "127.0.0.1", port = NULL, launch.browser = interactive()) {
    require_dashboard_packages()
    if (is.null(cache)) {
        cache <- default_fit_cache()
    }

    app <- launcher_app(analysis_dir, assemblies_dir, cache)
    args <- list(appDir = app, host = host, launch.browser = launch.browser)
    if (!is.null(port)) {
        args$port <- port
    }
    do.call(shiny::runApp, args)
}
