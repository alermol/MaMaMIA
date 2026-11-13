## Server-side language switching without rebuilding stateful inputs.
##
## The launcher server owns the per-session language and the file input lives
## in a renderUI whose expression deliberately does not read language(). The
## test asserts the label flips while that HTML stays byte-identical, which is
## the observable proxy for "the uploaded BAM was not re-read or rebuilt".

test_that("language defaults to en and toggling updates a header label", {
    skip_if_not_installed("shiny")
    state <- bam_launcher_state(tiny_config())
    expect_true(state$ok)

    shiny::testServer(bam_launcher_server(state), {
        expect_match(
            output$sidebar_input_label$html, "Input", fixed = TRUE
        )
        expect_match(
            output$mamamia_header_note$html,
            "Breadth-aware introgression mapping",
            fixed = TRUE
        )

        upload_before <- output$bam_upload_ui$html
        log_before <- output$calculation_log

        session$setInputs(mamamia_language = "ru")

        expect_match(
            output$sidebar_input_label$html,
            "Входные данные", fixed = TRUE
        )
        expect_match(
            output$mamamia_header_note$html,
            "Картирование интрогрессий по широте покрытия",
            fixed = TRUE
        )

        ## The stateful upload control was not rebuilt by the language change.
        expect_identical(output$bam_upload_ui$html, upload_before)
        expect_identical(output$calculation_log, log_before)
    })
})

test_that("an unknown language input falls back to en", {
    skip_if_not_installed("shiny")
    state <- bam_launcher_state(tiny_config())
    shiny::testServer(bam_launcher_server(state), {
        session$setInputs(mamamia_language = "de")
        expect_match(
            output$sidebar_input_label$html, "Input", fixed = TRUE
        )
    })
})
test_that("pass_only is preserved across a language switch", {
    skip_if_not_installed("shiny")
    state <- bam_launcher_state(tiny_config())
    shiny::testServer(bam_launcher_server(state), {
        session$setInputs(`explorer-pass_only` = FALSE)
        session$setInputs(mamamia_language = "ru")
        html <- output$view_filters_ui$html
        expect_false(grepl(
            'id="explorer-pass_only"[^>]*checked', html, perl = TRUE
        ))
    })
})

test_that("launcher server constructs outside an active reactive context", {
    skip_if_not_installed("shiny")
    ## testServer supplies a reactive domain, so it cannot catch an eager
    ## reactiveVal() read during session construction. Call the real server
    ## constructor at top level with a mock session instead. moduleServer()
    ## is stubbed because it requires a full Shiny session; initialization
    ## (including the initial log) happens before it and is what this guards.
    testthat::local_mocked_bindings(
        moduleServer = function(...) invisible(NULL),
        .package = "shiny"
    )
    state <- bam_launcher_state(tiny_config())
    server <- bam_launcher_server(state)
    session <- shiny::MockShinySession$new()
    expect_silent(server(session$input, session$output, session))
})
