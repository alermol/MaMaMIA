test_that("dashboard helpers annotate segment coordinates", {
    isa <- seg_step()
    dash <- isa_dashboard_segments(isa)

    expect_equal(dash$.row_id, seq_len(nrow(isa$segments)))
    expect_true(all(c("start.don", "stop.don", "start.rec", "stop.rec") %in% names(dash)))
    expect_true(all(dash$start.don <= dash$stop.don))
    expect_true(all(dash$start.rec <= dash$stop.rec))
})

test_that("dashboard filtering respects chromosomes, pairs, and putative flag", {
    isa <- seg_step()
    dash <- isa_dashboard_segments(isa)

    expect_equal(nrow(filter_dashboard_segments(dash, "d1", "d1 -> r1", FALSE)), 0L)
    filtered <- filter_dashboard_segments(dash, c("d1", "r1"), "d1 -> r1", TRUE)
    expect_true(all(filtered$putative_introgression))
    expect_equal(nrow(filtered), sum(dash$putative_introgression))
})

test_that("selected_isa_for_bedpe() marks selected source rows only", {
    isa <- seg_step()
    selected <- selected_isa_for_bedpe(isa, c(1L, 3L))

    expect_s3_class(selected, "ISA")
    expect_equal(which(selected$segments$putative_introgression), intersect(c(1L, 3L), seq_len(nrow(isa$segments))))
    expect_false(identical(selected$segments$putative_introgression, isa$segments$putative_introgression))
})

test_that("filter_dashboard_isa() follows the table filters", {
    isa <- seg_step()
    filtered <- filter_dashboard_isa(
        isa,
        chr_ids = c("d1", "r1"),
        pair_ids = "d1 -> r1",
        min_width = 2L
    )

    expect_s3_class(filtered, "ISA")
    expect_true(all(filtered$out$chr_id.don == "d1"))
    expect_true(all(filtered$out$chr_id.rec == "r1"))
    expect_equal(nrow(filtered$segments), nrow(isa$segments))

    expect_null(filter_dashboard_isa(isa, chr_ids = "nope"))
    expect_null(filter_dashboard_isa(isa, pair_ids = "nope"))
    expect_null(filter_dashboard_isa(isa, chr_ids = character()))

    wide <- filter_dashboard_isa(isa, min_width = 1e9)
    expect_s3_class(wide, "ISA")
    expect_equal(nrow(wide$segments), 0L)
})

test_that("dashboard genome uses chromosome lengths from ISA windows", {
    isa <- seg_step()
    genome <- isa_dashboard_genome(isa)

    expect_named(genome, c("d1", "r1"), ignore.order = TRUE)
    expect_true(all(unlist(genome) > 0))
    expect_named(isa_dashboard_genome(isa, "d1"), "d1")
})

test_that("the results card exposes the plot controls", {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("bslib")
    testthat::skip_if_not_installed("DT")

    html <- as.character(explorer_results_ui("explorer"))
    expect_true(any(grepl("explorer-plot_type", html)))
    expect_true(any(grepl("explorer-isa_plot", html)))
})

test_that("the explorer renders mirror and difference plots", {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("DT")
    testthat::skip_if_not_installed("shinyWidgets")

    isa <- seg_step()
    shiny::testServer(
        explorer_server,
        args = list(ISA = shiny::reactive(isa)),
        {
            session$setInputs(
                chr_ids = c("d1", "r1"),
                pair_ids = "d1 -> r1",
                putative_only = TRUE,
                threshold_type = "median",
                threshold_value = 0,
                min_width = 2,
                plot_type = "mirror"
            )
            expect_match(output$isa_plot$src, "^data:image/png")

            session$setInputs(plot_type = "diff")
            expect_match(output$isa_plot$src, "^data:image/png")

            session$setInputs(plot_type = "unexpected")
            expect_match(output$isa_plot$src, "^data:image/png")
        }
    )
})

test_that("the explorer plot validates without data or filters", {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("DT")
    testthat::skip_if_not_installed("shinyWidgets")

    shiny::testServer(
        explorer_server,
        args = list(ISA = shiny::reactive(NULL)),
        {
            session$setInputs(plot_type = "mirror")
            expect_error(output$isa_plot, "Run the analysis to plot")
        }
    )

    isa <- seg_step()
    shiny::testServer(
        explorer_server,
        args = list(ISA = shiny::reactive(isa)),
        {
            session$setInputs(
                chr_ids = character(),
                pair_ids = character(),
                threshold_type = "median",
                threshold_value = 0,
                min_width = 2,
                plot_type = "mirror"
            )
            expect_error(
                output$isa_plot,
                "Select at least one chromosome and pair to plot"
            )
        }
    )
})

test_that("the explorer download carries the sample name", {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("DT")
    testthat::skip_if_not_installed("shinyWidgets")

    isa <- seg_step()
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)
    shiny::testServer(
        explorer_server,
        args = list(ISA = shiny::reactive(isa), sample = shiny::reactive("sample-a")),
        {
            session$setInputs(
                chr_ids = c("d1", "r1"),
                pair_ids = "d1 -> r1",
                putative_only = TRUE,
                threshold_type = "median",
                threshold_value = 0,
                min_width = 2,
                segments_rows_selected = 1L
            )
            ## The mock session removes its download files when it closes, so
            ## copy the export while the session is still alive.
            expect_true(file.copy(output$download_bedpe, target))
        }
    )

    header <- strsplit(readLines(target)[1L], "\t")[[1L]]
    expect_identical(header[[length(header)]], "sample")
    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(nrow(out), 1L)
    expect_identical(out$sample, "sample-a")
})

test_that("the explorer download keeps the legacy layout without a sample", {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("DT")
    testthat::skip_if_not_installed("shinyWidgets")

    isa <- seg_step()
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)
    shiny::testServer(
        explorer_server,
        args = list(ISA = shiny::reactive(isa)),
        {
            session$setInputs(
                chr_ids = c("d1", "r1"),
                pair_ids = "d1 -> r1",
                putative_only = TRUE,
                threshold_type = "median",
                threshold_value = 0,
                min_width = 2,
                segments_rows_selected = 1L
            )
            expect_true(file.copy(output$download_bedpe, target))
        }
    )

    header <- strsplit(readLines(target)[1L], "\t")[[1L]]
    expect_identical(length(header), 14L)
    expect_false("sample" %in% header)
})

test_that("the explorer preselects rows from a selection map", {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("DT")
    testthat::skip_if_not_installed("shinyWidgets")

    isa <- seg_step()
    shiny::testServer(
        explorer_server,
        args = list(
            ISA = shiny::reactive(isa),
            selection = shiny::reactive(c(1L, 2L))
        ),
        {
            session$setInputs(
                chr_ids = c("d1", "r1"),
                pair_ids = "d1 -> r1",
                putative_only = FALSE,
                threshold_type = "median",
                threshold_value = 0,
                min_width = 2
            )
            rendered <- jsonlite::fromJSON(
                output$segments,
                simplifyVector = FALSE
            )
            expect_identical(
                as.integer(unlist(rendered$x$selection$selected)),
                c(1L, 2L)
            )
        }
    )
})

test_that("the explorer reports selections with the current sample", {
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("DT")
    testthat::skip_if_not_installed("shinyWidgets")

    isa <- seg_step()
    seen <- list()
    shiny::testServer(
        explorer_server,
        args = list(
            ISA = shiny::reactive(isa),
            sample = shiny::reactive("sample-a"),
            selection = shiny::reactive(integer()),
            on_selection = function(rows, sample) {
                seen[[length(seen) + 1L]] <<- list(
                    rows = rows, sample = sample
                )
            }
        ),
        {
            session$setInputs(
                chr_ids = c("d1", "r1"),
                pair_ids = "d1 -> r1",
                putative_only = FALSE,
                threshold_type = "median",
                threshold_value = 0,
                min_width = 2,
                segments_rows_selected = c(1L, 2L)
            )
            picked <- Filter(function(x) {
                identical(x$sample, "sample-a") &&
                    identical(x$rows, c(1L, 2L))
            }, seen)
            expect_gt(length(picked), 0L)
        }
    )
})
