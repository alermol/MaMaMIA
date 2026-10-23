## T6.1-T6.4: BAM-first Shiny application.

skip_if_no_dashboard <- function() {
    for (package in c("shiny", "bslib", "DT", "gridlayout", "shinyWidgets")) {
        testthat::skip_if_not_installed(package)
    }
}

app_state <- function(fixture) {
    state <- bam_launcher_state(fixture$config_path)
    expect_true(state$ok)
    state
}

app_inputs <- function(...) {
    defaults <- list(
        resolution = "100",
        reference_sample = ""
    )
    utils::modifyList(defaults, list(...))
}

set_app_inputs <- function(session, ...) {
    do.call(session$setInputs, app_inputs(...))
}

explorer_inputs <- function(...) {
    defaults <- list(
        `explorer-chr_ids` = c("d1", "r1"),
        `explorer-pair_ids` = "d1 -> r1",
        `explorer-threshold_type` = "median",
        `explorer-threshold_value` = 0,
        `explorer-min_width` = 2,
        `explorer-putative_only` = TRUE,
        `explorer-plot_type` = "mirror"
    )
    utils::modifyList(defaults, list(...))
}

set_explorer_inputs <- function(session, ...) {
    do.call(session$setInputs, explorer_inputs(...))
}

## `renderUI()` values arrive as a list of tags plus a trailing dependency
## marker in `testServer`; flatten them for regex assertions.
output_html <- function(output) {
    paste(as.character(output), collapse = "")
}

## Simulate the data frame Shiny produces for `fileInput("bams")`.
upload_input <- function(paths) {
    paths <- as.character(paths)
    data.frame(
        datapath = paths,
        name = basename(paths),
        stringsAsFactors = FALSE
    )
}

fake_bam_result <- function(samples = c("A", "B"), isas = NULL,
                            counts = NULL) {
    if (is.null(isas)) {
        isas <- stats::setNames(rep(list(seg_step()), length(samples)), samples)
    }
    if (is.null(counts)) {
        side <- matrix(
            10, nrow = 60L, ncol = length(samples),
            dimnames = list(NULL, samples)
        )
        counts <- list(donor = side, recipient = side)
    }
    mode <- if (length(samples) > 1L) "median" else "none"
    paths <- if (length(samples) == 0L) {
        character()
    } else {
        paste0(samples, ".bam")
    }
    structure(
        list(
            assembly = NULL,
            config = NULL,
            bams = data.frame(
                path = paths,
                sample = samples,
                index_path = rep("", length(samples)),
                n_targets = rep(2L, length(samples)),
                stringsAsFactors = FALSE
            ),
            samples = data.frame(
                sample = samples,
                path = paths,
                assembly = rep("variable", length(samples)),
                index_path = rep("", length(samples)),
                reference = rep(FALSE, length(samples)),
                normalization = rep(mode, length(samples)),
                n_windows = rep(60L, length(samples)),
                n_segments = vapply(
                    isas, function(x) nrow(x$segments), integer(1)
                ),
                n_putative_introgressions = vapply(
                    isas,
                    function(x) {
                        sum(x$segments$putative_introgression, na.rm = TRUE)
                    },
                    integer(1)
                ),
                stringsAsFactors = FALSE
            ),
            windows = list(
                table = data.frame(
                    chr_id = rep(c("d1", "r1"), each = 30L),
                    start = 0L,
                    stop = 1L,
                    role = rep(c("donor", "recipient"), each = 30L),
                    chr_name = "",
                    subgenome = "",
                    stringsAsFactors = FALSE
                ),
                n_per_chromosome = 30L
            ),
            gc = data.frame(
                chr_id = character(), start = integer(),
                stop = integer(), gc = numeric()
            ),
            counts = counts,
            rca = list(),
            isa = isas,
            normalization = list(
                mode = mode, reference = NULL,
                by_side = TRUE, pseudocount = 0.5
            ),
            provenance = list()
        ),
        class = c("BamAnalysisResult", "list")
    )
}

test_that("a missing configuration renders an error card", {
    skip_if_no_dashboard()
    old <- Sys.getenv("MAMAMIA_CONFIG", unset = NA_character_)
    on.exit(
        if (is.na(old)) {
            Sys.unsetenv("MAMAMIA_CONFIG")
        } else {
            Sys.setenv(MAMAMIA_CONFIG = old)
        },
        add = TRUE
    )
    Sys.unsetenv("MAMAMIA_CONFIG")

    state <- bam_launcher_state(NULL)
    expect_false(state$ok)
    expect_identical(state$error$stage, "config")
    expect_match(state$error$message, "MAMAMIA_CONFIG")

    html <- paste(as.character(bam_launcher_ui(state)), collapse = "")
    expect_match(html, "Configuration error")
    expect_match(html, "MAMAMIA_CONFIG")
})

test_that("an invalid configuration card names the registry file", {
    skip_if_no_dashboard()
    path <- tempfile(fileext = ".yml")
    on.exit(unlink(path), add = TRUE)
    writeLines(c("version: 1", "assemblies: {}"), path)

    state <- bam_launcher_state(path)
    expect_false(state$ok)
    expect_identical(state$error$stage, "config")
    expect_identical(state$error$files, path)

    html <- paste(as.character(bam_launcher_ui(state)), collapse = "")
    expect_match(html, "Configuration error")
    expect_match(html, "File:")
    expect_match(html, basename(path), fixed = TRUE)
})

test_that("MAMAMIA_CONFIG supplies the configuration when config is NULL", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    old <- Sys.getenv("MAMAMIA_CONFIG", unset = NA_character_)
    on.exit(
        if (is.na(old)) {
            Sys.unsetenv("MAMAMIA_CONFIG")
        } else {
            Sys.setenv(MAMAMIA_CONFIG = old)
        },
        add = TRUE
    )
    Sys.setenv(MAMAMIA_CONFIG = fixture$config_path)

    state <- bam_launcher_state(NULL)
    expect_true(state$ok)
    expect_identical(state$source, "environment")
    expect_identical(state$config$assemblies$variable$label,
        "Variable-GC synthetic assembly"
    )
})

test_that("resolution choices follow the configuration", {
    defaults <- .config_defaults()
    choices <- bam_resolution_choices(defaults)
    expect_identical(unname(choices), c("100", "1000", "10000"))
    expect_match(names(choices)[[1L]], "Low")
    expect_match(names(choices)[[2L]], "Medium")
    expect_match(names(choices)[[3L]], "High")
    expect_identical(bam_default_resolution(defaults), "1000")

    limited <- list(
        allowed_windows = c(250L, 500L),
        windows_per_chromosome = 250L
    )
    expect_identical(bam_default_resolution(limited), "250")
    expect_identical(bam_resolution_value("500", limited), 500L)
    expect_error(bam_resolution_value("1000", limited), "one of",
        class = "MaMaMIAWindowError"
    )
    expect_error(bam_resolution_value(NULL, limited), "resolution",
        class = "MaMaMIAWindowError"
    )
})

test_that("the assembly card counts windows at every resolution", {
    state <- bam_launcher_state(tiny_config_path())
    expect_true(state$ok)
    assembly <- state$config$assemblies$tiny

    info <- bam_assembly_info(assembly, c(100L, 1000L))
    expect_identical(info$key, "tiny")
    expect_identical(info$donor, "d1 (1D)")
    expect_identical(info$recipient, "r1 (1R)")
    expect_identical(info$pairs, "d1 -> r1")
    expect_identical(info$windows$n_windows, c(200L, 2000L))

    short <- list(
        key = "short",
        label = "Short",
        description = "",
        pairs = list(),
        contigs = list(
            recipient = data.frame(
                id = "r1", name = "1R", subgenome = "R", length = 500,
                stringsAsFactors = FALSE
            ),
            donor = data.frame(
                id = "d1", name = "1D", subgenome = "D", length = 2000,
                stringsAsFactors = FALSE
            )
        )
    )
    expect_identical(
        bam_assembly_windows(short, c(100L, 1000L)),
        c(200L, NA_integer_)
    )
    expect_match(
        paste(
            as.character(bam_assembly_summary_ui(short, c(100L, 1000L))),
            collapse = ""
        ),
        "too short"
    )
})

test_that("the preflight card lists file, sample, targets and index", {
    table <- data.frame(
        path = "/tmp/session/001-sample.bam",
        source = "/srv/mamamia/bams/sample.bam",
        sample = "sample",
        assembly = "tiny",
        n_targets = 3L,
        n_assembly_targets = 2L,
        index_action = "reuse",
        warnings = "",
        stringsAsFactors = FALSE
    )
    display <- bam_preflight_display(table)
    expect_identical(
        names(display),
        c("file", "sample", "assembly", "targets", "index", "warnings")
    )
    expect_identical(display$file, "sample.bam")
    expect_identical(display$targets, "3/2")
    expect_identical(display$index, "reuse")
    expect_identical(nrow(bam_preflight_display(NULL)), 0L)
})

test_that("errors carry a stage and an actionable suggestion", {
    cases <- list(
        config = quote(.config_error("bad config")),
        BAM = quote(.mamamia_error("MaMaMIABamHeaderError", "bad bam")),
        assembly = quote(.mamamia_error(
            "MaMaMIAAssemblyMatchError", "no match"
        )),
        index = quote(.mamamia_error("MaMaMIAIndexError", "no index")),
        windows = quote(.mamamia_error("MaMaMIAWindowError", "too short")),
        GC = quote(.mamamia_error("MaMaMIAGcrankError", "bad gc")),
        counts = quote(.mamamia_error("MaMaMIACountError", "bad counts")),
        correction = quote(.mamamia_error(
            "MaMaMIASampleError", "Sample 's' (GC correction): too few"
        )),
        normalization = quote(.mamamia_error(
            "MaMaMIAAnalysisError", "Normalization failed: x"
        )),
        segmentation = quote(.mamamia_error(
            "MaMaMIASampleError", "Sample 's' (segmentation): boom"
        )),
        limits = quote(.mamamia_error("MaMaMIALimitError", "too many"))
    )
    for (stage in names(cases)) {
        error <- tryCatch(eval(cases[[stage]]), error = function(e) e)
        info <- bam_error_info(error)
        expect_identical(info$stage, stage, info = stage)
        expect_true(nzchar(info$suggestion), info = stage)
        expect_match(
            bam_error_text(info),
            paste0("\\[", stage, "\\]"),
            info = stage
        )
    }

    plain <- tryCatch(stop("plain failure", call. = FALSE), error = function(e) e)
    expect_identical(bam_error_info(plain)$stage, "analysis")

    table <- data.frame(
        sample = c("A", "B"),
        source = c("/srv/bams/a.bam", "/srv/bams/b.bam"),
        stringsAsFactors = FALSE
    )
    sample_error <- tryCatch(
        .mamamia_error(
            "MaMaMIASampleError",
            "Sample 'B' (GC correction): boom"
        ),
        error = function(e) e
    )
    info <- bam_error_info(
        sample_error, files = "/tmp/fallback.bam", table = table
    )
    expect_identical(info$files, "/srv/bams/b.bam")
    expect_match(bam_error_text(info), "File: /srv/bams/b.bam", fixed = TRUE)
    expect_match(bam_error_text(info), "Suggested fix:")

    fallback <- bam_error_info(plain, files = c("/tmp/a.bam", "/tmp/b.bam"))
    expect_match(
        bam_error_text(fallback), "Files: /tmp/a.bam, /tmp/b.bam",
        fixed = TRUE
    )
    expect_null(bam_error_file_label(character()))
    expect_match(
        bam_error_file_label(paste0("/tmp/", 1:7, ".bam")),
        "\\(\\+2 more\\)"
    )
})

test_that("a contig shorter than the resolution blocks the run", {
    short <- list(
        key = "short",
        label = "Short",
        description = "",
        pairs = list(),
        contigs = list(
            recipient = data.frame(
                id = "r1", name = "1R", subgenome = "R", length = 500,
                stringsAsFactors = FALSE
            ),
            donor = data.frame(
                id = "d1", name = "1D", subgenome = "D", length = 2000,
                stringsAsFactors = FALSE
            )
        )
    )
    expect_null(bam_resolution_problem(short, 100L))
    problem <- bam_resolution_problem(short, 1000L)
    expect_match(problem, "Contig 'r1' \\(500 bp\\) is shorter")
    expect_match(problem, "choose a lower resolution")
    expect_null(bam_resolution_problem(NULL, 1000L))
})

test_that("progress stages follow the cohort size and normalization", {
    expect_identical(
        bam_progress_stages(1L, "none"),
        c(
            "preflight", "index", "windows", "gc", "counts",
            "correction", "segmentation"
        )
    )
    cohort <- bam_progress_stages(3L, "median")
    expect_identical(sum(cohort == "counts"), 3L)
    expect_identical(sum(cohort == "correction"), 3L)
    expect_identical(sum(cohort == "segmentation"), 3L)
    expect_identical(sum(cohort == "normalization"), 1L)
    expect_identical(length(bam_progress_stages(3L, "none")), 13L)
    expect_identical(
        bam_progress_label("counts", list(sample = "s1")),
        "Counting reads for s1"
    )
    expect_identical(bam_progress_label("preflight"), "Checking BAM headers")
    expect_identical(
        bam_progress_label(
            "preflight",
            list(n_bams = 2L, n_samples = 3L, assembly = "tiny")
        ),
        "Checking BAM headers (2 BAM files, 3 samples)"
    )
    expect_identical(
        bam_progress_label("index", list(n_bams = 2L, staged = 1L)),
        "Ensuring CSI indexes (1 of 2 staged)"
    )
    expect_identical(
        bam_progress_label(
            "windows", list(n_windows = 2000L, n_per_chromosome = 1000L)
        ),
        "Generating analysis windows (2,000 windows, 1,000 per chromosome)"
    )
    expect_identical(
        bam_progress_label("gc", list(n_windows = 2000L)),
        "Reading GC content for 2,000 windows"
    )
    expect_identical(
        bam_progress_label("counts", list(sample = "s1", total = 1234)),
        "Counting reads for s1 (1,234 reads)"
    )
    expect_identical(
        bam_progress_label(
            "correction", list(sample = "s1", n_ideal = 8000L)
        ),
        "Correcting GC bias for s1 (8,000 ideal windows)"
    )
    expect_identical(
        bam_progress_label(
            "normalization", list(n_samples = 3L, mode = "median")
        ),
        "Normalizing samples (3 samples, median)"
    )
    expect_identical(
        bam_progress_label(
            "normalization",
            list(n_samples = 3L, mode = "reference", reference = "A")
        ),
        "Normalizing samples (3 samples, reference A)"
    )
    expect_identical(
        bam_progress_label(
            "segmentation", list(sample = "s1", n_segments = 12L)
        ),
        "Segmenting chromosome pairs for s1 (12 segments)"
    )
})

test_that("run parameters use the registry defaults as-is", {
    config <- list(defaults = .config_defaults())
    base <- list(resolution = "1000", reference_sample = "")

    single <- bam_run_parameters(base, config, "s1")
    expect_identical(single$resolution, 1000L)
    expect_identical(single$normalization, "none")
    expect_null(single$reference_sample)
    expect_identical(single$mapq, as.integer(config$defaults$mapq))
    expect_identical(
        single$exclude_flags, as.integer(config$defaults$exclude_flags)
    )
    expect_identical(single$max_fit_windows, config$defaults$max_fit_windows)
    expect_identical(single$cores, min(30L, bam_max_cores()))

    cohort <- bam_run_parameters(base, config, c("s1", "s2"))
    expect_identical(cohort$normalization, "median")
    expect_null(cohort$reference_sample)

    reference <- bam_run_parameters(
        utils::modifyList(base, list(reference_sample = "s2")),
        config, c("s1", "s2")
    )
    expect_identical(reference$normalization, "reference")
    expect_identical(reference$reference_sample, "s2")

    expect_error(
        bam_run_parameters(
            utils::modifyList(base, list(reference_sample = "nope")),
            config, c("s1", "s2")
        ),
        "not one of"
    )

    bad_mapq <- config
    bad_mapq$defaults$mapq <- 500
    expect_error(
        bam_run_parameters(base, bad_mapq, "s1"),
        "MAPQ",
        class = "MaMaMIAAnalysisError"
    )
    bad_flags <- config
    bad_flags$defaults$exclude_flags <- 9000
    expect_error(
        bam_run_parameters(base, bad_flags, "s1"),
        "Exclude flags",
        class = "MaMaMIAAnalysisError"
    )
    bad_fit <- config
    bad_fit$defaults$max_fit_windows <- 0
    expect_error(
        bam_run_parameters(base, bad_fit, "s1"),
        "Max fit windows",
        class = "MaMaMIAAnalysisError"
    )
})

test_that("the cohort summary reports reads, calls and normalization", {
    result <- fake_bam_result(c("A", "B"))
    summary <- bam_cohort_summary(result)

    expect_identical(summary$sample, c("A", "B"))
    expect_identical(summary$reads_donor, c(600, 600))
    expect_identical(summary$reads_recipient, c(600, 600))
    expect_identical(summary$reads_total, c(1200, 1200))
    expect_identical(
        summary$n_segments,
        c(nrow(result$isa$A$segments), nrow(result$isa$B$segments))
    )
    expect_identical(summary$normalization, c("median", "median"))
    expect_identical(summary$reference, c("", ""))

    empty <- bam_cohort_summary(fake_bam_result(character()))
    expect_identical(nrow(empty), 0L)
    expect_identical(bam_sample_choices(NULL), character())
})

test_that("coverage warnings find zero, near-zero and all-zero references", {
    zero <- fake_bam_result(
        "A",
        counts = list(
            donor = matrix(0, 60L, 1L, dimnames = list(NULL, "A")),
            recipient = matrix(10, 60L, 1L, dimnames = list(NULL, "A"))
        )
    )
    warnings <- bam_coverage_warnings(zero)
    expect_true(any(grepl("no reads on the donor", warnings)))
    expect_false(any(grepl("recipient", warnings)))

    near <- fake_bam_result(
        "A",
        counts = list(
            donor = matrix(0.1, 60L, 1L, dimnames = list(NULL, "A")),
            recipient = matrix(10, 60L, 1L, dimnames = list(NULL, "A"))
        )
    )
    expect_true(any(grepl("near-zero coverage on the donor", bam_coverage_warnings(near))))

    healthy <- fake_bam_result("A")
    expect_identical(bam_coverage_warnings(healthy), character())

    reference <- fake_bam_result(c("A", "B"))
    reference$samples$reference[[1L]] <- TRUE
    reference$normalization$mode <- "reference"
    reference$counts$donor[1L, "A"] <- 0
    reference$counts$recipient[1L, "A"] <- 0
    reference_warnings <- bam_coverage_warnings(reference)
    expect_true(any(grepl("Reference sample 'A' has 1 donor and 1 recipient", reference_warnings)))
    expect_true(any(grepl("pseudocount", reference_warnings)))
})

test_that("selection maps keep rows per sample", {
    result <- fake_bam_result(c("A", "B"))
    map <- bam_selection_map(result)
    expect_identical(map, list(A = integer(), B = integer()))

    map <- bam_selection_update(map, "A", c(3, 1, 3, NA))
    expect_identical(map$A, c(1L, 3L))
    expect_identical(map$B, integer())

    map <- bam_selection_update(map, NULL, 2L)
    expect_identical(map$A, c(1L, 3L))
    map <- bam_selection_update(map, "B", integer())
    expect_identical(map$B, integer())
})

test_that("the app shows preflight, assembly and a disabled run button", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        expect_match(output_html(output$run_ui), "disabled")
        expect_match(output_html(output$assembly_summary), "Upload BAM files")
        expect_match(output$preflight_table, "Upload BAM files")

        session$setInputs(bams = upload_input(fixture$bam_path))
        pf <- preflight_state()
        expect_true(pf$ok)
        expect_identical(pf$table$sample, "synthetic")
        expect_match(output$preflight_table, "synthetic")
        expect_match(output$preflight_table, "synthetic.bam")
        expect_match(output_html(output$assembly_summary), "Variable-GC")
        expect_match(output_html(output$assembly_summary), "d1 \\(1D\\)")
        expect_match(
            paste(bam_preflight_log(pf), collapse = "\n"),
            "Preflight OK"
        )
        expect_false(grepl("disabled", output_html(output$run_ui)))
    })
})

test_that("the initial calculation log is built without a reactive context", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    time <- as.POSIXct("2026-01-02 03:04:05", tz = "UTC")
    entries <- bam_initial_log(state$config, time = time)
    expect_length(entries, 3L)
    expect_match(
        entries[[1L]],
        "^\\[03:04:05\\] MaMaMIA session started\\.$"
    )
    expect_match(entries[[2L]], "Configuration validated")
    expect_match(
        entries[[3L]],
        "Samples are detected from the BAM header",
        fixed = TRUE
    )
})

test_that("the launcher uses the grid layout with header, sidebar, main and footer", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    html <- paste(as.character(bam_launcher_ui(state)), collapse = "")
    layout <- bam_dashboard_layout()
    expect_match(html, 'id="mamamia-grid"', fixed = TRUE)
    expect_match(html, 'data-gridlayout-area="header"', fixed = TRUE)
    expect_match(html, 'data-gridlayout-area="sidebar"', fixed = TRUE)
    expect_match(html, 'data-gridlayout-area="main"', fixed = TRUE)
    expect_match(html, 'data-gridlayout-area="footer"', fixed = TRUE)
    expect_match(layout, "|0.5rem|240px   |1fr      |", fixed = TRUE)
    expect_match(layout, "|75px  |header  |header   |", fixed = TRUE)
    expect_match(layout, "|1fr   |sidebar |main     |", fixed = TRUE)
    expect_match(layout, "|100px |footer  |footer   |", fixed = TRUE)
    expect_match(html, "mamamia-header", fixed = TRUE)
    expect_match(html, "mamamia-sidebar", fixed = TRUE)
    expect_match(html, "mamamia-main", fixed = TRUE)
    expect_match(html, "mamamia-footer", fixed = TRUE)
    expect_match(html, "MaMaMIA</h1>", fixed = TRUE)
    expect_match(html, "Example data", fixed = TRUE)
    expect_match(html, 'id="explorer-isa_plot"', fixed = TRUE)
})

test_that("one upload control auto-detects the sample count", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    multi_bam <- add_multi_sample_gc_bam(
        fixture, c("sample-a", "sample-b"), steps = c(4L, 6L)
    )
    state <- app_state(fixture)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        upload <- output_html(output$bam_upload_ui)
        expect_match(upload, "BAM file")
        expect_false(grepl("multiple=\"multiple\"", upload, fixed = TRUE))
        expect_match(
            output_html(output$bam_mode_note),
            "samples are detected from the @RG SM"
        )

        session$setInputs(bams = upload_input(fixture$bam_path))
        expect_true(preflight_state()$ok)
        expect_match(
            output_html(output$bam_mode_note),
            "1 sample detected from @RG SM"
        )
        expect_null(output$reference_ui)

        session$setInputs(bams = upload_input(multi_bam))
        expect_true(preflight_state()$ok)
        expect_match(
            output_html(output$bam_mode_note),
            "2 samples detected from @RG SM"
        )
        expect_match(output_html(output$reference_ui), "Reference sample")
    })
})

test_that("the setup appears as single-open spoilers after preflight", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        expect_false(sidebar_unlocked())
        expect_null(output$sidebar_body)

        set_app_inputs(session, bams = upload_input(fixture$bam_path))
        expect_true(preflight_state()$ok)
        expect_true(sidebar_unlocked())

        sidebar <- output_html(output$sidebar_body)
        expect_match(sidebar, "accordion", fixed = TRUE)
        expect_match(sidebar, "Analysis")
        expect_match(sidebar, "Assembly")
        expect_false(grepl("Advanced", sidebar, fixed = TRUE))
        expect_false(grepl("Filters", sidebar, fixed = TRUE))
        expect_false(grepl("Export", sidebar, fixed = TRUE))
        expect_false(grepl("Minimum MAPQ", sidebar, fixed = TRUE))
        expect_false(grepl("Max fit windows", sidebar, fixed = TRUE))
        expect_match(sidebar, "mamamia-run-panel", fixed = TRUE)
        expect_gt(
            regexpr("mamamia-run-panel", sidebar, fixed = TRUE),
            regexpr("accordion", sidebar, fixed = TRUE)
        )
        expect_match(output_html(output$run_ui), "Run analysis")
    })
})

test_that("filters and exports live in the main result area", {
    skip_if_no_dashboard()

    sidebar <- paste(
        as.character(bam_sidebar_details_ui(list(defaults = .config_defaults()))),
        collapse = ""
    )
    main <- paste(as.character(bam_results_ui()), collapse = "")

    expect_false(grepl("Filters", sidebar, fixed = TRUE))
    expect_false(grepl("Export", sidebar, fixed = TRUE))
    expect_match(main, "mamamia-filter-bar", fixed = TRUE)
    expect_match(main, "sample_selector_ui", fixed = TRUE)
    expect_lt(
        regexpr("sample_selector_ui", main, fixed = TRUE),
        regexpr("Export", main, fixed = TRUE)
    )
    expect_match(main, "explorer-chr_ids", fixed = TRUE)
    expect_match(main, "explorer-pair_ids", fixed = TRUE)
    expect_match(main, "explorer-putative_only", fixed = TRUE)
    expect_match(main, "explorer-threshold_type", fixed = TRUE)
    expect_match(main, "explorer-min_width", fixed = TRUE)
    expect_match(main, "explorer-select_all", fixed = TRUE)
    expect_match(main, "explorer-download_bedpe", fixed = TRUE)
    expect_match(main, "cohort_ui", fixed = TRUE)
    expect_match(main, "Export")
})

test_that("the swap button exchanges donor and recipient species", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture(species = TRUE)
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    html <- paste(
        as.character(bam_sidebar_species_ui(state$config)),
        collapse = ""
    )
    expect_match(html, "Swap donor and recipient")
    expect_false(grepl("Detect the assembly from the BAM header", html))

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(
            recipient_species = "recipient", donor_species = "donor"
        )
        expect_identical(
            assembly_selection()$assembly$pairs, c(d1 = "r1")
        )

        session$setInputs(swap_species = 1L)
        expect_identical(recipient_species(), "donor")
        expect_identical(donor_species(), "recipient")
        expect_identical(
            assembly_selection()$key, "donor__D__recipient__R"
        )
        expect_identical(
            assembly_selection()$assembly$pairs, c(r1 = "d1")
        )

        session$setInputs(swap_species = 2L)
        expect_identical(recipient_species(), "recipient")
        expect_identical(donor_species(), "donor")
        expect_identical(
            assembly_selection()$assembly$pairs, c(d1 = "r1")
        )
    })
})

test_that("the assembly always follows the species selectors", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture(species = TRUE)
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)
    resolved <- resolveSpeciesAssembly(
        state$config,
        donor = "donor", recipient = "recipient",
        key = "preset_pair"
    )
    state$config <- resolved$config

    html <- paste(
        as.character(bam_sidebar_species_ui(state$config)),
        collapse = ""
    )
    expect_false(grepl("Detect the assembly from the BAM header", html))
    expect_match(html, "Recipient species")
    expect_match(html, "Donor species")

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(
            recipient_species = "recipient",
            donor_species = "donor",
            bams = upload_input(fixture$bam_path),
            resolution = "100"
        )
        pf <- preflight_state()
        expect_true(pf$ok)
        expect_identical(pf$assembly, "recipient__R__donor__D")
    })
})

test_that("a species-only registry exposes species selectors and runs", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture(species = TRUE)
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    html <- paste(
        as.character(bam_sidebar_species_ui(state$config)),
        collapse = ""
    )
    expect_match(html, "Recipient species")
    expect_match(html, "Donor species")
    expect_false(grepl("Detect the assembly from the BAM header", html))
    expect_false(grepl("Predefined assembly", html, fixed = TRUE))

    ## Even when the registry also declares a preset, the UI offers no
    ## predefined-assembly selector.
    resolved <- resolveSpeciesAssembly(
        state$config, donor = "donor", recipient = "recipient"
    )
    with_preset <- state
    with_preset$config <- resolved$config
    preset_html <- paste(
        as.character(bam_sidebar_species_ui(with_preset$config)),
        collapse = ""
    )
    expect_false(grepl("Predefined assembly", preset_html, fixed = TRUE))
    expect_false(grepl("Detect the assembly from the BAM header", preset_html))

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(
            recipient_species = "recipient",
            donor_species = "donor",
            bams = upload_input(fixture$bam_path),
            resolution = "100",
            mapq = 20,
            exclude_flags = 2308,
            cores = 1,
            max_fit_windows = 1000,
            reference_sample = ""
        )
        pf <- preflight_state()
        expect_true(pf$ok)
        expect_identical(pf$assembly, "recipient__R__donor__D")
        expect_identical(pf$assembly_config$pairs, c(d1 = "r1"))
        expect_match(output_html(output$assembly_summary), "Variable-GC recipient")
        expect_match(output_html(output$assembly_summary), "d1 \\(1D\\)")
        expect_false(grepl("disabled", output_html(output$run_ui)))

        suppressWarnings(session$setInputs(run = 1L))
        result <- result()
        expect_s3_class(result, "BamAnalysisResult")
        expect_identical(result$assembly$key, "recipient__R__donor__D")
        expect_identical(result$bams$sample, "synthetic")
        expect_match(output$cohort_table, "synthetic")
    })
})

test_that("an invalid species pair is shown before any BAM is read", {
    skip_if_no_dashboard()
    fixture <- variable_gc_fixture(species = TRUE)
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(
            recipient_species = "recipient",
            donor_species = "recipient"
        )
        pf <- preflight_state()
        expect_false(pf$ok)
        expect_match(pf$error$message, "two different species")
        log <- paste(bam_preflight_log(pf), collapse = "\n")
        expect_match(log, "two different species")
        expect_match(output_html(output$run_ui), "disabled")
    })
})

test_that("uploaded BAMs are staged before preflight", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    upload <- data.frame(
        name = basename(fixture$bam_path),
        size = file.info(fixture$bam_path)$size,
        type = "application/octet-stream",
        datapath = fixture$bam_path,
        stringsAsFactors = FALSE
    )
    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(bams = upload)
        pf <- preflight_state()
        expect_true(pf$ok)
        expect_false(identical(pf$staged, fixture$bam_path))
        expect_true(all(file.exists(pf$staged)))
        expect_true(all(dirname(pf$staged) != dirname(fixture$bam_path)))
    })
})

test_that("an assembly mismatch shows candidate scores", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    bad <- file.path(fixture$directory, "bad.bam")
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:r1\tLN:20000",
            "@SQ\tSN:d1\tLN:19000",
            "@RG\tID:bad\tSM:bad-sample"
        ),
        bad,
        index = FALSE
    )
    state <- app_state(fixture)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(bams = upload_input(bad))
        pf <- preflight_state()
        expect_false(pf$ok)
        expect_identical(pf$error$stage, "assembly")
        expect_true(is.data.frame(pf$candidates))
        expect_true("variable" %in% pf$candidates$key)
        expect_true(any(pf$candidates$n_length_mismatches > 0L))

        error_log <- paste(bam_preflight_log(pf), collapse = "\n")
        expect_match(error_log, "Preflight failed at stage 'assembly'")
        expect_match(error_log, "Suggested fix")
        expect_match(error_log, "File:")
        expect_match(error_log, "bad\\.bam")
        expect_match(error_log, "Assembly candidates")
        expect_match(error_log, "variable")
        expect_match(error_log, "50.0%")
        expect_match(output_html(output$run_ui), "disabled")
    })
})

test_that("missing BAMs and server limits fail without crashing the session", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    second <- add_variable_gc_bam(fixture, "synthetic-2")
    state <- app_state(fixture)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(
            bams = upload_input(file.path(fixture$directory, "missing.bam"))
        )
        pf <- preflight_state()
        expect_false(pf$ok)
        expect_identical(pf$error$stage, "BAM")
        expect_match(
            paste(bam_preflight_log(pf), collapse = "\n"),
            "Suggested fix"
        )
        expect_match(output_html(output$run_ui), "disabled")
    })

    config <- yaml::read_yaml(fixture$config_path)
    config$server$max_bams <- 1
    limited_path <- file.path(fixture$directory, "limited.yml")
    writeLines(yaml::as.yaml(config), limited_path)
    limited <- bam_launcher_state(limited_path)
    expect_true(limited$ok)

    shiny::testServer(bam_launcher_server(limited, cache = FALSE), {
        session$setInputs(bams = upload_input(fixture$bam_path))
        pf <- preflight_state()
        expect_true(pf$ok)
    })
    expect_error(
        preflightBams(
            c(fixture$bam_path, second), limited$config
        ),
        "Too many BAM files",
        class = "MaMaMIALimitError"
    )

    config$server$max_bams <- 4
    config$server$max_total_bytes <- 1
    byte_limit_path <- file.path(fixture$directory, "byte-limit.yml")
    writeLines(yaml::as.yaml(config), byte_limit_path)
    byte_limited <- bam_launcher_state(byte_limit_path)
    expect_true(byte_limited$ok)

    shiny::testServer(bam_launcher_server(byte_limited, cache = FALSE), {
        session$setInputs(bams = upload_input(fixture$bam_path))
        pf <- preflight_state()
        expect_false(pf$ok)
        expect_identical(pf$error$stage, "limits")
        expect_match(
            paste(bam_preflight_log(pf), collapse = "\n"),
            "exceeds"
        )
        expect_match(output_html(output$run_ui), "disabled")
    })
})

test_that("a single BAM runs end to end and cleans up the session", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)
    before <- list.files(tempdir(), pattern = "^mamamia-session-")
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        set_app_inputs(session, bams = upload_input(fixture$bam_path))
        expect_true(preflight_state()$ok)
        expect_false(grepl("result_sample", output_html(output$sample_selector_ui)))

        suppressWarnings(session$setInputs(run = 1))
        expect_match(output$cohort_table, "synthetic")
        expect_match(output$cohort_table, "none")
        expect_false(grepl("result_sample", output_html(output$sample_selector_ui)))
        expect_false(grepl("result_sample", output_html(output$cohort_ui)))

        set_explorer_inputs(session)
        expect_match(output[["explorer-isa_plot"]]$src, "^data:image/png")

        circos <- jsonlite::fromJSON(
            output[["explorer-circos"]], simplifyVector = FALSE
        )
        modules <- vapply(
            circos$x$moduleList_NGCircos,
            function(module) module[[1L]], character(1)
        )
        expect_true("CNV_difference" %in% modules)
        expect_true("ARC_segment_spans" %in% modules)
        expect_true("LINK_introgressions" %in% modules)

        session$setInputs(`explorer-segments_rows_selected` = 1L)
        expect_identical(selection_map()[["synthetic"]], 1L)
        expect_true(file.copy(output$download_cohort, target))
    })

    expect_false(file.exists(paste0(fixture$bam_path, ".csi")))
    expect_identical(
        setdiff(list.files(tempdir(), pattern = "^mamamia-session-"), before),
        character()
    )
    out <- utils::read.table(
        target,
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
    expect_identical(nrow(out), 1L)
    expect_identical(out$sample, "synthetic")
    expect_identical(names(out)[[ncol(out)]], "sample")
})

test_that("per-sample BEDPE helpers fall back to putative calls", {
    expect_null(bam_sample_selection(list(), "A"))
    expect_null(bam_sample_selection(list(A = integer()), "A"))
    expect_null(bam_sample_selection(list(A = 1L), NULL))
    expect_identical(
        bam_sample_selection(list(A = c(2, 1), B = 3L), "A"),
        list(A = c(1L, 2L))
    )

    expect_identical(
        bam_bedpe_filename("sample-a"),
        "mamamia-sample-a-selected-introgressions.bedpe"
    )
    expect_identical(
        bam_bedpe_filename("sample/a b"),
        "mamamia-sample-a-b-selected-introgressions.bedpe"
    )
    expect_identical(
        bam_bedpe_filename(NULL),
        "mamamia-selected-introgressions.bedpe"
    )
})

test_that("a cohort shows a selector and exports a combined BEDPE", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    bam <- add_multi_sample_gc_bam(
        fixture, c("sample-a", "sample-b"), steps = c(4L, 6L)
    )
    state <- app_state(fixture)
    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)
    sample_target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(sample_target), add = TRUE)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        set_app_inputs(
            session,
            bams = upload_input(bam)
        )
        expect_true(preflight_state()$ok)
        expect_match(
            output_html(output$reference_ui),
            "None \\(median normalization\\)"
        )
        before_run <- output_html(output$sample_selector_ui)
        expect_match(before_run, "result_sample")
        expect_match(before_run, "sample-a")
        expect_match(before_run, "sample-b")

        session$setInputs(reference_sample = "")
        suppressWarnings(session$setInputs(run = 1))
        expect_match(output$cohort_table, "sample-a")
        expect_match(output$cohort_table, "sample-b")
        expect_match(output$cohort_table, "median")
        selector <- output_html(output$sample_selector_ui)
        expect_match(selector, "result_sample")
        expect_match(selector, "sample-a")
        expect_match(selector, "sample-b")
        expect_false(grepl("result_sample", output_html(output$cohort_ui)))
        expect_match(
            output_html(output$cohort_ui),
            "Download selected BEDPE \\(all samples\\)"
        )
        expect_match(
            output_html(output$cohort_ui),
            "Download BEDPE \\(current sample\\)"
        )

        set_explorer_inputs(session)
        session$setInputs(`explorer-segments_rows_selected` = 1L)
        expect_identical(selection_map()[["sample-a"]], 1L)

        session$setInputs(
            result_sample = "sample-b",
            `explorer-segments_rows_selected` = NULL
        )
        session$setInputs(`explorer-segments_rows_selected` = 1L)
        expect_identical(selection_map()[["sample-b"]], 1L)
        expect_identical(selection_map()[["sample-a"]], 1L)

        ## Switching back restores the rows chosen for the first sample.
        session$setInputs(result_sample = "sample-a")
        rendered <- jsonlite::fromJSON(
            output[["explorer-segments"]],
            simplifyVector = FALSE
        )
        expect_identical(
            as.integer(unlist(rendered$x$selection$selected)),
            1L
        )

        expect_true(file.copy(output$download_cohort, target))

        ## The per-sample export follows the selector and the stored rows.
        expect_true(file.copy(output$download_sample, sample_target))
    })

    out <- utils::read.table(
        target,
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
    expect_identical(nrow(out), 2L)
    expect_identical(sort(unique(out$sample)), c("sample-a", "sample-b"))
    expect_identical(names(out)[[ncol(out)]], "sample")

    per_sample <- utils::read.table(
        sample_target,
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
    expect_identical(nrow(per_sample), 1L)
    expect_identical(per_sample$sample, "sample-a")
    expect_identical(names(per_sample)[[ncol(per_sample)]], "sample")
})

test_that("a three-sample cohort exports each sample and the cohort BEDPE", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    samples <- c("sample-a", "sample-b", "sample-c")
    bam <- add_multi_sample_gc_bam(
        fixture, samples, steps = c(4L, 6L, 8L)
    )
    state <- app_state(fixture)
    per_sample <- stats::setNames(
        lapply(samples, function(sample) tempfile(fileext = ".bedpe")),
        samples
    )
    cohort_target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(c(per_sample, cohort_target)), add = TRUE)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        set_app_inputs(
            session,
            bams = upload_input(bam)
        )
        expect_true(preflight_state()$ok)

        session$setInputs(reference_sample = "")
        suppressWarnings(session$setInputs(run = 1))
        for (sample in samples) {
            expect_match(output$cohort_table, sample)
        }
        expect_match(output$cohort_table, "median")
        expect_match(output_html(output$sample_selector_ui), "result_sample")
        expect_false(grepl("result_sample", output_html(output$cohort_ui)))

        set_explorer_inputs(session)
        for (sample in samples) {
            session$setInputs(
                result_sample = sample,
                `explorer-segments_rows_selected` = NULL
            )
            session$setInputs(`explorer-segments_rows_selected` = 1L)
            expect_identical(selection_map()[[sample]], 1L)
            expect_true(file.copy(output$download_sample, per_sample[[sample]]))
        }
        expect_identical(
            vapply(selection_map(), length, integer(1L)),
            c(`sample-a` = 1L, `sample-b` = 1L, `sample-c` = 1L)
        )
        expect_true(file.copy(output$download_cohort, cohort_target))
    })

    for (sample in samples) {
        out <- utils::read.table(
            per_sample[[sample]],
            header = TRUE, sep = "\t", stringsAsFactors = FALSE
        )
        expect_identical(nrow(out), 1L)
        expect_identical(out$sample, sample)
        expect_identical(names(out)[[ncol(out)]], "sample")
    }

    cohort <- utils::read.table(
        cohort_target,
        header = TRUE, sep = "\t", stringsAsFactors = FALSE
    )
    expect_identical(nrow(cohort), 3L)
    expect_identical(sort(unique(cohort$sample)), samples)
})

test_that("the cohort card shows zero-coverage warnings", {
    skip_if_no_dashboard()
    skip_if_no_tiny_fixtures()
    state <- bam_launcher_state(tiny_config_path())
    expect_true(state$ok)
    fake <- fake_bam_result(
        "A",
        counts = list(
            donor = matrix(0, 60L, 1L, dimnames = list(NULL, "A")),
            recipient = matrix(10, 60L, 1L, dimnames = list(NULL, "A"))
        )
    )

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        result(fake)
        expect_match(output$cohort_table, "A")
        expect_match(
            output_html(output$cohort_ui),
            "no reads on the donor"
        )
    })
})

test_that("a failing run reports the stage and keeps the session alive", {
    skip_if_no_dashboard()
    skip_if_no_tiny_fixtures()
    state <- bam_launcher_state(tiny_config_path())
    expect_true(state$ok)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(
            bams = upload_input(tiny_bam_path()),
            resolution = "100",
            mapq = 20,
            exclude_flags = 2308,
            cores = 1,
            max_fit_windows = 1000,
            reference_sample = ""
        )
        expect_true(preflight_state()$ok)
        expect_null(result())

        ## The committed tiny fixture has a constant GC profile, so the fit
        ## fails inside the correction stage. The handler must catch it.
        suppressWarnings(session$setInputs(run = 1L))
        expect_null(result())
        expect_match(
            output_html(output$cohort_ui),
            "Run the analysis to enable the exports"
        )
        info <- last_error()
        expect_false(is.null(info))
        expect_identical(info$stage, "correction")
        expect_true(length(info$files) > 0L)
        expect_match(
            output$calculation_log,
            "Run failed at stage 'correction'",
            fixed = TRUE
        )
        expect_true(preflight_state()$ok)
        expect_match(
            paste(bam_preflight_log(preflight_state()), collapse = "\n"),
            "Preflight OK"
        )
    })
})

test_that("the run button is disabled when a contig is too short", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture(
        contig_length = 5000L,
        allowed_windows = c(100L, 1000L)
    )
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)
    ## A config loaded at startup always fits the largest allowed resolution;
    ## this simulates a registry whose allowed set changed without a re-read.
    state$config$defaults$allowed_windows <- c(100L, 1000L, 10000L)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        set_app_inputs(
            session,
            bams = upload_input(fixture$bam_path),
            resolution = "10000"
        )
        expect_true(preflight_state()$ok)
        expect_match(output_html(output$run_ui), "shorter than")
        expect_match(output_html(output$run_ui), "disabled")

        session$setInputs(resolution = "1000")
        expect_false(grepl("disabled", output_html(output$run_ui)))
    })
})

test_that("every error category is reported without escaping the server", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    state <- app_state(fixture)

    cases <- list(
        config = quote(.config_error("bad config")),
        BAM = quote(.mamamia_error("MaMaMIABamHeaderError", "bad bam")),
        assembly = quote(.mamamia_error(
            "MaMaMIAAssemblyMatchError", "no match"
        )),
        index = quote(.mamamia_error("MaMaMIAIndexError", "no index")),
        windows = quote(.mamamia_error("MaMaMIAWindowError", "too short")),
        GC = quote(.mamamia_error("MaMaMIAGcrankError", "bad gc")),
        counts = quote(.mamamia_error("MaMaMIACountError", "bad counts")),
        correction = quote(.mamamia_error(
            "MaMaMIASampleError",
            "Sample 'synthetic' (GC correction): too few"
        )),
        normalization = quote(.mamamia_error(
            "MaMaMIAAnalysisError", "Normalization failed: x"
        )),
        segmentation = quote(.mamamia_error(
            "MaMaMIASampleError",
            "Sample 'synthetic' (segmentation): boom"
        )),
        limits = quote(.mamamia_error("MaMaMIALimitError", "too many")),
        analysis = quote(stop("plain failure", call. = FALSE))
    )

    check_case <- function(stage) {
        expression <- cases[[stage]]
        testthat::local_mocked_bindings(
            analyzeBamInput = function(...) {
                error <- tryCatch(
                    eval(expression), error = function(e) e
                )
                stop(error)
            },
            .package = "MaMaMIA"
        )
        shiny::testServer(bam_launcher_server(state, cache = FALSE), {
            set_app_inputs(session, bams = upload_input(fixture$bam_path))
            expect_true(preflight_state()$ok)
            previous <- fake_bam_result("previous")
            result(previous)

            suppressWarnings(session$setInputs(run = 1L))
            info <- last_error()
            expect_identical(info$stage, stage, info = stage)
            expect_true(nzchar(info$suggestion), info = stage)
            expect_true(length(info$files) > 0L, info = stage)
            expect_match(
                output$calculation_log,
                paste0("stage '", stage, "'"),
                fixed = TRUE,
                info = stage
            )
            expect_identical(result(), previous, info = stage)
            expect_true(preflight_state()$ok, info = stage)
        })
    }
    for (stage in names(cases)) {
        check_case(stage)
    }
})

test_that("invalid reference selection is reported before the analysis", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    bam <- add_multi_sample_gc_bam(fixture, c("sample-a", "sample-b"))
    state <- app_state(fixture)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        set_app_inputs(
            session,
            bams = upload_input(bam),
            reference_sample = "nope"
        )
        expect_true(preflight_state()$ok)

        suppressWarnings(session$setInputs(run = 1L))
        info <- last_error()
        expect_identical(info$stage, "normalization")
        expect_match(info$message, "Reference sample 'nope'", fixed = TRUE)
        expect_match(
            output$calculation_log,
            "Run parameters rejected at stage 'normalization'",
            fixed = TRUE
        )
        expect_null(result())
    })
})

test_that("the entrypoint builds the BAM app and validates upload size", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)

    app <- mamamia_app(config = fixture$config_path, cache = FALSE)
    expect_s3_class(app, "shiny.appobj")
    expect_identical(attr(app, "mamamia_input"), "bam")

    expect_error(runMaMaMIA(max_request_size = 0), "positive")
})

test_that("a multi-sample BAM is detected automatically", {
    skip_if_no_dashboard()
    skip_if_no_samtools()
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    bam <- add_multi_sample_gc_bam(
        fixture, c("sample-a", "sample-b"), steps = c(4L, 6L)
    )
    state <- app_state(fixture)

    shiny::testServer(bam_launcher_server(state, cache = FALSE), {
        session$setInputs(bams = upload_input(bam))
        pf <- preflight_state()
        expect_true(pf$ok)
        expect_identical(nrow(pf$table), 2L)
        expect_match(
            output_html(output$bam_mode_note),
            "2 samples detected from @RG SM",
            fixed = TRUE
        )
        expect_match(output_html(output$reference_ui), "Reference sample")
        expect_true(sidebar_unlocked())

        status <- output_html(output$footer_status)
        expect_match(status, "1 BAM file", fixed = TRUE)
        expect_false(grepl("sample-a", status, fixed = TRUE))
        expect_false(grepl("recipient__R", status, fixed = TRUE))
    })
})
