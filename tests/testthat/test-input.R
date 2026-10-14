## Fixture for the on-disk layout: one sample, one donor reference,
## two donor chromosomes (1At, 2At) and two recipient chromosomes (1A, 2A).
write_input_layout <- function(n = 60L) {
    root <- tempfile("mamamia-input")
    analysis <- file.path(root, "analysis")
    assemblies <- file.path(root, "assemblies")
    sample_dir <- file.path(analysis, "sample1")
    reference_dir <- file.path(assemblies, "ref1")
    dir.create(sample_dir, recursive = TRUE)
    dir.create(reference_dir, recursive = TRUE)

    windows <- function(chr, values) {
        start <- seq(0, by = 100, length.out = length(values))
        data.frame(chr, start, start + 100, values)
    }
    write_table <- function(x, path, append = FALSE) {
        utils::write.table(x, path,
            sep = "\t", row.names = FALSE,
            col.names = FALSE, quote = FALSE, append = append
        )
    }
    write_gc <- function(x, path) {
        writeLines("#1_usercol\t2_usercol\t3_usercol\t5_pct_gc", path)
        write_table(x, path, append = TRUE)
    }

    write_table(
        rbind(
            windows("d1", c(rep(8, n / 2), rep(16, n / 2))),
            windows("d2", rep(8, n))
        ),
        file.path(sample_dir, "ref1.w100.dcov.bed")
    )
    write_table(
        rbind(windows("r1", rep(10, n)), windows("r2", rep(10, n))),
        file.path(sample_dir, "ref1.w100.rcov.bed")
    )
    write_gc(
        rbind(
            windows("d1", seq(0.35, 0.55, length.out = n)),
            windows("d2", seq(0.35, 0.55, length.out = n))
        ),
        file.path(reference_dir, "donor.w100.gc.bed")
    )
    write_gc(
        rbind(
            windows("r1", seq(0.35, 0.55, length.out = n)),
            windows("r2", seq(0.35, 0.55, length.out = n))
        ),
        file.path(reference_dir, "recipient.w100.gc.bed")
    )
    ## The last row belongs to no coverage table and must be dropped on read.
    write_table(
        data.frame(
            c("d1", "d2", "r1", "r2", "unused1"),
            c("1At", "2At", "1A", "2A", "7Z"),
            c("At", "At", "A", "A", "Z")
        ),
        file.path(reference_dir, "chromosomes.tsv")
    )

    list(
        root = root,
        analysis = analysis,
        assemblies = assemblies,
        sample_dir = sample_dir,
        reference_dir = reference_dir
    )
}

test_that("readMaMaMIAInput() reads the on-disk layout into RCA-ready tables", {
    layout <- write_input_layout()
    input <- readMaMaMIAInput(layout$analysis, layout$assemblies, "sample1", "ref1")

    expect_named(input, c(
        "dcov", "dgc", "rcov", "rgc", "meta",
        "don_chr_ids", "rec_chr_ids", "sample", "reference", "window"
    ))
    expect_equal(input$sample, "sample1")
    expect_equal(input$reference, "ref1")
    expect_equal(input$window, "100")
    expect_equal(names(input$dcov), c("chr_id", "start", "stop", "cov"))
    expect_equal(names(input$dgc), c("chr_id", "start", "stop", "gc"))
    expect_equal(names(input$rcov), c("chr_id", "start", "stop", "cov"))
    expect_equal(names(input$rgc), c("chr_id", "start", "stop", "gc"))
    expect_equal(nrow(input$dcov), 120L)
    expect_equal(nrow(input$rcov), 120L)
    expect_setequal(input$don_chr_ids, c("d1", "d2"))
    expect_setequal(input$rec_chr_ids, c("r1", "r2"))
    ## chromosomes.tsv keeps only the chromosomes present in the coverage tables
    expect_setequal(input$meta$chr_id, c("d1", "d2", "r1", "r2"))

    rca <- RCA(input$dcov, input$dgc, input$rcov, input$rgc, input$meta)
    expect_s3_class(rca, "RCA")
    expect_equal(length(unique(rca$data$chr_id)), 4L)
})

test_that("readMaMaMIAInput() resolves and validates the windowing scheme", {
    layout <- write_input_layout()
    for (kind in c("dcov", "rcov")) {
        file.copy(
            file.path(layout$sample_dir, sprintf("ref1.w100.%s.bed", kind)),
            file.path(layout$sample_dir, sprintf("ref1.w200.%s.bed", kind))
        )
    }
    file.copy(
        file.path(layout$reference_dir, "donor.w100.gc.bed"),
        file.path(layout$reference_dir, "donor.w200.gc.bed")
    )
    file.copy(
        file.path(layout$reference_dir, "recipient.w100.gc.bed"),
        file.path(layout$reference_dir, "recipient.w200.gc.bed")
    )

    expect_setequal(
        list_window_labels(layout$analysis, layout$assemblies, "sample1", "ref1"),
        c(100L, 200L)
    )
    expect_error(
        readMaMaMIAInput(layout$analysis, layout$assemblies, "sample1", "ref1"),
        "Multiple windowing schemes"
    )
    expect_equal(
        readMaMaMIAInput(layout$analysis, layout$assemblies, "sample1", "ref1", window = 200)$window,
        "200"
    )
    expect_error(
        readMaMaMIAInput(layout$analysis, layout$assemblies, "sample1", "ref1", window = 999),
        "not available"
    )
    ## Coverage without assembly GC files is not offered as a scheme
    file.remove(file.path(layout$reference_dir, "donor.w200.gc.bed"))
    expect_equal(
        list_window_labels(layout$analysis, layout$assemblies, "sample1", "ref1"),
        100L
    )
})

test_that("readMaMaMIAInput() reports unknown samples, references and metadata", {
    layout <- write_input_layout()
    expect_error(
        readMaMaMIAInput(layout$analysis, layout$assemblies, "nope", "ref1"),
        "Sample directory does not exist"
    )
    expect_error(
        readMaMaMIAInput(layout$analysis, layout$assemblies, "sample1", "nope"),
        "Reference directory does not exist"
    )

    meta_path <- file.path(layout$reference_dir, "chromosomes.tsv")
    writeLines("d1\t1At\tAt\nd2\t2At\tAt", meta_path)
    expect_error(
        readMaMaMIAInput(layout$analysis, layout$assemblies, "sample1", "ref1"),
        "chromosomes.tsv.*missing"
    )
})

test_that("read_window_table() and read_chromosome_meta() reject malformed tables", {
    path <- tempfile()
    writeLines(c("a\t0\t10", "b\t10\t20"), path)
    expect_error(read_window_table(path, "cov", "coverage"), "at least 4 columns")

    writeLines(c("a\t0\t10\t1", "b\t10\t20"), path)
    expect_error(read_window_table(path, "cov", "coverage"), "non-numeric or missing")

    writeLines(c("a\t0\t10\tx"), path)
    expect_error(read_window_table(path, "cov", "coverage"), "non-numeric")

    writeLines(character(), path)
    expect_error(read_window_table(path, "cov", "coverage"), "no data rows")

    expect_error(read_window_table(tempfile(), "cov", "coverage"), "Missing coverage file")

    writeLines(c("a\t0\t10\t1", "a\t10\t20\t2"), path)
    expect_silent(read_window_table(path, "cov", "coverage"))
    ## a bare "#" header line is skipped, as written by bedtools nuc
    writeLines(c("#1_usercol\t2_usercol\t3_usercol\t5_pct_gc", "a\t0\t10\t0.4"), path)
    gc <- read_window_table(path, "gc", "GC")
    expect_equal(gc$gc, 0.4)

    writeLines("a\t1A\tA", path)
    expect_equal(names(read_chromosome_meta(path)), c("chr_id", "chr_name", "subgenome"))
    writeLines(c("a\t1A\tA", "a\t1B\tB"), path)
    expect_error(read_chromosome_meta(path), "duplicated chromosome IDs")
})

test_that("discovery helpers list usable entries only", {
    layout <- write_input_layout()
    expect_equal(listMaMaMIASamples(layout$analysis), "sample1")
    expect_equal(listMaMaMIAReferences(layout$assemblies), "ref1")
    expect_equal(listMaMaMIASamples(layout$analysis), "sample1")

    dir.create(file.path(layout$analysis, "empty"))
    dir.create(file.path(layout$assemblies, "no_meta"))
    expect_equal(listMaMaMIASamples(layout$analysis), "sample1")
    expect_equal(listMaMaMIAReferences(layout$assemblies), "ref1")

    expect_error(listMaMaMIASamples(tempfile()), "does not exist")
    expect_error(listMaMaMIAReferences(tempfile()), "does not exist")
    expect_error(listMaMaMIASamples(NULL), "single non-empty string")
})

test_that("defaultTargetPairs() matches chromosome numbers within subgenomes", {
    meta <- make_meta(
        c("d1", "d2", "d3", "r1", "r2", "r3"),
        c("1At", "2At", "3G", "1A", "2A", "1B"),
        c("At", "At", "G", "A", "A", "B")
    )
    pairs <- defaultTargetPairs(
        meta, c("d1", "d2", "d3"), c("r1", "r2", "r3"),
        don_subgenomes = "At", rec_subgenomes = "A"
    )
    expect_equal(pairs, c(d1 = "r1", d2 = "r2"))

    ## A second recipient subgenome makes every number ambiguous
    expect_error(
        defaultTargetPairs(meta, "d1", c("r1", "r2", "r3")),
        "No donor->recipient pairs"
    )
    expect_error(
        defaultTargetPairs(meta, "d1", "r1", don_subgenomes = "Z"),
        "No donor or recipient chromosomes match"
    )
})

test_that("analyze_input_data() runs the pipeline on the on-disk layout", {
    skip_on_cran()
    layout <- write_input_layout()
    input <- readMaMaMIAInput(layout$analysis, layout$assemblies, "sample1", "ref1")

    isa <- suppressWarnings(analyze_input_data(
        input,
        cores = 1L,
        don_subgenomes = "At",
        rec_subgenomes = "A",
        cache = file.path(layout$root, "fit-cache.rds"),
        verbose = FALSE
    ))

    expect_s3_class(isa, "ISA")
    expect_setequal(unique(isa$segments$chr_id.don), c("d1", "d2"))
    expect_setequal(unique(isa$segments$chr_id.rec), c("r1", "r2"))
    expect_true(all(is.finite(isa$segments$seg.median)))
    expect_true(any(isa$segments$putative_introgression))
})

test_that("launcher_app() scans, loads and runs the analysis", {
    skip_on_cran()
    testthat::skip_if_not_installed("shiny")
    testthat::skip_if_not_installed("bslib")
    testthat::skip_if_not_installed("DT")
    testthat::skip_if_not_installed("shinyWidgets")

    layout <- write_input_layout()
    app <- launcher_app(
        layout$analysis, layout$assemblies,
        cache = file.path(layout$root, "app-cache.rds")
    )

    shiny::testServer(app, {
        expect_match(output$input_summary, "Pick an input source")

        session$setInputs(
            source = "files",
            analysis_dir = layout$analysis,
            assemblies_dir = layout$assemblies,
            sample = "sample1",
            reference = "ref1",
            window = "100",
            don_subgenome = "At",
            rec_subgenome = "A"
        )
        session$flushReact()
        expect_match(
            output$input_summary,
            "2 donor-recipient pair\\(s\\) will be segmented \\(At -> A\\)"
        )

        suppressWarnings({
            session$setInputs(cores = 1L, run = 1)
            session$flushReact()
        })
        expect_match(output[["explorer-selection_status"]], "0 segment\\(s\\) selected")
        expect_true(nzchar(as.character(output[["explorer-segments"]])))
    })
})
