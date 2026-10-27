## T2.3: gcrank R wrapper, cache and real-data numerical regression.

read_gc_bed <- function(path, n = 3L) {
    table <- utils::read.delim(
        path,
        header = FALSE,
        comment.char = "#",
        col.names = c("chr_id", "start", "stop", "gc"),
        colClasses = c("character", "numeric", "numeric", "numeric"),
        stringsAsFactors = FALSE
    )
    utils::head(table, n)
}

copy_gcrank_index <- function() {
    target <- tempfile("gcrank-copy-")
    dir.create(target)
    files <- list.files(tiny_gcrank_index(), full.names = TRUE)
    copied <- file.copy(files, target, recursive = TRUE)
    if (!all(copied)) {
        stop("could not copy the tiny gcrank index")
    }
    target
}

test_that("gcrankContigs() reads the tiny index metadata", {
    skip_if_no_tiny_fixtures()
    on.exit(gcrankClearCache(), add = TRUE)

    contigs <- gcrankContigs(tiny_gcrank_index())
    expect_named(contigs, c("name", "length", "offset"))
    expect_identical(contigs$name, c("d1", "r1"))
    expect_identical(contigs$length, c(20000, 20000))
    expect_identical(contigs$offset, c(0, 20000))
})

test_that("gcrankGc() is vectorized over contigs and windows", {
    skip_if_no_tiny_fixtures()
    on.exit(gcrankClearCache(), add = TRUE)
    index <- tiny_gcrank_index()

    start <- c(0, 100, 19900)
    expect_equal(
        gcrankGc(index, "d1", start, start + 100),
        rep(0.4, 3L),
        tolerance = 1e-6
    )
    expect_equal(
        gcrankGc(index, "r1", start, start + 100),
        rep(0.6, 3L),
        tolerance = 1e-6
    )
    expect_equal(
        gcrankGc(index, c("d1", "r1", "d1"), start, start + 100),
        c(0.4, 0.6, 0.4),
        tolerance = 1e-6
    )
    expect_equal(
        gcrankGc(index, c("d1", "r1"), 0, 100),
        c(0.4, 0.6),
        tolerance = 1e-6
    )
    expect_equal(
        gcrankGc(index, "d1", 0, 20000, remove_ambiguous = TRUE),
        0.4,
        tolerance = 1e-6
    )
    expect_identical(gcrankGc(index, character(), numeric(), numeric()),
        numeric(0)
    )
})

test_that("gcrankGc() validates its inputs", {
    skip_if_no_tiny_fixtures()
    on.exit(gcrankClearCache(), add = TRUE)
    index <- tiny_gcrank_index()

    expect_error(
        gcrankGc(index, "missing", 0, 100),
        "gcrank query failed for missing:0-100",
        class = "MaMaMIAGcrankError"
    )
    expect_error(
        gcrankGc(index, "d1", 19999, 20001),
        "gcrank query failed for d1:19999-20001",
        class = "MaMaMIAGcrankError"
    )
    expect_error(
        gcrankGc(index, "d1", 100, 100),
        "end must be greater than start",
        class = "MaMaMIAGcrankError"
    )
    expect_error(
        gcrankGc(index, "d1", 0.5, 100),
        "whole numbers",
        class = "MaMaMIAGcrankError"
    )
    expect_error(
        gcrankGc(index, "d1", 0, Inf),
        "finite",
        class = "MaMaMIAGcrankError"
    )
    expect_error(
        gcrankGc(index, NA_character_, 0, 100),
        "without NA",
        class = "MaMaMIAGcrankError"
    )
    expect_error(
        gcrankGc(index, c("d1", "r1"), c(0, 0, 0), 100),
        "length 1 or 3",
        class = "MaMaMIAGcrankError"
    )
    expect_error(
        gcrankGc(index, "d1", 0, 100, remove_ambiguous = NA),
        "true or false",
        class = "MaMaMIAGcrankError"
    )
})

test_that("gcrankGc() rejects missing and non-index directories", {
    expect_error(
        gcrankGc(tempfile("missing-"), "d1", 0, 100),
        "does not exist",
        class = "MaMaMIAGcrankError"
    )
    empty <- tempfile("empty-index-")
    dir.create(empty)
    on.exit(unlink(empty, recursive = TRUE), add = TRUE)
    expect_error(
        gcrankGc(empty, "d1", 0, 100),
        "not a gcrank index",
        class = "MaMaMIAGcrankError"
    )
})

test_that("the gcrank cache reuses handles and keeps at most two indexes", {
    skip_if_no_tiny_fixtures()
    on.exit(gcrankClearCache(), add = TRUE)
    gcrankClearCache()
    index <- tiny_gcrank_index()

    gcrankContigs(index)
    expect_identical(
        length(ls(MaMaMIA:::.MaMaMIAgcrankCache)), 1L
    )
    gcrankGc(index, "d1", 0, 100)
    expect_identical(
        length(ls(MaMaMIA:::.MaMaMIAgcrankCache)), 1L
    )

    first <- copy_gcrank_index()
    second <- copy_gcrank_index()
    on.exit(unlink(c(first, second), recursive = TRUE), add = TRUE)
    gcrankContigs(first)
    gcrankContigs(second)
    expect_identical(
        length(ls(MaMaMIA:::.MaMaMIAgcrankCache)), 2L
    )
    expect_false(exists(
        normalizePath(index, winslash = "/", mustWork = TRUE),
        envir = MaMaMIA:::.MaMaMIAgcrankCache,
        inherits = FALSE
    ))
    expect_identical(
        gcrankGc(second, "r1", 0, 100),
        0.6,
        tolerance = 1e-6
    )
})

test_that("real wheat and rye GC values match the reference BED tables", {
    root <- skip_if_no_test_data()
    on.exit(gcrankClearCache(), add = TRUE)

    check <- function(index_dir, bed_path) {
        index_dir <- file.path(root, "gcidx", index_dir, "index.gcrank")
        bed_path <- file.path(root, "assemblies", "cs_t2t_se_cereale", bed_path)
        if (!dir.exists(index_dir) || !file.exists(bed_path)) {
            testthat::skip("real gcrank indexes or reference GC tables are missing")
        }
        expected <- read_gc_bed(bed_path, 3L)
        gc <- gcrankGc(
            index_dir, expected$chr_id, expected$start, expected$stop
        )
        expect_identical(sprintf("%.6f", gc), sprintf("%.6f", expected$gc))
    }

    check("triticum-aestivum", "recipient.w1000.gc.bed")
    check("secale-cereale", "donor.w1000.gc.bed")

    wheat <- file.path(root, "gcidx", "triticum-aestivum", "index.gcrank")
    contigs <- gcrankContigs(wheat)
    expect_identical(nrow(contigs), 21L)
    expect_identical(
        contigs$length[contigs$name == "CP155614.1"],
        603818960
    )
})
