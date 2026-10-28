## T2.1 and T2.2: contigs.tsv parsing and exact-N window generation.

synthetic_assembly <- function(donor_length = 20000, recipient_length = donor_length,
                               key = "synthetic") {
    list(
        key = key,
        label = "Synthetic assembly",
        description = "",
        gc_indexes = character(0),
        contigs = list(
            recipient = data.frame(
                id = "r1", name = "1R", component = "R",
                length = recipient_length, index_path = NA_character_,
                stringsAsFactors = FALSE
            ),
            donor = data.frame(
                id = "d1", name = "1D", component = "D",
                length = donor_length, index_path = NA_character_,
                stringsAsFactors = FALSE
            )
        ),
        pairs = c(d1 = "r1")
    )
}

tiny_contigs <- function() {
    readContigsTsv(file.path(tiny_gcrank_index(), "contigs.tsv"))
}

expect_exact_tiling <- function(table, chr_id, length_bp, n) {
    idx <- which(table$chr_id == chr_id)
    expect_identical(length(idx), n)
    start <- table$start[idx]
    stop <- table$stop[idx]
    expect_identical(start[[1L]], 0L)
    expect_identical(stop[[n]], as.integer(length_bp))
    expect_identical(start[-1L], stop[-n])
    expect_true(all(stop > start))
    expect_identical(
        start,
        as.integer(floor((seq_len(n) - 1) * length_bp / n))
    )
    expect_identical(
        stop,
        as.integer(floor(seq_len(n) * length_bp / n))
    )
}

test_that("readContigsTsv() reads the tiny gcrank contig table", {
    skip_if_no_tiny_fixtures()

    contigs <- tiny_contigs()
    expect_named(contigs, c("name", "length", "offset"))
    expect_identical(contigs$name, c("d1", "r1"))
    expect_identical(contigs$length, c(20000, 20000))
    expect_identical(contigs$offset, c(0, 20000))
})

test_that("readContigsTsv() rejects malformed contig tables", {
    write_table <- function(lines) {
        path <- tempfile(fileext = ".tsv")
        writeLines(lines, path)
        path
    }

    duplicate <- write_table(c("d1\t100\t0", "d1\t100\t100"))
    expect_error(
        readContigsTsv(duplicate),
        "Duplicate contig 'd1'",
        class = "MaMaMIAWindowError"
    )

    non_numeric <- write_table(c("d1\tnot-a-length\t0"))
    expect_error(
        readContigsTsv(non_numeric),
        "length.*not numeric",
        class = "MaMaMIAWindowError"
    )

    wrong_columns <- write_table(c("d1\t100"))
    expect_error(
        readContigsTsv(wrong_columns),
        "expected exactly 3",
        class = "MaMaMIAWindowError"
    )

    non_positive <- write_table(c("d1\t0\t0"))
    expect_error(
        readContigsTsv(non_positive),
        "positive integers",
        class = "MaMaMIAWindowError"
    )

    non_monotone <- write_table(c("d1\t100\t0", "d2\t100\t50"))
    expect_error(
        readContigsTsv(non_monotone),
        "cumulative lengths of the preceding contigs",
        class = "MaMaMIAWindowError"
    )

    empty <- tempfile(fileext = ".tsv")
    file.create(empty)
    expect_error(
        readContigsTsv(empty),
        "No contigs found",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        readContigsTsv(tempfile(fileext = ".tsv")),
        "does not exist",
        class = "MaMaMIAWindowError"
    )
})

test_that("readContigsTsv() reads the real wheat and rye metadata", {
    root <- skip_if_no_test_data()
    wheat <- file.path(
        root, "gcidx", "triticum-aestivum", "index.gcrank", "contigs.tsv"
    )
    rye <- file.path(
        root, "gcidx", "secale-cereale", "index.gcrank", "contigs.tsv"
    )
    if (!file.exists(wheat) || !file.exists(rye)) {
        testthat::skip("real gcrank contig tables are missing")
    }

    wheat_contigs <- readContigsTsv(wheat)
    rye_contigs <- readContigsTsv(rye)
    expect_identical(nrow(wheat_contigs), 21L)
    expect_identical(nrow(rye_contigs), 242L)
    expect_identical(
        wheat_contigs$length[wheat_contigs$name == "CP155614.1"],
        603818960
    )
    expect_identical(
        wheat_contigs$offset[[1L]], 0
    )
})

test_that("makeAnalysisWindows() builds both representations in target order", {
    skip_if_no_tiny_fixtures()
    config <- readMaMaMIAConfig(tiny_config_path())
    windows <- makeAnalysisWindows(
        tiny_assembly(config), tiny_contigs(), 1000L
    )

    expect_named(windows, c("gr", "table", "n_per_chromosome"))
    expect_identical(windows$n_per_chromosome, 1000L)
    table <- windows$table
    expect_named(table, c(
        "chr_id", "start", "stop", "iid", "role", "chr_name", "component",
        "group", "copy_number"
    ))
    expect_identical(nrow(table), 2000L)
    expect_identical(unique(table$chr_id), c("r1", "d1"))
    expect_identical(unique(table$role), c("recipient", "donor"))
    expect_identical(unique(table$chr_name), c("1R", "1D"))
    expect_identical(unique(table$component), c("R", "D"))
    expect_exact_tiling(table, "r1", 20000, 1000L)
    expect_exact_tiling(table, "d1", 20000, 1000L)

    gr <- windows$gr
    expect_s4_class(gr, "GRanges")
    expect_identical(length(gr), nrow(table))
    expect_identical(
        as.character(GenomicRanges::seqnames(gr)),
        table$chr_id
    )
    expect_identical(IRanges::start(gr), table$start + 1L)
    expect_identical(IRanges::end(gr), table$stop)
    expect_identical(IRanges::width(gr), table$stop - table$start)
    info <- GenomicRanges::seqinfo(gr)
    expect_identical(
        as.numeric(methods::slot(info, "seqlengths")),
        c(20000, 20000)
    )
})

test_that("makeAnalysisWindows() supports 100, 1000 and 10000 windows", {
    skip_if_no_tiny_fixtures()
    config <- readMaMaMIAConfig(tiny_config_path())

    for (n in c(100L, 1000L, 10000L)) {
        windows <- makeAnalysisWindows(
            tiny_assembly(config), tiny_contigs(), n
        )
        expect_identical(nrow(windows$table), 2L * n)
        expect_exact_tiling(windows$table, "r1", 20000, n)
        expect_exact_tiling(windows$table, "d1", 20000, n)
    }
})

test_that("makeAnalysisWindows() handles lengths not divisible by N", {
    assembly <- synthetic_assembly(donor_length = 10000)
    contigs <- data.frame(
        name = c("d1", "r1"), length = c(10000, 10000),
        stringsAsFactors = FALSE
    )
    windows <- makeAnalysisWindows(assembly, contigs, 333L)
    expect_identical(nrow(windows$table), 666L)
    expect_exact_tiling(windows$table, "r1", 10000, 333L)
    expect_exact_tiling(windows$table, "d1", 10000, 333L)
})

test_that("makeAnalysisWindows() rejects a contig shorter than N", {
    assembly <- synthetic_assembly(donor_length = 500)
    contigs <- data.frame(
        name = c("d1", "r1"), length = c(500, 500),
        stringsAsFactors = FALSE
    )

    expect_error(
        makeAnalysisWindows(assembly, contigs, 1000L),
        "contig 'r1' is 500 bp, shorter than the requested 1000 windows",
        class = "MaMaMIAWindowError"
    )
})

test_that("makeAnalysisWindows() validates the contig lookup", {
    assembly <- synthetic_assembly()

    missing <- data.frame(
        name = "d1", length = 20000, stringsAsFactors = FALSE
    )
    expect_error(
        makeAnalysisWindows(assembly, missing, 100L),
        "target contig\\(s\\) not found.*r1",
        class = "MaMaMIAWindowError"
    )

    mismatch <- data.frame(
        name = c("d1", "r1"), length = c(19999, 20000),
        stringsAsFactors = FALSE
    )
    expect_error(
        makeAnalysisWindows(assembly, mismatch, 100L),
        "contig 'd1' is 19999 bp.*20000 bp",
        class = "MaMaMIAWindowError"
    )

    duplicate <- data.frame(
        name = c("d1", "d1", "r1"), length = c(20000, 20000, 20000),
        stringsAsFactors = FALSE
    )
    expect_error(
        makeAnalysisWindows(assembly, duplicate, 100L),
        "duplicate contig 'd1'",
        class = "MaMaMIAWindowError"
    )

    expect_error(
        makeAnalysisWindows(assembly, tempfile(fileext = ".tsv"), 100L),
        "does not exist",
        class = "MaMaMIAWindowError"
    )
})

test_that("makeAnalysisWindows() validates N and accepts a contigs path", {
    skip_if_no_tiny_fixtures()
    config <- readMaMaMIAConfig(tiny_config_path())
    path <- file.path(tiny_gcrank_index(), "contigs.tsv")

    expect_error(
        makeAnalysisWindows(tiny_assembly(config), path, 0L),
        "n_per_chromosome",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        makeAnalysisWindows(tiny_assembly(config), path, 2.5),
        "n_per_chromosome",
        class = "MaMaMIAWindowError"
    )
    windows <- makeAnalysisWindows(tiny_assembly(config), path, 100L)
    expect_identical(nrow(windows$table), 200L)
})

test_that("makeAnalysisWindows() tiles random contig lengths exactly", {
    set.seed(42)
    resolutions <- c(100L, 1000L, 10000L)
    for (i in seq_len(30L)) {
        length_bp <- sample(10000:60000, 1L)
        n <- sample(resolutions[resolutions <= length_bp], 1L)
        assembly <- synthetic_assembly(donor_length = length_bp)
        contigs <- data.frame(
            name = c("d1", "r1"), length = c(length_bp, length_bp),
            stringsAsFactors = FALSE
        )
        windows <- makeAnalysisWindows(assembly, contigs, n)
        expect_identical(nrow(windows$table), 2L * n)
        expect_exact_tiling(windows$table, "r1", length_bp, n)
        expect_exact_tiling(windows$table, "d1", length_bp, n)
        expect_identical(length(windows$gr), 2L * n)
    }
})
