## T0.4: the tiny committed fixtures are self-consistent and usable on CI.

test_that("tiny fixtures are present and complete", {
    skip_if_no_tiny_fixtures()

    expect_true(file.exists(tiny_fixture("tiny.fa")))
    expect_gt(file.info(tiny_fixture("tiny.bam"))$size, 0)
    expect_gt(file.info(tiny_fixture("tiny.bam.csi"))$size, 0)
    expect_true(dir.exists(tiny_gcrank_index()))
})

test_that("the tiny assembly config loads and points at the fixture index", {
    skip_if_no_tiny_fixtures()

    config <- readMaMaMIAConfig(tiny_config_path())

    expect_named(config$assemblies, "tiny")
    assembly <- config$assemblies$tiny
    expect_identical(assembly$pairs, c(d1 = "r1"))
    expect_identical(assembly$contigs$donor$id, "d1")
    expect_identical(assembly$contigs$recipient$id, "r1")
    expect_identical(assembly$contigs$donor$length, 20000)
    expect_identical(assembly$contigs$recipient$length, 20000)
    expect_identical(
        assembly$contigs$donor$index_path,
        normalizePath(tiny_gcrank_index(), winslash = "/", mustWork = TRUE)
    )
    expect_true(config$defaults$normalization_by_side)
})

test_that("the tiny FASTA has the expected GC content", {
    skip_if_no_tiny_fixtures()

    lines <- readLines(tiny_fixture("tiny.fa"), warn = FALSE)
    headers <- grep("^>", lines)
    expect_identical(lines[headers], c(">d1", ">r1"))
    sequences <- lapply(seq_along(headers), function(i) {
        start <- headers[[i]] + 1L
        end <- if (i < length(headers)) headers[[i + 1L]] - 1L else length(lines)
        paste(lines[start:end], collapse = "")
    })
    expect_identical(nchar(sequences[[1L]]), 20000L)
    expect_identical(nchar(sequences[[2L]]), 20000L)
    gc_fraction <- function(x) {
        bases <- strsplit(toupper(x), "", fixed = TRUE)[[1L]]
        sum(bases %in% c("G", "C")) / length(bases)
    }
    expect_equal(gc_fraction(sequences[[1L]]), 0.4, tolerance = 1e-12)
    expect_equal(gc_fraction(sequences[[2L]]), 0.6, tolerance = 1e-12)
})

test_that("the tiny gcrank index answers contig-relative queries", {
    skip_if_no_tiny_fixtures()

    index <- MaMaMIA:::gcrank_open_cpp(
        tiny_gcrank_index(), mmap = TRUE, validate_crc = FALSE
    )
    on.exit(MaMaMIA:::gcrank_close_cpp(index), add = TRUE)

    contigs <- MaMaMIA:::gcrank_contigs_cpp(index)
    expect_named(contigs, c("name", "length", "offset"))
    expect_identical(contigs$name, c("d1", "r1"))
    expect_identical(contigs$length, c(20000, 20000))
    expect_identical(contigs$offset, c(0, 20000))

    starts <- seq(0, 19900, by = 100)
    expect_equal(
        MaMaMIA:::gcrank_gc_cpp(index, "d1", starts, starts + 100),
        rep(0.4, length(starts)),
        tolerance = 1e-6
    )
    expect_equal(
        MaMaMIA:::gcrank_gc_cpp(index, "r1", starts, starts + 100),
        rep(0.6, length(starts)),
        tolerance = 1e-6
    )
})

test_that("the tiny BAM header is readable with Rhtslib", {
    skip_if_no_tiny_fixtures()

    header <- MaMaMIA:::bam_header_info_cpp(tiny_bam_path())

    expect_identical(header$format, "BAM")
    expect_identical(header$n_targets, 2L)
    expect_named(header$targets, c("name", "length"))
    expect_identical(header$targets$name, c("d1", "r1"))
    expect_identical(header$targets$length, c(20000, 20000))
    expect_match(header$header_text, "SM:tiny-sample")
    expect_match(header$header_text, "^@HD")
})

test_that("the tiny BAM has exact counts for MAPQ/flag assertions", {
    skip_if_no_tiny_fixtures()

    gr <- GenomicRanges::GRanges(
        c("d1", "r1"),
        IRanges::IRanges(c(1L, 1L), c(20000L, 20000L))
    )
    count <- function(mapq, filtered_flag) {
        as.integer(bamsignals::bamCount(
            tiny_bam_path(), gr,
            mapqual = mapq, filteredFlag = filtered_flag, verbose = FALSE
        ))
    }

    ## bamsignals filters a read only when *all* bits in filteredFlag are set,
    ## unlike samtools -F, which filters when any bit is set. The fixture has
    ## one secondary, one supplementary, one duplicate and one low-MAPQ read so
    ## T3.1 can verify a samtools-compatible countReads() wrapper.
    expect_identical(count(20, -1), c(8L, 5L))   # MAPQ only
    expect_identical(count(20, 1024), c(8L, 4L)) # duplicate
    expect_identical(count(20, 256), c(7L, 5L))  # secondary
    expect_identical(count(20, 2048), c(7L, 5L)) # supplementary
    expect_identical(count(20, 2308), c(8L, 5L)) # no read has all three bits
    expect_identical(count(0, 1024), c(9L, 4L))
})

test_that("a tiny fixture can run the end-to-end BAM pipeline", {
    ## The committed tiny fixture has constant GC (40% and 60%), which cannot
    ## support the ZINB smooth. The variable-GC fixture exercises the same
    ## contract with a GC profile the correction model can fit.
    fixture <- variable_gc_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)

    result <- suppressWarnings(analyzeBamInput(
        fixture$bam_path, fixture$config_path,
        normalization = "none", windows_per_chromosome = 100L,
        cores = 1L
    ))

    expect_named(result, c(
        "assembly", "config", "bams", "samples", "windows", "gc", "counts",
        "rca", "isa", "normalization", "provenance"
    ))
    expect_identical(result$bams$sample, "synthetic")
    expect_identical(result$samples$sample, "synthetic")
    expect_identical(result$windows$n_per_chromosome, 100L)
    expect_named(result$rca, "synthetic")
    expect_named(result$isa, "synthetic")
    expect_identical(result$normalization$mode, "none")
})
