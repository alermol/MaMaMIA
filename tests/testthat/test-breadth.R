## P1: htslib breadth backend parity and edge cases.

breadth_windows <- function(contig_length = 20000L, n = 100L) {
    config <- readMaMaMIAConfig(tiny_config_path())
    makeAnalysisWindows(tiny_assembly(config), tiny_contigs(), n)
}

samtools_depth_coverage <- function(bam, contigs, n, mapq = 20L,
                                    flags = "0x704") {
    skip_if_no_samtools()
    output <- system2(
        "samtools",
        c(
            "depth", "-a", "-Q", as.character(mapq), "-G", flags, bam
        ),
        stdout = TRUE, stderr = FALSE
    )
    if (length(output) == 0L) {
        return(NULL)
    }
    fields <- strsplit(output, "\t", fixed = TRUE)
    data.frame(
        chr = vapply(fields, `[[`, character(1), 1L),
        pos = as.integer(vapply(fields, `[[`, character(1), 2L)),
        depth = as.integer(vapply(fields, `[[`, character(1), 3L)),
        stringsAsFactors = FALSE
    )
}

test_that("countWindowMetrics() matches samtools depth on the tiny fixture", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()
    windows <- breadth_windows()
    header <- bamHeaderInfo(tiny_bam_path())
    observed <- countWindowMetrics(
        tiny_bam_path(), windows, "tiny-sample", header$read_groups
    )
    depth <- samtools_depth_coverage(
        tiny_bam_path(), windows$table, windows$n_per_chromosome
    )
    for (i in seq_len(nrow(windows$table))) {
        row <- windows$table[i, ]
        bases <- depth$depth[
            depth$chr == row$chr_id &
                depth$pos > row$start & depth$pos <= row$stop
        ]
        expect_identical(
            unname(observed$covered[1L, i]),
            as.integer(sum(bases > 0)),
            info = sprintf("%s window %d", row$chr_id, row$iid)
        )
    }
})

test_that("countWindowMetrics() matches bedtools coverage counts", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()
    if (!nzchar(Sys.which("bedtools"))) {
        testthat::skip("bedtools is not available")
    }
    windows <- breadth_windows()
    header <- bamHeaderInfo(tiny_bam_path())
    observed <- countWindowMetrics(
        tiny_bam_path(), windows, "tiny-sample", header$read_groups
    )
    bed <- tempfile(fileext = ".bed")
    filtered <- tempfile(fileext = ".bam")
    reads <- tempfile(fileext = ".bed")
    on.exit(unlink(c(bed, filtered, reads)), add = TRUE)
    utils::write.table(
        data.frame(
            windows$table$chr_id, windows$table$start, windows$table$stop
        ),
        bed, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    status <- system2(
        "samtools",
        c(
            "view", "-b", "-F", "0x704", "-q", "20",
            "-o", filtered, tiny_bam_path()
        ),
        stdout = FALSE, stderr = FALSE
    )
    skip_if_not(identical(status, 0L), "samtools could not filter the BAM")
    status <- system2(
        "bedtools", c("bamtobed", "-i", filtered),
        stdout = reads, stderr = FALSE
    )
    skip_if_not(identical(status, 0L), "bedtools bamtobed failed")
    output <- system2(
        "bedtools",
        c("coverage", "-a", bed, "-b", reads, "-counts"),
        stdout = TRUE, stderr = FALSE
    )
    counts <- as.integer(sub("^.*\t", "", output))
    expect_identical(unname(as.integer(observed$counts[1L, ])), counts)
})

test_that("a boundary read contributes to every window it overlaps", {
    skip_if_no_samtools()
    directory <- tempfile("breadth-boundary-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- file.path(directory, "boundary.bam")
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:200",
            "inside\t0\tchr1\t11\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            ## 0-based 90..119 overlaps windows 9 and 10 (0-based 90-99,
            ## 100-109) when the contig is split into 20 windows of 10 bp.
            "boundary\t0\tchr1\t91\t60\t30M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
        ),
        bam,
        index = TRUE
    )
    table <- data.frame(
        chr_id = rep("chr1", 20L),
        start = seq(0L, 190L, 10L),
        stop = seq(10L, 200L, 10L),
        role = "donor",
        chr_name = "1D",
        component = "D",
        group = "1",
        copy_number = 2,
        stringsAsFactors = FALSE
    )
    windows <- list(
        gr = GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(seq(1L, 191L, 10L), seq(10L, 200L, 10L)),
            seqlengths = c(chr1 = 200)
        ),
        table = table,
        n_per_chromosome = 20L
    )
    observed <- countWindowMetrics(
        bam, windows, "s",
        data.frame(id = "rg", sample = "s"), mapq = 0L, exclude_flags = 0L
    )
    expect_identical(as.integer(observed$counts[1L, 1:12]),
        c(0L, 1L, rep(0L, 7L), 1L, 1L, 1L))
    expect_identical(as.integer(observed$covered[1L, 9:12]),
        c(0L, 10L, 10L, 10L))
})

test_that("countWindowMetrics() separates samples by RG", {
    skip_if_no_samtools()
    directory <- tempfile("breadth-rg-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- file.path(directory, "two-samples.bam")
    write_test_bam(
        multi_sample_variable_gc_sam(c("a", "b"), contig_length = 20000L, steps = c(4L, 7L)),
        bam,
        index = TRUE
    )
    header <- bamHeaderInfo(bam)
    config <- readMaMaMIAConfig(tiny_config_path())
    contigs <- data.frame(name = c("d1", "r1"), length = c(20000, 20000))
    windows <- makeAnalysisWindows(
        tiny_assembly(config), contigs, n_per_chromosome = 20L
    )
    observed <- countWindowMetrics(
        bam, windows, c("a", "b"), header$read_groups
    )
    expect_identical(dim(observed$counts), c(2L, 40L))
    expect_identical(rownames(observed$counts), c("a", "b"))
    ## Sample b has reads at the same positions but a different step; the
    ## totals must differ and both must be positive.
    expect_gt(sum(observed$counts["a", ]), 0)
    expect_gt(sum(observed$counts["b", ]), 0)
    expect_false(identical(
        as.integer(observed$counts["a", ]),
        as.integer(observed$counts["b", ])
    ))
})

test_that("flag and MAPQ filters are applied", {
    skip_if_no_samtools()
    directory <- tempfile("breadth-flags-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- file.path(directory, "flags.bam")
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:100",
            "keep\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "lowmq\t0\tchr1\t21\t5\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "secondary\t256\tchr1\t31\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "duplicate\t1024\tchr1\t41\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            "supplementary\t2048\tchr1\t51\t60\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII"
        ),
        bam,
        index = TRUE
    )
    table <- data.frame(
        chr_id = "chr1", start = 0L, stop = 100L, role = "donor",
        chr_name = "1D", component = "D", group = "1", copy_number = 2,
        stringsAsFactors = FALSE
    )
    windows <- list(
        gr = GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(1L, 100L), seqlengths = c(chr1 = 100)
        ),
        table = table, n_per_chromosome = 1L
    )
    observed <- countWindowMetrics(
        bam, windows, "s",
        data.frame(id = "rg", sample = "s"), mapq = 20L, exclude_flags = 1796L
    )
    ## keep + supplementary (1796 does not exclude 2048).
    expect_identical(unname(as.integer(observed$counts[1L, 1L])), 2L)
    expect_identical(unname(as.integer(observed$covered[1L, 1L])), 20L)
})

test_that("CIGAR deletions, skips and clips are handled", {
    skip_if_no_samtools()
    directory <- tempfile("breadth-cigar-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- file.path(directory, "cigar.bam")
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:100",
            ## 5M2D5M covers 12 reference bases.
            "del\t0\tchr1\t1\t60\t5M2D5M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            ## 5M10N5M covers 10 bases in two blocks.
            "skip\t0\tchr1\t21\t60\t5M10N5M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII",
            ## Soft clips do not cover reference bases.
            "clip\t0\tchr1\t41\t60\t5S10M5S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\tIIIIIIIIIIIIIIIIIIII"
        ),
        bam,
        index = TRUE
    )
    table <- data.frame(
        chr_id = "chr1", start = 0L, stop = 100L, role = "donor",
        chr_name = "1D", component = "D", group = "1", copy_number = 2,
        stringsAsFactors = FALSE
    )
    windows <- list(
        gr = GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(1L, 100L), seqlengths = c(chr1 = 100)
        ),
        table = table, n_per_chromosome = 1L
    )
    with_deletions <- countWindowMetrics(
        bam, windows, "s", data.frame(id = "rg", sample = "s"),
        mapq = 0L, exclude_flags = 0L, count_deletions = TRUE
    )
    without_deletions <- countWindowMetrics(
        bam, windows, "s", data.frame(id = "rg", sample = "s"),
        mapq = 0L, exclude_flags = 0L, count_deletions = FALSE
    )
    expect_identical(unname(as.integer(with_deletions$covered[1L, 1L])), 32L)
    expect_identical(unname(as.integer(without_deletions$covered[1L, 1L])), 30L)
    expect_identical(unname(as.integer(with_deletions$counts[1L, 1L])), 3L)
    expect_identical(with_deletions$read_length, 10)
})

test_that("read length is inferred and empty BAMs return NA", {
    skip_if_no_samtools()
    directory <- tempfile("breadth-empty-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- file.path(directory, "empty.bam")
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:100"
        ),
        bam,
        index = TRUE
    )
    table <- data.frame(
        chr_id = "chr1", start = 0L, stop = 100L, role = "donor",
        chr_name = "1D", component = "D", group = "1", copy_number = 2,
        stringsAsFactors = FALSE
    )
    windows <- list(
        gr = GenomicRanges::GRanges(
            "chr1", IRanges::IRanges(1L, 100L), seqlengths = c(chr1 = 100)
        ),
        table = table, n_per_chromosome = 1L
    )
    observed <- countWindowMetrics(
        bam, windows, "s", data.frame(id = "rg", sample = "s"),
        mapq = 0L, exclude_flags = 0L
    )
    expect_identical(unname(as.integer(observed$counts[1L, 1L])), 0L)
    expect_identical(unname(as.integer(observed$covered[1L, 1L])), 0L)
    expect_true(is.na(observed$read_length))
})

test_that("countWindowMetrics() aligns per-sample metrics with the matrices", {
    skip_if_no_samtools()
    directory <- tempfile("breadth-align-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- file.path(directory, "two-samples.bam")
    write_test_bam(
        multi_sample_variable_gc_sam(
            c("a", "b"), contig_length = 20000L, steps = c(4L, 7L)
        ),
        bam,
        index = TRUE
    )
    header <- bamHeaderInfo(bam)
    config <- readMaMaMIAConfig(tiny_config_path())
    contigs <- data.frame(name = c("d1", "r1"), length = c(20000, 20000))
    windows <- makeAnalysisWindows(
        tiny_assembly(config), contigs, n_per_chromosome = 20L
    )
    observed <- countWindowMetrics(
        bam, windows, c("a", "b"), header$read_groups
    )
    metrics <- observed$metrics
    for (sample in c("a", "b")) {
        rows <- metrics[metrics$sample == sample, , drop = FALSE]
        expect_identical(
            rows$count,
            unname(as.integer(observed$counts[sample, ])),
            info = paste("count for", sample)
        )
        expect_identical(
            rows$covered,
            unname(as.integer(observed$covered[sample, ])),
            info = paste("covered for", sample)
        )
    }
})
