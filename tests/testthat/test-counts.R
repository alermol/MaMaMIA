## T3.1: samtools-compatible counting through bamsignals.

count_test_bam <- function() {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("counts-bam-")
    dir.create(directory)
    path <- file.path(directory, "tiny.bam")
    if (!file.copy(tiny_bam_path(), path)) {
        stop("could not copy the tiny BAM")
    }
    list(path = path, directory = directory)
}

tiny_count_windows <- function() {
    skip_if_no_tiny_fixtures()
    config <- readMaMaMIAConfig(tiny_config_path())
    windows <- makeAnalysisWindows(config$assemblies$tiny, tiny_contigs(), 1L)
    ## makeAnalysisWindows() orders recipient (r1) first, then donor (d1).
    expect_identical(windows$table$chr_id, c("r1", "d1"))
    windows$gr
}

test_that("countReads() matches samtools -F semantics on the tiny fixture", {
    fixture <- count_test_bam()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    gr <- tiny_count_windows()

    ## Window order is r1, d1. The fixture has four normal and one duplicate
    ## read on r1, and six normal, one low-MAPQ, one secondary and one
    ## supplementary read on d1.
    expect_identical(countReads(fixture$path, gr, exclude_flags = 0L),
        c(5L, 8L))
    expect_identical(countReads(fixture$path, gr, exclude_flags = 1024L),
        c(4L, 8L))
    expect_identical(countReads(fixture$path, gr, exclude_flags = 256L),
        c(5L, 7L))
    expect_identical(countReads(fixture$path, gr, exclude_flags = 2048L),
        c(5L, 7L))
    expect_identical(countReads(fixture$path, gr, exclude_flags = 4L),
        c(5L, 8L))

    ## Correct samtools semantics remove the secondary and supplementary
    ## reads. A direct filteredFlag = 2308 call would keep both (c(5, 8)),
    ## because bamsignals filters only when all bits are set.
    expect_identical(countReads(fixture$path, gr, exclude_flags = 2308L),
        c(5L, 6L))
    expect_identical(
        countReads(fixture$path, gr, mapq = 0L, exclude_flags = 1024L),
        c(4L, 9L)
    )
    expect_identical(
        countReads(fixture$path, gr, mapq = 255L, exclude_flags = 0L),
        c(0L, 0L)
    )
})

test_that("countReads() preserves GRanges order and creates the CSI index", {
    fixture <- count_test_bam()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    index <- paste0(fixture$path, ".csi")
    expect_false(file.exists(index))

    gr <- GenomicRanges::GRanges(
        c("d1", "r1"),
        IRanges::IRanges(c(1L, 1L), c(20000L, 20000L)),
        seqlengths = c(d1 = 20000, r1 = 20000)
    )
    counts <- countReads(fixture$path, gr, exclude_flags = 2308L)
    expect_identical(counts, c(6L, 5L))
    expect_true(file.exists(index))
    expect_true(MaMaMIA:::.bam_index_magic(index))
})

test_that("countReads() assigns forward reads by their 5' end", {
    skip_if_no_samtools()
    directory <- tempfile("boundary-bam-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    path <- file.path(directory, "boundary.bam")
    sequence <- strrep("A", 100L)
    quality <- strrep("I", 100L)
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:3000",
            sprintf("crossing\t0\tchr1\t960\t60\t100M\t*\t0\t0\t%s\t%s",
                sequence, quality),
            sprintf("inside\t0\tchr1\t1001\t60\t100M\t*\t0\t0\t%s\t%s",
                sequence, quality)
        ),
        path
    )

    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(c(1L, 1001L), c(1000L, 2000L)),
        seqlengths = c(chr1 = 3000)
    )
    ## The forward read starting at 960 crosses into the second window but its
    ## 5' end is the alignment start, so it belongs to the first window.
    expect_identical(countReads(path, gr, exclude_flags = 2308L), c(1L, 1L))
})

test_that("countReads() assigns reverse-strand reads by their 5' end", {
    skip_if_no_samtools()
    directory <- tempfile("fiveprime-bam-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    path <- file.path(directory, "fiveprime.bam")
    sequence <- strrep("A", 100L)
    quality <- strrep("I", 100L)
    make_read <- function(name, flag, position) {
        sprintf(
            "%s\t%d\tchr1\t%d\t60\t100M\t*\t0\t0\t%s\t%s",
            name, flag, position, sequence, quality
        )
    }
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:chr1\tLN:3000",
            ## Forward read starting in window 1: 5' end 960 -> window 1.
            make_read("forward_cross", 0L, 960L),
            ## Reverse read starting in window 1, ending at 1069: 5' end ->
            ## window 2. bedtools would count it in both windows.
            make_read("reverse_cross", 16L, 970L),
            ## Reverse read fully inside window 2: 5' end 1599 -> window 2.
            make_read("reverse_inside", 16L, 1500L),
            ## Forward read inside window 3.
            make_read("forward_inside", 0L, 2100L)
        ),
        path
    )

    gr <- GenomicRanges::GRanges(
        "chr1",
        IRanges::IRanges(
            c(1L, 1001L, 2001L),
            c(1000L, 2000L, 3000L)
        ),
        seqlengths = c(chr1 = 3000)
    )
    expect_identical(countReads(path, gr, exclude_flags = 0L), c(1L, 2L, 1L))
})

test_that("countReads() supports optional chunking without changing counts", {
    fixture <- count_test_bam()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    gr <- tiny_count_windows()

    old_options <- options(MaMaMIA.count_chunk_size = 1L)
    on.exit(options(old_options), add = TRUE)
    expect_identical(
        countReads(fixture$path, gr, exclude_flags = 2308L),
        c(5L, 6L)
    )
})

test_that("countReads() returns integer(0) for empty GRanges", {
    fixture <- count_test_bam()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    empty <- GenomicRanges::GRanges()
    expect_identical(countReads(fixture$path, empty), integer(0))
})

test_that("countReads() validates its inputs", {
    fixture <- count_test_bam()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    gr <- tiny_count_windows()

    expect_error(
        countReads(tempfile(fileext = ".bam"), gr),
        "does not exist",
        class = "MaMaMIABamHeaderError"
    )
    expect_error(
        countReads(fixture$path, data.frame(chr = "r1")),
        "GRanges",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(
            fixture$path,
            GenomicRanges::GRanges(
                "r1", IRanges::IRanges(1L, 1000L)
            )
        ),
        "positive seqlength",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(
            fixture$path,
            suppressWarnings(GenomicRanges::GRanges(
                "r1", IRanges::IRanges(1L, 1000L), seqlengths = c(r1 = 500)
            ))
        ),
        "exceeds the seqlength",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(
            fixture$path,
            GenomicRanges::GRanges(
                "r1", IRanges::IRanges(10L, 9L), seqlengths = c(r1 = 1000)
            )
        ),
        "empty window",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(fixture$path, gr, mapq = -1L),
        "mapq",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(fixture$path, gr, mapq = 256L),
        "mapq",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(fixture$path, gr, exclude_flags = -1L),
        "exclude_flags",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(fixture$path, gr, exclude_flags = 4096L),
        "exclude_flags",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(fixture$path, gr, exclude_flags = 4095L),
        "supports at most 6 bits",
        class = "MaMaMIACountError"
    )
    expect_error(
        countReads(fixture$path, gr, verbose = NA),
        "verbose",
        class = "MaMaMIACountError"
    )
})

test_that("countReads() refuses a corrupt index in a read-only directory", {
    fixture <- count_test_bam()
    if (.Platform$OS.type == "windows") {
        testthat::skip("POSIX permission bits are required")
    }
    on.exit({
        Sys.chmod(fixture$directory, mode = "0755")
        unlink(fixture$directory, recursive = TRUE)
    }, add = TRUE)
    writeBin(charToRaw("not a csi"), paste0(fixture$path, ".csi"))
    Sys.chmod(fixture$directory, mode = "0555")

    expect_error(
        countReads(fixture$path, tiny_count_windows()),
        "not writable.*stage the BAM in session temp",
        class = "MaMaMIAIndexError"
    )
})

test_that("countReads() matches legacy tables within boundary effects", {
    root <- skip_if_no_test_data()
    bam_path <- file.path(root, "bam_files", "Sq54", "cs_t2t_se_cereale.bam")
    legacy_dir <- file.path(root, "mamamia_analysis", "arsenal", "Sq54")
    recipient_bed <- file.path(
        legacy_dir, "cs_t2t_se_cereale.w1000.rcov.bed"
    )
    donor_bed <- file.path(
        legacy_dir, "cs_t2t_se_cereale.w1000.dcov.bed"
    )
    if (!all(file.exists(c(bam_path, recipient_bed, donor_bed)))) {
        testthat::skip("Sq54 BAM or legacy count tables are missing")
    }

    read_counts <- function(path) {
        utils::read.delim(
            path,
            header = FALSE,
            col.names = c("chr_id", "start", "stop", "count"),
            colClasses = c("character", "numeric", "numeric", "numeric"),
            stringsAsFactors = FALSE
        )
    }
    recipient <- read_counts(recipient_bed)
    donor <- read_counts(donor_bed)
    legacy <- rbind(recipient, donor)

    ## Contig lengths come from the gcrank index metadata, not from the BED.
    wheat <- readContigsTsv(file.path(
        root, "gcidx", "triticum-aestivum", "index.gcrank", "contigs.tsv"
    ))
    rye <- readContigsTsv(file.path(
        root, "gcidx", "secale-cereale", "index.gcrank", "contigs.tsv"
    ))
    contig_lengths <- stats::setNames(
        c(wheat$length, rye$length), c(wheat$name, rye$name)
    )
    missing <- setdiff(unique(legacy$chr_id), names(contig_lengths))
    expect_length(missing, 0L)
    seqlengths <- contig_lengths[unique(legacy$chr_id)]

    gr <- GenomicRanges::GRanges(
        legacy$chr_id,
        IRanges::IRanges(legacy$start + 1, legacy$stop),
        seqlengths = seqlengths
    )

    session <- tempfile("count-session-")
    dir.create(session)
    on.exit(unlink(session, recursive = TRUE), add = TRUE)
    staged <- stageBam(bam_path, session)

    counts <- countReads(
        staged, gr, mapq = 20L, exclude_flags = 2308L
    )
    expected <- as.integer(legacy$count)
    delta <- counts - expected

    expect_identical(length(counts), nrow(legacy))
    expect_true(all(delta <= 0L))
    expect_identical(delta[seq_len(3L)], c(0L, 0L, 0L))
    expect_lte(mean(delta != 0L), 0.02)
    expect_lte(max(abs(delta)), 20L)
    expect_lt(sum(abs(delta)) / sum(expected), 0.001)
})
