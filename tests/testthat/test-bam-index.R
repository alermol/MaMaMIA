## T1.2: CSI index creation, reuse, staging and failure modes.

copy_tiny_bam <- function(directory) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    path <- file.path(directory, "tiny.bam")
    expect_true(file.copy(tiny_bam_path(), path))
    path
}

count_tiny_bam <- function(path) {
    gr <- GenomicRanges::GRanges(
        c("d1", "r1"),
        IRanges::IRanges(c(1L, 1L), c(20000L, 20000L)),
        seqlengths = c(d1 = 20000L, r1 = 20000L)
    )
    countReads(path, gr, mapq = 20L, exclude_flags = 0L)
}

test_that("ensureBamIndex() builds a CSI index and countReads() reads it", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("bam-index-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- copy_tiny_bam(directory)

    index <- ensureBamIndex(bam)
    expect_identical(index, paste0(bam, ".csi"))
    expect_true(file.exists(index))
    expect_gt(file.info(index)$size, 0)
    expect_true(MaMaMIA:::.bam_index_magic(index))
    expect_identical(count_tiny_bam(bam), c(8L, 5L))
})

test_that("ensureBamIndex() reuses a current index and rebuilds a stale one", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("bam-index-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- copy_tiny_bam(directory)
    index <- ensureBamIndex(bam)

    first_mtime <- file.info(index)$mtime
    expect_identical(ensureBamIndex(bam), index)
    expect_identical(file.info(index)$mtime, first_mtime)

    Sys.setFileTime(index, Sys.time() - 300)
    stale_mtime <- file.info(index)$mtime
    expect_identical(ensureBamIndex(bam), index)
    expect_gt(as.numeric(file.info(index)$mtime), as.numeric(stale_mtime))
    expect_identical(count_tiny_bam(bam), c(8L, 5L))
})

test_that("ensureBamIndex() replaces corrupt or empty CSI files", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("bam-index-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- copy_tiny_bam(directory)
    index <- paste0(bam, ".csi")

    writeBin(charToRaw("not a csi"), index)
    expect_identical(ensureBamIndex(bam), index)
    expect_true(MaMaMIA:::.bam_index_magic(index))

    file.create(index)
    expect_identical(ensureBamIndex(bam), index)
    expect_true(MaMaMIA:::.bam_index_magic(index))
    expect_identical(count_tiny_bam(bam), c(8L, 5L))
})

test_that("ensureBamIndex() refuses a read-only BAM directory", {
    skip_if_no_tiny_fixtures()
    if (.Platform$OS.type == "windows") {
        testthat::skip("POSIX permission bits are required")
    }
    directory <- tempfile("bam-ro-")
    on.exit({
        Sys.chmod(directory, mode = "0755")
        unlink(directory, recursive = TRUE)
    }, add = TRUE)
    bam <- copy_tiny_bam(directory)
    Sys.chmod(directory, mode = "0555")

    expect_error(
        ensureBamIndex(bam),
        "not writable.*stage the BAM in session temp",
        class = "MaMaMIAIndexError"
    )
})

test_that("stageBam() indexes a symlinked BAM in session temp", {
    skip_if_no_tiny_fixtures()
    source_dir <- tempfile("bam-source-")
    session_dir <- tempfile("bam-session-")
    on.exit(unlink(c(source_dir, session_dir), recursive = TRUE), add = TRUE)
    bam <- copy_tiny_bam(source_dir)
    dir.create(session_dir)

    staged <- stageBam(bam, session_dir)
    on.exit(unlink(dirname(staged), recursive = TRUE), add = TRUE)
    expect_true(file.exists(staged))
    expect_true(startsWith(
        normalizePath(dirname(staged), winslash = "/"),
        paste0(normalizePath(session_dir, winslash = "/"), "/")
    ))
    if (.Platform$OS.type == "unix") {
        expect_true(nzchar(Sys.readlink(staged)))
    }

    index <- ensureBamIndex(staged)
    expect_identical(index, paste0(staged, ".csi"))
    expect_true(file.exists(index))
    expect_false(file.exists(paste0(bam, ".csi")))
    expect_identical(count_tiny_bam(staged), c(8L, 5L))
})

test_that("ensureBamIndex() ignores a BAI sibling and only returns CSI", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("bam-index-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- copy_tiny_bam(directory)
    bai <- paste0(bam, ".bai")
    file.create(bai)

    index <- ensureBamIndex(bam)
    expect_identical(index, paste0(bam, ".csi"))
    expect_true(file.exists(bai))
    expect_true(MaMaMIA:::.bam_index_magic(index))
})

test_that("ensureBamIndex() validates min_shift and threads", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("bam-index-")
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- copy_tiny_bam(directory)

    expect_error(
        ensureBamIndex(bam, min_shift = 0L),
        "min_shift",
        class = "MaMaMIAIndexError"
    )
    expect_error(
        ensureBamIndex(bam, threads = 0L),
        "threads",
        class = "MaMaMIAIndexError"
    )
})

test_that("ensureBamIndex() builds a CSI for the real long-contig BAM", {
    root <- skip_if_no_test_data()
    bam <- file.path(root, "bam_files", "Sq54", "cs_t2t_se_cereale.bam")
    if (!file.exists(bam)) {
        testthat::skip("Sq54 BAM is missing")
    }
    session_dir <- tempfile("bam-session-")
    on.exit(unlink(session_dir, recursive = TRUE), add = TRUE)
    dir.create(session_dir)

    staged <- stageBam(bam, session_dir)
    index <- ensureBamIndex(staged)
    expect_true(file.exists(index))
    expect_true(MaMaMIA:::.bam_index_magic(index))

    header <- bamHeaderInfo(staged)
    contig <- "CP155614.1"
    contig_length <- header$targets$length[
        match(contig, header$targets$name)
    ]
    gr <- GenomicRanges::GRanges(
        contig, IRanges::IRanges(1L, 1000L),
        seqlengths = stats::setNames(contig_length, contig)
    )
    counts <- countReads(
        staged, gr, mapq = 20L, exclude_flags = 0L
    )
    expect_length(counts, 1L)
    expect_gte(counts[[1L]], 0L)
})
