## T1.1: Rhtslib header parsing and D2 sample resolution.

test_that("bamHeaderInfo() parses the tiny fixture header", {
    skip_if_no_tiny_fixtures()

    header <- bamHeaderInfo(tiny_bam_path())

    expect_identical(header$format, "BAM")
    expect_identical(header$sort_order, "coordinate")
    expect_identical(header$n_targets, 2L)
    expect_identical(header$targets$name, c("d1", "r1"))
    expect_identical(header$targets$length, c(20000, 20000))
    expect_identical(header$samples, "tiny-sample")
    expect_identical(header$sample, "tiny-sample")
    expect_identical(header$sample_source, "rg")
    expect_identical(header$read_groups$id, "tiny")
    expect_identical(header$read_groups$sample, "tiny-sample")
    expect_identical(header$programs$id, "make-test-fixtures")
    expect_identical(header$programs$name, "make-test-fixtures")
    expect_identical(header$programs$version, "0.0.1")
    expect_true(is.na(header$programs$command))
    expect_match(header$header_text, "^@HD")
    expect_identical(header$path, normalizePath(tiny_bam_path(), winslash = "/"))
})

test_that("bamHeaderInfo() falls back to the basename without @RG SM", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()

    path <- tempfile("fallback-", fileext = ".bam")
    on.exit(unlink(path), add = TRUE)
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:d1\tLN:20000",
            "@SQ\tSN:r1\tLN:20000"
        ),
        path
    )

    header <- bamHeaderInfo(path)
    fallback <- sub("\\.bam$", "", basename(path), ignore.case = TRUE)
    expect_identical(header$samples, fallback)
    expect_identical(header$sample, fallback)
    expect_identical(header$sample_source, "basename")
    expect_identical(header$read_groups$id, "<basename>")
    expect_identical(header$read_groups$sample, fallback)
    expect_false(header$multi_sample)
})

test_that("bamHeaderInfo() returns every sample from several @RG SM values", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()

    path <- tempfile("multi-sm-", fileext = ".bam")
    on.exit(unlink(path), add = TRUE)
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:d1\tLN:20000",
            "@RG\tID:a\tSM:one",
            "@RG\tID:b\tSM:two"
        ),
        path
    )

    header <- bamHeaderInfo(path)
    expect_identical(header$samples, c("one", "two"))
    expect_identical(header$read_groups$id, c("a", "b"))
    expect_identical(header$read_groups$sample, c("one", "two"))
    expect_true(is.na(header$sample))
    expect_true(header$multi_sample)
    expect_identical(header$sample_source, "rg")
})

test_that("bamHeaderInfo() rejects incomplete @RG lines", {
    skip_if_no_samtools()

    missing_sm <- tempfile("missing-sm-", fileext = ".bam")
    missing_id <- tempfile("missing-id-", fileext = ".bam")
    duplicate_id <- tempfile("duplicate-id-", fileext = ".bam")
    on.exit(unlink(c(missing_sm, missing_id, duplicate_id)), add = TRUE)
    header <- c(
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:d1\tLN:20000"
    )
    write_test_bam(c(header, "@RG\tID:a"), missing_sm)
    write_test_bam(c(header, "@RG\tSM:one"), missing_id)
    write_test_bam(c(header, "@RG\tID:a\tSM:one", "@RG\tID:a\tSM:two"), duplicate_id)

    expect_error(
        bamHeaderInfo(missing_sm),
        "without a sample ID.*@RG SM",
        class = "MaMaMIASampleError"
    )
    expect_error(
        bamHeaderInfo(missing_id),
        "without a read group ID.*@RG ID",
        class = "MaMaMIASampleError"
    )
    expect_error(
        bamHeaderInfo(duplicate_id),
        "duplicate @RG ID.*'a'",
        class = "MaMaMIASampleError"
    )
})

test_that("bamHeaderInfo() rejects non-coordinate BAMs", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()

    paths <- c(
        unsorted = tempfile("unsorted-", fileext = ".bam"),
        unknown = tempfile("unknown-", fileext = ".bam")
    )
    on.exit(unlink(paths), add = TRUE)
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:unsorted",
            "@SQ\tSN:d1\tLN:20000"
        ),
        paths[["unsorted"]]
    )
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:unknown",
            "@SQ\tSN:d1\tLN:20000"
        ),
        paths[["unknown"]]
    )

    expect_error(
        bamHeaderInfo(paths[["unsorted"]]),
        "not coordinate-sorted .*SO:unsorted",
        class = "MaMaMIABamHeaderError"
    )
    expect_error(
        bamHeaderInfo(paths[["unknown"]]),
        "not coordinate-sorted .*SO:unknown",
        class = "MaMaMIABamHeaderError"
    )
})

test_that("bamHeaderInfo() rejects BAMs without targets", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()

    path <- tempfile("no-targets-", fileext = ".bam")
    on.exit(unlink(path), add = TRUE)
    write_test_bam("@HD\tVN:1.6\tSO:coordinate", path)

    expect_error(
        bamHeaderInfo(path),
        "contains no target sequences",
        class = "MaMaMIABamHeaderError"
    )
})

test_that("bamHeaderInfo() rejects missing, empty and truncated BAMs", {
    missing <- tempfile("missing-", fileext = ".bam")
    expect_error(
        bamHeaderInfo(missing),
        "does not exist",
        class = "MaMaMIABamHeaderError"
    )

    empty <- tempfile("empty-", fileext = ".bam")
    file.create(empty)
    on.exit(unlink(empty), add = TRUE)
    expect_error(
        bamHeaderInfo(empty),
        "is empty",
        class = "MaMaMIABamHeaderError"
    )

    skip_if_no_tiny_fixtures()
    truncated <- tempfile("truncated-", fileext = ".bam")
    on.exit(unlink(truncated), add = TRUE)
    connection <- file(tiny_bam_path(), "rb")
    bytes <- readBin(connection, what = "raw", n = 20L)
    close(connection)
    writeBin(bytes, truncated)
    expect_error(
        bamHeaderInfo(truncated),
        "Could not",
        class = "MaMaMIABamHeaderError"
    )
})

test_that("bamHeaderInfo() rejects SAM and other non-BAM formats", {
    sam <- tempfile("sample-", fileext = ".sam")
    on.exit(unlink(sam), add = TRUE)
    writeLines(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:d1\tLN:20000"
        ),
        sam
    )

    expect_error(
        bamHeaderInfo(sam),
        "Unsupported input format 'SAM'",
        class = "MaMaMIAUnsupportedFormatError"
    )
})

test_that("bamHeaderInfo() reads the real Sq54 BAM header without an index", {
    root <- skip_if_no_test_data()
    bam_path <- file.path(root, "bam_files", "Sq54", "cs_t2t_se_cereale.bam")
    if (!file.exists(bam_path)) {
        testthat::skip("Sq54 BAM is missing")
    }

    header <- bamHeaderInfo(bam_path)
    expect_identical(header$sample, "Sq54")
    expect_identical(header$sample_source, "rg")
    expect_identical(header$n_targets, 263L)
    expect_identical(
        unname(header$targets$length[header$targets$name == "CP155614.1"]),
        603818960
    )
})
