## T1.4: multi-BAM preflight, limits, sample resolution and containment.

preflight_tiny_config <- function() {
    skip_if_no_tiny_fixtures()
    readMaMaMIAConfig(tiny_config_path())
}

write_header_bam <- function(path, targets = c("d1", "r1"),
                             rg = "@RG\tID:test\tSM:test-sample") {
    lines <- c(
        "@HD\tVN:1.6\tSO:coordinate",
        sprintf("@SQ\tSN:%s\tLN:20000", targets)
    )
    if (!is.null(rg)) {
        lines <- c(lines, rg)
    }
    write_test_bam(lines, path)
}

test_that("preflightBams() describes the tiny fixture", {
    skip_if_no_tiny_fixtures()
    config <- readMaMaMIAConfig(tiny_config_path())

    table <- preflightBams(tiny_bam_path(), config)

    expect_s3_class(table, "data.frame")
    expect_identical(nrow(table), 1L)
    expect_named(table, c(
        "path", "file_size", "sample", "sample_source", "n_samples",
        "assembly", "n_targets", "n_assembly_targets", "index_path",
        "index_status", "index_action", "warnings", "valid"
    ))
    expect_identical(
        table$path,
        normalizePath(tiny_bam_path(), winslash = "/", mustWork = TRUE)
    )
    expect_equal(
        table$file_size,
        as.numeric(file.info(tiny_bam_path())$size)
    )
    expect_identical(table$sample, "tiny-sample")
    expect_identical(table$sample_source, "rg")
    expect_identical(table$n_samples, 1L)
    expect_identical(table$assembly, "tiny")
    expect_identical(table$n_targets, 2L)
    expect_identical(table$n_assembly_targets, 2L)
    expect_identical(table$index_path, paste0(table$path, ".csi"))
    expect_identical(table$index_status, "valid")
    expect_identical(table$index_action, "reuse")
    expect_identical(table$warnings, "")
    expect_true(table$valid)
})

test_that("preflightBams() accepts a config path and several samples", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()
    second <- tempfile("second-sample-", fileext = ".bam")
    on.exit(unlink(second), add = TRUE)
    write_header_bam(second, rg = "@RG\tID:second\tSM:second-sample")

    table <- preflightBams(
        c(tiny_bam_path(), second), tiny_config_path()
    )
    expect_identical(table$sample, c("tiny-sample", "second-sample"))
    expect_identical(table$assembly, c("tiny", "tiny"))
    expect_identical(table$index_status, c("valid", "missing"))
    expect_identical(table$index_action, c("reuse", "create"))
})

test_that("preflightBams() falls back to the file name and warns", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()
    path <- tempfile("fallback-", fileext = ".bam")
    on.exit(unlink(path), add = TRUE)
    write_header_bam(path, rg = NULL)

    table <- preflightBams(path, preflight_tiny_config())
    expect_identical(
        table$sample,
        sub("\\.bam$", "", basename(path), ignore.case = TRUE)
    )
    expect_identical(table$sample_source, "basename")
    expect_match(table$warnings, "derived from the file name", fixed = TRUE)
    expect_identical(table$index_status, "missing")
    expect_identical(table$index_action, "create")
    expect_match(table$warnings, "CSI index will be created", fixed = TRUE)
})

test_that("preflightBams() rejects duplicate sample names", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("duplicate-samples-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    first <- file.path(directory, "first.bam")
    second <- file.path(directory, "second.bam")
    expect_true(file.copy(tiny_bam_path(), first))
    expect_true(file.copy(tiny_bam_path(), second))

    expect_error(
        preflightBams(c(first, second), preflight_tiny_config()),
        "Duplicate sample name 'tiny-sample'.*first[.]bam.*second[.]bam",
        class = "MaMaMIASampleError"
    )
})

test_that("preflightBams() rejects a basename-fallback collision", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()
    first_dir <- tempfile("collision-a-")
    second_dir <- tempfile("collision-b-")
    dir.create(first_dir)
    dir.create(second_dir)
    on.exit(unlink(c(first_dir, second_dir), recursive = TRUE), add = TRUE)
    first <- file.path(first_dir, "sample.bam")
    second <- file.path(second_dir, "sample.bam")
    write_header_bam(first, rg = NULL)
    write_header_bam(second, rg = NULL)

    expect_error(
        preflightBams(c(first, second), preflight_tiny_config()),
        "Duplicate sample name 'sample'",
        class = "MaMaMIASampleError"
    )
})

test_that("preflightBams() expands a multi-sample BAM into one row per sample", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()
    path <- tempfile("multi-sm-", fileext = ".bam")
    on.exit(unlink(path), add = TRUE)
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:d1\tLN:20000",
            "@SQ\tSN:r1\tLN:20000",
            "@RG\tID:a\tSM:one",
            "@RG\tID:b\tSM:two"
        ),
        path
    )

    table <- preflightBams(path, preflight_tiny_config())
    expect_identical(nrow(table), 2L)
    expect_identical(table$sample, c("one", "two"))
    expect_identical(table$n_samples, c(2L, 2L))
    expect_identical(table$path, rep(
        normalizePath(path, winslash = "/", mustWork = TRUE), 2L
    ))
    expect_identical(table$index_status, c("missing", "missing"))
})

test_that("preflightBams() rejects BAMs matching different assemblies", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()
    config <- preflight_tiny_config()
    config$assemblies$other <- list(
        key = "other",
        label = "Other assembly",
        description = "",
        gc_indexes = character(0),
        contigs = list(
            recipient = data.frame(
                id = "x1", name = "1X", subgenome = "X", length = 20000,
                index_path = NA_character_, stringsAsFactors = FALSE
            ),
            donor = data.frame(
                id = "x2", name = "2X", subgenome = "X", length = 20000,
                index_path = NA_character_, stringsAsFactors = FALSE
            )
        ),
        pairs = c(x2 = "x1")
    )
    path <- tempfile("other-assembly-", fileext = ".bam")
    on.exit(unlink(path), add = TRUE)
    write_header_bam(
        path, targets = c("x1", "x2"),
        rg = "@RG\tID:other\tSM:other-sample"
    )

    expect_error(
        preflightBams(c(tiny_bam_path(), path), config),
        "match different assemblies",
        class = "MaMaMIAAssemblyMatchError"
    )
})

test_that("preflightBams() enforces max_bams and max_total_bytes", {
    skip_if_no_tiny_fixtures()
    config <- preflight_tiny_config()
    config$server$max_bams <- 1L
    expect_error(
        preflightBams(rep(tiny_bam_path(), 2L), config),
        "Too many BAM files: 2 \\(limit 1\\)",
        class = "MaMaMIALimitError"
    )

    config <- preflight_tiny_config()
    config$server$max_total_bytes <- file.info(tiny_bam_path())$size - 1
    expect_error(
        preflightBams(tiny_bam_path(), config),
        "Total BAM size .* exceeds",
        class = "MaMaMIALimitError"
    )
})

test_that("preflightBams() accepts uploads without a server BAM root", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("bam-upload-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    upload <- file.path(directory, "upload.bam")
    expect_true(file.copy(tiny_bam_path(), upload))

    config <- preflight_tiny_config()
    ## Old registries may still carry the removed key; it must not restrict
    ## uploads or fail validation.
    config$server$bam_root <- tempfile("legacy-root-")
    expect_no_error(preflightBams(upload, config))
})

test_that("preflightBams() reports CSI status and planned action", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("bam-status-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- file.path(directory, "tiny.bam")
    expect_true(file.copy(tiny_bam_path(), bam))
    config <- preflight_tiny_config()

    missing <- preflightBams(bam, config)
    expect_identical(missing$index_status, "missing")
    expect_identical(missing$index_action, "create")

    index <- paste0(bam, ".csi")
    expect_true(file.copy(paste0(tiny_bam_path(), ".csi"), index))
    valid <- preflightBams(bam, config)
    expect_identical(valid$index_status, "valid")
    expect_identical(valid$index_action, "reuse")

    Sys.setFileTime(index, Sys.time() - 300)
    stale <- preflightBams(bam, config)
    expect_identical(stale$index_status, "stale")
    expect_identical(stale$index_action, "rebuild")

    writeBin(charToRaw("not a csi"), index)
    corrupt <- preflightBams(bam, config)
    expect_identical(corrupt$index_status, "corrupt")
    expect_identical(corrupt$index_action, "rebuild")

    file.create(paste0(bam, ".bai"))
    with_bai <- preflightBams(bam, config)
    expect_match(
        with_bai$warnings, "BAI index is present but ignored", fixed = TRUE
    )
})

test_that("preflightBams() plans staging for read-only directories", {
    skip_if_no_tiny_fixtures()
    if (.Platform$OS.type == "windows") {
        testthat::skip("POSIX permission bits are required")
    }
    directory <- tempfile("bam-readonly-")
    dir.create(directory)
    on.exit({
        Sys.chmod(directory, mode = "0755")
        unlink(directory, recursive = TRUE)
    }, add = TRUE)
    bam <- file.path(directory, "tiny.bam")
    expect_true(file.copy(tiny_bam_path(), bam))
    Sys.chmod(directory, mode = "0555")

    table <- preflightBams(bam, preflight_tiny_config())
    expect_identical(table$index_status, "missing")
    expect_identical(table$index_action, "stage")
    expect_match(table$warnings, "staged in session temp", fixed = TRUE)
})

test_that("preflightBams() warns about extra header targets", {
    skip_if_no_tiny_fixtures()
    skip_if_no_samtools()
    path <- tempfile("extra-targets-", fileext = ".bam")
    on.exit(unlink(path), add = TRUE)
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:d1\tLN:20000",
            "@SQ\tSN:r1\tLN:20000",
            "@SQ\tSN:x1\tLN:5000",
            "@RG\tID:test\tSM:extra-sample"
        ),
        path
    )

    table <- preflightBams(path, preflight_tiny_config())
    expect_match(
        table$warnings,
        "1 target\\(s\\) not in assembly 'Tiny synthetic assembly'",
        fixed = FALSE
    )
})

test_that("preflightBams() rejects unreadable BAM files", {
    skip_if_no_tiny_fixtures()
    if (.Platform$OS.type == "windows") {
        testthat::skip("POSIX permission bits are required")
    }
    directory <- tempfile("bam-unreadable-")
    dir.create(directory)
    on.exit({
        Sys.chmod(file.path(directory, "tiny.bam"), mode = "0644")
        unlink(directory, recursive = TRUE)
    }, add = TRUE)
    bam <- file.path(directory, "tiny.bam")
    expect_true(file.copy(tiny_bam_path(), bam))
    Sys.chmod(bam, mode = "0000")

    expect_error(
        preflightBams(bam, preflight_tiny_config()),
        "not readable",
        class = "MaMaMIABamHeaderError"
    )
})

test_that("preflightBams() validates its inputs", {
    skip_if_no_tiny_fixtures()
    config <- preflight_tiny_config()
    expect_error(
        preflightBams(character(0), config),
        "non-empty character vector",
        class = "MaMaMIAPreflightError"
    )
    expect_error(
        preflightBams(tempfile(fileext = ".bam"), config),
        "does not exist",
        class = "MaMaMIABamHeaderError"
    )
    expect_error(
        preflightBams(tiny_bam_path(), list(version = 1L)),
        "validated MaMaMIA configuration",
        class = "MaMaMIAConfigError"
    )
    expect_error(
        preflightBams(tiny_bam_path(), config, assembly = "unknown"),
        "Unknown assembly",
        class = "MaMaMIAAssemblyMatchError"
    )
})
