## T3.3: one BAM in, one RCA + one ISA out.

pipeline_fixture <- function() {
    variable_gc_fixture()
}

run_pipeline <- function(fixture, bams = fixture$bam_path, ...) {
    stages <- list()
    progress <- function(stage, detail) {
        stages[[length(stages) + 1L]] <<- list(
            stage = stage, detail = detail
        )
    }
    result <- suppressWarnings(analyzeBamInput(
        bams,
        fixture$config_path,
        normalization = "none",
        windows_per_chromosome = 100L,
        cores = 1L,
        max_fit_windows = Inf,
        progress = progress,
        ...
    ))
    list(result = result, stages = stages)
}

test_that("analyzeBamInput() returns a valid BamAnalysisResult", {
    fixture <- pipeline_fixture()
    run <- run_pipeline(fixture)
    result <- run$result

    expect_s3_class(result, "BamAnalysisResult")
    expect_named(result, c(
        "assembly", "config", "bams", "windows", "gc", "counts",
        "rca", "isa", "normalization", "provenance"
    ))

    expect_identical(result$bams$sample, "synthetic")
    expect_identical(result$bams$n_targets, 2L)
    expect_true(file.exists(result$bams$index_path))
    expect_identical(result$windows$n_per_chromosome, 100L)
    expect_identical(nrow(result$windows$table), 200L)
    expect_identical(nrow(result$gc), 200L)
    expect_true(all(result$gc$gc >= 0 & result$gc$gc <= 1))

    expect_identical(dim(result$counts$donor), c(100L, 1L))
    expect_identical(dim(result$counts$recipient), c(100L, 1L))
    expect_identical(colnames(result$counts$donor), "synthetic")
    expect_gt(sum(result$counts$donor), 0)
    expect_gt(sum(result$counts$recipient), 0)

    expect_named(result$rca, "synthetic")
    expect_s3_class(result$rca$synthetic, "RCA")
    expect_true(result$rca$synthetic$corrected)

    expect_named(result$isa, "synthetic")
    expect_s3_class(result$isa$synthetic, "ISA")
    expect_identical(
        result$isa$synthetic$meta$target_pairs,
        data.frame(
            chr_id.don = "d1", chr_id.rec = "r1",
            stringsAsFactors = FALSE
        )
    )
    expect_gt(nrow(result$isa$synthetic$segments), 0L)
    expect_true(all(c(
        "chr_id.don", "chr_id.rec", "loc.start", "loc.end",
        "seg.median", "putative_introgression"
    ) %in% names(result$isa$synthetic$segments)))

    expect_identical(result$normalization$mode, "none")
    expect_null(result$normalization$reference)
    expect_true(result$normalization$by_side)
    expect_identical(result$normalization$pseudocount, 0.5)
})

test_that("analyzeBamInput() reports every progress stage in order", {
    fixture <- pipeline_fixture()
    run <- run_pipeline(fixture)
    stages <- vapply(run$stages, `[[`, character(1), "stage")

    expect_identical(stages, c(
        "preflight", "index", "windows", "gc", "counts",
        "correction", "segmentation"
    ))

    details <- lapply(run$stages, `[[`, "detail")
    expect_identical(details[[1L]]$n_bams, 1L)
    expect_identical(details[[1L]]$assembly, "variable")
    expect_identical(details[[3L]]$n_windows, 200L)
    expect_identical(details[[4L]]$n_windows, 200L)
    expect_identical(details[[5L]]$sample, "synthetic")
    expect_gt(details[[5L]]$total, 0)
    expect_gt(details[[6L]]$n_ideal, 0L)
    expect_gt(details[[7L]]$n_segments, 0L)
    expect_identical(details[[7L]]$sample, "synthetic")
})

test_that("analyzeBamInput() records provenance and writes the fit cache", {
    fixture <- pipeline_fixture()
    run <- run_pipeline(fixture)
    provenance <- run$result$provenance

    expect_identical(provenance$package, "MaMaMIA")
    expect_match(provenance$package_version, "^[0-9.]+$")
    expect_match(provenance$r_version, "^[0-9.]+$")
    expect_true(all(c(
        "bamsignals", "Rhtslib", "GenomicRanges", "DNAcopy", "glmmTMB"
    ) %in% names(provenance$dependency_versions)))
    expect_match(provenance$htslib_version, "^[0-9]+\\.[0-9]+\\.[0-9]+$")
    expect_true(all(c("uvr_version", "uvr_lock", "uvr_lock_hash") %in% names(provenance)))
    expect_identical(provenance$gcrank_commit, MaMaMIA:::.MaMaMIAgcrankCommit)
    expect_identical(provenance$assembly, "variable")
    expect_identical(provenance$window_rule, "exact_equal_windows")
    expect_identical(provenance$windows_per_chromosome, 100L)
    expect_identical(provenance$n_windows, 200L)
    expect_identical(provenance$mapq, 20L)
    expect_identical(provenance$exclude_flags, 2308L)
    expect_identical(provenance$max_fit_windows, Inf)
    expect_false(provenance$max_fit_windows_enforced)
    expect_identical(provenance$count_semantics, "bamsignals_start_based")
    expect_identical(provenance$gc_source, "gcrank")
    expect_identical(provenance$normalization_mode, "none")
    expect_identical(provenance$normalization_requested, "none")
    expect_true(file.exists(provenance$config_path))
    expect_match(provenance$config_hash, "^[0-9a-f]{32}$")
    expect_true(file.exists(provenance$cache))
    expect_true(all(c(
        "path", "sample", "sample_source", "file_size", "mtime"
    ) %in% names(provenance$bams)))
    unlink(provenance$cache)
})

test_that("the gcrank commit matches the vendored source record", {
    vendor <- testthat::test_path("..", "..", "tools", "vendor-gcrank.txt")
    if (!file.exists(vendor)) {
        testthat::skip("tools/vendor-gcrank.txt is not available")
    }
    commit <- grep("^commit = ", readLines(vendor, warn = FALSE), value = TRUE)
    if (length(commit) != 1L) {
        testthat::skip("tools/vendor-gcrank.txt has no commit record")
    }
    expect_identical(
        MaMaMIA:::.MaMaMIAgcrankCommit,
        sub("^commit = ", "", commit)
    )
})

test_that("several BAMs work with normalization = 'none'", {
    fixture <- pipeline_fixture()
    second <- add_variable_gc_bam(fixture, "synthetic-2")
    run <- run_pipeline(fixture, bams = c(fixture$bam_path, second))
    result <- run$result

    expect_identical(result$bams$sample, c("synthetic", "synthetic-2"))
    expect_named(result$rca, c("synthetic", "synthetic-2"))
    expect_named(result$isa, c("synthetic", "synthetic-2"))
    expect_identical(dim(result$counts$donor), c(100L, 2L))
    expect_identical(colnames(result$counts$donor), c("synthetic", "synthetic-2"))
    expect_identical(result$normalization$mode, "none")

    stages <- vapply(run$stages, `[[`, character(1), "stage")
    expect_identical(stages[1:4], c("preflight", "index", "windows", "gc"))
    expect_identical(sum(stages == "counts"), 2L)
    expect_identical(sum(stages == "correction"), 2L)
    expect_identical(sum(stages == "segmentation"), 2L)
})

test_that("multi-sample normalization is rejected until T4.2", {
    fixture <- pipeline_fixture()
    second <- add_variable_gc_bam(fixture, "synthetic-2")

    expect_error(
        analyzeBamInput(
            c(fixture$bam_path, second),
            fixture$config_path,
            normalization = "median"
        ),
        "after T3\\.3",
        class = "MaMaMIAAnalysisError"
    )
})

test_that("a constant-GC sample fails inside the correction stage", {
    skip_if_no_tiny_fixtures()
    directory <- tempfile("constant-gc-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    bam <- file.path(directory, "tiny.bam")
    expect_true(file.copy(tiny_bam_path(), bam))

    ## The committed tiny fixture has exactly 40% and 60% GC, so the ZINB
    ## smooth cannot be constructed. The pipeline must surface that as a
    ## sample-scoped correction error, not as a raw glmmTMB failure.
    expect_error(
        suppressWarnings(analyzeBamInput(
            bam,
            tiny_config_path(),
            normalization = "none",
            cores = 1L
        )),
        "GC correction",
        class = "MaMaMIASampleError"
    )
})

test_that("a zero-coverage sample fails before the GC fit", {
    fixture <- pipeline_fixture()
    zero_bam <- file.path(fixture$directory, "zero.bam")
    write_test_bam(variable_gc_sam("zero")[1:4], zero_bam)

    expect_error(
        analyzeBamInput(
            zero_bam,
            fixture$config_path,
            normalization = "none",
            windows_per_chromosome = 100L
        ),
        "zero reads",
        class = "MaMaMIASampleError"
    )
})

test_that("a BAM covering only one side is rejected by preflight", {
    fixture <- pipeline_fixture()
    one_side <- file.path(fixture$directory, "one-side.bam")
    sequence <- strrep("A", 50L)
    quality <- strrep("I", 50L)
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:d1\tLN:20000",
            "@RG\tID:one\tSM:one-side",
            sprintf(
                "read\t0\td1\t1\t60\t50M\t*\t0\t0\t%s\t%s",
                sequence, quality
            )
        ),
        one_side
    )

    expect_error(
        analyzeBamInput(
            one_side, fixture$config_path, assembly = "variable"
        ),
        "covers only one side",
        class = "MaMaMIAAssemblyMatchError"
    )
})

test_that("a corrupted BAM is rejected by preflight", {
    fixture <- pipeline_fixture()
    corrupted <- file.path(fixture$directory, "corrupted.bam")
    writeLines("this is not a BAM file", corrupted)

    expect_error(
        analyzeBamInput(corrupted, fixture$config_path),
        "BAM",
        class = "MaMaMIABamHeaderError"
    )
})

test_that("analyzeBamInput() validates its arguments", {
    fixture <- pipeline_fixture()

    expect_error(
        analyzeBamInput(character(0), fixture$config_path),
        "bam_paths",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        analyzeBamInput(
            fixture$bam_path, fixture$config_path,
            windows_per_chromosome = 0L
        ),
        "windows_per_chromosome",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        analyzeBamInput(
            fixture$bam_path, fixture$config_path,
            mapq = 256L
        ),
        "mapq",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        analyzeBamInput(
            fixture$bam_path, fixture$config_path,
            progress = 1L
        ),
        "progress",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        analyzeBamInput(
            fixture$bam_path, fixture$config_path,
            normalization = "invalid"
        ),
        "should be one of",
        class = "MaMaMIAAnalysisError"
    )
})

test_that("a read-only BAM directory is staged in session temp", {
    if (.Platform$OS.type == "windows") {
        testthat::skip("POSIX permission bits are required")
    }
    fixture <- pipeline_fixture()
    read_only <- file.path(fixture$directory, "read-only")
    dir.create(read_only)
    source <- file.path(read_only, "sample.bam")
    expect_true(file.copy(fixture$bam_path, source))
    Sys.chmod(read_only, mode = "0555")
    on.exit(Sys.chmod(read_only, mode = "0755"), add = TRUE)

    result <- suppressWarnings(analyzeBamInput(
        source,
        fixture$config_path,
        normalization = "none",
        windows_per_chromosome = 100L,
        cores = 1L
    ))

    expect_true(result$provenance$staged)
    expect_false(file.exists(paste0(source, ".csi")))
    expect_true(file.exists(result$bams$index_path))
    expect_true(startsWith(
        result$bams$index_path,
        normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
    ))
})
