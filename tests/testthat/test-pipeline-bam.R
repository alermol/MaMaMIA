## T3.3: one BAM in, one RCA + one ISA out.

pipeline_fixture <- function() {
    variable_gc_fixture()
}

run_pipeline <- function(fixture, bams = fixture$bam_path,
                         normalization = "none", reference_sample = NULL,
                         max_fit_windows = Inf, windows = 100L, ...) {
    stages <- list()
    progress <- function(stage, detail) {
        stages[[length(stages) + 1L]] <<- list(
            stage = stage, detail = detail
        )
    }
    result <- suppressWarnings(analyzeBamInput(
        bams,
        fixture$config_path,
        normalization = normalization,
        reference_sample = reference_sample,
        windows_per_chromosome = windows,
        cores = 1L,
        max_fit_windows = max_fit_windows,
        progress = progress,
        ...
    ))
    list(result = result, stages = stages)
}

test_that("a species-derived assembly runs in either donor/recipient direction", {
    fixture <- variable_gc_fixture(species = TRUE)
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    config <- readMaMaMIAConfig(fixture$config_path)

    forward <- resolveSpeciesAssembly(config, donor = "donor", recipient = "recipient")
    expect_identical(forward$assembly$pairs, c(d1 = "r1"))
    result <- suppressWarnings(analyzeBamInput(
        fixture$bam_path,
        forward$config,
        assembly = forward$key,
        windows_per_chromosome = 100L,
        cores = 1L,
        max_fit_windows = Inf
    ))
    expect_s3_class(result, "BamAnalysisResult")
    expect_identical(result$assembly$key, "recipient__R__donor__D")
    expect_identical(result$assembly$pairs, c(d1 = "r1"))
    expect_identical(result$bams$sample, "synthetic")

    reversed <- resolveSpeciesAssembly(config, donor = "recipient", recipient = "donor")
    expect_identical(reversed$assembly$pairs, c(r1 = "d1"))
    expect_identical(reversed$assembly$key, "donor__D__recipient__R")
    reversed_result <- suppressWarnings(analyzeBamInput(
        fixture$bam_path,
        reversed$config,
        assembly = reversed$key,
        windows_per_chromosome = 100L,
        cores = 1L,
        max_fit_windows = Inf
    ))
    expect_s3_class(reversed_result, "BamAnalysisResult")
    expect_identical(reversed_result$assembly$pairs, c(r1 = "d1"))
})

test_that("analyzeBamInput() returns a valid BamAnalysisResult", {
    fixture <- pipeline_fixture()
    run <- run_pipeline(fixture)
    result <- run$result

    expect_s3_class(result, "BamAnalysisResult")
    expect_named(result, c(
        "assembly", "config", "bams", "samples", "windows", "gc", "counts",
        "rca", "isa", "normalization", "provenance"
    ))

    expect_identical(result$bams$sample, "synthetic")
    expect_identical(result$samples$sample, "synthetic")
    expect_false(result$samples$reference)
    expect_identical(result$samples$normalization, "none")
    expect_identical(result$samples$n_segments, nrow(result$isa$synthetic$segments))
    expect_identical(
        result$samples$n_putative_introgressions,
        sum(result$isa$synthetic$segments$putative_introgression)
    )
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
    expect_identical(result$isa$synthetic$meta$value_col, "cor.gc")
    expect_true(all(c(
        "chr_id.don", "chr_id.rec", "loc.start", "loc.end",
        "seg.median", "putative_introgression"
    ) %in% names(result$isa$synthetic$segments)))

    expect_identical(result$normalization$mode, "none")
    expect_null(result$normalization$reference)
    expect_true(result$normalization$by_side)
    expect_identical(result$normalization$pseudocount, 0.5)
    expect_false(result$rca$synthetic$fit_sampling$sampled)
    expect_identical(result$rca$synthetic$fit_sampling$max_fit_windows, Inf)
    expect_identical(
        result$rca$synthetic$fit_sampling$n_fitted,
        result$rca$synthetic$fit_sampling$n_ideal
    )
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
    expect_identical(details[[7L]]$value_col, "cor.gc")
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
    expect_true(provenance$max_fit_windows_enforced)
    expect_identical(provenance$n_samples, 1L)
    expect_identical(provenance$total_windows, 200L)
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

cohort_fixture <- function(steps = c(4L, 6L, 8L)) {
    fixture <- pipeline_fixture()
    samples <- paste0("sample-", seq_along(steps))
    bams <- vapply(seq_along(steps), function(i) {
        add_variable_gc_bam(fixture, samples[[i]], step = steps[[i]])
    }, character(1L))
    list(fixture = fixture, bams = bams, samples = samples)
}

test_that("several BAMs work with normalization = 'none'", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    second <- add_variable_gc_bam(fixture, "synthetic-2")
    run <- run_pipeline(fixture, bams = c(fixture$bam_path, second))
    result <- run$result

    expect_identical(result$bams$sample, c("synthetic", "synthetic-2"))
    expect_named(result$rca, c("synthetic", "synthetic-2"))
    expect_named(result$isa, c("synthetic", "synthetic-2"))
    expect_identical(dim(result$counts$donor), c(100L, 2L))
    expect_identical(colnames(result$counts$donor), c("synthetic", "synthetic-2"))
    expect_identical(result$normalization$mode, "none")
    expect_identical(result$samples$normalization, rep("none", 2L))
    expect_false(any(result$samples$reference))
    for (sample in names(result$isa)) {
        expect_identical(result$isa[[sample]]$meta$value_col, "cor.gc")
    }

    stages <- vapply(run$stages, `[[`, character(1), "stage")
    expect_identical(stages[1:4], c("preflight", "index", "windows", "gc"))
    expect_identical(sum(stages == "counts"), 2L)
    expect_identical(sum(stages == "correction"), 2L)
    expect_identical(sum(stages == "normalization"), 0L)
    expect_identical(sum(stages == "segmentation"), 2L)
})

test_that("median normalization runs after all corrections and before CBS", {
    cohort <- cohort_fixture()
    on.exit(unlink(cohort$fixture$directory, recursive = TRUE), add = TRUE)
    run <- run_pipeline(
        cohort$fixture,
        bams = cohort$bams,
        normalization = "median"
    )
    result <- run$result

    expect_identical(result$normalization$mode, "median")
    expect_null(result$normalization$reference)
    expect_true(result$normalization$by_side)
    expect_identical(result$normalization$pseudocount, 0.5)
    expect_identical(result$samples$sample, cohort$samples)
    expect_identical(result$samples$normalization, rep("median", 3L))
    expect_false(any(result$samples$reference))
    expect_identical(dim(result$counts$donor), c(100L, 3L))
    expect_named(result$rca, cohort$samples)
    expect_named(result$isa, cohort$samples)

    for (sample in cohort$samples) {
        expect_true("dose" %in% names(result$rca[[sample]]$data))
        expect_true(all(is.finite(result$rca[[sample]]$data$dose)))
        expect_identical(result$rca[[sample]]$normalization$mode, "median")
        expect_identical(result$isa[[sample]]$meta$value_col, "dose")
        expect_gt(nrow(result$isa[[sample]]$segments), 0L)
        expect_true(all(is.finite(result$isa[[sample]]$segments$seg.median)))
    }

    stages <- vapply(run$stages, `[[`, character(1), "stage")
    expect_identical(stages[1:4], c("preflight", "index", "windows", "gc"))
    expect_identical(sum(stages == "counts"), 3L)
    expect_identical(sum(stages == "correction"), 3L)
    expect_identical(sum(stages == "normalization"), 1L)
    expect_identical(sum(stages == "segmentation"), 3L)
    normalization_stage <- which(stages == "normalization")
    expect_true(all(which(stages == "correction") < normalization_stage))
    expect_true(all(which(stages == "segmentation") > normalization_stage))
    normalization_detail <- run$stages[[normalization_stage]]$detail
    expect_identical(normalization_detail$mode, "median")
    expect_identical(normalization_detail$n_samples, 3L)
})

test_that("reference normalization makes the reference sample's dose one", {
    cohort <- cohort_fixture()
    on.exit(unlink(cohort$fixture$directory, recursive = TRUE), add = TRUE)
    run <- run_pipeline(
        cohort$fixture,
        bams = cohort$bams,
        normalization = "reference",
        reference_sample = "sample-2"
    )
    result <- run$result

    expect_identical(result$normalization$mode, "reference")
    expect_identical(result$normalization$reference, "sample-2")
    expect_identical(result$samples$reference, c(FALSE, TRUE, FALSE))
    expect_true(all(result$rca$`sample-2`$data$dose == 1))
    expect_identical(result$rca$`sample-2`$normalization$reference, "sample-2")
    for (sample in cohort$samples) {
        expect_identical(result$isa[[sample]]$meta$value_col, "dose")
        expect_true(all(is.finite(result$isa[[sample]]$segments$seg.median)))
    }
})

test_that("a cohort result exports a combined BEDPE", {
    cohort <- cohort_fixture(c(4L, 6L))
    on.exit(unlink(cohort$fixture$directory, recursive = TRUE), add = TRUE)
    result <- run_pipeline(
        cohort$fixture,
        bams = cohort$bams,
        normalization = "median"
    )$result

    target <- tempfile(fileext = ".bedpe")
    on.exit(unlink(target), add = TRUE)
    selection <- stats::setNames(
        lapply(result$isa, function(isa) 1L),
        names(result$isa)
    )
    suppressMessages(writeCohortToBedpe(result, target, selected = selection))

    out <- utils::read.table(target, header = TRUE, sep = "\t")
    expect_identical(nrow(out), 2L)
    expect_identical(out$chrom1, rep("d1", 2L))
    expect_identical(out$chrom2, rep("r1", 2L))
    expect_identical(out$sample, cohort$samples)
    expect_identical(names(out)[[ncol(out)]], "sample")
})

test_that("reference normalization validates the reference sample", {
    cohort <- cohort_fixture(c(4L, 6L))
    on.exit(unlink(cohort$fixture$directory, recursive = TRUE), add = TRUE)

    expect_error(
        run_pipeline(
            cohort$fixture,
            bams = cohort$bams,
            normalization = "reference"
        ),
        "requires `reference_sample`",
        class = "MaMaMIAAnalysisError"
    )
    expect_error(
        run_pipeline(
            cohort$fixture,
            bams = cohort$bams,
            normalization = "reference",
            reference_sample = "nope"
        ),
        "not one of the input samples",
        class = "MaMaMIAAnalysisError"
    )
})

test_that("duplicate sample names are rejected before normalization", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    duplicate <- file.path(fixture$directory, "duplicate.bam")
    expect_true(file.copy(fixture$bam_path, duplicate))

    expect_error(
        run_pipeline(
            fixture,
            bams = c(fixture$bam_path, duplicate),
            normalization = "median"
        ),
        "Duplicate sample name",
        class = "MaMaMIASampleError"
    )
})

test_that("a nearly empty sample fails with a clear sample-scoped error", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    sparse <- file.path(fixture$directory, "sparse.bam")
    sequence <- strrep("A", 50L)
    quality <- strrep("I", 50L)
    write_test_bam(
        c(
            "@HD\tVN:1.6\tSO:coordinate",
            "@SQ\tSN:d1\tLN:20000",
            "@SQ\tSN:r1\tLN:20000",
            "@RG\tID:sparse\tSM:sparse",
            sprintf("r1\t0\td1\t1\t60\t50M\t*\t0\t0\t%s\t%s", sequence, quality),
            sprintf("r2\t0\tr1\t1\t60\t50M\t*\t0\t0\t%s\t%s", sequence, quality)
        ),
        sparse
    )

    expect_error(
        run_pipeline(
            fixture,
            bams = c(fixture$bam_path, sparse),
            normalization = "median"
        ),
        "GC correction",
        class = "MaMaMIASampleError"
    )
})

test_that("the cohort performance guard warns above max_fit_windows", {
    cohort <- cohort_fixture(c(4L, 6L))
    on.exit(unlink(cohort$fixture$directory, recursive = TRUE), add = TRUE)

    warnings <- character()
    result <- withCallingHandlers(
        analyzeBamInput(
            cohort$bams,
            cohort$fixture$config_path,
            normalization = "none",
            windows_per_chromosome = 100L,
            cores = 1L,
            max_fit_windows = 300L,
            cache = FALSE
        ),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )

    expect_true(any(grepl("Total analysis windows", warnings)))
    expect_identical(result$provenance$max_fit_windows, 300L)
    expect_identical(result$provenance$total_windows, 400L)
})

test_that("the pipeline applies the fit-window cap per sample", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    run <- run_pipeline(
        fixture,
        max_fit_windows = 500L,
        windows = 500L
    )
    sampling <- run$result$rca$synthetic$fit_sampling

    expect_true(sampling$sampled)
    expect_identical(sampling$max_fit_windows, 500)
    expect_identical(sampling$n_fitted, 500L)
    expect_gt(sampling$n_ideal, sampling$n_fitted)
    expect_identical(sum(sampling$strata$n_fitted), 500L)
    expect_true(all(
        sampling$strata$n_fitted >= pmin(sampling$strata$n_ideal, 10L)
    ))
    expect_true(all(is.finite(run$result$rca$synthetic$data$cor.gc)))
    expect_true(run$result$provenance$max_fit_windows_enforced)
    expect_identical(run$result$provenance$max_fit_windows, 500L)
    expect_gte(nrow(sampling$strata), 2L)
    expect_lte(nrow(sampling$strata), 20L)
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

test_that("a zero-coverage sample in a cohort fails before the GC fit", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    zero_bam <- file.path(fixture$directory, "zero.bam")
    write_test_bam(variable_gc_sam("zero")[1:4], zero_bam)

    expect_error(
        run_pipeline(
            fixture,
            bams = c(zero_bam, fixture$bam_path),
            normalization = "median"
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
