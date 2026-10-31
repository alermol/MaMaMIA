## The pipeline on a synthetic BAM with a known event.

event_sam <- function(sample, donor_event = TRUE, contig_length = 20000L,
                      read_length = 50L, step = 200L,
                      event_from = 6000L, event_to = 10000L,
                      event_step = 20L) {
    sequence <- strrep("A", read_length)
    quality <- strrep("I", read_length)
    make_reads <- function(contig, positions) {
        sprintf(
            "%s_%s_%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\tRG:Z:%s",
            sample, contig, positions, contig, positions, read_length,
            sequence, quality, sample
        )
    }
    base_positions <- seq(1L, contig_length - read_length + 1L, by = step)
    donor_positions <- if (isTRUE(donor_event)) {
        sort(unique(c(
            base_positions,
            seq(event_from, event_to, by = event_step)
        )))
    } else {
        base_positions
    }
    c(
        "@HD\tVN:1.6\tSO:coordinate",
        sprintf("@SQ\tSN:d1\tLN:%d", contig_length),
        sprintf("@SQ\tSN:r1\tLN:%d", contig_length),
        sprintf("@RG\tID:%s\tSM:%s", sample, sample),
        make_reads("d1", donor_positions),
        make_reads("r1", base_positions)
    )
}

event_bam <- function(fixture, sample, donor_event = TRUE, filename = NULL) {
    if (is.null(filename)) {
        filename <- paste0(sample, ".bam")
    }
    path <- file.path(fixture$directory, filename)
    write_test_bam(
        event_sam(sample, donor_event = donor_event),
        path,
        index = TRUE
    )
    path
}

pipeline_fixture <- function() {
    fixture <- variable_gc_fixture(sample = "P", windows = 100L)
    fixture$bam_path <- event_bam(fixture, "P", donor_event = TRUE)
    fixture$control_path <- event_bam(
        fixture, "C", donor_event = FALSE, filename = "control.bam"
    )
    fixture
}

test_that("analyzeBamInput() accepts a runtime donor aggregation", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    result <- analyzeBamInput(
        fixture$bam_path, fixture$config,
        assembly = fixture$assembly,
        windows_per_chromosome = 100L,
        donor_aggregation = "sum",
        cores = 1L
    )
    expect_identical(result$provenance$donor_aggregation, "sum")
    expect_error(
        analyzeBamInput(
            fixture$bam_path, fixture$config,
            assembly = fixture$assembly,
            windows_per_chromosome = 100L,
            donor_aggregation = "median",
            cores = 1L
        ),
        "donor_aggregation"
    )
})

test_that("the pipeline returns the BamAnalysisResult contract", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    progress <- character()
    result <- analyzeBamInput(
        fixture$bam_path, fixture$config,
        assembly = fixture$assembly,
        windows_per_chromosome = 100L,
        calibration = list(
            min_mb = 0.001,
            length_bins_mb = c(0.001, 0.01, Inf),
            floors = c(0.3, 0.3),
            fdr = 0.01
        ),
        progress = function(stage, detail) {
            progress <<- c(progress, stage)
        },
        cores = 1L
    )
    expect_s3_class(result, "BamAnalysisResult")
    expect_named(result, c(
        "assembly", "config", "bams", "samples", "windows", "gc", "metrics",
        "scores", "masks", "calibration", "calls", "provenance"
    ))
    expect_identical(result$samples$sample, "P")
    expect_identical(result$samples$mode, "single")
    expect_true(all(c(
        "metrics", "score", "mask", "segmentation", "calibration"
    ) %in% progress))
    ## Every stage the pipeline emits must be a planned stage: an unplanned
    ## one advances a progress bar that has no room for it, and the app's
    ## label table used to abort the analysis on an unknown name.
    expect_true(all(
        unique(progress) %in%
            MaMaMIA:::bam_progress_stages(length(unique(result$samples$sample)))
    ))
    expect_identical(result$provenance$pipeline_version, "0.1.0")
    expect_identical(result$provenance$breadth_backend, "htslib")
    expect_identical(result$provenance$gc_source, "gcrank")
    expect_equal(sum(result$masks$mask), 20)
    expect_gte(sum(result$calls$pass), 1L)
    expect_true(any(result$calls$class == "addition"))
    expect_silent(validateBamAnalysisResult(result))
})

test_that("a cohort analysis keeps the event in the positive sample", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    result <- analyzeBamInput(
        c(fixture$bam_path, fixture$control_path),
        fixture$config,
        assembly = fixture$assembly,
        windows_per_chromosome = 100L,
        calibration = list(
            min_mb = 0.001,
            length_bins_mb = c(0.001, 0.01, Inf),
            floors = c(0.25, 0.25),
            fdr = 0.01
        ),
        cores = 1L
    )
    expect_identical(sort(result$samples$sample), c("C", "P"))
    expect_identical(unique(result$samples$mode), "multi")
    expect_gte(sum(result$calls$pass[result$calls$sample == "P"]), 1L)
    expect_equal(sum(result$calls$pass[result$calls$sample == "C"]), 0L)
    expect_true(all(result$scores$mode == "multi"))
})

test_that("writeCallsToBedpe() writes donor and recipient intervals", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    result <- bam_example_result(1000)
    output <- tempfile(fileext = ".bedpe")
    on.exit(unlink(output), add = TRUE)
    writeCallsToBedpe(result, output)
    bedpe <- utils::read.delim(output, stringsAsFactors = FALSE)
    expect_true(nrow(bedpe) >= 1L)
    expect_true(all(c(
        "sample", "chr_id.don", "start.don", "stop.don",
        "chr_id.rec", "start.rec", "stop.rec", "score_mean", "class",
        "class_group", "class_evidence", "class_confidence", "class_status",
        "copy_state", "total_evidence", "component_status",
        "mappability_status", "confidence", "method", "mode",
        "recipient_component",
        "recipient_chr_name", "recipient_group", "recipient_assignment",
        "recipient_cross_group", "recipient_expected_breadth",
        "recipient_deletion_width_bp", "recipient_observed_breadth",
        "recipient_copy_number", "affected_fraction", "affected_copies",
        "dosage_confidence", "recipient_depletion_se",
        "recipient_width_overlap", "recipient_mapq_fraction",
        "donor_copies", "total_mean", "donor_component", "donor_chr_name"
    ) %in% names(bedpe)))
    hosted <- bedpe$recipient_assignment != "no_recipient_loss"
    expect_true(any(hosted))
    expect_true(all(!is.na(bedpe$recipient_component[hosted])))
    expect_true(all(!is.na(bedpe$recipient_chr_name[hosted])))
    expect_true(all(bedpe$start.don >= 0))
    expect_true(all(bedpe$stop.don > bedpe$start.don))
    expect_true(all(bedpe$start.rec[hosted] >= 0))
    expect_true(all(bedpe$stop.rec[hosted] > bedpe$start.rec[hosted]))
    expect_true(all(is.na(bedpe$start.rec[!hosted])))
})

test_that("donor_only scoring runs through the pipeline", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    result <- analyzeBamInput(
        fixture$bam_path, fixture$config,
        assembly = fixture$assembly,
        windows_per_chromosome = 100L,
        score = "donor_only",
        calibration = list(
            min_mb = 0.001,
            length_bins_mb = c(0.001, 0.01, Inf),
            floors = c(0.3, 0.3)
        ),
        cores = 1L
    )
    expect_identical(unique(result$scores$score_type), "donor_only")
    expect_true(all(result$calls$method == "donor_only"))
})

test_that("the pipeline rejects a version 1 config", {
    fixture <- pipeline_fixture()
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    lines <- readLines(fixture$config_path, warn = FALSE)
    lines <- sub("^version:.*$", "version: 1", lines)
    writeLines(lines, fixture$config_path)
    expect_error(
        analyzeBamInput(
            fixture$bam_path, fixture$config_path,
            assembly = fixture$assembly, windows_per_chromosome = 100L
        ),
        "version 1"
    )
})
