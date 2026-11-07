## Brief call table and downloadable result bundle.

test_that("the brief call table keeps the inspection columns", {
    many <- data.frame(
        sample = c("A", "A", "B", "B"),
        stringsAsFactors = FALSE
    )
    expect_identical(
        bam_brief_call_columns(many),
        c(
            "sample", "group", "class", "copy_state", "class_status",
            "recipient_chr_name", "donor_chr_name", "recipient_assignment",
            "start_bp", "end_bp", "width_bp", "pass"
        )
    )
    one <- data.frame(sample = c("A", "A"), stringsAsFactors = FALSE)
    expect_identical(
        bam_brief_call_columns(one),
        c(
            "group", "class", "copy_state", "class_status",
            "recipient_chr_name", "donor_chr_name", "recipient_assignment",
            "start_bp", "end_bp", "width_bp", "pass"
        )
    )
    labels <- bam_call_column_labels(bam_brief_call_columns(one))
    expect_identical(labels[[1L]], "Group")
    expect_identical(labels[[5L]], "Recipient")
    expect_identical(labels[[6L]], "Donor")
    expect_identical(labels[[11L]], "Pass")
})

test_that("the bundle column dictionary is complete and unique", {
    dictionary <- bam_bundle_column_dictionary()
    expect_identical(nrow(dictionary), 69L)
    expect_identical(anyDuplicated(dictionary$old), 0L)
    expect_identical(anyDuplicated(dictionary$name), 0L)
    expect_true(all(c(
        "call_class", "recipient_mean_depletion",
        "donor_to_recipient_width_ratio", "passes_calibration"
    ) %in% dictionary$name))
    renamed <- bam_bundle_descriptive_calls(data.frame(
        group = "1", class = "addition", stringsAsFactors = FALSE
    ))
    expect_named(renamed, c("chromosome_group", "call_class"))
})

test_that("every call column has a label and a description", {
    metadata <- bam_call_column_metadata()
    expect_identical(nrow(metadata), 69L)
    expect_identical(anyDuplicated(metadata$column), 0L)
    expect_identical(anyDuplicated(metadata$label), 0L)
    expect_false(anyNA(metadata$label))
    expect_false(any(!nzchar(metadata$label)))
    expect_false(anyNA(metadata$description))
    expect_false(any(!nzchar(metadata$description)))
    expect_identical(
        metadata$description[metadata$column == "donor_chr_name"],
        "Display name of the donor contig."
    )
})

test_that("the extra-column picker keeps defaults out of its choices", {
    choices <- bam_call_extra_column_choices()
    expect_true(all(nzchar(names(choices))))
    expect_true("Mean depletion" %in% names(choices))
    expect_true("recipient_mean_depletion" %in% unname(choices))
    expect_false(any(unname(choices) %in% bam_brief_call_columns()))
    expect_false("sample" %in% unname(choices))
})

test_that("the bundle records every effective profile parameter", {
    provenance <- list(
        analysis_profile = "custom",
        analysis_profile_base = "sensitive",
        analysis_profile_parameters = list(
            windows_per_chromosome = 1000L,
            metric = "combined",
            score = "paired_log_ratio",
            donor_aggregation = "max",
            mapq = 30L,
            calibration = list(
                fdr = 0.01,
                floors = c(1.5, 1.2, 1.0, 0.8, 0.4),
                min_mb = 5,
                merge_gap = 100
            )
        )
    )
    table <- bam_bundle_profile_parameters(provenance)
    expect_true(all(c(
        "windows_per_chromosome", "metric", "score",
        "donor_aggregation", "mapq", "calibration_fdr",
        "calibration_floor_1", "calibration_floor_5",
        "calibration_min_mb", "calibration_merge_gap"
    ) %in% table$parameter))
    expect_identical(unique(table$profile), "custom")
    expect_identical(unique(table$base_profile), "sensitive")
    expect_true(all(table$changed[table$parameter == "windows_per_chromosome"]))
    expect_true(all(table$changed[table$parameter == "mapq"]))
    expect_false(any(table$changed[table$parameter == "metric"]))
    expect_equal(
        table$default[table$parameter == "windows_per_chromosome"],
        "100"
    )
    expect_equal(table$default[table$parameter == "mapq"], "20")
})

test_that("the bundle filename is fixed and sample folders are safe", {
    expect_identical(
        bam_bundle_filename("Sq 74/7D"), "mamamia-results.zip"
    )
    expect_identical(bam_bundle_filename(NULL), "mamamia-results.zip")
    expect_identical(bam_bundle_filename(""), "mamamia-results.zip")
    expect_identical(
        bam_cohort_bundle_filename(), "mamamia-results.zip"
    )
    expect_identical(
        unname(bam_bundle_sample_folders(c("qc", "run", "A", "a"))),
        c("qc-2", "run-2", "A", "a-2")
    )
})

test_that("the static karyotype renders calls and states", {
    result <- bam_example_result(1000)
    sample <- unique(as.character(result$calls$sample))[[1L]]
    payload <- karyotype_payload(result, sample, "calls")
    svg <- bam_karyotype_static_svg(
        payload,
        subtitle = paste0("Sample: ", sample, " \u00b7 All calls"),
        footer = "MaMaMIA test"
    )
    expect_match(svg, "<svg", fixed = TRUE)
    expect_match(svg, 'preserveAspectRatio="xMidYMid meet"', fixed = TRUE)
    expect_match(svg, 'text-anchor="middle"', fixed = TRUE)
    expect_match(svg, "MaMaMIA test", fixed = TRUE)
    expect_match(svg, "Substitution", fixed = TRUE)
    expect_match(svg, "Deleted recipient", fixed = TRUE)
    expect_match(svg, "Added donor", fixed = TRUE)
    expect_match(svg, sample, fixed = TRUE)
    empty <- bam_karyotype_static_svg(
        list(chromosomes = list(), sample = "empty", palette = NULL)
    )
    expect_match(empty, "No chromosomes", fixed = TRUE)
})

test_that("the confident set is calibrated and reviewed as confident", {
    result <- bam_example_result(1000)
    confident <- bam_confident_calls(result$calls)
    expect_gt(nrow(confident), 0L)
    expect_lt(nrow(confident), nrow(result$calls))
    expect_true(all(isTRUE_vec(confident$pass)))
    expect_true(all(as.character(confident$class_status) == "confident"))
})

test_that("the result bundle contains the table, BEDPE and karyotypes", {
    skip_if(inherits(try(bam_bundle_zip_tool(), silent = TRUE),
        "try-error"), "the 7z command is unavailable")
    result <- bam_example_result(1000)
    sample <- unique(as.character(result$calls$sample))[[1L]]
    expected <- result$calls[
        as.character(result$calls$sample) == sample, ,
        drop = FALSE
    ]
    confident_calls <- bam_confident_calls(expected)
    current <- karyotype_export_payload(
        result, sample = sample, mode = "calls",
        calls = expected[isTRUE_vec(expected$pass), , drop = FALSE],
        selected_groups = NULL
    )
    full <- bam_full_karyotype_payload(current)
    expect_null(current$full_calls)
    expect_true(length(current$full_chromosomes) > 0L)
    confident_payload <- karyotype_payload(
        result, sample = sample, mode = "calls",
        calls = confident_calls, selected_groups = NULL
    )
    expect_true(all(vapply(
        full$chromosomes, function(x) isTRUE(x$selected), logical(1)
    )))
    expect_null(full$full_calls)
    expect_null(full$full_chromosomes)
    expect_identical(
        sum(vapply(full$chromosomes, function(x) length(x$calls), integer(1))),
        sum(vapply(
            current$full_chromosomes,
            function(x) length(x$calls), integer(1)
        ))
    )
    archive <- tempfile(fileext = ".zip")
    on.exit(unlink(archive), add = TRUE)
    bam_write_result_bundle(
        result, archive,
        sample = sample,
        current_payload = current,
        full_payload = full,
        confident_payload = confident_payload
    )
    expect_true(file.exists(archive))
    listing <- utils::unzip(archive, list = TRUE)
    names <- basename(listing$Name)
    expect_true(all(c(
        "README.txt", "calls-full.csv", "calls-confident.csv",
        "calls-all.bedpe", "calls-passing.bedpe",
        "calls-confident.bedpe", "karyotype-current.svg",
        "karyotype-full.svg", "karyotype-confident.svg",
        "report.html", "manifest.json", "checksums.md5",
        "calls-passing.csv", "calls-review.csv", "calls-summary.csv",
        "calls-by-group.csv", "calls-by-status.csv",
        "calls-by-copy-state.csv", "assignment-alternatives.csv",
        "component-evidence.csv"
    ) %in% names))
    directory <- tempfile("mamamia-bundle-test-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    utils::unzip(archive, exdir = directory)
    root <- file.path(directory, "mamamia-results")
    sample_folder <- unname(bam_bundle_sample_folders(sample))
    sample_root <- file.path(root, sample_folder)
    bams <- utils::read.csv(
        file.path(root, "run", "bams.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_false(any(grepl("Rtmp|/root/", bams$path)))
    indexes <- utils::read.csv(
        file.path(root, "run", "reference-indexes.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_false(any(grepl("Rtmp|/root/", indexes$index_path)))
    calls <- utils::read.csv(
        file.path(sample_root, "calls-full.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(calls), nrow(expected))
    expect_gt(ncol(calls), length(bam_brief_call_columns(expected)))
    dictionary <- bam_bundle_column_dictionary()
    expect_true(all(dictionary$name %in% names(calls)))
    expect_true(all(c(
        "recipient_chromosome", "recipient_mean_depletion",
        "affected_copy_number", "mean_introgression_score",
        "passes_calibration"
    ) %in% names(calls)))
    expect_false("recipient_chr_name" %in% names(calls))
    bedpe <- utils::read.delim(
        file.path(sample_root, "calls-all.bedpe"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true(all(c(
        "donor_contig_id", "donor_start_bp", "donor_end_bp",
        "recipient_contig_id", "recipient_start_bp", "recipient_end_bp",
        "call_class", "passes_calibration"
    ) %in% names(bedpe)))
    expect_true(all(names(bedpe) %in% dictionary$name))
    component_evidence <- utils::read.csv(
        file.path(sample_root, "qc", "component-evidence.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true(all(c(
        "call_number", "chromosome_group", "donor_component",
        "component_rank", "component_selected", "mean_component_score",
        "supported_fraction", "component_supported"
    ) %in% names(component_evidence)))
    confident_csv <- utils::read.csv(
        file.path(sample_root, "calls-confident.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(confident_csv), nrow(confident_calls))
    expect_true(all(confident_csv$passes_calibration))
    expect_true(all(confident_csv$review_status == "confident"))
    confident_bedpe <- utils::read.delim(
        file.path(sample_root, "calls-confident.bedpe"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(confident_bedpe), nrow(confident_calls))
    expect_true(all(confident_bedpe$passes_calibration))
    readme <- paste(
        readLines(file.path(sample_root, "README.txt")),
        collapse = "\n"
    )
    expect_match(readme, sample, fixed = TRUE)
    expect_match(readme, "calls-full.csv", fixed = TRUE)
    expect_match(readme, "Confident set", fixed = TRUE)
    expect_match(readme, "calls-confident.csv", fixed = TRUE)
    expect_match(readme, "Column descriptions", fixed = TRUE)
    expect_match(readme, "recipient_mean_depletion", fixed = TRUE)
    expect_match(readme, "donor_to_recipient_width_ratio", fixed = TRUE)
    expect_match(readme, "| Column | CSV | BEDPE | Description |", fixed = TRUE)
    expect_true(
        file.info(file.path(sample_root, "karyotype-full.svg"))$size > 0
    )
    expect_true(
        file.info(file.path(sample_root, "karyotype-confident.svg"))$size > 0
    )
    expect_true(
        file.info(file.path(sample_root, "calls-all.bedpe"))$size > 0
    )

    relative <- sub("^mamamia-results/", "", listing$Name)
    expect_true(all(c(
        "README.txt", "report.html", "cohort-summary.csv",
        "calls-all.csv", "calls-passing.csv", "calls-confident.csv",
        "calls-summary.csv", "calls-by-status.csv",
        "calls-by-copy-state.csv", "calls-all.bedpe",
        "calls-passing.bedpe", "calls-confident.bedpe",
        "manifest.json", "checksums.md5",
        "qc/masks-summary.csv", "qc/masks-windows.csv.gz",
        "qc/calibration.csv", "qc/plots/masked-fraction-by-group.png",
        "qc/plots/calls-by-class.png",
        "run/provenance.json", "run/config.yml", "run/parameters.csv",
        "run/profile-parameters.csv",
        "run/versions.csv", "run/session-info.txt", "run/assembly.csv",
        "run/reference-indexes.csv", "run/samples.csv", "run/bams.csv",
        "tracks/masks.bed"
    ) %in% relative))
    expect_true(all(paste0(sample_folder, "/", c(
        "README.txt", "report.html", "calls-full.csv",
        "calls-passing.csv", "calls-confident.csv", "calls-review.csv",
        "calls-summary.csv", "calls-by-group.csv", "calls-by-status.csv",
        "calls-by-copy-state.csv", "assignment-alternatives.csv",
        "calls-all.bedpe", "calls-passing.bedpe",
        "calls-confident.bedpe", "karyotype-current.svg",
        "karyotype-full.svg", "karyotype-confident.svg",
        "qc/metrics-summary.csv", "qc/contig-summary.csv",
        "qc/scores-summary.csv", "qc/metrics-windows.csv.gz",
        "qc/scores-windows.csv.gz", "qc/component-evidence.csv",
        "tracks/calls-recipient.bed", "tracks/calls-donor.bed",
        "tracks/breadth.bedGraph.gz"
    )) %in% relative))

    plot_files <- c(
        "qc/plots/coverage-by-group.png",
        "qc/plots/score-by-group.png",
        "qc/plots/score-profile.png",
        "qc/plots/depth-residuals.png",
        "qc/plots/calls-by-class.png"
    )
    expect_true(all(paste0(sample_folder, "/", plot_files) %in% relative))
    top_plot_files <- c(
        "qc/plots/masked-fraction-by-group.png",
        "qc/plots/calls-by-class.png"
    )
    expect_true(all(top_plot_files %in% relative))
    png_signature <- as.raw(c(
        0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a
    ))
    for (plot_file in plot_files) {
        plot_path <- file.path(sample_root, plot_file)
        expect_true(file.info(plot_path)$size > 0)
        expect_identical(
            readBin(plot_path, "raw", n = 8L), png_signature
        )
    }
    for (plot_file in top_plot_files) {
        plot_path <- file.path(root, plot_file)
        expect_true(file.info(plot_path)$size > 0)
        expect_identical(
            readBin(plot_path, "raw", n = 8L), png_signature
        )
    }

    summary <- utils::read.csv(
        file.path(sample_root, "calls-summary.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true(all(c(
        "call_class", "calls", "passing_calls", "confident_calls",
        "needs_review_calls", "total_width_bp",
        "mean_introgression_score"
    ) %in% names(summary)))
    expect_identical(sum(summary$calls), nrow(expected))
    by_group <- utils::read.csv(
        file.path(sample_root, "calls-by-group.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(sum(by_group$calls), nrow(expected))
    passing_csv <- utils::read.csv(
        file.path(sample_root, "calls-passing.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(passing_csv), sum(isTRUE_vec(expected$pass)))
    expect_true(all(passing_csv$passes_calibration))
    review_csv <- utils::read.csv(
        file.path(sample_root, "calls-review.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_false(any(
        review_csv$review_status == "confident" &
            review_csv$passes_calibration
    ))

    metrics <- utils::read.csv(
        gzfile(file.path(sample_root, "qc/metrics-windows.csv.gz")),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(
        nrow(metrics),
        sum(as.character(result$metrics$sample) == sample)
    )
    expect_true(all(c(
        "contig_id", "breadth", "gc_fraction", "high_mapq_read_fraction"
    ) %in% names(metrics)))
    scores <- utils::read.csv(
        gzfile(file.path(sample_root, "qc/scores-windows.csv.gz")),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(
        nrow(scores),
        sum(as.character(result$scores$sample) == sample)
    )
    expect_true(all(c(
        "introgression_score", "donor_residual",
        "recipient_residual", "total_depth_residual"
    ) %in% names(scores)))
    masks <- utils::read.csv(
        gzfile(file.path(root, "qc/masks-windows.csv.gz")),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(masks), nrow(result$masks))
    metrics_summary <- utils::read.csv(
        file.path(sample_root, "qc/metrics-summary.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true(all(c("donor", "recipient") %in% metrics_summary$role))
    expect_true(all(c(
        "mean_breadth", "mean_gc_fraction", "mean_expected_breadth"
    ) %in% names(metrics_summary)))
    scores_summary <- utils::read.csv(
        file.path(sample_root, "qc/scores-summary.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_gt(nrow(scores_summary), 0L)
    expect_true(all(c(
        "mean_introgression_score", "mean_donor_residual",
        "mean_recipient_residual", "mean_total_depth_residual"
    ) %in% names(scores_summary)))
    calibration <- utils::read.csv(
        file.path(root, "qc/calibration.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true("calibration_threshold" %in% names(calibration))

    config <- yaml::read_yaml(file.path(root, "run/config.yml"))
    expect_identical(config$version, result$config$version)
    provenance <- jsonlite::fromJSON(file.path(root, "run/provenance.json"))
    expect_identical(provenance$app, "MaMaMIA")
    expect_identical(provenance$config_hash, result$provenance$config_hash)
    profile_parameters <- utils::read.csv(
        file.path(root, "run/profile-parameters.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true(all(c(
        "profile", "base_profile", "parameter", "value",
        "default", "changed"
    ) %in% names(profile_parameters)))
    expect_true(all(c(
        "windows_per_chromosome", "metric", "score",
        "donor_aggregation", "mapq", "calibration_fdr",
        "calibration_floor_1", "calibration_floor_5",
        "calibration_min_mb", "calibration_merge_gap"
    ) %in% profile_parameters$parameter))
    assembly <- utils::read.csv(
        file.path(root, "run/assembly.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(assembly), 16L)
    expect_true(all(c("recipient", "donor") %in% assembly$role))
    versions <- utils::read.csv(
        file.path(root, "run/versions.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true(all(c("MaMaMIA", "R", "htslib") %in% versions$component))
    bams <- utils::read.csv(
        file.path(root, "run/bams.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true(all(c(
        "file_size", "mtime", "n_targets", "sample_source"
    ) %in% names(bams)))
    indexes <- utils::read.csv(
        file.path(root, "run/reference-indexes.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_gt(nrow(indexes), 0L)
    expect_true(all(c(
        "index_path", "exists", "bytes", "modified"
    ) %in% names(indexes)))

    report <- paste(
        readLines(file.path(sample_root, "report.html")),
        collapse = "\n"
    )
    expect_match(report, "<!DOCTYPE html>", fixed = TRUE)
    expect_match(report, "Run overview", fixed = TRUE)
    expect_match(report, "Calls by class", fixed = TRUE)
    expect_match(report, "Coverage and GC by group", fixed = TRUE)
    expect_match(report, "Score summary by group", fixed = TRUE)
    expect_match(report, "Score profile", fixed = TRUE)
    expect_match(report, "Depth residuals", fixed = TRUE)
    expect_match(report, "Coverage breadth by chromosome group", fixed = TRUE)
    expect_match(report, "data:image/png;base64,", fixed = TRUE)
    expect_match(report, "manifest.json", fixed = TRUE)
    expect_match(report, "qc/metrics-summary.csv", fixed = TRUE)

    top_report <- paste(
        readLines(file.path(root, "report.html")),
        collapse = "\n"
    )
    expect_match(top_report, "Run overview", fixed = TRUE)
    expect_match(top_report, "Sample summary", fixed = TRUE)
    expect_match(top_report, sample_folder, fixed = TRUE)
    combined <- utils::read.csv(
        file.path(root, "calls-all.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(combined), nrow(expected))
    cohort <- utils::read.csv(
        file.path(root, "cohort-summary.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(cohort), 1L)
    expect_identical(as.character(cohort$sample), sample)

    manifest <- jsonlite::fromJSON(file.path(root, "manifest.json"))
    expect_identical(manifest$sample, sample)
    expect_identical(manifest$samples, sample)
    expect_gt(nrow(manifest$files), 25L)
    expect_true(all(file.exists(
        file.path(root, manifest$files$path)
    )))
    expect_identical(
        unname(tools::md5sum(file.path(root, manifest$files$path))),
        manifest$files$md5
    )
    checksums <- readLines(file.path(root, "checksums.md5"))
    expect_true(any(grepl("manifest.json", checksums, fixed = TRUE)))
    expect_true(any(grepl("report.html", checksums, fixed = TRUE)))
    first_checksum <- strsplit(checksums[[1L]], "  ", fixed = TRUE)[[1L]]
    expect_identical(
        unname(tools::md5sum(file.path(root, first_checksum[[2L]]))),
        first_checksum[[1L]]
    )

    recipient_bed <- readLines(
        file.path(sample_root, "tracks/calls-recipient.bed")
    )
    expect_match(recipient_bed[[1L]], "track name=", fixed = TRUE)
    expect_gt(length(recipient_bed), 1L)
    breadth_lines <- readLines(
        gzfile(file.path(sample_root, "tracks/breadth.bedGraph.gz"))
    )
    expect_match(breadth_lines[[1L]], "track type=bedGraph", fixed = TRUE)
    expect_gt(length(breadth_lines), 1L)

    expect_match(readme, "report.html", fixed = TRUE)
    expect_match(readme, "qc/metrics-windows.csv.gz", fixed = TRUE)
    expect_match(readme, "qc/plots/", fixed = TRUE)
    expect_match(readme, "manifest.json", fixed = TRUE)
    expect_match(readme, "Diagnostic conventions", fixed = TRUE)
})

test_that("diagnostic plots handle empty inputs and multi-sample cohorts", {
    expect_null(bam_bundle_plot_coverage(NULL))
    expect_null(bam_bundle_plot_scores(data.frame()))
    expect_null(bam_bundle_plot_score_profile(data.frame()))
    expect_null(bam_bundle_plot_residuals(data.frame()))
    expect_null(bam_bundle_plot_masks(data.frame()))
    expect_null(bam_bundle_plot_calls(data.frame()))

    calls <- data.frame(
        call_class = c("addition", "deletion"),
        calls = c(2L, 1L),
        passing_calls = c(1L, 1L),
        confident_calls = c(1L, 0L),
        stringsAsFactors = FALSE
    )
    expect_s3_class(bam_bundle_plot_calls(calls), "ggplot")

    cohort <- data.frame(
        sample = c("A", "B"),
        n_all_calls = c(3L, 5L),
        n_class_confident = c(1L, 2L),
        stringsAsFactors = FALSE
    )
    expect_s3_class(bam_bundle_plot_cohort(cohort), "ggplot")
    expect_null(bam_bundle_plot_cohort(cohort[1L, , drop = FALSE]))
})

test_that("SVG numbers ignore the locale decimal mark", {
    old <- options(OutDec = ",")
    on.exit(options(old), add = TRUE)
    expect_identical(bam_svg_number(1.5), "1.50")
})

test_that("a relative bundle path is written to the caller's directory", {
    skip_if_not(bam_example_available(), "bundled example is not installed")
    result <- bam_example_result(1000)
    sample <- unique(as.character(result$calls$sample))[[1L]]
    directory <- tempfile("mamamia-relative-")
    dir.create(directory)
    old <- setwd(directory)
    on.exit({
        setwd(old)
        unlink(directory, recursive = TRUE)
    }, add = TRUE)
    bam_write_result_bundle(result, "relative-bundle.zip", sample = sample)
    expect_true(file.exists(file.path(directory, "relative-bundle.zip")))
    unpacked <- file.path(directory, "unpacked")
    dir.create(unpacked)
    utils::unzip(
        file.path(directory, "relative-bundle.zip"), exdir = unpacked
    )
    root <- file.path(unpacked, "mamamia-results")
    sample_root <- file.path(
        root, unname(bam_bundle_sample_folders(sample))
    )
    readme <- readLines(file.path(sample_root, "README.txt"))
    expect_true(any(grepl("karyotype-full.svg", readme, fixed = TRUE)))
    expect_true(any(grepl("report.html", readme, fixed = TRUE)))
})

test_that("a cohort bundle has per-sample folders and top-level summaries", {
    skip_if(inherits(try(bam_bundle_zip_tool(), silent = TRUE),
        "try-error"), "the 7z command is unavailable")
    result <- bam_example_result(1000)
    samples <- c("cohort A", "cohort B")
    multi <- result
    multi$calls <- do.call(rbind, lapply(samples, function(sample) {
        calls <- result$calls
        calls$sample <- sample
        calls
    }))
    multi$metrics <- do.call(rbind, lapply(samples, function(sample) {
        metrics <- result$metrics
        metrics$sample <- sample
        metrics
    }))
    multi$scores <- do.call(rbind, lapply(samples, function(sample) {
        scores <- result$scores
        scores$sample <- sample
        scores
    }))
    multi$samples <- do.call(rbind, lapply(samples, function(sample) {
        meta <- result$samples[1L, , drop = FALSE]
        meta$sample <- sample
        meta$path <- paste0(sample, ".bam")
        meta
    }))
    multi$bams <- do.call(rbind, lapply(samples, function(sample) {
        bam <- result$bams[1L, , drop = FALSE]
        bam$sample <- sample
        bam$path <- paste0(sample, ".bam")
        bam
    }))
    class(multi) <- class(result)

    archive <- tempfile(fileext = ".zip")
    on.exit(unlink(archive), add = TRUE)
    bam_write_cohort_bundle(
        multi, archive,
        generated = as.POSIXct("2026-01-01 00:00:00", tz = "UTC")
    )
    expect_true(file.exists(archive))
    listing <- utils::unzip(archive, list = TRUE)
    relative <- sub("^mamamia-results/", "", listing$Name)
    expect_true(all(c(
        "README.txt", "report.html", "cohort-summary.csv",
        "calls-all.csv", "calls-passing.csv", "calls-confident.csv",
        "calls-summary.csv", "calls-by-status.csv",
        "calls-by-copy-state.csv", "calls-all.bedpe",
        "calls-passing.bedpe", "calls-confident.bedpe",
        "manifest.json", "checksums.md5", "qc/calibration.csv",
        "qc/masks-summary.csv", "qc/masks-windows.csv.gz",
        "qc/plots/cohort-calls.png", "qc/plots/calls-by-class.png",
        "qc/plots/masked-fraction-by-group.png", "tracks/masks.bed",
        "run/provenance.json", "run/samples.csv", "run/bams.csv"
    ) %in% relative))

    folders <- unname(bam_bundle_sample_folders(samples))
    sample_files <- c(
        "calls-full.csv", "calls-passing.csv", "calls-confident.csv",
        "calls-review.csv", "calls-summary.csv", "calls-all.bedpe",
        "calls-passing.bedpe", "calls-confident.bedpe",
        "karyotype-current.svg", "karyotype-full.svg",
        "karyotype-confident.svg", "report.html", "README.txt",
        "qc/metrics-summary.csv", "qc/scores-summary.csv",
        "qc/metrics-windows.csv.gz", "qc/scores-windows.csv.gz",
        "qc/component-evidence.csv", "tracks/calls-recipient.bed",
        "tracks/calls-donor.bed", "tracks/breadth.bedGraph.gz"
    )
    for (i in seq_along(samples)) {
        prefix <- paste0(folders[[i]], "/")
        expect_true(all(paste0(prefix, sample_files) %in% relative))
        expect_false(any(grepl(
            paste0("^", prefix, "(run/|qc/calibration\\.csv|qc/masks-)"),
            relative
        )))
    }

    directory <- tempfile("mamamia-cohort-test-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    utils::unzip(archive, exdir = directory)
    root <- file.path(directory, "mamamia-results")
    cohort <- utils::read.csv(
        file.path(root, "cohort-summary.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(cohort), length(samples))
    expect_setequal(as.character(cohort$sample), samples)
    calls <- utils::read.csv(
        file.path(root, "calls-all.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(calls), nrow(result$calls) * length(samples))
    expect_setequal(unique(as.character(calls$sample_name)), samples)
    bedpe <- utils::read.delim(
        file.path(root, "calls-all.bedpe"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(bedpe), nrow(result$calls) * length(samples))
    expect_true("sample_name" %in% names(bedpe))

    for (i in seq_along(samples)) {
        sample_calls <- utils::read.csv(
            file.path(root, folders[[i]], "calls-full.csv"),
            check.names = FALSE, stringsAsFactors = FALSE
        )
        expect_identical(nrow(sample_calls), nrow(result$calls))
        expect_true(all(as.character(sample_calls$sample_name) == samples[[i]]))
        sample_readme <- paste(
            readLines(file.path(
                root, folders[[i]], "README.txt"
            )),
            collapse = "\n"
        )
        expect_match(sample_readme, "Column descriptions", fixed = TRUE)
        expect_false(grepl("qc/calibration.csv", sample_readme, fixed = TRUE))
        expect_false(grepl("run/provenance.json", sample_readme, fixed = TRUE))
    }

    readme <- paste(
        readLines(file.path(root, "README.txt")),
        collapse = "\n"
    )
    expect_match(readme, "Per-sample folders", fixed = TRUE)
    expect_match(readme, "cohort-summary.csv", fixed = TRUE)
    expect_match(readme, "cohort-A/", fixed = TRUE)
    report <- paste(
        readLines(file.path(root, "report.html")),
        collapse = "\n"
    )
    expect_match(report, "Cohort results", fixed = TRUE)
    expect_match(report, "Per-sample bundles", fixed = TRUE)
    expect_match(report, "cohort-A/report.html", fixed = TRUE)

    manifest <- jsonlite::fromJSON(file.path(root, "manifest.json"))
    expect_identical(manifest$sample, "cohort")
    expect_identical(manifest$n_samples, length(samples))
    expect_setequal(manifest$samples, samples)
    expect_true(all(file.exists(file.path(root, manifest$files$path))))
})

test_that("bundle workers follow the explicit option first", {
    withr::local_options(
        MaMaMIA.bundle_workers = 3L,
        MaMaMIA.cores = 2L
    )
    expect_identical(bam_bundle_workers(), 3L)

    withr::local_options(MaMaMIA.bundle_workers = NULL)
    expect_identical(bam_bundle_workers(), 2L)
})

test_that("bundle workers default to the CPU and memory budget", {
    withr::local_options(
        MaMaMIA.bundle_workers = NULL,
        MaMaMIA.cores = NULL
    )
    expect_identical(
        bam_bundle_workers(),
        as.integer(min(
            mamamia_default_cores(),
            mamamia_memory_worker_cap()
        ))
    )
})

test_that("the archive uses 7-Zip at maximum compression and the worker budget", {
    args <- bam_bundle_archive_args("bundle.zip", "results", workers = 11L)
    expect_true("-tzip" %in% args)
    expect_true("-mx=9" %in% args)
    expect_true("-mmt=11" %in% args)
})
