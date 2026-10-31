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
            "recipient_chr_name", "recipient_assignment",
            "start_bp", "end_bp", "width_bp", "pass"
        )
    )
    one <- data.frame(sample = c("A", "A"), stringsAsFactors = FALSE)
    expect_identical(
        bam_brief_call_columns(one),
        c(
            "group", "class", "copy_state", "class_status",
            "recipient_chr_name", "recipient_assignment",
            "start_bp", "end_bp", "width_bp", "pass"
        )
    )
    labels <- bam_brief_call_labels(bam_brief_call_columns(one))
    expect_identical(labels[[1L]], "Group")
    expect_identical(labels[[5L]], "Recipient")
    expect_identical(labels[[10L]], "Pass")
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

test_that("the bundle filename is safe", {
    expect_identical(
        bam_bundle_filename("Sq 74/7D"), "mamamia-Sq-74-7D-results.zip"
    )
    expect_identical(bam_bundle_filename(NULL), "mamamia-sample-results.zip")
    expect_identical(bam_bundle_filename(""), "mamamia-sample-results.zip")
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
    skip_if_not(nzchar(Sys.which("zip")), "the zip command is unavailable")
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
        "karyotype-full.svg", "karyotype-confident.svg"
    ) %in% names))
    directory <- tempfile("mamamia-bundle-test-")
    dir.create(directory)
    on.exit(unlink(directory, recursive = TRUE), add = TRUE)
    utils::unzip(archive, exdir = directory)
    root <- file.path(
        directory, paste0("mamamia-", sample, "-results")
    )
    calls <- utils::read.csv(
        file.path(root, "calls-full.csv"),
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
        file.path(root, "calls-all.bedpe"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_true(all(c(
        "donor_contig_id", "donor_start_bp", "donor_end_bp",
        "recipient_contig_id", "recipient_start_bp", "recipient_end_bp",
        "call_class", "passes_calibration"
    ) %in% names(bedpe)))
    expect_true(all(names(bedpe) %in% dictionary$name))
    confident_csv <- utils::read.csv(
        file.path(root, "calls-confident.csv"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(confident_csv), nrow(confident_calls))
    expect_true(all(confident_csv$passes_calibration))
    expect_true(all(confident_csv$review_status == "confident"))
    confident_bedpe <- utils::read.delim(
        file.path(root, "calls-confident.bedpe"),
        check.names = FALSE, stringsAsFactors = FALSE
    )
    expect_identical(nrow(confident_bedpe), nrow(confident_calls))
    expect_true(all(confident_bedpe$passes_calibration))
    readme <- paste(
        readLines(file.path(root, "README.txt")),
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
    expect_true(file.info(file.path(root, "karyotype-full.svg"))$size > 0)
    expect_true(
        file.info(file.path(root, "karyotype-confident.svg"))$size > 0
    )
    expect_true(file.info(file.path(root, "calls-all.bedpe"))$size > 0)
})
