## Interactive karyotype payload and widget.

test_that("karyotype_payload() lays out diploid pairs in config order", {
    result <- bam_example_result(100)
    payload <- karyotype_payload(result, "example-sample", "heatmap")
    expect_identical(payload$mode, "heatmap")
    expect_identical(payload$sample, "example-sample")
    expect_true(is.finite(payload$score_limit))
    expect_length(payload$chromosomes, 2L)
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$label, character(1)),
        c("1A", "1A")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copy, integer(1)),
        c(1L, 2L)
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copies, integer(1)),
        c(2L, 2L)
    )
    chromosome <- payload$chromosomes[[1L]]
    expect_identical(chromosome$component, "A")
    expect_false(chromosome$deleted)
    expect_false(chromosome$added)
    expect_identical(chromosome$length, 20000)
    expect_length(chromosome$scores, 100L)
    hosted <- sum(result$calls$recipient_assignment != "no_recipient_loss",
        na.rm = TRUE
    )
    expect_length(chromosome$calls, hosted)
    expect_length(payload$chromosomes[[2L]]$calls, hosted)
    if (hosted > 0L) {
        expect_true(all(c(
            "affected_copies", "recipient_copy_number",
            "affected_fraction", "dosage_confidence",
            "class_group", "donor_mean", "recipient_mean",
            "donor_start_bp", "donor_end_bp",
            "recipient_run_start_bp", "recipient_run_end_bp",
            "recipient_deletion_width_bp", "recipient_margin",
            "recipient_fraction_ratio"
        ) %in% names(chromosome$calls[[1L]])))
        expect_true(is.finite(chromosome$calls[[1L]]$donor_start_bp))
        expect_true(is.finite(chromosome$calls[[1L]]$donor_end_bp))
    }
})

test_that("a full donor chromosome is added with its diploid copies", {
    result <- bam_example_result(100)
    addition <- result$calls[2, , drop = FALSE]
    addition$start_iid <- 0L
    addition$end_iid <- 99L
    addition$start_bp <- 0
    addition$end_bp <- 20000
    addition$width_bp <- 20000
    addition$score_mean <- 2
    addition$donor_mean <- 2
    addition$recipient_mean <- 0
    addition$class <- "addition"
    addition$confidence <- "high"
    addition$pass <- TRUE
    addition$n_windows <- 100L
    result$calls <- rbind(result$calls, addition)

    payload <- karyotype_payload(result, "example-sample", "calls")
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "1R", "1R")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$label, character(1)),
        c("1A", "1A", "1R", "1R")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$added, logical(1)),
        c(FALSE, FALSE, TRUE, TRUE)
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copy, integer(1)),
        c(1L, 2L, 1L, 2L)
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copies, integer(1)),
        c(2L, 2L, 2L, 2L)
    )
})

test_that("a partial donor call is not a whole-chromosome addition", {
    result <- bam_example_result(100)
    addition <- result$calls[2, , drop = FALSE]
    addition$start_iid <- 0L
    addition$end_iid <- 50L
    addition$start_bp <- 0
    addition$end_bp <- 10000
    addition$width_bp <- 10000
    addition$score_mean <- 2
    addition$donor_mean <- 2
    addition$recipient_mean <- 0
    addition$class <- "addition"
    addition$confidence <- "high"
    addition$pass <- TRUE
    addition$n_windows <- 51L
    result$calls <- rbind(result$calls, addition)

    payload <- karyotype_payload(result, "example-sample", "calls")
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A")
    )
})

test_that("a deleted recipient chromosome stays in the karyotype", {
    result <- bam_example_result(100)
    extra <- result$assembly$contigs$recipient[1, , drop = FALSE]
    extra$id <- "r2"
    extra$name <- "1B"
    extra$component <- "B"
    result$assembly$contigs$recipient <- rbind(
        result$assembly$contigs$recipient, extra
    )
    metric <- result$metrics[result$metrics$role == "recipient", ,
        drop = FALSE
    ][1, , drop = FALSE]
    metric$chr_id <- "r2"
    metric$chr_name <- "1B"
    metric$component <- "B"
    metric$breadth <- 0.02
    result$metrics <- rbind(result$metrics, metric)
    result$metrics$breadth[result$metrics$chr_id == "r1"] <- 0.0009

    payload <- karyotype_payload(result, "example-sample", "heatmap")
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "1B", "1B")
    )
    expect_true(payload$chromosomes[[1L]]$deleted)
    expect_true(payload$chromosomes[[2L]]$deleted)
    expect_false(payload$chromosomes[[3L]]$deleted)
    expect_false(payload$chromosomes[[4L]]$deleted)
})

test_that("a cross-group assignment is painted on the recipient chromosome", {
    result <- bam_example_result(100)
    extra <- result$assembly$contigs$recipient[1, , drop = FALSE]
    extra$id <- "r2"
    extra$name <- "2A"
    extra$component <- "A"
    extra$group <- "2"
    extra$length <- 40000
    result$assembly$contigs$recipient <- rbind(
        result$assembly$contigs$recipient, extra
    )
    call <- result$calls[2, , drop = FALSE]
    call$group <- "1"
    call$recipient_group <- "2"
    call$recipient_cross_group <- TRUE
    call$recipient_chr_id <- "r2"
    call$recipient_chr_name <- "2A"
    call$recipient_start_bp <- 0
    call$recipient_end_bp <- 40000
    result$calls <- call
    payload <- karyotype_payload(result, "example-sample", "calls")
    chromosomes <- payload$chromosomes
    expect_identical(
        vapply(chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "2A", "2A")
    )
    expect_length(chromosomes[[1L]]$calls, 0L)
    expect_length(chromosomes[[2L]]$calls, 0L)
    expect_length(chromosomes[[3L]]$calls, 1L)
    expect_length(chromosomes[[4L]]$calls, 1L)
    moved <- chromosomes[[3L]]$calls[[1L]]
    expect_true(moved$cross_group)
    expect_identical(moved$donor_group, "1")
    expect_identical(moved$recipient_group, "2")
    expect_identical(moved$recipient_chr_name, "2A")
    expect_equal(moved$start, 0)
    expect_equal(moved$end, 1)
})

test_that("the karyotype includes every recipient component", {
    fixture <- multi_component_fixture(samples = "sample1")
    on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
    config <- readMaMaMIAConfig(fixture$config_path)
    resolved <- bam_resolve_species_pair(config, "recipient", "donor")
    result <- analyzeBamInput(
        fixture$bam_path, resolved$config,
        assembly = resolved$key,
        windows_per_chromosome = 100L,
        cores = 1L
    )
    payload <- karyotype_payload(result, "sample1", "heatmap")
    expect_identical(
        unique(vapply(payload$chromosomes, function(x) x$component,
            character(1)
        )),
        c("A", "B")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "1B", "1B")
    )
})

test_that("mamamia_karyotype() returns a widget with its stylesheet", {
    skip_if_not_installed("htmlwidgets")
    skip_if_not_installed("htmltools")
    result <- bam_example_result(100)
    widget <- mamamia_karyotype(karyotype_payload(result))
    expect_s3_class(widget, "mamamiaKaryotypeCard")
    dependencies <- htmltools::resolveDependencies(widget$dependencies)
    expect_silent(dependencies)
    expect_s3_class(widget, "htmlwidget")
    binding <- Filter(
        function(dependency) {
            identical(
                dependency$name, "mamamiaKaryotypeCard-binding"
            )
        },
        dependencies
    )
    expect_length(binding, 1L)
    expect_identical(binding[[1L]]$version, "0.99.2.2")
    expect_true(any(vapply(
        widget$dependencies,
        function(dependency) {
            identical(dependency$name, "mamamiaKaryotypeCard-css") &&
                identical(dependency$stylesheet, "mamamiaKaryotypeCard.css")
        },
        logical(1)
    )))
})

test_that("the explorer exposes the karyotype tab and mode switch", {
    skip_if_not_installed("htmlwidgets")
    skip_if_not_installed("htmltools")
    html <- as.character(explorer_results_ui("explorer"))
    expect_match(html, "Karyotype", fixed = TRUE)
    expect_match(html, "explorer-karyotype_mode", fixed = TRUE)
    expect_match(html, "explorer-karyotype", fixed = TRUE)
})

test_that("added donor chromosomes preserve registry order", {
    result <- bam_example_result(100)
    donor <- result$assembly$contigs$donor
    extra <- donor[1, , drop = FALSE]
    extra$id <- "d2"
    extra$name <- "2R"
    extra$group <- "2"
    extra$component <- "R"
    extra$length <- 20000
    extra$copy_number <- 2
    result$assembly$contigs$donor <- rbind(extra, donor)

    scores2 <- result$scores[
        as.character(result$scores$group) == "1", , drop = FALSE
    ]
    scores2$group <- "2"
    scores2$chr_name <- "2R"
    result$scores <- rbind(result$scores, scores2)

    base <- result$calls[1, , drop = FALSE]
    first <- base
    first$group <- "1"
    first$start_iid <- 0L
    first$end_iid <- 99L
    first$width_bp <- 20000
    first$donor_mean <- 2
    first$pass <- TRUE
    first$donor_chr_id <- donor$id[[1L]]
    first$donor_chr_name <- donor$name[[1L]]
    second <- first
    second$group <- "2"
    second$donor_chr_id <- "d2"
    second$donor_chr_name <- "2R"
    result$calls <- rbind(first, second)

    payload <- karyotype_payload(result, "example-sample", "calls")
    added <- vapply(payload$chromosomes, function(x) isTRUE(x$added), logical(1))
    expect_identical(
        vapply(payload$chromosomes[added], function(x) x$name, character(1)),
        c("2R", "2R", "1R", "1R")
    )
})
