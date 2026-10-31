## Interactive karyotype payload and widget.

test_that("karyotype_payload() lays out diploid pairs in config order", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "heatmap")
    expect_identical(payload$mode, "heatmap")
    expect_identical(payload$sample, "example-sample")
    expect_true(is.finite(payload$score_limit))
    expect_length(payload$chromosomes, 8L)
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "2A", "2A", "3A", "3A", "4A", "4A")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$label, character(1)),
        c("1A", "1A", "2A", "2A", "3A", "3A", "4A", "4A")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copy, integer(1)),
        c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L)
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copies, integer(1)),
        rep(2L, 8L)
    )
    chromosome <- payload$chromosomes[[1L]]
    expect_identical(chromosome$component, "A")
    expect_false(chromosome$deleted)
    expect_false(chromosome$added)
    expect_identical(chromosome$length, 160000)
    expect_length(chromosome$scores, 1000L)
    hosted <- sum(
        result$calls$recipient_assignment != "no_recipient_loss" &
            as.character(result$calls$group) == "1",
        na.rm = TRUE
    )
    expect_length(chromosome$calls, hosted)
    expect_length(payload$chromosomes[[2L]]$calls, hosted)
    expect_equal(hosted, 1L)
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
    result <- bam_example_result(1000)
    addition <- result$calls[result$calls$class == "substitution" &
        result$calls$pass %in% TRUE, , drop = FALSE][1L, , drop = FALSE]
    addition$start_iid <- 0L
    addition$end_iid <- 999L
    addition$start_bp <- 0
    addition$end_bp <- 160000
    addition$width_bp <- 160000
    addition$score_mean <- 2
    addition$donor_mean <- 2
    addition$recipient_mean <- 0
    addition$class <- "addition"
    addition$confidence <- "high"
    addition$pass <- TRUE
    addition$n_windows <- 1000L
    result$calls <- rbind(result$calls, addition)

    payload <- karyotype_payload(result, "example-sample", "calls")
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "2A", "2A", "3A", "3A", "4A", "4A",
            "1R", "1R")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$label, character(1)),
        c("1A", "1A", "2A", "2A", "3A", "3A", "4A", "4A",
            "1R", "1R")
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$added, logical(1)),
        c(rep(FALSE, 8L), TRUE, TRUE)
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copy, integer(1)),
        c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L)
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copies, integer(1)),
        rep(2L, 10L)
    )
})

test_that("a partial donor call is not a whole-chromosome addition", {
    result <- bam_example_result(1000)
    addition <- result$calls[result$calls$class == "substitution" &
        result$calls$pass %in% TRUE, , drop = FALSE][1L, , drop = FALSE]
    addition$start_iid <- 0L
    addition$end_iid <- 500L
    addition$start_bp <- 0
    addition$end_bp <- 80000
    addition$width_bp <- 80000
    addition$score_mean <- 2
    addition$donor_mean <- 2
    addition$recipient_mean <- 0
    addition$class <- "addition"
    addition$confidence <- "high"
    addition$pass <- TRUE
    addition$n_windows <- 501L
    result$calls <- rbind(result$calls, addition)

    payload <- karyotype_payload(result, "example-sample", "calls")
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "2A", "2A", "3A", "3A", "4A", "4A")
    )
})

test_that("a deleted recipient chromosome stays in the karyotype", {
    result <- bam_example_result(1000)
    extra <- result$assembly$contigs$recipient[1, , drop = FALSE]
    extra$id <- "r5"
    extra$name <- "1B"
    extra$component <- "B"
    result$assembly$contigs$recipient <- rbind(
        result$assembly$contigs$recipient, extra
    )
    metric <- result$metrics[result$metrics$role == "recipient", ,
        drop = FALSE
    ][1, , drop = FALSE]
    metric$chr_id <- "r5"
    metric$chr_name <- "1B"
    metric$component <- "B"
    metric$breadth <- 0.02
    result$metrics <- rbind(result$metrics, metric)
    result$metrics$breadth[result$metrics$chr_id == "r1"] <- 0.0009

    payload <- karyotype_payload(result, "example-sample", "heatmap")
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "2A", "2A", "3A", "3A", "4A", "4A", "1B", "1B")
    )
    expect_true(payload$chromosomes[[1L]]$deleted)
    expect_true(payload$chromosomes[[2L]]$deleted)
    expect_false(payload$chromosomes[[3L]]$deleted)
    expect_false(payload$chromosomes[[4L]]$deleted)
    expect_false(payload$chromosomes[[9L]]$deleted)
    expect_false(payload$chromosomes[[10L]]$deleted)
})

test_that("a cross-group assignment is painted on the recipient chromosome", {
    result <- bam_example_result(1000)
    call <- result$calls[result$calls$class == "substitution" &
        result$calls$pass %in% TRUE, , drop = FALSE][1L, , drop = FALSE]
    call$group <- "1"
    call$recipient_group <- "2"
    call$recipient_cross_group <- TRUE
    call$recipient_chr_id <- "r2"
    call$recipient_chr_name <- "2A"
    call$recipient_start_bp <- 0
    call$recipient_end_bp <- 220000
    result$calls <- call
    payload <- karyotype_payload(result, "example-sample", "calls")
    chromosomes <- payload$chromosomes
    expect_identical(
        vapply(chromosomes, function(x) x$name, character(1)),
        c("1A", "1A", "2A", "2A", "3A", "3A", "4A", "4A")
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
    result <- bam_example_result(1000)
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
    expect_identical(binding[[1L]]$version, "0.1.0")
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
    result <- bam_example_result(1000)
    base <- result$calls[result$calls$class == "substitution" &
        result$calls$pass %in% TRUE, , drop = FALSE][1L, , drop = FALSE]
    make_addition <- function(group, donor_id, donor_name, length_bp) {
        call <- base
        call$group <- as.character(group)
        call$start_iid <- 0L
        call$end_iid <- 999L
        call$start_bp <- 0
        call$end_bp <- length_bp
        call$width_bp <- length_bp
        call$donor_mean <- 2
        call$class <- "addition"
        call$pass <- TRUE
        call$n_windows <- 1000L
        call$donor_chr_id <- donor_id
        call$donor_chr_name <- donor_name
        call
    }
    result$calls <- rbind(
        make_addition(1L, "d1", "1R", 160000),
        make_addition(3L, "d3", "3R", 180000)
    )

    payload <- karyotype_payload(result, "example-sample", "calls")
    added <- vapply(payload$chromosomes, function(x) isTRUE(x$added), logical(1))
    expect_identical(
        vapply(payload$chromosomes[added], function(x) x$name, character(1)),
        c("1R", "1R", "3R", "3R")
    )
})
