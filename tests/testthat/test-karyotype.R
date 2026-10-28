## Interactive karyotype payload and widget.

test_that("karyotype_payload() lays out chromosomes by subgenome", {
    result <- bam_example_result(100)
    payload <- karyotype_payload(result, "example-sample", "heatmap")
    expect_identical(payload$mode, "heatmap")
    expect_identical(payload$sample, "example-sample")
    expect_true(is.finite(payload$score_limit))
    expect_length(payload$subgenomes, 1L)
    subgenome <- payload$subgenomes[[1L]]
    expect_identical(subgenome$subgenome, "A")
    expect_identical(subgenome$label, "Subgenome A")
    expect_true(subgenome$selected)
    expect_length(subgenome$chromosomes, 1L)
    chromosome <- subgenome$chromosomes[[1L]]
    expect_identical(chromosome$name, "1A")
    expect_identical(chromosome$subgenome, "A")
    expect_false(chromosome$deleted)
    expect_false(chromosome$added)
    expect_identical(chromosome$length, 20000)
    expect_length(chromosome$scores, 100L)
    hosted <- sum(result$calls$recipient_assignment != "no_recipient_loss",
        na.rm = TRUE
    )
    expect_length(chromosome$calls, hosted)
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
        vapply(payload$subgenomes, function(x) x$subgenome, character(1)),
        c("A", "R")
    )
    added <- payload$subgenomes[[2L]]$chromosomes
    expect_length(added, 2L)
    expect_identical(
        vapply(added, function(x) x$name, character(1)),
        c("1R", "1R")
    )
    expect_identical(
        vapply(added, function(x) x$label, character(1)),
        c("1R 1", "1R 2")
    )
    expect_identical(
        vapply(added, function(x) x$copy, integer(1)),
        1:2
    )
    expect_identical(
        vapply(added, function(x) x$copies, integer(1)),
        c(2L, 2L)
    )
    expect_true(all(vapply(added, function(x) x$added, logical(1))))
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
        vapply(payload$subgenomes, function(x) x$subgenome, character(1)),
        "A"
    )
})

test_that("a deleted recipient chromosome stays in the karyotype", {
    result <- bam_example_result(100)
    extra <- result$assembly$contigs$recipient[1, , drop = FALSE]
    extra$id <- "r2"
    extra$name <- "1B"
    extra$subgenome <- "B"
    result$assembly$contigs$recipient <- rbind(
        result$assembly$contigs$recipient, extra
    )
    metric <- result$metrics[result$metrics$role == "recipient", ,
        drop = FALSE
    ][1, , drop = FALSE]
    metric$chr_id <- "r2"
    metric$chr_name <- "1B"
    metric$subgenome <- "B"
    metric$breadth <- 0.02
    result$metrics <- rbind(result$metrics, metric)
    result$metrics$breadth[result$metrics$chr_id == "r1"] <- 0.0009

    payload <- karyotype_payload(result, "example-sample", "heatmap")
    expect_identical(
        vapply(payload$subgenomes, function(x) x$subgenome, character(1)),
        c("A", "B")
    )
    recipient <- payload$subgenomes[[1L]]$chromosomes
    expect_length(recipient, 1L)
    expect_identical(recipient[[1L]]$name, "1A")
    expect_true(recipient[[1L]]$deleted)
    expect_false(recipient[[1L]]$added)
})

test_that("a cross-group assignment is painted on the recipient chromosome", {
    result <- bam_example_result(100)
    extra <- result$assembly$contigs$recipient[1, , drop = FALSE]
    extra$id <- "r2"
    extra$name <- "2A"
    extra$subgenome <- "A"
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
    expect_length(payload$subgenomes, 1L)
    chromosomes <- payload$subgenomes[[1L]]$chromosomes
    expect_identical(
        vapply(chromosomes, function(x) x$name, character(1)),
        c("1A", "2A")
    )
    expect_length(chromosomes[[1L]]$calls, 0L)
    expect_length(chromosomes[[2L]]$calls, 1L)
    moved <- chromosomes[[2L]]$calls[[1L]]
    expect_true(moved$cross_group)
    expect_identical(moved$donor_group, "1")
    expect_identical(moved$recipient_group, "2")
    expect_identical(moved$recipient_chr_name, "2A")
    expect_equal(moved$start, 0)
    expect_equal(moved$end, 1)
})

test_that("the karyotype includes every recipient subgenome", {
    fixture <- multi_subgenome_fixture(samples = "sample1")
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
        vapply(payload$subgenomes, function(x) x$subgenome, character(1)),
        c("A", "B")
    )
})

test_that("mamamia_karyotype() returns a widget with its stylesheet", {
    skip_if_not_installed("htmlwidgets")
    skip_if_not_installed("htmltools")
    result <- bam_example_result(100)
    widget <- mamamia_karyotype(karyotype_payload(result))
    expect_s3_class(widget, "mamamiaKaryotype")
    expect_s3_class(widget, "htmlwidget")
    expect_true(any(vapply(
        widget$dependencies,
        function(dependency) {
            identical(dependency$name, "mamamiaKaryotype-css") &&
                identical(dependency$stylesheet, "mamamiaKaryotype.css")
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
