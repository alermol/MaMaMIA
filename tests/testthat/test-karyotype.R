## Interactive karyotype payload and widget.

test_that("karyotype_payload() lays out one homolog per chromosome", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "heatmap")
    expect_identical(payload$mode, "heatmap")
    expect_identical(payload$sample, "example-sample")
    expect_true(is.finite(payload$score_limit))
    expected_names <- c(
        "1A", "2A", "2R", "3A", "4A", "5A", "5R",
        "6A", "7A", "8A", "8R"
    )
    expect_length(payload$chromosomes, 11L)
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        expected_names
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$label, character(1)),
        expected_names
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copy, integer(1)),
        rep(1L, 11L)
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copies, integer(1)),
        rep(1L, 11L)
    )
    expect_true(any(vapply(
        payload$chromosomes, function(x) isTRUE(x$added), logical(1)
    )))
    expect_true(any(vapply(
        payload$chromosomes, function(x) isTRUE(x$deleted), logical(1)
    )))
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
    expect_equal(hosted, 2L)
    expect_setequal(
        vapply(chromosome$calls, function(call) call$class, character(1)),
        c("substitution", "addition")
    )
    expect_identical(payload$chromosomes[[2L]]$name, "2A")
    expect_length(payload$chromosomes[[2L]]$calls, 0L)
    expect_identical(payload$chromosomes[[3L]]$name, "2R")
    expect_true(isTRUE(payload$chromosomes[[3L]]$added_partial))
    if (hosted > 0L) {
        expect_true(all(c(
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

test_that("the bundled example covers every karyotype rendering state", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "calls")
    chromosomes <- payload$chromosomes
    added <- Filter(function(x) isTRUE(x$added), chromosomes)
    recipient <- Filter(function(x) !isTRUE(x$added), chromosomes)
    added_classes <- function(x) {
        vapply(x$calls, function(call) call$class, character(1))
    }
    recipient_classes <- unlist(lapply(recipient, added_classes))

    expect_true(any(vapply(
        chromosomes, function(x) isTRUE(x$deleted), logical(1)
    )))
    expect_true(any(recipient_classes == "deletion"))
    expect_true(any(vapply(added, function(x) {
        !isTRUE(x$added_partial) && "addition" %in% added_classes(x)
    }, logical(1))))
    expect_true(any(vapply(added, function(x) {
        isTRUE(x$added_partial) && "addition" %in% added_classes(x)
    }, logical(1))))
    expect_true(any(vapply(added, function(x) {
        !isTRUE(x$added_partial) && "substitution" %in% added_classes(x)
    }, logical(1))))
    expect_true(any(recipient_classes == "addition"))
    expect_true(any(recipient_classes == "substitution"))
    expect_true(any(recipient_classes == "complex"))
    expect_true(any(recipient_classes == "other"))
})

test_that("a full donor chromosome is added with one homolog bar", {
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
    expected_names <- c(
        "1A", "1R", "2A", "2R", "3A", "4A", "5A", "5R",
        "6A", "7A", "8A", "8R"
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        expected_names
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$label, character(1)),
        expected_names
    )
    added <- vapply(
        payload$chromosomes, function(x) isTRUE(x$added), logical(1)
    )
    expect_identical(which(added), c(2L, 4L, 8L, 12L))
    expect_false(isTRUE(payload$chromosomes[[2L]]$added_partial))
    expect_true(isTRUE(payload$chromosomes[[4L]]$added_partial))
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copy, integer(1)),
        rep(1L, 12L)
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$copies, integer(1)),
        rep(1L, 12L)
    )
})

test_that("a partial addition below 20% stays on the recipient bar", {
    result <- bam_example_result(1000)
    addition <- result$calls[result$calls$class == "substitution" &
        result$calls$pass %in% TRUE, , drop = FALSE][1L, , drop = FALSE]
    addition$start_iid <- 0L
    addition$end_iid <- 500L
    addition$start_bp <- 0
    addition$end_bp <- 24000
    addition$width_bp <- 24000
    ## The example already has a short 1A addition; stay below 20% after
    ## combining the two donor intervals.
    addition$donor_start_bp <- 30000
    addition$donor_end_bp <- 40000
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
        c("1A", "2A", "2R", "3A", "4A", "5A", "5R",
            "6A", "7A", "8A", "8R")
    )
    recipient <- Filter(
        function(chromosome) identical(chromosome$name, "1A"),
        payload$chromosomes
    )[[1L]]
    expect_true(any(vapply(
        recipient$calls, function(call) identical(call$class, "addition"),
        logical(1)
    )))
})

test_that("a partial addition at 20% gets a partial donor bar", {
    result <- bam_example_result(1000)
    addition <- result$calls[result$calls$class == "substitution" &
        result$calls$pass %in% TRUE, , drop = FALSE][1L, , drop = FALSE]
    addition$start_iid <- 0L
    addition$end_iid <- 500L
    addition$start_bp <- 0
    addition$end_bp <- 40000
    addition$width_bp <- 40000
    addition$donor_start_bp <- 0
    addition$donor_end_bp <- 40000
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
        c("1A", "1R", "2A", "2R", "3A", "4A", "5A", "5R",
            "6A", "7A", "8A", "8R")
    )
    donor <- Filter(
        function(chromosome) identical(chromosome$name, "1R"),
        payload$chromosomes
    )[[1L]]
    expect_true(donor$added)
    expect_true(donor$added_partial)
    expect_true(any(vapply(
        donor$calls, function(call) identical(call$class, "addition"),
        logical(1)
    )))
    recipient <- Filter(
        function(chromosome) identical(chromosome$name, "1A"),
        payload$chromosomes
    )[[1L]]
    expect_false(any(vapply(
        recipient$calls, function(call) identical(call$class, "addition"),
        logical(1)
    )))

    svg <- bam_karyotype_static_svg(payload)
    expect_match(svg, 'stroke-dasharray="4 3"', fixed = TRUE)
    expect_match(svg, 'fill="none" stroke="', fixed = TRUE)
    expect_match(svg, 'fill="#1f78b4"', fixed = TRUE)
})

test_that("split additions covering a donor chromosome add it to the plot", {
    result <- bam_example_result(1000)
    donor <- result$assembly$contigs$donor[1L, , drop = FALSE]
    template <- result$calls[1L, , drop = FALSE]
    make_segment <- function(start_bp, end_bp) {
        segment <- template
        segment$class <- "addition"
        segment$class_evidence <- "donor"
        segment$pass <- TRUE
        segment$donor_mean <- 2
        segment$recipient_mean <- 0
        segment$score_mean <- 2
        segment$start_iid <- 10L
        segment$end_iid <- 20L
        segment$width_bp <- 1
        segment$donor_chr_id <- donor$id
        segment$donor_chr_name <- donor$name
        segment$donor_component <- donor$component
        segment$group <- donor$group
        segment$donor_start_bp <- start_bp
        segment$donor_end_bp <- end_bp
        segment$recipient_chr_id <- NA_character_
        segment$recipient_chr_name <- NA_character_
        segment$recipient_component <- NA_character_
        segment$recipient_start_bp <- NA_real_
        segment$recipient_end_bp <- NA_real_
        segment$recipient_assignment <- "no_recipient_loss"
        segment
    }
    result$calls <- rbind(
        make_segment(0, 80000),
        make_segment(80001, as.numeric(donor$length))
    )

    payload <- karyotype_payload(result, "example-sample", "calls")
    added <- Filter(
        function(chromosome) isTRUE(chromosome$added),
        payload$chromosomes
    )
    expect_identical(
        vapply(added, `[[`, "", "name"),
        as.character(donor$name)
    )
    expect_true(all(vapply(
        added, function(chromosome) length(chromosome$calls), integer(1)
    ) == 2L))
    expect_false(any(vapply(
        added, function(chromosome) isTRUE(chromosome$added_partial),
        logical(1)
    )))

    ## One large partial segment still has no recipient bar; show it on the
    ## donor chromosome as a partial addition.
    result$calls <- result$calls[1L, , drop = FALSE]
    payload <- karyotype_payload(result, "example-sample", "calls")
    added <- Filter(
        function(chromosome) isTRUE(chromosome$added),
        payload$chromosomes
    )
    expect_identical(
        vapply(added, `[[`, "", "name"),
        as.character(donor$name)
    )
    expect_true(all(vapply(
        added, function(chromosome) length(chromosome$calls), integer(1)
    ) == 1L))
    expect_true(all(vapply(
        added, function(chromosome) isTRUE(chromosome$added_partial),
        logical(1)
    )))

    ## A tiny partial segment does not justify a whole donor bar.
    result$calls$donor_end_bp <- 8000
    payload <- karyotype_payload(result, "example-sample", "calls")
    added <- Filter(
        function(chromosome) isTRUE(chromosome$added),
        payload$chromosomes
    )
    expect_false(as.character(donor$name) %in% vapply(
        added, `[[`, "", "name"
    ))
})

test_that("a deleted recipient chromosome stays in the karyotype", {
    result <- bam_example_result(1000)
    extra <- result$assembly$contigs$recipient[1, , drop = FALSE]
    extra$id <- "r9"
    extra$name <- "1B"
    extra$component <- "B"
    result$assembly$contigs$recipient <- rbind(
        result$assembly$contigs$recipient, extra
    )
    metric <- result$metrics[result$metrics$role == "recipient", ,
        drop = FALSE
    ][1, , drop = FALSE]
    metric$chr_id <- "r9"
    metric$chr_name <- "1B"
    metric$component <- "B"
    metric$breadth <- 0.02
    result$metrics <- rbind(result$metrics, metric)
    result$metrics$breadth[result$metrics$chr_id == "r1"] <- 0.0009

    payload <- karyotype_payload(result, "example-sample", "heatmap")
    expected_names <- c(
        "1A", "1B", "2A", "2R", "3A", "4A", "5A", "5R",
        "6A", "7A", "8A", "8R"
    )
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        expected_names
    )
    deleted <- vapply(
        payload$chromosomes, function(x) isTRUE(x$deleted), logical(1)
    )
    expect_identical(which(deleted), c(1L, 9L, 11L))
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
        paste0(1:8, "A")
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
    expect_identical(binding[[1L]]$version, "0.1.4")
    expect_true(any(vapply(
        widget$dependencies,
        function(dependency) {
            identical(dependency$name, "mamamiaKaryotypeCard-css") &&
                identical(dependency$stylesheet, "mamamiaKaryotypeCard.css")
        },
        logical(1)
    )))
})

test_that("the karyotype marks donor support for the tooltip", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "calls")
    calls <- unlist(
        lapply(payload$chromosomes, `[[`, "calls"),
        recursive = FALSE
    )
    classes <- vapply(calls, function(call) call$class, character(1))
    donor_supported <- vapply(
        calls, function(call) isTRUE(call$donor_supported), logical(1)
    )
    recipient_supported <- vapply(
        calls, function(call) isTRUE(call$recipient_supported), logical(1)
    )
    expect_true(any(classes == "deletion"))
    expect_false(any(donor_supported[classes == "deletion"]))
    expect_true(any(recipient_supported[classes == "deletion"]))
    expect_true(any(classes == "addition"))
    expect_true(any(donor_supported[classes == "addition"]))
    expect_false(any(recipient_supported[classes == "addition"]))
    expect_true(any(classes == "substitution"))
    expect_true(any(donor_supported[classes == "substitution"]))
    expect_true(any(recipient_supported[classes == "substitution"]))

    added <- Filter(
        function(chromosome) isTRUE(chromosome$added),
        payload$chromosomes
    )
    added_classes <- unlist(lapply(added, function(chromosome) {
        vapply(chromosome$calls, function(call) call$class, character(1))
    }))
    expect_true("addition" %in% added_classes)
    expect_false("deletion" %in% added_classes)

    recipient_bars <- Filter(
        function(chromosome) !isTRUE(chromosome$added),
        payload$chromosomes
    )
    recipient_classes <- unlist(lapply(
        recipient_bars,
        function(chromosome) {
            vapply(chromosome$calls, function(call) call$class, character(1))
        }
    ))
    ## The short 1A addition stays on its recipient bar, while the longer 2A
    ## addition is drawn on the 2R donor bar.
    expect_true("addition" %in% recipient_classes)
})

test_that("the Shiny karyotype output loads htmlwidgets before its binding", {
    skip_if_not_installed("htmlwidgets")
    skip_if_not_installed("htmltools")
    tag <- mamamia_karyotype_output("karyotype")
    names <- vapply(
        htmltools::htmlDependencies(tag),
        function(dependency) dependency$name,
        character(1)
    )
    expect_true("htmlwidgets" %in% names)
    expect_true("mamamiaKaryotypeCard-binding" %in% names)
    expect_lt(
        match("htmlwidgets", names),
        match("mamamiaKaryotypeCard-binding", names)
    )
})

test_that("the explorer exposes the karyotype tab and mode switch", {
    skip_if_not_installed("htmlwidgets")
    skip_if_not_installed("htmltools")
    html <- as.character(explorer_results_ui("explorer"))
    expect_match(html, "Karyotype", fixed = TRUE)
    expect_match(html, "explorer-karyotype_mode", fixed = TRUE)
    expect_match(html, "explorer-hide_deletions", fixed = TRUE)
    expect_match(html, "explorer-deletion_info", fixed = TRUE)
    expect_match(html, "Hide deletions", fixed = TRUE)
    expect_match(html, "explorer-hide_unplaced", fixed = TRUE)
    expect_match(html, "explorer-unplaced_info", fixed = TRUE)
    expect_match(html, "Hide unplaced", fixed = TRUE)
    expect_match(html, "explorer-karyotype", fixed = TRUE)
    expect_match(html, "explorer-tooltip_help", fixed = TRUE)
    expect_match(html, "Tooltip help", fixed = TRUE)
    expect_match(html, "explorer-download_karyotype_current", fixed = TRUE)
    expect_match(html, "explorer-download_karyotype_full", fixed = TRUE)
    expect_match(html, "Current view", fixed = TRUE)
    expect_match(html, "Full view", fixed = TRUE)
    expect_match(html, "mamamiaDownloadKaryotype", fixed = TRUE)
    expect_match(html, "mamamia-analysis-loading", fixed = TRUE)
    expect_match(html, "mamamia-analysis-stage", fixed = TRUE)
    expect_match(html, "Preparing analysis", fixed = TRUE)
    expect_match(html, "mamamia-analysis-timer", fixed = TRUE)
    expect_match(html, "mamamia-analysis-timer-value", fixed = TRUE)
    expect_match(html, "Processing time", fixed = TRUE)
    expect_match(html, "00:00", fixed = TRUE)
})

test_that("the karyotype export payload carries every call", {
    result <- bam_example_result(1000)
    visible <- result$calls[
        result$calls$class == "substitution" &
            result$calls$pass %in% TRUE, , drop = FALSE
    ]
    payload <- karyotype_export_payload(
        result, "example-sample", "calls", calls = visible
    )
    expect_null(payload$full_calls)
    expect_null(payload$full_chromosomes)
    expect_true(is.list(payload$views$full))
    expect_true(
        length(payload$views$full$chromosomes) >=
            length(payload$chromosomes)
    )
    current <- sum(vapply(
        payload$chromosomes,
        function(chromosome) length(chromosome$calls),
        integer(1)
    ))
    full <- sum(vapply(
        payload$views$full$chromosomes,
        function(chromosome) length(chromosome$calls),
        integer(1)
    ))
    expect_gt(current, 0L)
    expect_gt(full, current)
})

test_that("the tooltip help explains the call card fields", {
    html <- as.character(explorer_tooltip_help_ui())
    expect_match(html, "Quick interpretation", fixed = TRUE)
    expect_match(html, "Placement", fixed = TRUE)
    expect_match(html, "Recipient", fixed = TRUE)
    expect_match(html, "Score", fixed = TRUE)
    expect_match(html, "Status chips", fixed = TRUE)
    expect_match(html, "Copy state", fixed = TRUE)
    expect_match(html, "Width ratio", fixed = TRUE)
    expect_match(html, "No-call card", fixed = TRUE)
    expect_match(html, "Caveats", fixed = TRUE)
    expect_match(html, "total depth rises", fixed = TRUE)
})

test_that("a full substitution splits into a lost recipient and added donor", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "calls")
    recipient <- Filter(
        function(x) identical(x$name, "8A"), payload$chromosomes
    )
    donor <- Filter(
        function(x) identical(x$name, "8R"), payload$chromosomes
    )
    expect_length(recipient, 1L)
    expect_length(donor, 1L)
    expect_true(all(vapply(
        recipient, function(x) isTRUE(x$deleted), logical(1)
    )))
    expect_true(all(vapply(
        recipient, function(x) length(x$calls) == 0L, logical(1)
    )))
    expect_false(any(vapply(
        donor, function(x) isTRUE(x$deleted), logical(1)
    )))
    expect_true(all(vapply(
        donor, function(x) isTRUE(x$added), logical(1)
    )))
    expect_true(all(vapply(donor, function(x) {
        any(vapply(x$calls, function(call) {
            identical(call$class, "substitution")
        }, logical(1)))
    }, logical(1))))
    scene_donor <- Filter(
        function(x) identical(x$name, "8R"), payload$scene$chromosomes
    )[[1L]]
    expect_identical(scene_donor$base$fill, payload$scene$palette$added)
    expect_identical(
        scene_donor$calls[[1L]]$fill, payload$scene$palette$substitution
    )
    heatmap_payload <- karyotype_payload(
        result, "example-sample", "heatmap"
    )
    heatmap_donor <- Filter(
        function(x) identical(x$name, "8R"),
        heatmap_payload$scene$chromosomes
    )[[1L]]
    expect_identical(
        heatmap_donor$base$fill, heatmap_payload$scene$palette$zero
    )
})

test_that("hiding deletions removes whole chromosomes and partial segments", {
    result <- bam_example_result(1000)
    shown <- karyotype_payload(result, "example-sample", "calls")
    hidden <- karyotype_payload(
        result, "example-sample", "calls", hide_deletions = TRUE
    )
    labels <- vapply(hidden$chromosomes, function(x) x$name, character(1))
    expect_identical(
        labels,
        c("1A", "2A", "2R", "3A", "4A", "5A", "5R", "7A", "8R")
    )
    expect_true(hidden$hide_deletions)
    expect_identical(hidden$hidden_deletions$segments, 3L)
    expect_gt(hidden$hidden_deletions$bases, 0)
    expect_false(any(vapply(
        hidden$chromosomes, function(x) isTRUE(x$deleted), logical(1)
    )))

    three <- Filter(
        function(x) identical(x$name, "3A"), hidden$chromosomes
    )[[1L]]
    expect_lt(three$display_length, three$length)
    expect_equal(three$display_length, three$length * (1 - 0.022))
    expect_length(three$deleted_intervals, 1L)
    expect_identical(three$deleted_intervals[[1L]]$start, 0.5)
    expect_false(any(vapply(
        three$calls, function(call) identical(call$class, "deletion"),
        logical(1)
    )))

    expect_true(any(vapply(
        shown$chromosomes, function(x) identical(x$name, "6A"), logical(1)
    )))
    expect_identical(shown$hidden_deletions$segments, 0L)
})

test_that("hiding deletions compresses heatmap scores with the chromosome", {
    result <- bam_example_result(1000)
    shown <- karyotype_payload(result, "example-sample", "heatmap")
    hidden <- karyotype_payload(
        result, "example-sample", "heatmap", hide_deletions = TRUE
    )
    shown_three <- Filter(
        function(x) identical(x$name, "3A"), shown$chromosomes
    )[[1L]]
    hidden_three <- Filter(
        function(x) identical(x$name, "3A"), hidden$chromosomes
    )[[1L]]
    expect_lt(length(hidden_three$scores), length(shown_three$scores))
    expect_equal(
        length(hidden_three$scores),
        max(2L, as.integer(round(length(shown_three$scores) * 0.978)))
    )
})

test_that("deletion intervals are merged and remapped after splicing", {
    intervals <- karyotype_merge_intervals(list(
        c(0.20, 0.30), c(0.25, 0.35), c(0.60, 0.70)
    ))
    expect_length(intervals, 2L)
    expect_equal(intervals[[1L]], c(0.20, 0.35))
    expect_equal(intervals[[2L]], c(0.60, 0.70))
    surviving <- 1 - karyotype_interval_width(intervals)
    expect_equal(surviving, 0.75)
    expect_equal(
        karyotype_compress_fraction(c(0.10, 0.50, 0.80), intervals, surviving),
        c(0.10 / 0.75, 0.35 / 0.75, 0.55 / 0.75)
    )
    expect_equal(
        vapply(
            c(0, 0.4, 1), karyotype_display_to_source, numeric(1),
            intervals = intervals, surviving = surviving
        ),
        c(0, 0.45, 1)
    )
})

test_that("the full export keeps deletions when the live view hides them", {
    result <- bam_example_result(1000)
    visible <- result$calls[result$calls$pass %in% TRUE, , drop = FALSE]
    live <- karyotype_export_payload(
        result, "example-sample", "calls",
        calls = visible, hide_deletions = TRUE
    )
    expect_true(live$hide_deletions)
    expect_false(any(vapply(
        live$chromosomes, function(x) identical(x$name, "8A"), logical(1)
    )))
    expect_true(any(vapply(
        live$views$full$chromosomes,
        function(x) identical(x$name, "8A"),
        logical(1)
    )))
    full <- bam_full_karyotype_payload(live)
    expect_false(full$hide_deletions)
    expect_identical(full$hidden_deletions$segments, 0L)
    expect_true(any(vapply(
        full$chromosomes, function(x) identical(x$name, "8A"), logical(1)
    )))
    svg <- bam_karyotype_static_svg(live)
    expect_match(svg, "3 segments omitted", fixed = TRUE)
    expect_false(grepl("Deleted recipient", svg, fixed = TRUE))
})

test_that("a full substitution follows the donor scale and adds only its donor", {
    result <- bam_example_result(1000)
    ## Donor group 8 is 200 kb, while the recipient group maximum is longer.
    ## The call width comes from the donor windows, so the full-chromosome test
    ## must not divide it by the recipient length.
    result$scores$length[as.character(result$scores$group) == "8"] <- 250000
    extra_donor <- result$assembly$contigs$donor[
        result$assembly$contigs$donor$id == "d8", , drop = FALSE
    ]
    extra_donor$id <- "d8b"
    extra_donor$name <- "8S"
    extra_donor$component <- "S"
    result$assembly$contigs$donor <- rbind(
        result$assembly$contigs$donor, extra_donor
    )

    payload <- karyotype_payload(result, "example-sample", "calls")
    labels <- vapply(payload$chromosomes, function(x) x$name, character(1))
    added <- vapply(
        payload$chromosomes, function(x) isTRUE(x$added), logical(1)
    )
    expect_true(all(c("5R", "8R") %in% labels[added]))
    expect_false("8S" %in% labels)
})

test_that("a full substitution only deletes its assigned recipient", {
    result <- bam_example_result(1000)
    extra <- result$assembly$contigs$recipient[
        result$assembly$contigs$recipient$id == "r8", , drop = FALSE
    ]
    extra$id <- "r8b"
    extra$name <- "8B"
    extra$component <- "B"
    result$assembly$contigs$recipient <- rbind(
        result$assembly$contigs$recipient, extra
    )
    normal <- result$metrics[
        result$metrics$role == "recipient" &
            result$metrics$chr_id == "r7", , drop = FALSE
    ]
    normal$chr_id <- "r8b"
    normal$chr_name <- "8B"
    normal$component <- "B"
    normal$group <- "8"
    result$metrics <- rbind(result$metrics, normal)

    payload <- karyotype_payload(result, "example-sample", "calls")
    deleted <- vapply(
        payload$chromosomes, function(x) isTRUE(x$deleted), logical(1)
    )
    labels <- vapply(
        payload$chromosomes, function(x) x$name, character(1)
    )
    expect_true(all(deleted[labels == "8A"]))
    expect_false(any(deleted[labels == "8B"]))

    hidden <- karyotype_payload(
        result, "example-sample", "calls", hide_deletions = TRUE
    )
    hidden_labels <- vapply(
        hidden$chromosomes, function(x) x$name, character(1)
    )
    expect_false("8A" %in% hidden_labels)
    expect_true(all(c("8B", "8R") %in% hidden_labels))
    expect_lt(
        match("8R", hidden_labels), match("8B", hidden_labels)
    )
})

test_that("added donor chromosomes follow their recipient group", {
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
    expect_identical(
        vapply(payload$chromosomes, function(x) x$name, character(1)),
        c(
            "1A", "1R", "2A", "3A", "3R", "4A",
            "5A", "6A", "7A", "8A"
        )
    )
    added <- vapply(payload$chromosomes, function(x) isTRUE(x$added), logical(1))
    expect_identical(which(added), c(2L, 5L))
})
