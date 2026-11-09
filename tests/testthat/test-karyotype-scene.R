## Shared karyotype scene and self-contained export views.

scene_legend_keys <- function(scene) {
    vapply(scene$legend$entries, function(entry) entry$key, character(1))
}

scene_bin_count <- function(scene) {
    sum(vapply(
        scene$chromosomes,
        function(chromosome) length(chromosome$heatmap_bins),
        integer(1)
    ))
}

## A synthetic St02-like donor: the call is a full substitution judged by the
## recipient loss, but the donor segment stops at 92.8% of the donor
## chromosome. This is the shape that exposed the class-coloured base bug.
st02_like_result <- function() {
    result <- bam_example_result(1000)
    donor_name <- "8R"
    recipient_name <- "8A"
    donor_row <- which(
        as.character(result$assembly$contigs$donor$name) == donor_name
    )
    recipient_row <- which(
        as.character(result$assembly$contigs$recipient$name) == recipient_name
    )
    result$assembly$contigs$donor$name[donor_row] <- "St02"
    result$assembly$contigs$donor$component[donor_row] <- "St"
    result$assembly$contigs$donor$length[donor_row] <- 495766686
    result$assembly$contigs$recipient$name[recipient_row] <- "2A"
    result$assembly$contigs$recipient$length[recipient_row] <- 788135348

    windows <- result$windows$table
    donor_rows <- which(
        as.character(windows$role) == "donor" &
            as.character(windows$group) == "8"
    )
    windows$stop[donor_rows[which.max(windows$stop[donor_rows])]] <- 495766686
    result$windows$table <- windows

    index <- which(
        as.character(result$calls$donor_chr_name) == donor_name &
            as.character(result$calls$class) == "substitution"
    )[[1L]]
    calls <- result$calls
    calls$start_iid[index] <- 0L
    calls$end_iid[index] <- 999L
    calls$start_bp[index] <- 0
    calls$end_bp[index] <- 664187716
    calls$width_bp[index] <- 664187716
    calls$donor_chr_name[index] <- "St02"
    calls$donor_component[index] <- "St"
    calls$donor_start_bp[index] <- 0
    calls$donor_end_bp[index] <- 460071484
    calls$recipient_chr_name[index] <- "2A"
    calls$recipient_start_bp[index] <- 0
    calls$recipient_end_bp[index] <- 731389602
    calls$recipient_run_start_bp[index] <- 0
    calls$recipient_run_end_bp[index] <- 773948911
    calls$recipient_deletion_width_bp[index] <- 773948911
    calls$recipient_assignment[index] <- "assigned"
    calls$class_status[index] <- "confident"
    calls$confidence[index] <- "high"
    calls$pass[index] <- TRUE
    result$calls <- calls
    result
}

svg_groups <- function(svg) {
    regmatches(svg, gregexpr("<g>.*?</g>", svg, perl = TRUE))[[1L]]
}

## One full-span neutral call with no recipient host, on a donor contig that
## is not otherwise added to the example karyotype.
unplaced_result <- function() {
    result <- bam_example_result(1000)
    donor <- result$assembly$contigs$donor[3L, , drop = FALSE]
    call <- result$calls[1L, , drop = FALSE]
    call$sample <- "example-sample"
    call$group <- as.character(donor$group)
    call$start_iid <- 0L
    call$end_iid <- 999L
    call$start_bp <- 0
    call$end_bp <- as.numeric(donor$length)
    call$width_bp <- as.numeric(donor$length)
    call$score_mean <- 0.1
    call$donor_mean <- 0.1
    call$recipient_mean <- 0
    call$total_mean <- 0
    call$class <- "other"
    call$class_evidence <- "none"
    call$class_confidence <- "low"
    call$class_status <- "needs_review"
    call$confidence <- "low"
    call$copy_state <- "neutral"
    call$total_evidence <- "neutral"
    call$pass <- FALSE
    call$recipient_assignment <- "no_recipient_loss"
    call$recipient_chr_id <- NA_character_
    call$recipient_chr_name <- NA_character_
    call$recipient_component <- NA_character_
    call$recipient_group <- NA_character_
    call$recipient_start_bp <- NA_real_
    call$recipient_end_bp <- NA_real_
    call$recipient_run_start_bp <- NA_real_
    call$recipient_run_end_bp <- NA_real_
    call$recipient_deletion_width_bp <- NA_real_
    call$donor_chr_id <- as.character(donor$id)
    call$donor_chr_name <- as.character(donor$name)
    call$donor_component <- as.character(donor$component)
    call$donor_start_bp <- 0
    call$donor_end_bp <- as.numeric(donor$length)
    call$n_windows <- 1000L
    result$calls <- rbind(result$calls, call)
    result
}

test_that("the scene resolves legend and call styles for the static SVG", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "calls")
    scene <- karyotype_scene(payload)
    svg <- bam_karyotype_static_svg(payload)

    expect_identical(scene$legend$kind, "calls")
    expect_identical(
        scene_legend_keys(scene),
        c(
            "substitution", "addition", "deletion", "complex",
            "other", "deleted", "added", "unplaced"
        )
    )
    for (entry in scene$legend$entries) {
        expect_match(svg, entry$label, fixed = TRUE)
        expect_match(
            svg, sprintf('fill="%s"', entry$fill), fixed = TRUE
        )
    }

    chromosome_index <- which(vapply(
        scene$chromosomes,
        function(chromosome) length(chromosome$calls) > 0L,
        logical(1)
    ))[[1L]]
    calls <- scene$chromosomes[[chromosome_index]]$calls
    expect_gt(length(calls), 0L)
    for (call in calls) {
        expect_match(svg, sprintf('fill="%s"', call$fill), fixed = TRUE)
        expect_match(
            svg,
            sprintf(
                'opacity="%s"',
                formatC(call$opacity, format = "f", digits = 2)
            ),
            fixed = TRUE
        )
    }
})

test_that("added donors use zero in heatmap and state colour in calls", {
    result <- st02_like_result()
    payload <- karyotype_payload(result, "example-sample", "heatmap")
    palette <- payload$scene$palette
    st02 <- Filter(
        function(chromosome) identical(chromosome$name, "St02"),
        payload$scene$chromosomes
    )[[1L]]

    expect_true(st02$added)
    expect_false(st02$added_partial)
    expect_identical(st02$base$fill, palette$zero)
    expect_identical(st02$base$stroke, palette$added_text)
    expect_identical(st02$caption$fill, palette$added_text)
    expect_length(st02$calls, 1L)
    expect_identical(st02$calls[[1L]]$class, "substitution")
    expect_identical(st02$calls[[1L]]$fill, palette$substitution)
    expect_lt(st02$calls[[1L]]$end, 1)
    expect_lt(tail(st02$heatmap_bins, 1L)[[1L]]$end, 1)
    expect_gt(tail(st02$heatmap_bins, 1L)[[1L]]$end, 0.9)

    svg <- bam_karyotype_static_svg(payload)
    groups <- svg_groups(svg)
    st02_group <- groups[grepl("St02", groups, fixed = TRUE)]
    expect_length(st02_group, 1L)
    expect_match(
        st02_group, sprintf('fill="%s"', palette$zero), fixed = TRUE
    )
    ## Heatmap mode draws bins over the call and no call rectangle, so neither
    ## the added-donor fill nor the class colour appears in the group.
    expect_false(grepl(
        sprintf('fill="%s"', palette$added),
        st02_group,
        fixed = TRUE
    ))
    expect_false(grepl(
        sprintf('fill="%s"', palette$substitution),
        st02_group,
        fixed = TRUE
    ))

    calls_payload <- karyotype_payload(result, "example-sample", "calls")
    calls_st02 <- Filter(
        function(chromosome) identical(chromosome$name, "St02"),
        calls_payload$scene$chromosomes
    )[[1L]]
    expect_identical(calls_st02$base$fill, palette$added)
    calls_svg <- bam_karyotype_static_svg(calls_payload)
    calls_group <- svg_groups(calls_svg)
    calls_st02_group <- calls_group[grepl("St02", calls_group, fixed = TRUE)]
    expect_length(calls_st02_group, 1L)
    expect_match(
        calls_st02_group, sprintf('fill="%s"', palette$added), fixed = TRUE
    )
    expect_match(
        calls_st02_group, sprintf('fill="%s"', palette$substitution),
        fixed = TRUE
    )
})

test_that("full additions use the state base and partial additions stay dashed", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "calls")
    palette <- payload$scene$palette
    added <- Filter(
        function(chromosome) isTRUE(chromosome$added),
        payload$scene$chromosomes
    )
    full_addition <- Filter(
        function(chromosome) {
            !isTRUE(chromosome$added_partial) &&
                length(chromosome$calls) == 1L &&
                identical(chromosome$calls[[1L]]$class, "addition")
        },
        added
    )[[1L]]
    expect_identical(full_addition$base$fill, palette$added)
    expect_identical(full_addition$base$stroke, palette$added_text)
    expect_identical(full_addition$caption$fill, palette$added_text)
    expect_identical(full_addition$calls[[1L]]$fill, palette$addition)

    partial <- Filter(
        function(chromosome) isTRUE(chromosome$added_partial),
        added
    )[[1L]]
    expect_identical(partial$base$fill, "none")
    expect_identical(partial$base$stroke, palette$negative)
    expect_identical(partial$calls[[1L]]$fill, palette$addition)

    heatmap_payload <- karyotype_payload(
        result, "example-sample", "heatmap"
    )
    heatmap_added <- Filter(
        function(chromosome) isTRUE(chromosome$added),
        heatmap_payload$scene$chromosomes
    )
    heatmap_full <- Filter(
        function(chromosome) {
            !isTRUE(chromosome$added_partial) &&
                length(chromosome$calls) == 1L &&
                identical(chromosome$calls[[1L]]$class, "addition")
        },
        heatmap_added
    )[[1L]]
    heatmap_partial <- Filter(
        function(chromosome) isTRUE(chromosome$added_partial),
        heatmap_added
    )[[1L]]
    expect_identical(heatmap_full$base$fill, palette$zero)
    expect_identical(heatmap_partial$base$fill, palette$zero)
    expect_identical(heatmap_partial$base$stroke, palette$negative)
})

test_that("unplaced calls render neutral bars when not hidden", {
    result <- unplaced_result()
    fake_id <- nrow(result$calls)

    hidden <- karyotype_payload(
        result, "example-sample", "calls", hide_unplaced = TRUE
    )
    expect_false(any(vapply(
        hidden$chromosomes,
        function(chromosome) isTRUE(chromosome$unplaced),
        logical(1)
    )))
    expect_false(fake_id %in% hidden$rendered_call_ids)
    expect_true(fake_id %in% hidden$unplaced_call_ids)

    shown <- karyotype_payload(
        result, "example-sample", "calls",
        hide_unplaced = FALSE, show_unplaced = TRUE
    )
    bars <- Filter(
        function(chromosome) isTRUE(chromosome$unplaced),
        shown$chromosomes
    )
    expect_gte(length(bars), 1L)
    bar <- Filter(
        function(chromosome) identical(chromosome$name, "3R"),
        bars
    )[[1L]]
    expect_identical(bar$name, "3R")
    expect_identical(bar$calls[[1L]]$class, "other")
    expect_true(fake_id %in% shown$rendered_call_ids)
    expect_true(fake_id %in% shown$unplaced_call_ids)
    scene_bar <- Filter(
        function(chromosome) {
            isTRUE(chromosome$unplaced) &&
                identical(chromosome$name, "3R")
        },
        shown$scene$chromosomes
    )[[1L]]
    expect_identical(scene_bar$base$fill, shown$scene$palette$unplaced)
    expect_identical(
        scene_bar$base$stroke, shown$scene$palette$unplaced_text
    )
    expect_identical(scene_bar$caption$text, "Unplaced")

    heatmap <- karyotype_payload(
        result, "example-sample", "heatmap",
        hide_unplaced = FALSE, show_unplaced = TRUE
    )
    heatmap_bar <- Filter(
        function(chromosome) {
            isTRUE(chromosome$unplaced) &&
                identical(chromosome$name, "3R")
        },
        heatmap$scene$chromosomes
    )[[1L]]
    expect_identical(heatmap_bar$base$fill, heatmap$scene$palette$zero)
    expect_gt(length(heatmap_bar$heatmap_bins), 0L)

    calls <- result$calls[
        as.character(result$calls$sample) == "example-sample", ,
        drop = FALSE
    ]
    calls$.call_id <- seq_len(nrow(calls))
    rendered <- explorer_rendered_calls(calls, shown$rendered_call_ids)
    expect_true(fake_id %in% rendered$.call_id)
    expect_false(any(
        !is.na(rendered$recipient_assignment) &
            rendered$recipient_assignment == "no_recipient_loss" &
            !rendered$.call_id %in% shown$rendered_call_ids
    ))

    export <- karyotype_export_payload(
        result, "example-sample", "calls", hide_unplaced = TRUE
    )
    expect_false(any(vapply(
        export$views$full$chromosomes,
        function(chromosome) isTRUE(chromosome$unplaced),
        logical(1)
    )))
})

test_that("donor-only additions stay placed when unplaced calls are hidden", {
    result <- bam_example_result(1000)
    calls <- result$calls[
        as.character(result$calls$sample) == "example-sample", ,
        drop = FALSE
    ]
    calls$.call_id <- seq_len(nrow(calls))
    call_id <- which(
        as.character(calls$donor_chr_name) == "5R" &
            as.character(calls$class) == "addition"
    )[[1L]]

    payload <- karyotype_payload(
        result, "example-sample", "calls",
        calls = calls, hide_unplaced = TRUE
    )
    bar <- Filter(
        function(chromosome) identical(chromosome$name, "5R"),
        payload$chromosomes
    )
    expect_length(bar, 1L)
    expect_true(bar[[1L]]$added)
    expect_true(call_id %in% payload$rendered_call_ids)
    expect_false(call_id %in% payload$unplaced_call_ids)
})

test_that("the scene respects calls, heatmap and hidden deletions", {
    result <- bam_example_result(1000)
    calls_payload <- karyotype_payload(result, "example-sample", "calls")
    calls_scene <- karyotype_scene(calls_payload, heatmap_bins = 80L)
    expect_identical(calls_scene$mode, "calls")
    expect_identical(scene_bin_count(calls_scene), 0L)

    heatmap_payload <- karyotype_payload(
        result, "example-sample", "heatmap"
    )
    heatmap_scene <- heatmap_payload$scene
    expect_identical(heatmap_scene$legend$kind, "heatmap")
    expect_gt(scene_bin_count(heatmap_scene), 0L)
    heatmap_svg <- bam_karyotype_static_svg(heatmap_payload)
    expect_match(heatmap_svg, "Score", fixed = TRUE)
    expect_match(heatmap_svg, "linearGradient", fixed = TRUE)
    with_bins <- Filter(
        function(chromosome) length(chromosome$heatmap_bins) > 0L,
        heatmap_scene$chromosomes
    )[[1L]]
    expect_match(
        heatmap_svg,
        sprintf('fill="%s"', with_bins$heatmap_bins[[1L]]$fill),
        fixed = TRUE
    )

    hidden_payload <- karyotype_payload(
        result, "example-sample", "calls", hide_deletions = TRUE
    )
    hidden_scene <- karyotype_scene(hidden_payload)
    expect_false(
        any(scene_legend_keys(hidden_scene) %in% c("deletion", "deleted"))
    )
    expect_identical(hidden_scene$omitted_note$segments, 3L)
    expect_identical(
        hidden_scene$omitted_note$text, "3 segments omitted"
    )
    expect_false(any(vapply(
        hidden_scene$chromosomes,
        function(chromosome) isTRUE(chromosome$deleted),
        logical(1)
    )))
    hidden_svg <- bam_karyotype_static_svg(hidden_payload)
    expect_match(hidden_svg, "3 segments omitted", fixed = TRUE)
    expect_false(grepl("Deleted recipient", hidden_svg, fixed = TRUE))
})

test_that("the export payload carries complete self-contained views", {
    result <- bam_example_result(1000)
    sample <- "example-sample"
    visible <- result$calls[
        as.character(result$calls$class) == "substitution" &
            isTRUE_vec(result$calls$pass), ,
        drop = FALSE
    ]
    payload <- karyotype_export_payload(
        result, sample, "calls", calls = visible, hide_deletions = TRUE
    )
    expect_null(payload$full_calls)
    expect_null(payload$full_chromosomes)
    expect_true(is.list(payload$views$full))
    expect_false(payload$views$full$hide_deletions)

    self_contained <- function(chromosomes) {
        length(chromosomes) > 0L && all(vapply(
            chromosomes,
            function(chromosome) {
                is.list(chromosome$calls) &&
                    is.numeric(chromosome$scores)
            },
            logical(1)
        ))
    }
    expect_true(self_contained(payload$chromosomes))
    expect_true(self_contained(payload$views$full$chromosomes))

    expected <- karyotype_payload(
        result, sample, "calls", calls = NULL, hide_deletions = FALSE
    )
    count_calls <- function(chromosomes) {
        sum(vapply(
            chromosomes,
            function(chromosome) length(chromosome$calls),
            integer(1)
        ))
    }
    expect_identical(
        count_calls(payload$views$full$chromosomes),
        count_calls(expected$chromosomes)
    )
    expect_gt(
        count_calls(payload$views$full$chromosomes),
        count_calls(payload$chromosomes)
    )
    expect_true(any(vapply(
        payload$views$full$chromosomes,
        function(chromosome) identical(chromosome$name, "8A"),
        logical(1)
    )))
})

test_that("the full view does not depend on positional export side channels", {
    result <- bam_example_result(1000)
    payload <- karyotype_export_payload(
        result, "example-sample", "calls",
        calls = result$calls[
            isTRUE_vec(result$calls$pass), , drop = FALSE
        ]
    )
    view_names <- vapply(
        payload$views$full$chromosomes,
        function(chromosome) as.character(chromosome$name),
        character(1)
    )
    ## A stale side channel must not override the self-contained view.
    payload$full_chromosomes <- list(list(
        name = "BOGUS", calls = list(), selected = FALSE
    ))
    payload$full_calls <- list(list(list(class = "addition")))
    full <- bam_full_karyotype_payload(payload)
    expect_identical(
        vapply(
            full$chromosomes,
            function(chromosome) as.character(chromosome$name),
            character(1)
        ),
        view_names
    )
    expect_false("BOGUS" %in% vapply(
        full$chromosomes,
        function(chromosome) as.character(chromosome$name),
        character(1)
    ))

    ## A payload without the full view has no full export. The old positional
    ## side channels must not reconstruct one.
    legacy <- list(
        chromosomes = list(list(
            id = "r1", name = "1A", component = "A", group = "1",
            length = 100, display_length = 100, copy = 1L, copies = 1L,
            calls = list(), scores = numeric(), selected = FALSE
        )),
        full_calls = list(list(list(
            class = "addition", start = 0, end = 1
        ))),
        full_chromosomes = list(list(
            id = "r1", name = "1A", calls = list(), selected = FALSE
        ))
    )
    expect_null(bam_full_karyotype_payload(legacy))
})

test_that("the live scene carries the shared heatmap bins", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "heatmap")
    scene <- payload$scene
    policy <- karyotype_scene_heatmap_bin_count()

    expect_identical(policy, 120L)
    expect_gt(scene_bin_count(scene), 0L)
    with_bins <- Filter(
        function(chromosome) length(chromosome$heatmap_bins) > 0L,
        scene$chromosomes
    )
    for (chromosome in with_bins) {
        expect_lte(length(chromosome$heatmap_bins), policy)
        expect_true(all(vapply(
            chromosome$heatmap_bins,
            function(bin) bin$start >= 0 && bin$end <= 1,
            logical(1)
        )))
    }

    recomputed <- karyotype_scene(payload, heatmap_bins = policy)
    expect_identical(
        vapply(scene$chromosomes, function(chromosome) {
            length(chromosome$heatmap_bins)
        }, integer(1)),
        vapply(recomputed$chromosomes, function(chromosome) {
            length(chromosome$heatmap_bins)
        }, integer(1))
    )
})

test_that("the static SVG draws the scene bins, not a plot-height policy", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "heatmap")
    svg <- bam_karyotype_static_svg(payload, width = 1600, height = 900)
    with_bins <- Filter(
        function(chromosome) length(chromosome$heatmap_bins) > 0L,
        payload$scene$chromosomes
    )
    expect_gt(length(with_bins), 0L)
    for (bin in with_bins[[1L]]$heatmap_bins) {
        expect_match(svg, sprintf('fill="%s"', bin$fill), fixed = TRUE)
    }

    ## A plot-height-derived resolution would use 345 bins here, so its first
    ## bin is a different width. The scene policy governs the rendered bin.
    plot_height_bins <- max(40L, as.integer(floor((900 - 118 - 92) / 2)))
    expect_gt(plot_height_bins, karyotype_scene_heatmap_bin_count())
})

test_that("hidden deletions keep heatmap bins aligned with the spliced bar", {
    result <- bam_example_result(1000)
    hidden <- karyotype_payload(
        result, "example-sample", "heatmap", hide_deletions = TRUE
    )
    ## 3A's only calls are the partial deletion that was spliced out, so it
    ## has no surviving call interval and therefore no heatmap bins.
    three <- Filter(
        function(chromosome) identical(chromosome$name, "3A"),
        hidden$chromosomes
    )[[1L]]
    three_scene <- Filter(
        function(chromosome) identical(chromosome$name, "3A"),
        hidden$scene$chromosomes
    )[[1L]]
    expect_gt(length(three$deleted_intervals), 0L)
    expect_length(three_scene$heatmap_bins, 0L)

    with_bins <- which(vapply(
        hidden$scene$chromosomes,
        function(chromosome) length(chromosome$heatmap_bins) > 0L,
        logical(1)
    ))
    expect_gt(length(with_bins), 0L)
    for (i in with_bins) {
        bins <- hidden$scene$chromosomes[[i]]$heatmap_bins
        expect_true(all(vapply(
            bins,
            function(bin) bin$start >= 0 && bin$end <= 1,
            logical(1)
        )))
        ## The hidden scene bins are built from the same compressed scores
        ## and calls that the browser and static SVG draw.
        expect_identical(
            length(bins),
            length(karyotype_scene_bins(
                hidden$chromosomes[[i]]$scores,
                hidden$chromosomes[[i]]$calls,
                karyotype_scene_heatmap_bin_count(),
                hidden$scene$palette,
                hidden$scene$score_limit
            ))
        )
    }
})

test_that("the browser uses the scene, its bins and the self-contained full view", {
    js <- readLines(
        mamamia_path("www", "htmlwidgets", "mamamiaKaryotypeCard.js"),
        warn = FALSE
    )
    expect_true(any(grepl("payload.scene", js, fixed = TRUE)))
    expect_true(any(grepl("heatmap_bins", js, fixed = TRUE)))
    expect_false(any(grepl("bar.height / 2", js, fixed = TRUE)))
    expect_true(any(grepl(
        "payload.views && payload.views.full", js, fixed = TRUE
    )))
    expect_false(any(grepl("full_chromosomes", js, fixed = TRUE)))
    expect_false(any(grepl("full_calls", js, fixed = TRUE)))
})
