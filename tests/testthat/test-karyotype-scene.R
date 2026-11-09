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
            "other", "deleted", "added"
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
