## Resolved tooltip projection for karyotype calls.

payload_calls <- function(payload) {
    unlist(
        lapply(payload$chromosomes, `[[`, "calls"),
        recursive = FALSE
    )
}

tooltip_section_labels <- function(tooltip) {
    vapply(tooltip$sections, function(section) section$label, character(1))
}

tooltip_row_labels <- function(tooltip) {
    unlist(lapply(tooltip$sections, function(section) {
        if (is.null(section$rows)) {
            return(character())
        }
        vapply(section$rows, function(row) row$label, character(1))
    }))
}

test_that("every new payload call carries a resolved tooltip", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "calls")
    calls <- payload_calls(payload)
    expect_gt(length(calls), 0L)
    expect_true(all(vapply(
        calls,
        function(call) is.list(call$tooltip),
        logical(1)
    )))
    for (call in calls) {
        tooltip <- call$tooltip
        expect_true(nzchar(tooltip$title))
        expect_true(is.list(tooltip$status))
        expect_true(nzchar(tooltip$status$text))
        expect_true(grepl("^mamamia-tip-status-", tooltip$status$class))
        expect_true(nzchar(tooltip$subtitle))
        expect_identical(
            tooltip_section_labels(tooltip)[[1L]], "Placement"
        )
        expect_true("Evidence" %in% tooltip_section_labels(tooltip))
    }
})

test_that("tooltip labels come from the call schema", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "calls")
    calls <- payload_calls(payload)
    labels <- call_schema_labels()

    expect_identical(
        call_tooltip_label("recipient_chr_name"), labels[["recipient_chr_name"]]
    )
    expect_identical(
        call_tooltip_metric_label("donor_mean"), "Donor"
    )
    expect_identical(
        call_tooltip_metric_label("score_mean"), labels[["score_mean"]]
    )
    expect_identical(
        call_tooltip_metric_label("recipient_mean"), "Recipient"
    )
    expect_identical(
        call_tooltip_metric_label("total_mean"), "Total"
    )

    substitution <- Filter(
        function(call) identical(call$class, "substitution"), calls
    )[[1L]]
    tooltip <- substitution$tooltip
    metrics <- Filter(
        function(section) identical(section$label, "Evidence"),
        tooltip$sections
    )[[1L]]$metrics
    expect_identical(
        vapply(metrics, function(metric) metric$label, character(1)),
        c("Score", "Donor", "Recipient", "Total")
    )
    expect_true(all(vapply(metrics, function(metric) {
        is.character(metric$value) && nzchar(metric$value) &&
            metric$class %in% c("pos", "neg", "muted")
    }, logical(1))))
    expect_true(any(grepl(
        "Donor", tooltip_row_labels(tooltip), fixed = TRUE
    )))
    expect_true(any(grepl(
        "Recipient", tooltip_row_labels(tooltip), fixed = TRUE
    )))
})

test_that("the matched-call browser card only reads the tooltip model", {
    js <- readLines(
        mamamia_path("www", "htmlwidgets", "mamamiaKaryotypeCard.js"),
        warn = FALSE
    )
    matched_refs <- unique(unlist(regmatches(
        js,
        gregexpr("matched\\.[A-Za-z_]+", js)
    )))
    expect_gt(length(matched_refs), 0L)
    expect_identical(matched_refs, "matched.tooltip")

    ## Canonical call fields are no longer inspected in the matched-call card.
    canonical <- call_schema_names()
    matched_lines <- grep("matched", js, value = TRUE)
    for (field in canonical) {
        expect_false(any(grepl(
            paste0("matched.", field), matched_lines, fixed = TRUE
        )))
    }
})

test_that("an extra schema-shaped field needs no browser change", {
    result <- bam_example_result(1000)
    payload <- karyotype_payload(result, "example-sample", "calls")
    call <- payload_calls(payload)[[1L]]
    call$future_schema_field <- "ignored by the renderer"

    tooltip <- call_tooltip(
        call,
        sample = payload$sample,
        chromosome_group = call$donor_group
    )
    expect_identical(tooltip_section_labels(tooltip)[[1L]], "Placement")
    expect_true(all(vapply(
        c("title", "status", "subtitle", "sections", "accent"),
        function(name) name %in% names(tooltip),
        logical(1)
    )))
})
