## Canonical call schema and shared call-status rules.

test_that("the call schema is the single declaration of the call record", {
    schema <- call_schema()
    expect_identical(nrow(schema), 64L)
    expect_identical(anyDuplicated(schema$name), 0L)
    expect_identical(anyDuplicated(schema$export), 0L)
    expect_false(anyNA(schema$label))
    expect_false(anyNA(schema$description))
    expect_true(all(nzchar(schema$label)))
    expect_true(all(nzchar(schema$description)))
    expect_true(all(schema$type %in% c(
        "character", "numeric", "integer", "logical"
    )))

    contracts <- bam_result_contracts()
    expect_identical(contracts$calls, call_schema_names())
    dictionary <- bam_bundle_column_dictionary()
    expect_identical(dictionary$old, call_schema_names())
    expect_identical(dictionary$name, unname(call_schema_export()))
    expect_identical(dictionary$description, schema$description)
    expect_identical(
        unname(bam_call_column_labels(call_schema_names())),
        schema$label
    )
})

test_that("the empty and non-empty BEDPE schemas are identical", {
    empty <- .bedpe_table(
        data.frame(sample = character(0), stringsAsFactors = FALSE),
        windows = list(table = data.frame()), pass_only = FALSE
    )
    expect_identical(names(empty), bam_bedpe_columns())
    expect_identical(anyDuplicated(names(empty)), 0L)
    expect_identical(sum(names(empty) == "mapq_fraction"), 1L)
    expect_true(all(c("start_iid", "end_iid") %in% names(empty)))
    expect_identical(empty$pass, logical(0))
    expect_identical(empty$mapq_fraction, numeric(0))

    result <- bam_example_result(1000)
    populated <- .bedpe_table(
        result$calls, windows = result$windows, pass_only = FALSE
    )
    expect_identical(names(populated), names(empty))
    expect_gt(nrow(populated), 0L)
})

test_that("the call record normalizer applies schema types and defaults", {
    record <- call_record(list(
        sample = "S1",
        group = 7,
        pass = "TRUE",
        extra_presentation = "kept"
    ))
    expect_identical(record$sample, "S1")
    expect_identical(record$group, "7")
    expect_true(record$pass)
    expect_identical(record$extra_presentation, "kept")
    expect_true(all(call_schema_names() %in% names(record)))
    expect_identical(names(record)[seq_along(call_schema_names())],
        call_schema_names()
    )

    partial <- call_record(
        list(sample = "S1", class = "addition"),
        include_missing = FALSE
    )
    expect_named(partial, c("sample", "class"))
    expect_identical(partial$class, "addition")

    normalized <- call_record_normalize(list(start_iid = "4"))
    expect_identical(normalized$start_iid, 4L)
})

test_that("component status rules agree across classification and assignment", {
    calls <- data.frame(
        donor_component_spread = c(0, 2.5, 0, 0, NA_real_),
        recipient_component_spread = c(0, 0, 2.1, 0, NA_real_),
        n_donor_supported = c(1, 1, 1, 2, NA_real_),
        stringsAsFactors = FALSE
    )
    expected <- c(
        "consistent", "heterogeneous", "heterogeneous",
        "heterogeneous", "consistent"
    )
    expect_identical(.call_component_status(calls), expected)
    expect_identical(call_component_status_from(calls), expected)
    for (i in seq_len(nrow(calls))) {
        expect_identical(
            .assignment_component_status(calls[i, , drop = FALSE]),
            expected[[i]]
        )
    }
})

test_that("mappability status rules agree across classification and assignment", {
    calls <- data.frame(
        mapq_fraction = c(0.1, 0.5, 0.9, NA_real_),
        stringsAsFactors = FALSE
    )
    expected <- c("low", "ok", "ok", NA_character_)
    expect_identical(.call_mappability_status(calls), expected)
    expect_identical(call_mappability_status_from(calls), expected)
    prepared <- .assignment_empty(calls)
    expect_identical(prepared$mappability_status, expected)

    missing <- data.frame(sample = "S1", stringsAsFactors = FALSE)
    expect_identical(
        .call_mappability_status(missing), NA_character_
    )
    expect_identical(
        .assignment_empty(missing)$mappability_status, NA_character_
    )
})

test_that("the class-confidence downgrade rule is shared", {
    expect_identical(
        call_status_downgrade(
            c("heterogeneous", "consistent", "consistent", NA_character_),
            c("ok", "low", "ok", "low"),
            c("high", "high", "moderate", "high")
        ),
        c(TRUE, TRUE, FALSE, TRUE)
    )
})
