## Global assignment of donor segments to recipient contigs.

assignment_fixture <- function(samples = paste0("S", 1:5), n_per = 20L) {
    contigs <- data.frame(
        role = c(rep("recipient", 4L), rep("donor", 2L)),
        id = c("r1a", "r1b", "r2a", "r2b", "d1", "d2"),
        name = c("1A", "1B", "2A", "2B", "1R", "2R"),
        component = c("A", "B", "A", "B", "R", "R"),
        group = c("1", "1", "2", "2", "1", "2"),
        length = 2000,
        copy_number = 2,
        stringsAsFactors = FALSE
    )
    table <- do.call(rbind, lapply(seq_len(nrow(contigs)), function(i) {
        start <- floor(seq(0L, 2000L, length.out = n_per + 1L))
        data.frame(
            chr_id = contigs$id[[i]],
            start = start[-length(start)],
            stop = start[-1L],
            iid = seq.int(0L, n_per - 1L),
            role = contigs$role[[i]],
            chr_name = contigs$name[[i]],
            component = contigs$component[[i]],
            group = contigs$group[[i]],
            copy_number = 2,
            stringsAsFactors = FALSE
        )
    }))
    metrics <- expand.grid(
        sample = samples,
        chr_id = contigs$id,
        iid = seq.int(0L, n_per - 1L),
        stringsAsFactors = FALSE
    )
    metrics <- merge(
        metrics, table,
        by = c("chr_id", "iid"), sort = FALSE
    )
    metrics$breadth <- ifelse(metrics$role == "recipient", 0.02, 0)
    metrics$breadth[
        metrics$role == "donor" &
            metrics$iid >= 5L & metrics$iid <= 9L
    ] <- 0.5
    metrics$lambda_breadth <- -log(1 - metrics$breadth)
    metrics$lambda_count <- metrics$lambda_breadth
    metrics$gc <- 0.5
    metrics$count <- metrics$breadth * 10
    metrics$covered <- metrics$breadth * 100
    metrics$length <- 100
    metrics$read_length <- 50

    scores <- expand.grid(
        sample = samples,
        group = c("1", "2"),
        iid = seq.int(0L, n_per - 1L),
        stringsAsFactors = FALSE
    )
    scores$mode <- "multi"
    scores$donor_component <- "R"
    scores$score <- 0
    scores$start <- table$start[match(scores$iid, table$iid)]
    scores$stop <- table$stop[match(scores$iid, table$iid)]

    calls <- data.frame(
        sample = samples[[1L]],
        group = c("1", "2"),
        start_iid = 5L,
        end_iid = 9L,
        start_bp = 500,
        end_bp = 1000,
        width_bp = 500,
        score_mean = 1,
        donor_mean = 1,
        recipient_mean = 0,
        class = "addition",
        confidence = "high",
        method = "paired_lr",
        mode = "multi",
        pass = TRUE,
        stringsAsFactors = FALSE
    )

    list(
        assembly = list(contigs = list(
            recipient = contigs[contigs$role == "recipient", ,
                drop = FALSE
            ],
            donor = contigs[contigs$role == "donor", , drop = FALSE]
        )),
        metrics = metrics,
        windows = list(table = table, n_per_chromosome = n_per),
        scores = scores,
        calls = calls
    )
}

assign_fixture <- function(fixture) {
    assignCallPositions(
        fixture$calls, fixture$metrics, fixture$windows,
        fixture$scores, fixture$assembly
    )
}

one_call <- function(fixture, group = "1") {
    fixture$calls <- fixture$calls[
        as.character(fixture$calls$group) == group, , drop = FALSE
    ]
    fixture
}

delete_interval <- function(fixture, sample, chr_id, from = 5L, to = 9L,
                            value = 0.0001) {
    rows <- as.character(fixture$metrics$sample) == sample &
        as.character(fixture$metrics$chr_id) == chr_id &
        fixture$metrics$iid >= from & fixture$metrics$iid <= to
    fixture$metrics$breadth[rows] <- value
    fixture$metrics$lambda_breadth[rows] <- -log(1 - value)
    fixture
}

test_that("a same-group deletion is assigned to the depleted contig", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r1a")
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$recipient_chr_name, "1A")
    expect_identical(call$recipient_group, "1")
    expect_identical(call$recipient_assignment, "assigned")
    expect_false(call$recipient_cross_group)
    expect_identical(call$class, "substitution")
    expect_identical(call$class_group, "addition")
    expect_gt(call$recipient_mean_depletion, 1.5)
    expect_identical(call$donor_chr_name, "1R")
})

test_that("a cross-group deletion is reported as a translocation candidate", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r2a")
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$recipient_chr_name, "2A")
    expect_identical(call$recipient_group, "2")
    expect_true(call$recipient_cross_group)
    ## Without a reliable candidate in the call's own group this is a
    ## translocation candidate, not a confirmed host. This stops a
    ## copy-neutral event from being assigned to an unrelated deletion.
    expect_identical(call$recipient_assignment, "ambiguous")
    expect_identical(call$class, "addition")
    expect_identical(call$class_group, "addition")
    expect_true(is.finite(call$recipient_margin))
})

test_that("an addition with no recipient loss is not forced onto a host", {
    fixture <- one_call(assignment_fixture())
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$recipient_assignment, "no_recipient_loss")
    expect_true(is.na(call$recipient_chr_id))
    expect_identical(call$class, "addition")
    expect_identical(call$class_group, "addition")
})

test_that("a copy-neutral call is not assigned to an unrelated deletion", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r2a")
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$recipient_chr_name, "2A")
    expect_true(call$recipient_cross_group)
    expect_identical(call$recipient_assignment, "ambiguous")
    expect_identical(call$class, "addition")
    expect_identical(call$class_group, "addition")
})

test_that("a constitutive low-coverage region is not a host", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r2a", value = 0.0002)
    fixture <- delete_interval(fixture, "S2", "r2a", value = 0.0002)
    fixture <- delete_interval(fixture, "S3", "r2a", value = 0.0002)
    fixture <- delete_interval(fixture, "S4", "r2a", value = 0.0002)
    fixture <- delete_interval(fixture, "S5", "r2a", value = 0.0002)
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_false(identical(call$recipient_chr_name, "2A"))
    expect_identical(call$recipient_assignment, "no_recipient_loss")
})

test_that("one call cannot take a recipient contig from a stronger call", {
    fixture <- assignment_fixture()
    fixture <- delete_interval(fixture, "S1", "r1a")
    calls <- assign_fixture(fixture)
    first <- calls[calls$group == "1", , drop = FALSE]
    second <- calls[calls$group == "2", , drop = FALSE]
    expect_identical(first$recipient_chr_name, "1A")
    expect_identical(second$recipient_assignment, "no_recipient_loss")
})

test_that("a whole-chromosome deletion does not capture a small segment", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r2a", from = 0L, to = 19L)
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_true(is.na(call$recipient_chr_id))
    expect_identical(call$recipient_assignment, "no_recipient_loss")
})

test_that("a weak cross-group candidate is not reported as assigned", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r2a", value = 0.01)
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_true(call$recipient_assignment %in% c("ambiguous", "no_recipient_loss"))
})

test_that("assignment reclassifies a deletion as a substitution", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r1a")
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$class_group, "addition")
    expect_identical(call$class, "substitution")
})

test_that("the assignment is stable when it is run twice", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r1a")
    first <- assign_fixture(fixture)
    second <- assignCallPositions(
        first, fixture$metrics, fixture$windows,
        fixture$scores, fixture$assembly
    )
    expect_identical(first$recipient_chr_name, second$recipient_chr_name)
    expect_identical(first$class, second$class)
    expect_identical(first$class_group, second$class_group)
})

test_that("a non-passing call cannot displace a passing call", {
    candidates <- data.frame(
        call = c(1L, 1L, 2L, 2L),
        chr_id = c("r1a", "r2a", "r1a", "r2a"),
        score = c(5, 4.5, 4.4, 3.5),
        run_start_iid = c(0L, 0L, 0L, 0L),
        run_end_iid = c(1L, 1L, 1L, 1L),
        reliable = c(TRUE, TRUE, TRUE, TRUE),
        pass = c(TRUE, TRUE, FALSE, FALSE),
        stringsAsFactors = FALSE
    )
    contigs <- data.frame(
        chr_id = c("r1a", "r2a"), stringsAsFactors = FALSE
    )
    chosen <- .assignment_match_sample(candidates, 2L, contigs)
    expect_identical(chosen[[1L]], 1L)
    expect_identical(chosen[[2L]], 4L)
})

test_that("dosage distinguishes one and two affected copies", {
    half_breadth <- 1 - sqrt(1 - 0.02)
    full <- .assignment_dosage("assigned", 1, 0.02, 0.0001, 2)
    expect_identical(full$affected_copies, 2L)
    expect_identical(full$confidence, "high")
    expect_gt(full$affected_fraction, 0.9)

    half <- .assignment_dosage("assigned", 1, 0.02, half_breadth, 2)
    expect_identical(half$affected_copies, 1L)
    expect_identical(half$confidence, "high")
    expect_equal(half$affected_fraction, 0.5, tolerance = 0.02)

    partial <- .assignment_dosage("assigned", 1, 0.02, 0.015, 2)
    expect_true(is.na(partial$affected_copies))
    expect_identical(partial$confidence, "low")

    expect_identical(
        .assignment_dosage("ambiguous", 1, 0.02, 0.0001, 2)$confidence,
        "unavailable"
    )
    expect_identical(
        .assignment_dosage("assigned", 0.1, 0.02, 0.0001, 2)$confidence,
        "unavailable"
    )
    expect_identical(
        .assignment_dosage("assigned", 1, 0.02, 0.02, 2)$confidence,
        "unavailable"
    )
})

test_that("assignment records recipient dosage", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r1a", value = 0.0001)
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$recipient_copy_number, 2)
    expect_identical(call$affected_copies, 2L)
    expect_identical(call$dosage_confidence, "high")
    expect_gt(call$affected_fraction, 0.9)
})

test_that("a deletion-only call is not classified as a substitution", {
    fixture <- one_call(assignment_fixture())
    fixture$calls$donor_mean <- 0.1
    fixture <- delete_interval(fixture, "S1", "r1a", value = 0.0001)
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$recipient_assignment, "assigned")
    expect_identical(call$class, "deletion")
    expect_identical(call$class_evidence, "recipient")
    expect_identical(call$class_group, "addition")
    expect_true(is.na(call$affected_copies))
})

test_that("assignment records donor and recipient evidence", {
    fixture <- one_call(assignment_fixture())
    fixture <- delete_interval(fixture, "S1", "r1a", value = 0.0001)
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$class, "substitution")
    expect_identical(call$class_evidence, "donor+recipient")
    expect_identical(call$class_confidence, "high")
})

test_that("a width-mismatched substitution is complex and needs review", {
    fixture <- one_call(assignment_fixture())
    fixture$calls$width_bp <- 200
    fixture$calls$end_bp <- 700
    fixture <- delete_interval(fixture, "S1", "r1a", value = 0.0001)
    calls <- assign_fixture(fixture)
    call <- calls[calls$group == "1", , drop = FALSE]
    expect_identical(call$class, "complex")
    expect_identical(call$class_evidence, "donor+recipient; width mismatch")
    expect_identical(call$class_status, "needs_review")
})
