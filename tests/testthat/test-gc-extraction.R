## T2.4: route windows to gcrank indexes and build the complete GC table.

## Replace the global gcrankGc() binding for the duration of one test. The
## source loader keeps application functions in the global environment, so a
## namespace-level mock would not reach them.
local_gcrank_mock <- function(fn) {
    env <- globalenv()
    had <- exists("gcrankGc", envir = env, inherits = FALSE)
    old <- if (had) get("gcrankGc", envir = env, inherits = FALSE) else NULL
    assign("gcrankGc", fn, envir = env)
    withr::defer({
        if (had) {
            assign("gcrankGc", old, envir = env)
        } else if (exists("gcrankGc", envir = env, inherits = FALSE)) {
            rm(list = "gcrankGc", envir = env)
        }
    }, envir = parent.frame())
    invisible(fn)
}

two_index_assembly <- function() {
    list(
        key = "two-index",
        label = "Two-index assembly",
        description = "",
        gc_indexes = c("idx-recipient", "idx-donor"),
        contigs = list(
            recipient = data.frame(
                id = "r1", name = "1R", component = "R", length = 20000,
                index_path = "idx-recipient", stringsAsFactors = FALSE
            ),
            donor = data.frame(
                id = "d1", name = "1D", component = "D", length = 20000,
                index_path = "idx-donor", stringsAsFactors = FALSE
            )
        ),
        pairs = c(d1 = "r1")
    )
}

test_that("extractWindowGc() fills a complete GC table in window order", {
    skip_if_no_tiny_fixtures()
    on.exit(gcrankClearCache(), add = TRUE)
    config <- readMaMaMIAConfig(tiny_config_path())
    windows <- makeAnalysisWindows(
        tiny_assembly(config), tiny_contigs(), 100L
    )

    gc <- extractWindowGc(windows, tiny_assembly(config))

    expect_named(gc, c("chr_id", "start", "stop", "gc"))
    expect_identical(nrow(gc), 200L)
    expect_identical(gc$chr_id, windows$table$chr_id)
    expect_identical(gc$start, windows$table$start)
    expect_identical(gc$stop, windows$table$stop)
    expected <- ifelse(gc$chr_id == "r1", 0.6, 0.4)
    expect_equal(gc$gc, expected, tolerance = 1e-6)
    expect_false(anyNA(gc$gc))
    expect_true(all(gc$gc >= 0 & gc$gc <= 1))
})

test_that("extractWindowGc() handles a real two-index assembly", {
    skip_if_no_tiny_fixtures()
    on.exit(gcrankClearCache(), add = TRUE)
    copy_index <- function() {
        target <- tempfile("gcrank-copy-")
        dir.create(target)
        copied <- file.copy(
            list.files(tiny_gcrank_index(), full.names = TRUE),
            target, recursive = TRUE
        )
        if (!all(copied)) {
            stop("could not copy the tiny gcrank index")
        }
        target
    }
    first <- copy_index()
    second <- copy_index()
    on.exit(unlink(c(first, second), recursive = TRUE), add = TRUE)

    assembly <- two_index_assembly()
    assembly$contigs$recipient$index_path <- first
    assembly$contigs$donor$index_path <- second
    table <- data.frame(
        chr_id = c("r1", "d1", "r1", "d1"),
        start = c(0, 0, 1000, 1000),
        stop = c(1000, 1000, 2000, 2000),
        stringsAsFactors = FALSE
    )

    gc <- extractWindowGc(list(table = table), assembly)
    expect_equal(gc$gc, c(0.6, 0.4, 0.6, 0.4), tolerance = 1e-6)
})

test_that("extractWindowGc() calls gcrank once per index", {
    assembly <- two_index_assembly()
    table <- data.frame(
        chr_id = c("r1", "d1", "r1", "d1"),
        start = c(0, 0, 1000, 1000),
        stop = c(1000, 1000, 2000, 2000),
        stringsAsFactors = FALSE
    )
    calls <- new.env(parent = emptyenv())
    calls$args <- list()

    local_gcrank_mock(function(index_dir, contig, start, end, ...) {
        calls$args[[length(calls$args) + 1L]] <- list(
            index = index_dir, contig = contig,
            start = start, end = end
        )
        value <- if (identical(index_dir, "idx-recipient")) 0.25 else 0.75
        rep(value, length(contig))
    })

    gc <- extractWindowGc(list(table = table), assembly)

    expect_identical(gc$gc, c(0.25, 0.75, 0.25, 0.75))
    expect_identical(gc$start, c(0L, 0L, 1000L, 1000L))
    expect_length(calls$args, 2L)
    indexes <- vapply(calls$args, `[[`, character(1), "index")
    expect_setequal(indexes, c("idx-recipient", "idx-donor"))
    for (call in calls$args) {
        if (identical(call$index, "idx-recipient")) {
            expect_identical(call$contig, c("r1", "r1"))
            expect_identical(call$start, c(0, 1000))
        } else {
            expect_identical(call$contig, c("d1", "d1"))
            expect_identical(call$start, c(0, 1000))
        }
    }
})

test_that("extractWindowGc() rejects duplicate target contigs", {
    assembly <- two_index_assembly()
    assembly$contigs$donor$id <- "r1"
    table <- data.frame(
        chr_id = "r1", start = 0, stop = 1000, stringsAsFactors = FALSE
    )

    expect_error(
        extractWindowGc(list(table = table), assembly),
        "more than one role",
        class = "MaMaMIAWindowError"
    )
})

test_that("extractWindowGc() requires an index for every contig", {
    assembly <- two_index_assembly()
    assembly$contigs$donor$index_path <- NA_character_
    table <- data.frame(
        chr_id = c("r1", "d1"), start = 0, stop = 1000,
        stringsAsFactors = FALSE
    )

    expect_error(
        extractWindowGc(list(table = table), assembly),
        "no gc_index for contig 'd1'",
        class = "MaMaMIAGcrankError"
    )
})

test_that("extractWindowGc() validates window coordinates", {
    assembly <- two_index_assembly()
    make_windows <- function(chr_id, start, stop) {
        list(table = data.frame(
            chr_id = chr_id, start = start, stop = stop,
            stringsAsFactors = FALSE
        ))
    }

    expect_error(
        extractWindowGc(
            make_windows("unknown", 0, 100), assembly
        ),
        "no target contig 'unknown'",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        extractWindowGc(
            make_windows("r1", -1, 100), assembly
        ),
        "start must be non-negative",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        extractWindowGc(
            make_windows("r1", 100, 100), assembly
        ),
        "smaller than stop",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        extractWindowGc(
            make_windows("r1", 0, 20001), assembly
        ),
        "exceeds the contig length",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        extractWindowGc(
            make_windows("r1", 0.5, 100), assembly
        ),
        "whole numbers",
        class = "MaMaMIAWindowError"
    )
})

test_that("extractWindowGc() rejects missing and out-of-range GC", {
    assembly <- two_index_assembly()
    table <- list(table = data.frame(
        chr_id = "r1", start = 0, stop = 1000, stringsAsFactors = FALSE
    ))

    local_gcrank_mock(function(...) NA_real_)
    expect_error(
        extractWindowGc(table, assembly),
        "GC content is missing",
        class = "MaMaMIAGcrankError"
    )

    local_gcrank_mock(function(...) 1.5)
    expect_error(
        extractWindowGc(table, assembly),
        "outside \\[0, 1\\]",
        class = "MaMaMIAGcrankError"
    )

    local_gcrank_mock(function(...) stop("backend failure"))
    expect_error(
        extractWindowGc(table, assembly),
        "GC extraction failed for index 'idx-recipient': backend failure",
        class = "MaMaMIAGcrankError"
    )
})

test_that("extractWindowGc() handles empty windows", {
    gc <- extractWindowGc(
        list(table = data.frame(
            chr_id = character(0), start = numeric(0), stop = numeric(0),
            stringsAsFactors = FALSE
        )),
        two_index_assembly()
    )
    expect_named(gc, c("chr_id", "start", "stop", "gc"))
    expect_identical(nrow(gc), 0L)
    expect_identical(gc$gc, numeric(0))
})
