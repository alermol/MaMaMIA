## T2.4: route windows to gcrank indexes and build the complete GC table.

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

    gc <- MaMaMIA:::extractWindowGc(windows, tiny_assembly(config))

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

    gc <- MaMaMIA:::extractWindowGc(list(table = table), assembly)
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

    testthat::local_mocked_bindings(
        gcrankGc = function(index_dir, contig, start, end, ...) {
            calls$args[[length(calls$args) + 1L]] <- list(
                index = index_dir, contig = contig,
                start = start, end = end
            )
            value <- if (identical(index_dir, "idx-recipient")) 0.25 else 0.75
            rep(value, length(contig))
        },
        .package = "MaMaMIA"
    )

    gc <- MaMaMIA:::extractWindowGc(list(table = table), assembly)

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
        MaMaMIA:::extractWindowGc(list(table = table), assembly),
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
        MaMaMIA:::extractWindowGc(list(table = table), assembly),
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
        MaMaMIA:::extractWindowGc(
            make_windows("unknown", 0, 100), assembly
        ),
        "no target contig 'unknown'",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        MaMaMIA:::extractWindowGc(
            make_windows("r1", -1, 100), assembly
        ),
        "start must be non-negative",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        MaMaMIA:::extractWindowGc(
            make_windows("r1", 100, 100), assembly
        ),
        "smaller than stop",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        MaMaMIA:::extractWindowGc(
            make_windows("r1", 0, 20001), assembly
        ),
        "exceeds the contig length",
        class = "MaMaMIAWindowError"
    )
    expect_error(
        MaMaMIA:::extractWindowGc(
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

    testthat::local_mocked_bindings(
        gcrankGc = function(...) NA_real_,
        .package = "MaMaMIA"
    )
    expect_error(
        MaMaMIA:::extractWindowGc(table, assembly),
        "GC content is missing",
        class = "MaMaMIAGcrankError"
    )

    testthat::local_mocked_bindings(
        gcrankGc = function(...) 1.5,
        .package = "MaMaMIA"
    )
    expect_error(
        MaMaMIA:::extractWindowGc(table, assembly),
        "outside \\[0, 1\\]",
        class = "MaMaMIAGcrankError"
    )

    testthat::local_mocked_bindings(
        gcrankGc = function(...) stop("backend failure"),
        .package = "MaMaMIA"
    )
    expect_error(
        MaMaMIA:::extractWindowGc(table, assembly),
        "GC extraction failed for index 'idx-recipient': backend failure",
        class = "MaMaMIAGcrankError"
    )
})

test_that("extractWindowGc() handles empty windows", {
    gc <- MaMaMIA:::extractWindowGc(
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

test_that("real GC windows match the reference BED after rounding", {
    root <- skip_if_no_test_data()
    recipient_bed <- file.path(
        root, "assemblies", "cs_t2t_se_cereale", "recipient.w1000.gc.bed"
    )
    donor_bed <- file.path(
        root, "assemblies", "cs_t2t_se_cereale", "donor.w1000.gc.bed"
    )
    wheat_index <- file.path(
        root, "gcidx", "triticum-aestivum", "index.gcrank"
    )
    rye_index <- file.path(root, "gcidx", "secale-cereale", "index.gcrank")
    if (!all(file.exists(c(recipient_bed, donor_bed))) ||
        !dir.exists(wheat_index) || !dir.exists(rye_index)) {
        testthat::skip("real GC indexes or reference GC tables are missing")
    }
    on.exit(gcrankClearCache(), add = TRUE)

    read_reference_gc_bed <- function(path) {
        utils::read.delim(
            path,
            header = FALSE,
            comment.char = "#",
            col.names = c("chr_id", "start", "stop", "gc"),
            colClasses = c("character", "numeric", "numeric", "numeric"),
            stringsAsFactors = FALSE
        )
    }
    reference <- rbind(
        read_reference_gc_bed(recipient_bed),
        read_reference_gc_bed(donor_bed)
    )

    chromosomes <- utils::read.delim(
        file.path(root, "assemblies", "cs_t2t_se_cereale", "chromosomes.tsv"),
        header = FALSE,
        col.names = c("id", "name", "component"),
        stringsAsFactors = FALSE
    )
    entries <- function(x) {
        lapply(seq_len(nrow(x)), function(i) {
            list(
                id = x$id[[i]], name = x$name[[i]],
                component = x$component[[i]]
            )
        })
    }
    config <- validateMaMaMIAConfig(list(
        version = 2L,
        species = list(
            triticum_aestivum = list(
                label = "Triticum aestivum CS T2T",
                gc_indexes = list(wheat_index),
                contigs = entries(
                    chromosomes[chromosomes$component != "R", ]
                )
            ),
            secale_cereale = list(
                label = "Secale cereale",
                gc_indexes = list(rye_index),
                contigs = entries(
                    chromosomes[chromosomes$component == "R", ]
                )
            )
        )
    ))
    assembly <- resolveSpeciesAssembly(
        config, donor = "secale_cereale", recipient = "triticum_aestivum",
        key = "cs_t2t_se_cereale"
    )$assembly

    ## The reference BED tables use their own tiling, so compare at those
    ## exact intervals. makeAnalysisWindows() distributes the remainder
    ## (T2.2), which differs from the reference width after the first window
    ## on most contigs.
    gc <- MaMaMIA:::extractWindowGc(
        list(table = reference[, c("chr_id", "start", "stop")]),
        assembly
    )

    expect_identical(gc$chr_id, reference$chr_id)
    expect_identical(gc$start, as.integer(reference$start))
    expect_identical(gc$stop, as.integer(reference$stop))
    expect_identical(
        sprintf("%.6f", gc$gc),
        sprintf("%.6f", reference$gc)
    )
})
