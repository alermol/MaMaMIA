## T3.2: windows + GC + counts -> SampleTables -> valid RCA.

sample_tables_fixture <- function(n = 10L) {
    skip_if_no_tiny_fixtures()
    config <- readMaMaMIAConfig(tiny_config_path())
    assembly <- config$assemblies$tiny
    windows <- makeAnalysisWindows(assembly, tiny_contigs(), n)
    gc <- MaMaMIA:::extractWindowGc(windows, assembly)
    ## makeAnalysisWindows() orders recipient (r1) first, then donor (d1).
    expect_identical(unique(windows$table$chr_id), c("r1", "d1"))
    list(
        assembly = assembly,
        windows = windows,
        gc = gc,
        counts = c(rep(5L, n), rep(7L, n))
    )
}

sample_tables_without_pairs <- function(assembly) {
    assembly$pairs <- stats::setNames(character(0), character(0))
    assembly
}

test_that("buildSampleTables() builds the documented SampleTables contract", {
    fixture <- sample_tables_fixture()
    on.exit(gcrankClearCache(), add = TRUE)

    tables <- buildSampleTables(
        fixture$counts, fixture$gc, fixture$windows,
        fixture$assembly, "sample-1"
    )

    expect_s3_class(tables, "SampleTables")
    expect_named(tables, c(
        "dcov", "dgc", "rcov", "rgc", "meta",
        "don_chr_ids", "rec_chr_ids", "sample", "reference", "window",
        "target_pairs"
    ))

    expect_named(tables$dcov, c("chr_id", "start", "stop", "cov"))
    expect_named(tables$dgc, c("chr_id", "start", "stop", "gc"))
    expect_named(tables$rcov, c("chr_id", "start", "stop", "cov"))
    expect_named(tables$rgc, c("chr_id", "start", "stop", "gc"))
    expect_identical(nrow(tables$dcov), 10L)
    expect_identical(nrow(tables$rcov), 10L)
    expect_identical(unique(tables$dcov$chr_id), "d1")
    expect_identical(unique(tables$rcov$chr_id), "r1")
    expect_equal(tables$dcov$cov, rep(7, 10L))
    expect_equal(tables$rcov$cov, rep(5, 10L))
    expect_equal(tables$dgc$gc, rep(0.4, 10L), tolerance = 1e-6)
    expect_equal(tables$rgc$gc, rep(0.6, 10L), tolerance = 1e-6)

    expect_identical(tables$meta$chr_id, c("r1", "d1"))
    expect_identical(tables$meta$chr_name, c("1R", "1D"))
    expect_identical(tables$meta$subgenome, c("R", "D"))
    expect_identical(tables$don_chr_ids, "d1")
    expect_identical(tables$rec_chr_ids, "r1")
    expect_identical(tables$sample, "sample-1")
    expect_identical(tables$reference, "tiny")
    expect_identical(tables$window, 10L)
    expect_identical(tables$target_pairs, c(d1 = "r1"))

    rca <- do.call(RCA, tables[c("dcov", "dgc", "rcov", "rgc", "meta")])
    expect_s3_class(rca, "RCA")
    expect_false(rca$corrected)
    iid <- split(rca$data$iid, rca$data$chr_id)
    expect_true(all(vapply(iid, function(x) identical(as.integer(x), 1:10), logical(1L))))
})

test_that("buildSampleTables() accepts donor/recipient and vector inputs", {
    fixture <- sample_tables_fixture()
    on.exit(gcrankClearCache(), add = TRUE)

    by_role <- buildSampleTables(
        list(
            donor = rep(7L, 10L),
            recipient = rep(5L, 10L)
        ),
        fixture$gc,
        fixture$windows,
        fixture$assembly,
        "sample-1"
    )
    by_vector <- buildSampleTables(
        fixture$counts,
        fixture$gc$gc,
        fixture$windows,
        fixture$assembly,
        "sample-1"
    )

    expect_identical(by_role$dcov, by_vector$dcov)
    expect_identical(by_role$rcov, by_vector$rcov)
    expect_identical(by_role$dgc, by_vector$dgc)
    expect_identical(by_role$rgc, by_vector$rgc)
})

test_that("explicit assembly pairs win over the chromosome-number rule", {
    fixture <- sample_tables_fixture()
    on.exit(gcrankClearCache(), add = TRUE)

    ## A display name that the number rule cannot pair with 1R.
    windows <- fixture$windows
    windows$table$chr_name[windows$table$chr_id == "d1"] <- "9D"

    explicit <- buildSampleTables(
        fixture$counts, fixture$gc, windows,
        fixture$assembly, "sample-1"
    )
    expect_identical(explicit$target_pairs, c(d1 = "r1"))

    expect_error(
        buildSampleTables(
            fixture$counts, fixture$gc, windows,
            sample_tables_without_pairs(fixture$assembly), "sample-1"
        ),
        "no explicit assembly pairs",
        class = "MaMaMIASampleError"
    )
})

test_that("buildSampleTables() derives pairs when the assembly has none", {
    fixture <- sample_tables_fixture()
    on.exit(gcrankClearCache(), add = TRUE)

    tables <- buildSampleTables(
        fixture$counts,
        fixture$gc,
        fixture$windows,
        sample_tables_without_pairs(fixture$assembly),
        "sample-1"
    )
    expect_identical(tables$target_pairs, c(d1 = "r1"))
})

test_that("buildSampleTables() rejects missing roles and misaligned values", {
    fixture <- sample_tables_fixture()
    on.exit(gcrankClearCache(), add = TRUE)

    recipient_only <- fixture$windows$table[
        fixture$windows$table$role == "recipient", ,
        drop = FALSE
    ]
    expect_error(
        buildSampleTables(
            rep(5L, 10L), fixture$gc$gc[
                fixture$windows$table$role == "recipient"
            ],
            recipient_only, fixture$assembly, "sample-1"
        ),
        "no donor windows",
        class = "MaMaMIASampleError"
    )

    expect_error(
        buildSampleTables(
            fixture$counts[-1L], fixture$gc, fixture$windows,
            fixture$assembly, "sample-1"
        ),
        "length 20",
        class = "MaMaMIASampleError"
    )

    bad_gc <- fixture$gc
    bad_gc$start[[1L]] <- bad_gc$start[[1L]] + 1L
    expect_error(
        buildSampleTables(
            fixture$counts, bad_gc, fixture$windows,
            fixture$assembly, "sample-1"
        ),
        "not aligned",
        class = "MaMaMIASampleError"
    )

    bad_counts <- fixture$counts
    bad_counts[[1L]] <- -1L
    expect_error(
        buildSampleTables(
            bad_counts, fixture$gc, fixture$windows,
            fixture$assembly, "sample-1"
        ),
        "non-negative",
        class = "MaMaMIASampleError"
    )
})

test_that("buildSampleTables() requires equal window counts", {
    fixture <- sample_tables_fixture()
    on.exit(gcrankClearCache(), add = TRUE)

    keep <- -which(fixture$windows$table$chr_id == "d1")[[1L]]
    expect_error(
        buildSampleTables(
            fixture$counts[keep],
            fixture$gc$gc[keep],
            fixture$windows$table[keep, , drop = FALSE],
            fixture$assembly, "sample-1"
        ),
        "same number of windows",
        class = "MaMaMIAWindowError"
    )
})

test_that("buildSampleTables() rejects inconsistent roles and pairs", {
    fixture <- sample_tables_fixture()
    on.exit(gcrankClearCache(), add = TRUE)

    mixed <- fixture$windows$table
    mixed$role[mixed$chr_id == "d1"][[1L]] <- "recipient"
    expect_error(
        buildSampleTables(
            fixture$counts, fixture$gc, mixed,
            fixture$assembly, "sample-1"
        ),
        "one role",
        class = "MaMaMIAWindowError"
    )

    bad_pairs <- fixture$assembly
    bad_pairs$pairs <- c(d1 = "r1", d1 = "r1")
    expect_error(
        buildSampleTables(
            fixture$counts, fixture$gc, fixture$windows,
            bad_pairs, "sample-1"
        ),
        "duplicate donor",
        class = "MaMaMIASampleError"
    )

    unknown_pairs <- fixture$assembly
    unknown_pairs$pairs <- c(d1 = "unknown")
    expect_error(
        buildSampleTables(
            fixture$counts, fixture$gc, fixture$windows,
            unknown_pairs, "sample-1"
        ),
        "unknown recipient",
        class = "MaMaMIASampleError"
    )
})

test_that("RCA validation failures carry the sample name", {
    fixture <- sample_tables_fixture()
    on.exit(gcrankClearCache(), add = TRUE)

    gapped <- fixture$windows$table
    idx <- which(gapped$chr_id == "d1")[[5L]]
    gapped$start[[idx]] <- gapped$start[[idx]] + 100L
    gapped$stop[[idx]] <- gapped$stop[[idx]] - 100L

    expect_error(
        buildSampleTables(
            fixture$counts, fixture$gc$gc, list(
                table = gapped,
                n_per_chromosome = 10L
            ),
            fixture$assembly, "broken-sample"
        ),
        "Sample 'broken-sample'.*gaps",
        class = "MaMaMIASampleError"
    )
})
