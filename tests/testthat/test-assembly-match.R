## T1.3: BAM header to assembly matching.

tiny_header <- function(names, lengths = rep(20000, length(names)),
                       path = "/tmp/test.bam") {
    list(
        path = path,
        format = "BAM",
        targets = data.frame(
            name = names,
            length = lengths,
            stringsAsFactors = FALSE
        )
    )
}

with_second_assembly <- function(config) {
    second <- tiny_assembly(config)
    second$key <- "tiny2"
    second$label <- "Tiny synthetic assembly with a second donor"
    second$contigs$donor <- transform(
        second$contigs$donor,
        id = "d2", name = "2D"
    )
    config$assemblies$tiny2 <- second
    config
}

with_superset_assembly <- function(config) {
    superset <- tiny_assembly(config)
    superset$key <- "tiny_full"
    superset$label <- "Tiny assembly with two donor contigs"
    superset$contigs$donor <- rbind(
        superset$contigs$donor,
        transform(superset$contigs$donor, id = "d2", name = "2D")
    )
    config$assemblies$tiny_full <- superset
    config
}

test_that("matchAssembly() selects an exact match and ignores extra contigs", {
    skip_if_no_tiny_fixtures()
    config <- tiny_resolved_config()

    exact <- matchAssembly(tiny_header(c("d1", "r1")), config)
    expect_identical(exact$assembly, "tiny")
    expect_true(exact$exact)
    expect_true(exact$matched)
    expect_identical(exact$score, 1)
    expect_length(exact$missing, 0L)
    expect_identical(nrow(exact$length_mismatches), 0L)
    expect_true(all(c(
        "key", "label", "n_targets", "n_matched", "n_missing",
        "n_length_mismatches", "score", "exact"
    ) %in% names(exact$candidates)))
    expect_match(exact$summary, "matches assembly 'tiny' exactly")

    extra <- matchAssembly(
        tiny_header(c("unplaced1", "d1", "r1", "unplaced2"),
            lengths = c(100, 20000, 20000, 200)
        ),
        config
    )
    expect_identical(extra$assembly, "tiny")
    expect_true(extra$exact)
})

test_that("matchAssembly() distinguishes assemblies that share a recipient", {
    skip_if_no_tiny_fixtures()
    config <- with_second_assembly(tiny_resolved_config())

    expect_identical(
        matchAssembly(tiny_header(c("d1", "r1")), config)$assembly,
        "tiny"
    )
    expect_identical(
        matchAssembly(tiny_header(c("d2", "r1")), config)$assembly,
        "tiny2"
    )
    expect_error(
        matchAssembly(tiny_header("r1"), config),
        "does not match any configured assembly",
        class = "MaMaMIAAssemblyMatchError"
    )
})

test_that("matchAssembly() prefers the most specific exact match", {
    skip_if_no_tiny_fixtures()
    config <- with_superset_assembly(tiny_resolved_config())

    match <- matchAssembly(tiny_header(c("d1", "d2", "r1")), config)
    expect_identical(match$assembly, "tiny_full")
    expect_identical(match$candidates$n_targets[1L], 3L)
})

test_that("matchAssembly() reports ambiguity between equal exact matches", {
    skip_if_no_tiny_fixtures()
    config <- with_second_assembly(tiny_resolved_config())

    expect_error(
        matchAssembly(tiny_header(c("d1", "d2", "r1")), config),
        "matches more than one assembly: tiny, tiny2",
        class = "MaMaMIAAssemblyMatchError"
    )
})

test_that("matchAssembly() validates an explicit assembly key", {
    skip_if_no_tiny_fixtures()
    config <- tiny_resolved_config()

    expect_identical(
        matchAssembly(tiny_header(c("d1", "r1")), config,
            assembly = "tiny"
        )$assembly,
        "tiny"
    )
    expect_error(
        matchAssembly(tiny_header(c("d1", "r1")), config,
            assembly = "unknown"
        ),
        "Unknown assembly 'unknown'",
        class = "MaMaMIAAssemblyMatchError"
    )
})

test_that("matchAssembly() reports length mismatches with contig context", {
    skip_if_no_tiny_fixtures()
    config <- tiny_resolved_config()
    header <- tiny_header(c("d1", "r1"), lengths = c(19999, 20000))

    error <- tryCatch(
        matchAssembly(header, config, assembly = "tiny"),
        error = identity
    )
    expect_s3_class(error, "MaMaMIAHeaderMismatchError")
    expect_match(
        conditionMessage(error),
        "target 'd1' has length 19999; assembly 'tiny' expects 20000"
    )

    auto_error <- tryCatch(matchAssembly(header, config), error = identity)
    expect_s3_class(auto_error, "MaMaMIAAssemblyMatchError")
    expect_match(
        conditionMessage(auto_error),
        "does not match any configured assembly; closest match: 'tiny'"
    )
    expect_identical(auto_error$data$candidates$n_length_mismatches[1L], 1L)
})

test_that("matchAssembly() reports a BAM that covers only one side", {
    skip_if_no_tiny_fixtures()
    config <- tiny_resolved_config()

    error <- tryCatch(
        matchAssembly(tiny_header("r1"), config, assembly = "tiny"),
        error = identity
    )
    expect_s3_class(error, "MaMaMIAAssemblyMatchError")
    expect_match(conditionMessage(error), "covers only one side of assembly 'tiny'")
})

test_that("matchAssembly() reports partial candidates when nothing matches", {
    skip_if_no_tiny_fixtures()
    config <- tiny_resolved_config()

    error <- tryCatch(
        matchAssembly(tiny_header("d1"), config),
        error = identity
    )
    expect_s3_class(error, "MaMaMIAAssemblyMatchError")
    candidates <- error$data$candidates
    expect_identical(candidates$key, "tiny")
    expect_identical(candidates$score, 0.5)
    expect_identical(candidates$n_matched, 1L)
    expect_identical(candidates$n_missing, 1L)
})

test_that("matchAssemblyCohort() requires one assembly for all BAMs", {
    skip_if_no_tiny_fixtures()
    config <- with_second_assembly(tiny_resolved_config())

    matches <- matchAssemblyCohort(
        list(
            tiny_header(c("d1", "r1"), path = "/tmp/a.bam"),
            tiny_header(c("d1", "r1"), path = "/tmp/b.bam")
        ),
        config
    )
    expect_length(matches, 2L)
    expect_identical(
        vapply(matches, `[[`, character(1), "assembly"),
        c("tiny", "tiny")
    )

    expect_error(
        matchAssemblyCohort(
            list(
                tiny_header(c("d1", "r1"), path = "/tmp/a.bam"),
                tiny_header(c("d2", "r1"), path = "/tmp/b.bam")
            ),
            config
        ),
        "match different assemblies.*a[.]bam.*tiny.*b[.]bam.*tiny2",
        class = "MaMaMIAAssemblyMatchError"
    )
})

test_that("matchAssembly() matches the real Sq54 BAM exactly", {
    root <- skip_if_no_test_data()
    bam_path <- file.path(root, "bam_files", "Sq54", "cs_t2t_se_cereale.bam")
    if (!file.exists(bam_path)) {
        testthat::skip("Sq54 BAM is missing")
    }

    chromosomes <- utils::read.delim(
        file.path(root, "assemblies", "cs_t2t_se_cereale", "chromosomes.tsv"),
        header = FALSE,
        col.names = c("id", "name", "component"),
        stringsAsFactors = FALSE
    )
    entries <- function(x) {
        lapply(seq_len(nrow(x)), function(i) {
            list(id = x$id[[i]], name = x$name[[i]],
                component = x$component[[i]]
            )
        })
    }
    config <- validateMaMaMIAConfig(list(
        version = 2L,
        species = list(
            triticum_aestivum = list(
                label = "Triticum aestivum CS T2T",
                gc_indexes = list(
                    file.path(root, "gcidx", "triticum-aestivum", "index.gcrank")
                ),
                contigs = entries(
                    chromosomes[chromosomes$component != "R", ]
                )
            ),
            secale_cereale = list(
                label = "Secale cereale",
                gc_indexes = list(
                    file.path(root, "gcidx", "secale-cereale", "index.gcrank")
                ),
                contigs = entries(
                    chromosomes[chromosomes$component == "R", ]
                )
            )
        )
    ))
    config <- resolveSpeciesAssembly(
        config, donor = "secale_cereale", recipient = "triticum_aestivum",
        key = "cs_t2t_se_cereale"
    )$config

    match <- matchAssembly(bamHeaderInfo(bam_path), config)
    expect_identical(match$assembly, "cs_t2t_se_cereale")
    expect_true(match$exact)
    expect_identical(match$score, 1)
    expect_identical(nrow(match$length_mismatches), 0L)
})
