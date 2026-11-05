## Configuration: ploidy, multi-component assemblies and migration errors.

tiny_species_config <- function() {
    config <- list(
        version = 2L,
        defaults = list(windows_per_chromosome = 100L),
        species = list(
            tiny_recipient = list(
                label = "Tiny recipient",
                default_component = "R",
                gc_indexes = "tiny.gcrank",
                contigs = list(list(
                    id = "r1", name = "1R", component = "R"
                ))
            ),
            tiny_donor = list(
                label = "Tiny donor",
                default_component = "D",
                gc_indexes = "tiny.gcrank",
                contigs = list(list(
                    id = "d1", name = "1D", component = "D"
                ))
            )
        )
    )
    attr(config, "config_path") <- tiny_config_path()
    validateMaMaMIAConfig(config)
}

test_that("the tiny fixture is a valid configuration", {
    config <- tiny_config()
    expect_identical(config$version, 2L)
    expect_identical(config$defaults$metric, "breadth")
    expect_identical(config$defaults$score, "paired_log_ratio")
    expect_identical(config$defaults$donor_aggregation, "max")
    expect_identical(config$defaults$exclude_flags, 1796L)
    expect_identical(
        config$defaults$calibration$fdr, 0.01
    )
    assembly <- tiny_assembly(config)
    expect_identical(assembly$recipient_components, "R")
    expect_identical(assembly$donor_components, "D")
    expect_equal(unname(assembly$sample_composition), 2)
})

test_that("a single-component species can omit component labels", {
    entries <- .species_contig_entries(
        list(
            list(id = "chr1", name = "1"),
            list(id = "chr2", name = "2")
        ),
        "diploid"
    )
    expect_identical(entries$component, c("default", "default"))
    expect_identical(entries$group, c("1", "2"))
    expect_identical(entries$copy_number, c(2, 2))

    named <- .species_contig_entries(
        list(list(id = "chr1", name = "1")),
        "diploid", default_component = "main"
    )
    expect_identical(named$component, "main")

    mixed <- .species_contig_entries(
        list(
            list(id = "chr1", name = "1", component = "A"),
            list(id = "chr2", name = "2")
        ),
        "diploid"
    )
    expect_identical(mixed$component, c("A", "A"))

    expect_error(
        .species_contig_entries(
            list(
                list(id = "chr1", name = "1", component = "A"),
                list(id = "chr2", name = "1", component = "B"),
                list(id = "chr3", name = "2")
            ),
            "polyploid"
        ),
        "every contig must declare"
    )
})

test_that("a diploid species without component labels resolves", {
    config <- tiny_species_config()
    config$species$tiny_recipient$default_component <- NULL
    config$species$tiny_recipient$contigs <- list(
        list(id = "r1", name = "1")
    )
    config$species$tiny_donor$default_component <- NULL
    config$species$tiny_donor$contigs <- list(
        list(id = "d1", name = "1")
    )
    config <- validateMaMaMIAConfig(config)
    expect_identical(
        unique(config$species$tiny_recipient$contigs$component),
        "default"
    )
    expect_identical(
        unique(config$species$tiny_recipient$contigs$group), "1"
    )
    expect_identical(
        config$species$tiny_recipient$default_component, "default"
    )
    pair <- resolveSpeciesAssembly(
        config, donor = "tiny_donor", recipient = "tiny_recipient"
    )
    expect_identical(pair$recipient_components, "default")
    expect_identical(pair$donor_components, "default")
    expect_equal(unname(pair$assembly$sample_composition), 2)
})

test_that("resolution accepts ten-window steps between the allowed ticks", {
    config <- tiny_species_config()
    config$defaults$allowed_windows <- c(100L, 1000L, 5000L)
    config$defaults$windows_per_chromosome <- 1010L
    validated <- validateMaMaMIAConfig(config)
    expect_identical(validated$defaults$windows_per_chromosome, 1010L)

    config$defaults$windows_per_chromosome <- 1005L
    expect_error(validateMaMaMIAConfig(config), "steps of 10")
    config$defaults$windows_per_chromosome <- 5001L
    expect_error(validateMaMaMIAConfig(config), "between")
})

test_that("non-integer versions are rejected", {
    config <- list(version = 2.9, species = list())
    expect_error(validateMaMaMIAConfig(config), "integer")
    config$version <- "2.9"
    expect_error(validateMaMaMIAConfig(config), "integer")
    config$version <- c(2, 3)
    expect_error(validateMaMaMIAConfig(config), "integer")
    config$version <- list(2)
    expect_error(validateMaMaMIAConfig(config), "integer")
})

test_that("a component must agree on its copy number", {
    config <- tiny_species_config()
    config$species$tiny_donor$contigs <- list(
        list(id = "d1", name = "1D", component = "D", copy_number = 2),
        list(id = "d2", name = "2D", component = "D", copy_number = 4)
    )
    expect_error(validateMaMaMIAConfig(config), "copy_number")
})

test_that("calibration bins must end with Inf", {
    config <- tiny_species_config()
    config$defaults$calibration <- list(length_bins_mb = c(5, 10, 25))
    expect_error(validateMaMaMIAConfig(config), "end with Inf")
})

test_that("single-sample floors match the length bins", {
    config <- tiny_species_config()
    config$defaults$calibration <- list(single_sample_floors = c(1, 1))
    expect_error(
        validateMaMaMIAConfig(config),
        "one value per length bin"
    )
})

test_that("allowed_windows supplies a compatible default resolution", {
    config <- tiny_species_config()
    config$defaults$windows_per_chromosome <- NULL
    config$defaults$allowed_windows <- c(200, 400)
    validated <- validateMaMaMIAConfig(config)
    expect_identical(validated$defaults$windows_per_chromosome, 400L)

    config$defaults$windows_per_chromosome <- 200L
    validated <- validateMaMaMIAConfig(config)
    expect_identical(validated$defaults$windows_per_chromosome, 200L)
})

test_that("version 1 is rejected with migration guidance", {
    config <- list(version = 1L, species = list())
    expect_error(
        validateMaMaMIAConfig(config),
        "version 1 is no longer supported"
    )
    expect_error(validateMaMaMIAConfig(config), "MIGRATION.md")
})

test_that("unknown keys are rejected at every level", {
    config <- tiny_species_config()
    config$unknown_top <- 1
    expect_error(validateMaMaMIAConfig(config), "Unknown top-level")

    config <- tiny_species_config()
    config$species$tiny_recipient$gc_index <- "tiny.gcrank"
    expect_error(validateMaMaMIAConfig(config), "unknown key")

    config <- tiny_species_config()
    contigs <- config$species$tiny_recipient$contigs
    contigs$copy_numbr <- 4
    config$species$tiny_recipient$contigs <- contigs
    expect_error(validateMaMaMIAConfig(config), "copy_numbr")

    config <- tiny_species_config()
    config$server <- list(max_bam = 1)
    expect_error(validateMaMaMIAConfig(config), "Unknown `server`")
})

test_that("malformed gcrank manifest metadata is rejected", {
    source <- testthat::test_path("fixtures", "tiny.gcrank")
    make_index <- function(edit) {
        target <- tempfile("gcrank-manifest-")
        dir.create(target)
        file.copy(
            list.files(source, full.names = TRUE), target, recursive = TRUE
        )
        path <- file.path(target, "manifest.json")
        manifest <- jsonlite::fromJSON(path, simplifyVector = TRUE)
        manifest <- edit(manifest)
        jsonlite::write_json(manifest, path, auto_unbox = TRUE, pretty = TRUE)
        target
    }

    bad_count <- make_index(function(manifest) {
        manifest$source$num_contigs <- "2"
        manifest
    })
    on.exit(unlink(bad_count, recursive = TRUE), add = TRUE)
    config <- tiny_species_config()
    config$species$tiny_recipient$gc_indexes <- bad_count
    config$species$tiny_donor$gc_indexes <- bad_count
    expect_error(
        validateMaMaMIAConfig(config), "source.num_contigs"
    )

    bad_source <- make_index(function(manifest) {
        manifest$source <- "not a mapping"
        manifest
    })
    on.exit(unlink(bad_source, recursive = TRUE), add = TRUE)
    config$species$tiny_recipient$gc_indexes <- bad_source
    config$species$tiny_donor$gc_indexes <- bad_source
    expect_error(validateMaMaMIAConfig(config), "`source` mapping")
})

test_that("unknown defaults keys are rejected", {
    config <- tiny_species_config()
    config$defaults$fitted_windows <- 1000L
    expect_error(validateMaMaMIAConfig(config), "fitted_windows")
    config$defaults$fitted_windows <- NULL
    config$defaults$normalize_sides <- TRUE
    expect_error(validateMaMaMIAConfig(config), "normalize_sides")
})

test_that("metric, score and aggregation defaults are validated", {
    config <- tiny_species_config()
    config$defaults$metric <- "depth"
    expect_error(validateMaMaMIAConfig(config), "metric")
    config$defaults$metric <- NULL
    config$defaults$score <- "difference"
    expect_error(validateMaMaMIAConfig(config), "score")
    config$defaults$score <- NULL
    config$defaults$donor_aggregation <- "median"
    expect_error(validateMaMaMIAConfig(config), "donor_aggregation")
    config$defaults$donor_aggregation <- NULL
    config$defaults$recipient_aggregation <- "max"
    expect_error(validateMaMaMIAConfig(config), "recipient_aggregation")
})

test_that("mask and calibration defaults are validated", {
    config <- tiny_species_config()
    config$defaults$mask <- list(recurrence_threshold = 0.4)
    expect_error(validateMaMaMIAConfig(config), "Unknown `defaults` key")
    config$defaults$mask <- NULL
    config$defaults$calibration <- list(fdr = 0)
    expect_error(validateMaMaMIAConfig(config), "fdr")
    config$defaults$calibration <- list(length_bins_mb = c(10, 5))
    expect_error(validateMaMaMIAConfig(config), "increasing")
    config$defaults$calibration <- list(single_sample_floors = c(-1))
    expect_error(validateMaMaMIAConfig(config), "floors")
})

test_that("listMaMaMIAAssemblies and listMaMaMIASpecies expose ploidy", {
    config <- tiny_resolved_config()
    assemblies <- listMaMaMIAAssemblies(config)
    expect_true(all(c(
        "recipient_components", "donor_components"
    ) %in% names(assemblies)))
    expect_identical(assemblies$recipient_components, "R")
    species <- listMaMaMIASpecies(config)
    expect_true(all(c("components", "default_component") %in% names(species)))
})

test_that("resolveSpeciesAssembly accepts singular and plural components", {
    config <- tiny_species_config()
    pair <- resolveSpeciesAssembly(
        config, donor = "tiny_donor", recipient = "tiny_recipient",
        recipient_components = "R", donor_components = "D"
    )
    expect_identical(pair$recipient_components, "R")
    expect_identical(pair$donor_components, "D")
    expect_identical(nrow(pair$assembly$contigs$recipient), 1L)
    expect_identical(nrow(pair$assembly$contigs$donor), 1L)
    expect_identical(pair$assembly$contigs$recipient$group, "1")
    expect_identical(pair$assembly$contigs$recipient$copy_number, 2)
})

test_that("unknown components and copy numbers are rejected", {
    config <- tiny_species_config()
    expect_error(
        resolveSpeciesAssembly(
            config, donor = "tiny_donor", recipient = "tiny_recipient",
            recipient_components = "Z"
        ),
        "no recipient component"
    )
    config$species$tiny_recipient$contigs$copy_number <- -1
    expect_error(validateMaMaMIAConfig(config), "copy_number")
})

test_that("default_sample_composition validates names and values", {
    config <- tiny_species_config()
    pair <- resolveSpeciesAssembly(
        config, donor = "tiny_donor", recipient = "tiny_recipient",
        default_sample_composition = c(R = 4)
    )
    expect_equal(unname(pair$assembly$sample_composition), 4)
    expect_error(
        resolveSpeciesAssembly(
            config, donor = "tiny_donor", recipient = "tiny_recipient",
            default_sample_composition = c(Z = 2)
        ),
        "unknown recipient component"
    )
})

test_that("the derived assembly carries group and copy_number", {
    config <- tiny_config()
    assembly <- tiny_assembly(config)
    expect_true(all(c("group", "copy_number") %in% names(
        assembly$contigs$recipient
    )))
    windows <- makeAnalysisWindows(
        assembly, tiny_contigs(), n_per_chromosome = 100L
    )
    expect_true(all(c("group", "copy_number") %in% names(windows$table)))
    expect_identical(unique(windows$table$group), "1")
})

test_that("the server accepts one BAM by default", {
    config <- tiny_species_config()
    expect_identical(config$server$max_bams, 1L)
})

test_that("the shipped registry example is role-free", {
    path <- mamamia_path("config", "assemblies.example.yml")
    skip_if_not(nzchar(path) && file.exists(path),
        "the example registry is not available")
    example <- yaml::read_yaml(path, handlers = list(int = as.numeric))
    expect_identical(as.integer(example$version), 2L)
    expect_identical(as.integer(example$server$max_bams), 1L)
    expect_null(example$assemblies)
    expect_gte(length(example$species), 2L)
    for (species in example$species) {
        expect_true(all(c("label", "gc_indexes", "contigs") %in%
            names(species)))
        expect_true(all(vapply(species$contigs, function(contig) {
            all(c("id", "name") %in% names(contig))
        }, logical(1))))
    }
    expect_identical(example$defaults$exclude_flags, 1796)
})

test_that("zero-padded chromosome groups normalise to bare numbers", {
    expect_identical(
        chr_number(c("J01", "St01", "1A", "X")),
        c("1", "1", "1", NA_character_)
    )
    entries <- .species_contig_entries(
        list(
            list(id = "d1", name = "J01", component = "D"),
            list(id = "d2", name = "2D", component = "D"),
            list(id = "d3", name = "St03", component = "D", group = "03")
        ),
        "donor"
    )
    expect_identical(entries$group, c("1", "2", "3"))
})

test_that("group-matched homoeologs are not reported unpaired", {
    derived <- .species_derived_pairs(
        data.frame(id = "d1", group = "1", stringsAsFactors = FALSE),
        data.frame(
            id = c("r1", "r2", "r3"),
            group = c("1", "1", "2"),
            stringsAsFactors = FALSE
        )
    )
    expect_identical(derived$pairs[["d1"]], "r1")
    expect_identical(derived$unpaired$contig, "r3")
    expect_identical(
        derived$unpaired$reason, "no donor contig in group 2"
    )
})

test_that("readMaMaMIAConfig caches by path and invalidates on change", {
    source <- dirname(tiny_config_path())
    target <- tempfile("config-cache-")
    dir.create(target)
    file.copy(
        list.files(source, full.names = TRUE), target, recursive = TRUE
    )
    on.exit(unlink(target, recursive = TRUE), add = TRUE)
    path <- file.path(target, "assemblies.yml")

    first <- readMaMaMIAConfig(path)
    first$defaults$mapq <- 1L
    second <- readMaMaMIAConfig(path)
    expect_identical(second$defaults$mapq, 20L)

    lines <- readLines(path)
    lines <- sub("mapq: 20", "mapq: 21", lines, fixed = TRUE)
    lines <- c(lines, "# cache invalidation")
    writeLines(lines, path)
    third <- readMaMaMIAConfig(path)
    expect_identical(third$defaults$mapq, 21L)
})

test_that("declared assemblies are rejected by the registry schema", {
    config <- tiny_species_config()
    config$assemblies <- list(tiny = list(label = "tiny"))
    expect_error(
        validateMaMaMIAConfig(config),
        "not part of the registry schema"
    )
})

test_that("a registry needs at least two species", {
    config <- tiny_species_config()
    config$species <- config$species["tiny_recipient"]
    expect_error(validateMaMaMIAConfig(config), "at least two")
})
