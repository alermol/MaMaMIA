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

test_that("version 1 is rejected with migration guidance", {
    config <- list(version = 1L, species = list())
    expect_error(
        validateMaMaMIAConfig(config),
        "version 1 is no longer supported"
    )
    expect_error(validateMaMaMIAConfig(config), "MIGRATION.md")
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
    config$defaults$mask <- list(recurrence_threshold = 2)
    expect_error(validateMaMaMIAConfig(config), "recurrence_threshold")
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

test_that("the shipped registry example is role-free", {
    path <- system.file("config", "assemblies.example.yml",
        package = "MaMaMIA")
    skip_if_not(nzchar(path) && file.exists(path),
        "the example registry is not installed")
    example <- yaml::read_yaml(path, handlers = list(int = as.numeric))
    expect_identical(as.integer(example$version), 2L)
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
