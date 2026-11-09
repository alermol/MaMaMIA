## Candidate 5: resolve a registry and a species pair into one runnable context.

test_that("a registry path plus species pair opens a runnable context", {
    skip_if_no_tiny_fixtures()
    context <- openAnalysisConfig(
        tiny_config_path(),
        donor = "tiny_donor",
        recipient = "tiny_recipient"
    )
    expect_s3_class(context, "MaMaMIAAnalysisContext")
    expect_named(context, c(
        "config", "assembly_key", "assembly", "donor_species",
        "recipient_species", "donor_components", "recipient_components",
        "unpaired"
    ))
    expect_identical(
        context$assembly_key, "tiny_recipient__R__tiny_donor__D"
    )
    expect_identical(
        context$assembly, context$config$assemblies[[context$assembly_key]]
    )
    expect_identical(context$recipient_species, "tiny_recipient")
    expect_identical(context$donor_species, "tiny_donor")
    expect_identical(context$recipient_components, "R")
    expect_identical(context$donor_components, "D")
    expect_identical(
        names(analysisContextConfig(context)), names(context$config)
    )
    expect_identical(
        analysisContextAssemblyKey(context), context$assembly_key
    )
})

test_that("a validated registry list opens the same context", {
    skip_if_no_tiny_fixtures()
    path_context <- openAnalysisConfig(
        tiny_config_path(),
        donor = "tiny_donor", recipient = "tiny_recipient"
    )
    list_context <- openAnalysisConfig(
        tiny_config(), donor = "tiny_donor", recipient = "tiny_recipient"
    )
    expect_identical(list_context$assembly_key, path_context$assembly_key)
    expect_identical(list_context$assembly, path_context$assembly)
})

test_that("a resolved config may select its only assembly", {
    skip_if_no_tiny_fixtures()
    resolved <- tiny_resolved_config()
    context <- openAnalysisConfig(resolved)
    expect_identical(context$assembly_key, "tiny")
    expect_identical(context$config, resolved)

    keyed <- openAnalysisConfig(resolved, assembly = "tiny")
    expect_identical(keyed$assembly_key, "tiny")
})

test_that("unknown species and assemblies name the available choices", {
    skip_if_no_tiny_fixtures()
    expect_error(
        openAnalysisConfig(
            tiny_config_path(),
            donor = "missing_donor", recipient = "tiny_recipient"
        ),
        "Unknown donor species 'missing_donor'.*tiny_donor",
        class = "MaMaMIAConfigError"
    )
    expect_error(
        openAnalysisConfig(
            tiny_config_path(),
            donor = "tiny_donor", recipient = "missing_recipient"
        ),
        "Unknown recipient species 'missing_recipient'.*tiny_recipient",
        class = "MaMaMIAConfigError"
    )
    resolved <- tiny_resolved_config()
    expect_error(
        openAnalysisConfig(resolved, assembly = "missing"),
        "Unknown assembly 'missing'.*tiny",
        class = "MaMaMIAConfigError"
    )
    expect_error(
        openAnalysisConfig(tiny_config(), donor = "tiny_donor"),
        "Pass both `donor` and `recipient`",
        class = "MaMaMIAConfigError"
    )
})

test_that("the pipeline fills missing defaults without re-validating assemblies", {
    skip_if_no_tiny_fixtures()
    config <- tiny_resolved_config()
    config$defaults <- NULL
    config$server <- NULL
    filled <- .pipeline_config(config)
    expect_identical(filled$defaults$metric, "breadth")
    expect_identical(filled$server$max_bams, 1L)
    expect_identical(names(filled$assemblies), "tiny")
})

test_that("preflight and the pipeline accept an analysis context", {
    skip_if_no_tiny_fixtures()
    context <- openAnalysisConfig(
        tiny_config_path(),
        donor = "tiny_donor", recipient = "tiny_recipient"
    )
    table <- preflightBams(tiny_bam_path(), context)
    expect_identical(table$assembly, context$assembly_key)

    result <- analyzeBamInput(
        tiny_bam_path(), context,
        windows_per_chromosome = 100L,
        calibration = list(
            length_bins_mb = c(0.001, 0.01, Inf),
            floors = c(0.1, 0.1),
            min_mb = 0.001
        ),
        cores = 1L
    )
    expect_s3_class(result, "BamAnalysisResult")
    expect_identical(result$samples$sample, "tiny-sample")
    expect_identical(
        names(result$config$assemblies), context$assembly_key
    )
    expect_identical(result$assembly$key, context$assembly_key)
})
