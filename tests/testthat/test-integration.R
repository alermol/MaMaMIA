## P9: real 24-BAM validation on the Aegilops speltoides panel.

study_contigs <- function(root) {
    path <- file.path(root, "introgression-study", "data", "windows_n1000",
        "contigs.tsv")
    if (!file.exists(path)) {
        testthat::skip("study contig table is missing")
    }
    contigs <- utils::read.delim(path, stringsAsFactors = FALSE)
    contigs$group <- as.character(contigs$group)
    contigs$component <- ifelse(
        contigs$role == "donor",
        "S",
        sub("^[0-9]+", "", contigs$group)
    )
    contigs$group <- ifelse(
        contigs$role == "donor",
        contigs$group,
        sub("[A-Za-z]+$", "", contigs$group)
    )
    contigs
}

real_config <- function(root) {
    contigs <- study_contigs(root)
    speltoides_index <- file.path(
        root, "gcidx", "speltoides-study.gcrank"
    )
    if (!dir.exists(speltoides_index)) {
        testthat::skip(
            "the CH1S..CH7S gcrank index is missing (build it from the study GC table)"
        )
    }
    species <- list(
        triticum_aestivum = list(
            label = "Triticum aestivum",
            gc_indexes = file.path(
                root, "gcidx", "triticum-aestivum", "index.gcrank"
            ),
            contigs = lapply(
                seq_len(nrow(contigs[contigs$role == "recipient", ])),
                function(i) {
                    x <- contigs[contigs$role == "recipient", ][i, ]
                    list(
                        id = x$chr, name = x$component,
                        component = x$component, group = x$group,
                        copy_number = 2
                    )
                }
            )
        ),
        aegilops_speltoides = list(
            label = "Aegilops speltoides",
            gc_indexes = speltoides_index,
            contigs = lapply(
                seq_len(nrow(contigs[contigs$role == "donor", ])),
                function(i) {
                    x <- contigs[contigs$role == "donor", ][i, ]
                    list(
                        id = x$chr, name = paste0(x$group, "S"),
                        component = "S", group = x$group, copy_number = 2
                    )
                }
            )
        )
    )
    config <- list(
        version = 2L,
        defaults = list(
            windows_per_chromosome = 1000L,
            allowed_windows = c(100L, 1000L, 10000L),
            mapq = 20L,
            exclude_flags = 1796L,
            mask = list(physical_fallback_mb = 1),
            calibration = list(
                length_bins_mb = c(5, 10, 25, 50, 100, Inf),
                fdr = 0.01,
                single_sample_floors = c(2.0, 1.5, 1.2, 1.0, 0.5)
            )
        ),
        server = list(max_bams = 30L, max_total_bytes = 20 * 1024^3),
        species = species
    )
    attr(config, "config_path") <- file.path(root, "assemblies-study.yml")
    config <- validateMaMaMIAConfig(config)
    resolveSpeciesAssembly(
        config, donor = "aegilops_speltoides",
        recipient = "triticum_aestivum",
        recipient_components = c("A", "B", "D"),
        donor_components = "S",
        default_sample_composition = c(A = 2, B = 2, D = 2),
        key = "wheat_speltoides"
    )$config
}

real_bams <- function(root) {
    directory <- file.path(root, "bam_files", "aegilops-speltoides")
    if (!dir.exists(directory)) {
        testthat::skip("Aegilops speltoides BAM directory is missing")
    }
    paths <- list.files(directory, pattern = "[.]bam$", full.names = TRUE)
    if (length(paths) < 20L) {
        testthat::skip("fewer than 20 Aegilops speltoides BAMs are available")
    }
    paths
}

reference_truth <- function() {
    data.frame(
        sample = c("Sq74", "Sq75", "Sq75", "Sq75", "Sq75", "Sq76"),
        group = c("7", "7", "2", "6", "5", "7"),
        stringsAsFactors = FALSE
    )
}

test_that("the real assembly is ploidy-aware", {
    root <- skip_if_no_test_data()
    config <- real_config(root)
    expect_identical(config$version, 2L)
    expect_identical(
        config$assemblies$wheat_speltoides$recipient_components,
        c("A", "B", "D")
    )
    expect_identical(
        config$assemblies$wheat_speltoides$donor_components, "S"
    )
    expect_identical(
        nrow(config$assemblies$wheat_speltoides$contigs$recipient), 21L
    )
    expect_identical(
        nrow(config$assemblies$wheat_speltoides$contigs$donor), 7L
    )
    expect_identical(
        sort(unique(config$assemblies$wheat_speltoides$contigs$recipient$group)),
        as.character(1:7)
    )
})

test_that("the 24-BAM validation recovers all reference events", {
    skip_if_not(
        identical(Sys.getenv("MAMAMIA_SLOW_TESTS"), "1"),
        "set MAMAMIA_SLOW_TESTS=1 to run the 24-BAM validation"
    )
    root <- skip_if_no_test_data()
    config <- real_config(root)
    bams <- real_bams(root)
    result <- analyzeBamInput(
        bams, config,
        assembly = "wheat_speltoides",
        windows_per_chromosome = 1000L,
        cores = min(30L, max(1L, parallel::detectCores()))
    )
    calls <- result$calls[result$calls$pass, , drop = FALSE]
    truth <- reference_truth()
    truth$detected <- vapply(seq_len(nrow(truth)), function(i) {
        any(
            calls$sample == truth$sample[[i]] &
                calls$group == truth$group[[i]]
        )
    }, logical(1))
    expect_equal(sum(truth$detected), 6L)
    expect_lte(
        sum(calls$width_bp >= 20e6 &
            !calls$sample %in% c("Sq74", "Sq75", "Sq76")),
        2L
    )
    sq76_7 <- calls[calls$sample == "Sq76" & calls$group == "7", ]
    expect_gte(nrow(sq76_7), 1L)
    expect_lte(abs(min(sq76_7$start_iid) - 750), 20L)

    panel <- buildPanelMask(result$metrics, z = 0.5,
        recurrence_threshold = 0.4)
    expect_true(mean(panel$mask) > 0.1 && mean(panel$mask) < 0.2)
    expect_true(any(panel$mask[panel$group == "7" & panel$iid <= 10]))
    sq76_5 <- calls[calls$sample == "Sq76" & calls$group == "5", ]
    expect_gte(nrow(sq76_5), 1L)
    expect_lt(max(sq76_5$width_bp), 639356089)
})

test_that("single-sample analysis recovers the reference events", {
    skip_if_not(
        identical(Sys.getenv("MAMAMIA_SLOW_TESTS"), "1"),
        "set MAMAMIA_SLOW_TESTS=1 to run the single-sample validation"
    )
    root <- skip_if_no_test_data()
    config <- real_config(root)
    samples <- c("Sq74", "Sq75", "Sq76")
    calls <- do.call(rbind, lapply(samples, function(sample) {
        bam <- file.path(
            root, "bam_files", "aegilops-speltoides", paste0(sample, ".bam")
        )
        result <- analyzeBamInput(
            bam, config, assembly = "wheat_speltoides",
            windows_per_chromosome = 1000L, cores = 30L
        )
        result$calls[result$calls$pass, , drop = FALSE]
    }))
    truth <- reference_truth()
    truth$detected <- vapply(seq_len(nrow(truth)), function(i) {
        any(
            calls$sample == truth$sample[[i]] &
                calls$group == truth$group[[i]]
        )
    }, logical(1))
    expect_equal(sum(truth$detected), 6L)
})
