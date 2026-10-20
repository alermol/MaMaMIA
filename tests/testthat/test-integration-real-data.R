## T0.4: gated integration checks against the real Sq54 data.
##
## These tests only read headers, manifests and small contig tables. The full
## BAM pipeline regression is added by T3.3/T3.4; this file proves that the
## real-data fixtures are discoverable and consistent with the registry format.

real_assembly_config <- function(root) {
    index_paths <- c(
        normalizePath(
            file.path(root, "gcidx", "triticum-aestivum", "index.gcrank"),
            winslash = "/", mustWork = TRUE
        ),
        normalizePath(
            file.path(root, "gcidx", "secale-cereale", "index.gcrank"),
            winslash = "/", mustWork = TRUE
        )
    )
    chromosomes <- utils::read.delim(
        file.path(root, "assemblies", "cs_t2t_se_cereale", "chromosomes.tsv"),
        header = FALSE,
        col.names = c("id", "name", "subgenome"),
        stringsAsFactors = FALSE
    )
    entries <- function(x) {
        lapply(seq_len(nrow(x)), function(i) {
            list(id = x$id[[i]], name = x$name[[i]], subgenome = x$subgenome[[i]])
        })
    }
    list(
        version = 1L,
        assemblies = list(cs_t2t_se_cereale = list(
            label = "Triticum aestivum CS T2T + Secale cereale",
            description = "Real-data integration fixture",
            gc_indexes = as.list(index_paths),
            contigs = list(
                recipient = entries(chromosomes[chromosomes$subgenome != "R", ]),
                donor = entries(chromosomes[chromosomes$subgenome == "R", ])
            ),
            pairs = list(
                OZ284382.1 = "CP155614.1",
                OZ284383.1 = "CP155611.1",
                OZ284384.1 = "CP155608.1",
                OZ284385.1 = "CP155605.1",
                OZ284386.1 = "CP155602.1",
                OZ284387.1 = "CP155599.1",
                OZ284388.1 = "CP155597.1"
            )
        ))
    )
}

test_that("the real Sq54 BAM and legacy tables are discoverable", {
    root <- skip_if_no_test_data()
    bam_path <- file.path(root, "bam_files", "Sq54", "cs_t2t_se_cereale.bam")
    legacy_rcov <- file.path(
        root, "mamamia_analysis", "arsenal", "Sq54",
        "cs_t2t_se_cereale.w1000.rcov.bed"
    )
    legacy_dcov <- file.path(
        root, "mamamia_analysis", "arsenal", "Sq54",
        "cs_t2t_se_cereale.w1000.dcov.bed"
    )
    if (!all(file.exists(c(bam_path, legacy_rcov, legacy_dcov)))) {
        testthat::skip("Sq54 BAM or legacy BED tables are missing")
    }

    expect_gt(file.info(bam_path)$size, 0)
    expect_gt(file.info(legacy_rcov)$size, 0)
    expect_gt(file.info(legacy_dcov)$size, 0)
})

test_that("the real assembly registry matches the Sq54 header", {
    root <- skip_if_no_test_data()
    bam_path <- file.path(root, "bam_files", "Sq54", "cs_t2t_se_cereale.bam")
    if (!file.exists(bam_path)) {
        testthat::skip("Sq54 BAM is missing")
    }

    config_path <- file.path(tempdir(), "assemblies-real.yml")
    writeLines(
        yaml::as.yaml(real_assembly_config(root)),
        config_path
    )
    config <- readMaMaMIAConfig(config_path)
    assembly <- config$assemblies$cs_t2t_se_cereale
    expect_identical(length(assembly$pairs), 7L)

    header <- MaMaMIA:::bam_header_info_cpp(bam_path)
    header_lengths <- stats::setNames(
        header$targets$length, header$targets$name
    )
    configured <- rbind(
        assembly$contigs$recipient, assembly$contigs$donor
    )
    expect_true(all(configured$id %in% names(header_lengths)))
    expect_equal(
        unname(header_lengths[configured$id]),
        configured$length
    )
    expect_identical(header$n_targets, length(header$targets$name))
})
