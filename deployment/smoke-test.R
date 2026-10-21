config <- MaMaMIA::readMaMaMIAConfig(Sys.getenv("MAMAMIA_CONFIG"))

## Species-first registries have no predefined assembly: resolve the first two
## species, or the pair named by MAMAMIA_SMOKE_RECIPIENT/MAMAMIA_SMOKE_DONOR.
analyze_config <- config
assembly_key <- NULL
if (length(config$species) >= 2L) {
    recipient_key <- Sys.getenv(
        "MAMAMIA_SMOKE_RECIPIENT", unset = names(config$species)[[1L]]
    )
    donor_key <- Sys.getenv(
        "MAMAMIA_SMOKE_DONOR", unset = names(config$species)[[2L]]
    )
    resolved <- MaMaMIA::resolveSpeciesAssembly(
        config, donor = donor_key, recipient = recipient_key
    )
    analyze_config <- resolved$config
    assembly_key <- resolved$key
    cat("species pair:", assembly_key, "\n")
}

# The production BAM root is read-only. Stage a symlink in session temp so
# htslib can write the CSI index next to the link, exactly as the app does.
source_bam <- Sys.getenv(
    "MAMAMIA_SMOKE_BAM",
    unset = "/srv/mamamia/bams/Sq54/cs_t2t_se_cereale.bam"
)
if (!file.exists(source_bam)) {
    stop("Smoke-test BAM does not exist: ", source_bam, call. = FALSE)
}
staged_bam <- file.path(tempdir(), basename(source_bam))
if (!file.symlink(source_bam, staged_bam)) {
    stop("Could not stage ", source_bam, call. = FALSE)
}

result <- MaMaMIA::analyzeBamInput(
    staged_bam,
    analyze_config,
    assembly = assembly_key,
    windows_per_chromosome = 1000L,
    normalization = "none",
    max_fit_windows = 5000L,
    cores = 1L,
    cache = NULL,
    progress = function(stage, detail = NULL) {
        message("stage: ", stage, if (is.null(detail)) "" else paste0(" (", detail, ")"))
    }
)
print(result$samples)
cat("segments:", sum(result$samples$n_segments), "\n")
cat("putative:", sum(result$samples$n_putative_introgressions), "\n")
unlink(c(staged_bam, paste0(staged_bam, ".csi")))
stopifnot(inherits(result, "BamAnalysisResult"))
stopifnot(sum(result$samples$n_segments) > 0)
