config <- MaMaMIA::readMaMaMIAConfig(Sys.getenv("MAMAMIA_CONFIG"))

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
    config,
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
