## Smoke test for the production image. It runs the full pipeline without
## starting Shiny, and renders both a cairo text probe and the real score
## profile.
##
## By default it uses the bundled demo BAM and registry, so the container needs
## no mounted BAM. Set MAMAMIA_SMOKE_BAM (and MAMAMIA_CONFIG) to exercise a
## production registry on the host instead.

font_probe <- tempfile(fileext = ".png")
grDevices::png(font_probe, width = 400, height = 300, type = "cairo")
graphics::plot(1:3, main = "fontconfig probe")
grDevices::dev.off()
if (!file.exists(font_probe) || !isTRUE(file.info(font_probe)$size > 0)) {
    stop(
        "Cairo could not render text; install fontconfig and fonts",
        call. = FALSE
    )
}
unlink(font_probe)

source_bam <- Sys.getenv("MAMAMIA_SMOKE_BAM", unset = "")
bundled <- !nzchar(source_bam)
if (nzchar(source_bam)) {
    config_path <- Sys.getenv("MAMAMIA_CONFIG", unset = "")
    if (!nzchar(config_path)) {
        stop(
            "MAMAMIA_CONFIG must be set when MAMAMIA_SMOKE_BAM is used",
            call. = FALSE
        )
    }
    config <- MaMaMIA::readMaMaMIAConfig(config_path)
    analyze_config <- config
    assembly_key <- NULL

    ## Species-first registries have no predefined assembly: resolve the first
    ## two species, or the pair named by MAMAMIA_SMOKE_RECIPIENT/_DONOR.
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
} else {
    demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
    if (!nzchar(demo) || !dir.exists(demo)) {
        stop(
            "The bundled bam-demo fixture is missing; ",
            "set MAMAMIA_SMOKE_BAM to use an external BAM",
            call. = FALSE
        )
    }
    config <- MaMaMIA::readMaMaMIAConfig(file.path(demo, "assemblies.yml"))
    resolved <- MaMaMIA::resolveSpeciesAssembly(
        config, donor = "secale_like", recipient = "triticum_like"
    )
    analyze_config <- resolved$config
    assembly_key <- resolved$key
    source_bam <- file.path(demo, "demo.bam")
    cat("bundled demo pair:", assembly_key, "\n")
}

if (!file.exists(source_bam)) {
    stop("Smoke-test BAM does not exist: ", source_bam, call. = FALSE)
}

## Stage a symlink in session temp so htslib can write the CSI index next to
## the link, exactly as the app does for uploaded BAMs.
staged_bam <- tempfile(pattern = "smoke-", fileext = ".bam")
if (!file.symlink(source_bam, staged_bam)) {
    stop("Could not stage ", source_bam, call. = FALSE)
}

## The bundled demo contigs are 20 kb, far below the 5 Mb production minimum
## call width, and its depth is uniform by construction (breadth saturates),
## so the demo runs at the resolution the README uses, on the count metric and
## with a toy-scale calibration. An external BAM keeps the production
## defaults and only has to prove that the pipeline completes.
if (bundled) {
    resolution <- 100L
    metric <- "count"
    calibration <- list(
        min_mb = 0.001,
        length_bins_mb = c(0.001, 0.01, Inf),
        floors = c(0.3, 0.3)
    )
} else {
    resolution <- 1000L
    metric <- NULL
    calibration <- NULL
}

result <- MaMaMIA::analyzeBamInput(
    staged_bam,
    analyze_config,
    assembly = assembly_key,
    windows_per_chromosome = resolution,
    metric = metric,
    calibration = calibration,
    cores = 1L,
    progress = function(stage, detail = list()) {
        message("stage: ", stage)
    }
)
print(result$samples)
cat("calibrated calls:", sum(result$calls$pass, na.rm = TRUE), "\n")

## The dashboard's renderPlot() prints a ggplot to a cairo device; exercise
## that path on the real result as well.
plot_file <- tempfile(fileext = ".png")
grDevices::png(plot_file, width = 1200, height = 400, type = "cairo")
print(
    ggplot2::ggplot(
        result$scores,
        ggplot2::aes(
            x = .data$start / 1e6, y = .data$score, colour = .data$mode
        )
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey70") +
        ggplot2::geom_line(linewidth = 0.3) +
        ggplot2::facet_wrap(~group, scales = "free_x", ncol = 1) +
        ggplot2::theme_minimal(base_size = 12)
)
grDevices::dev.off()
if (!file.exists(plot_file) || !isTRUE(file.info(plot_file)$size > 0)) {
    stop("Could not render the score profile", call. = FALSE)
}
cat("plot:", file.info(plot_file)$size, "bytes\n")
unlink(plot_file)

unlink(c(staged_bam, paste0(staged_bam, ".csi")))
stopifnot(inherits(result, "BamAnalysisResult"))
stopifnot(nrow(result$metrics) > 0L)
stopifnot(nrow(result$scores) > 0L)
stopifnot(nrow(result$calibration) > 0L)
if (bundled) {
    ## The demo fixture is built to contain a donor depth change, so a run
    ## that finds nothing means segmentation or calibration is broken.
    stopifnot(nrow(result$calls) > 0L)
}
