#!/usr/bin/env Rscript
## Regenerate inst/extdata/example/example.rds for the pipeline.
##
## The bundled example is a pre-computed BamAnalysisResult replayed by the
## Shiny app. It uses the synthetic 20 kb demo index and an in-memory BAM with
## one donor addition, so the example shows a real calibrated call.

suppressPackageStartupMessages(library(MaMaMIA))
root <- normalizePath(file.path(getwd()), winslash = "/")
example_dir <- file.path(root, "inst", "extdata", "example")
config_path <- file.path(example_dir, "assemblies.yml")
config <- readMaMaMIAConfig(config_path)
pair <- resolveSpeciesAssembly(
    config, donor = "secale_like", recipient = "triticum_like"
)

event_sam <- function(sample, contig_length = 20000L, read_length = 50L,
                      step = 200L, event_from = 6000L, event_to = 10000L,
                      event_step = 20L) {
    sequence <- strrep("A", read_length)
    quality <- strrep("I", read_length)
    make_reads <- function(contig, positions) {
        sprintf(
            "%s_%s_%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\tRG:Z:%s",
            sample, contig, positions, contig, positions, read_length,
            sequence, quality, sample
        )
    }
    base <- seq(1L, contig_length - read_length + 1L, by = step)
    donor <- sort(unique(c(base, seq(event_from, event_to, by = event_step))))
    c(
        "@HD\tVN:1.6\tSO:coordinate",
        sprintf("@SQ\tSN:d1\tLN:%d", contig_length),
        sprintf("@SQ\tSN:r1\tLN:%d", contig_length),
        sprintf("@RG\tID:%s\tSM:%s", sample, sample),
        make_reads("d1", donor),
        make_reads("r1", base)
    )
}

directory <- tempfile("mamamia-example-")
dir.create(directory)
on.exit(unlink(directory, recursive = TRUE), add = TRUE)
sam <- file.path(directory, "example.sam")
bam <- file.path(directory, "example.bam")
writeLines(event_sam("example-sample"), sam)
status <- system2(
    "samtools", c("view", "--no-PG", "-b", "-o", bam, sam),
    stdout = FALSE, stderr = FALSE
)
stopifnot(identical(status, 0L))
ensureBamIndex(bam)

preflight <- preflightBams(bam, pair$config, assembly = pair$key)
result <- analyzeBamInput(
    bam, pair$config, assembly = pair$key,
    windows_per_chromosome = 100L,
    calibration = list(
        min_mb = 0.001,
        length_bins_mb = c(0.001, 0.01, Inf),
        floors = c(0.3, 0.3)
    ),
    cores = 1L
)

data <- list(
    results = list("100" = result),
    preflight = preflight,
    species_selection = list(
        recipient = "triticum_like",
        recipient_subgenome = "A",
        donor = "secale_like",
        donor_subgenome = "R"
    )
)
saveRDS(data, file.path(example_dir, "example.rds"), version = 2)
cat(
    "wrote", file.path(example_dir, "example.rds"),
    "with", sum(result$calls$pass), "calibrated call(s)\n"
)
