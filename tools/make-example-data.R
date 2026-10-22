#!/usr/bin/env Rscript

## Generate the pre-computed data behind the "bundled example" in the Shiny
## app.
##
## Usage:
##   Rscript tools/make-example-data.R [--force]
##       [--out inst/extdata/example]
##
## The example is generated from the small documentation fixture in
## inst/extdata/bam-demo. Only the resulting tables are shipped in
## inst/extdata/example, so production containers do not need an example BAM
## (or a server-side BAM directory) mounted. The app replays the stored
## progress details, which keeps the user-facing stage messages identical to a
## real BAM run.
##
## The demo contigs are 20 kb, so the 10,000-windows-per-chromosome resolution
## has only 2 bp windows and too few distinct GC values for the ZINB fit. Those
## resolutions are skipped; the generated example supports every resolution
## that succeeds.

args <- commandArgs(trailingOnly = TRUE)

value_after <- function(args, flag) {
    i <- match(flag, args)
    if (is.na(i)) {
        return(NULL)
    }
    if (i == length(args)) {
        stop(flag, " requires a value", call. = FALSE)
    }
    args[[i + 1L]]
}

force <- "--force" %in% args
out_dir <- value_after(args, "--out")
if (is.null(out_dir)) {
    out_dir <- file.path("inst", "extdata", "example")
}
demo_dir <- file.path("inst", "extdata", "bam-demo")
if (!dir.exists(demo_dir)) {
    demo_dir <- system.file("extdata", "bam-demo", package = "MaMaMIA")
}
if (!nzchar(demo_dir) || !dir.exists(demo_dir)) {
    stop("bam-demo fixture not found; run from the package root", call. = FALSE)
}

if (!requireNamespace("MaMaMIA", quietly = TRUE)) {
    stop("MaMaMIA must be installed to generate the example data", call. = FALSE)
}

if (dir.exists(out_dir) && !force) {
    if (file.exists(file.path(out_dir, "example.rds"))) {
        stop(
            "example data already exists in ", out_dir,
            "; pass --force to regenerate it",
            call. = FALSE
        )
    }
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

demo_bam <- file.path(demo_dir, "demo.bam")
demo_config_path <- file.path(demo_dir, "assemblies.yml")
if (!file.exists(demo_bam) || !file.exists(demo_config_path)) {
    stop("bam-demo fixture is incomplete", call. = FALSE)
}

config <- MaMaMIA::readMaMaMIAConfig(demo_config_path)
pair <- MaMaMIA::resolveSpeciesAssembly(
    config, donor = "secale_like", recipient = "triticum_like"
)
allowed <- sort(unique(as.integer(config$defaults$allowed_windows)))

message("Generating example results for: ", paste(allowed, collapse = ", "))
results <- list()
for (resolution in allowed) {
    result <- tryCatch(
        suppressWarnings(MaMaMIA::analyzeBamInput(
            demo_bam,
            pair$config,
            assembly = pair$key,
            windows_per_chromosome = resolution,
            normalization = "none",
            cores = 1L,
            max_fit_windows = Inf,
            cache = FALSE
        )),
        error = function(e) {
            message(
                "  resolution ", resolution, " skipped: ", conditionMessage(e)
            )
            NULL
        }
    )
    if (is.null(result)) {
        next
    }
    results[[as.character(resolution)]] <- result
    message(
        "  resolution ", resolution, ": ",
        nrow(result$windows$table), " windows, ",
        result$samples$n_segments[[1L]], " segments"
    )
}
if (length(results) == 0L) {
    stop("No example resolution could be generated", call. = FALSE)
}

sanitize_result <- function(result) {
    result$config <- NULL
    result$bams$path <- "demo.bam"
    result$bams$index_path <- "demo.bam.csi"
    result$samples$path <- "demo.bam"
    result$samples$index_path <- "demo.bam.csi"
    result$assembly$gc_indexes <- "demo.gcrank"
    result$provenance$config_path <- NA_character_
    result$provenance$config_hash <- NA_character_
    result$provenance$bams$path <- "demo.bam"
    result$provenance$bams$mtime <- NA_character_
    result$provenance$cache <- NA_character_
    result$provenance$uvr_lock <- NA_character_
    result$provenance$uvr_lock_hash <- NA_character_
    result$provenance$started_at <- NA_character_
    result$provenance$duration_secs <- NA_real_
    result
}
results <- lapply(results, sanitize_result)

preflight <- MaMaMIA::preflightBams(
    demo_bam, pair$config, assembly = pair$key
)
preflight$path <- "demo.bam"
preflight$index_path <- "demo.bam.csi"
preflight$source <- "demo.bam"
preflight$warnings <- ""

## Rebuild the registry YAML from the demo text, dropping the server-side BAM
## root and limiting the offered resolutions to the ones that were generated.
yaml_lines <- readLines(demo_config_path, warn = FALSE)
yaml_lines <- yaml_lines[!grepl("^\\s*bam_root:", yaml_lines)]
start <- grep("^\\s*allowed_windows:\\s*$", yaml_lines)
if (length(start) == 1L) {
    end <- start
    while (end < length(yaml_lines) &&
        grepl("^\\s+-\\s*[0-9]+\\s*$", yaml_lines[[end + 1L]])) {
        end <- end + 1L
    }
    yaml_lines <- c(
        yaml_lines[seq_len(start)],
        sprintf("  - %d", sort(unique(as.integer(names(results))))),
        yaml_lines[seq.int(end + 1L, length(yaml_lines))]
    )
}
writeLines(yaml_lines, file.path(out_dir, "assemblies.yml"))

gcrank_source <- file.path(demo_dir, "demo.gcrank")
gcrank_target <- file.path(out_dir, "demo.gcrank")
if (dir.exists(gcrank_target)) {
    unlink(gcrank_target, recursive = TRUE, force = TRUE)
}
if (!file.copy(gcrank_source, out_dir, recursive = TRUE)) {
    stop("Could not copy demo.gcrank", call. = FALSE)
}

example <- list(
    version = 1L,
    label = "Bundled synthetic example",
    description = paste(
        "Two synthetic 20 kb contigs (Triticum-like recipient and",
        "Secale-like donor) with a donor excess in the first half."
    ),
    assembly = pair$key,
    species_selection = list(
        recipient = "triticum_like",
        recipient_subgenome = "A",
        donor = "secale_like",
        donor_subgenome = "R"
    ),
    preflight = preflight,
    results = results
)
saveRDS(
    example,
    file.path(out_dir, "example.rds"),
    version = 2,
    compress = "xz"
)

message(
    "Wrote example data to ",
    normalizePath(out_dir, winslash = "/"),
    " (", format(file.info(file.path(out_dir, "example.rds"))$size),
    " bytes)"
)
