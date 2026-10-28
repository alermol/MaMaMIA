#!/usr/bin/env Rscript
## Run the pipeline on the 24 Aegilops speltoides validation BAMs.
##
## Usage:
##   Rscript tools/validate-24.R [--bams Sq74,Sq75,Sq76] [--out FILE]
##
## The assembly is built from the study's 1000-window contig table, which maps
## every wheat and speltoides contig to its component and chromosome group.

suppressPackageStartupMessages(library(MaMaMIA))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
    position <- match(name, args)
    if (is.na(position)) {
        return(default)
    }
    if (position == length(args)) {
        stop("Missing value for ", name, call. = FALSE)
    }
    args[[position + 1L]]
}

root <- normalizePath(getwd(), winslash = "/")
study <- get_arg(
    "--study", file.path(root, "test_data", "introgression-study")
)
bam_dir <- get_arg(
    "--bam-dir", file.path(root, "test_data", "bam_files", "aegilops-speltoides")
)
selected <- get_arg("--bams", NULL)
out <- get_arg("--out", file.path(study, "results", "validation_24bam.rds"))
cores <- as.integer(get_arg("--cores", "30"))
score_type <- get_arg("--score", "paired_log_ratio")

contigs <- utils::read.delim(
    file.path(study, "data", "windows_n1000", "contigs.tsv"),
    stringsAsFactors = FALSE
)
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

species <- list(
    triticum_aestivum = list(
        label = "Triticum aestivum",
        gc_indexes = file.path(
            root, "test_data", "gcidx", "triticum-aestivum", "index.gcrank"
        ),
        contigs = lapply(
            seq_len(nrow(contigs[contigs$role == "recipient", ])),
            function(i) {
                x <- contigs[contigs$role == "recipient", ][i, ]
                list(
                    id = x$chr, name = x$component, component = x$component,
                    group = x$group, copy_number = 2
                )
            }
        )
    ),
    aegilops_speltoides = list(
        label = "Aegilops speltoides",
        gc_indexes = file.path(
            root, "test_data", "gcidx", "speltoides-study.gcrank"
        ),
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
        metric = "breadth",
        score = "paired_log_ratio",
        donor_aggregation = "mean",
        mask = list(
            physical_fallback_mb = 1,
            recurrence_threshold = 0.4
        ),
        calibration = list(
            length_bins_mb = c(5, 10, 25, 50, 100, Inf),
            fdr = 0.01,
            single_sample_floors = c(2.0, 1.5, 1.2, 1.0, 0.5)
        )
    ),
    server = list(max_bams = 30L, max_total_bytes = 20 * 1024^3),
    species = species
)
attr(config, "config_path") <- file.path(study, "assemblies-study.yml")
config <- validateMaMaMIAConfig(config)
config <- resolveSpeciesAssembly(
    config, donor = "aegilops_speltoides",
    recipient = "triticum_aestivum",
    recipient_components = c("A", "B", "D"),
    donor_components = "S",
    default_sample_composition = c(A = 2, B = 2, D = 2),
    key = "wheat_speltoides"
)$config

bams <- list.files(bam_dir, pattern = "[.]bam$", full.names = TRUE)
if (!is.null(selected)) {
    wanted <- trimws(strsplit(selected, ",", fixed = TRUE)[[1L]])
    bams <- bams[sub("[.]bam$", "", basename(bams)) %in% wanted]
}
if (length(bams) == 0L) {
    stop("No BAMs selected", call. = FALSE)
}
message(
    "Analyzing ", length(bams), " BAM(s) with ", cores, " core(s), score ",
    score_type
)

started <- Sys.time()
result <- analyzeBamInput(
    bams, config,
    assembly = "wheat_speltoides",
    windows_per_chromosome = 1000L,
    score = score_type,
    cores = cores
)
elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
message("Finished in ", round(elapsed, 1), " s")

calls <- result$calls[result$calls$pass, , drop = FALSE]
print(
    calls[order(calls$sample, calls$group, calls$start_iid),
        c("sample", "group", "start_iid", "end_iid", "width_bp",
          "score_mean", "donor_mean", "recipient_mean", "class",
          "confidence", "method")],
    row.names = FALSE
)

truth <- data.frame(
    sample = c("Sq74", "Sq75", "Sq75", "Sq75", "Sq75", "Sq76"),
    group = c("7", "7", "2", "6", "5", "7"),
    stringsAsFactors = FALSE
)
truth$detected <- vapply(seq_len(nrow(truth)), function(i) {
    any(
        calls$sample == truth$sample[[i]] &
            calls$group == truth$group[[i]]
    )
}, logical(1))
cat("\nReference events detected:", sum(truth$detected), "/", nrow(truth), "\n")
print(truth, row.names = FALSE)
cat(
    "Control calls >=20 Mb:",
    sum(calls$width_bp >= 20e6 & !calls$sample %in% c("Sq74", "Sq75", "Sq76")),
    "\n"
)
cat("Masked fraction:", round(mean(result$masks$mask), 3), "\n")
cat("Elapsed seconds:", round(elapsed, 1), "\n")

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(result, out)
message("Saved ", out)
