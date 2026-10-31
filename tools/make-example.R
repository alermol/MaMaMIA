#!/usr/bin/env Rscript
## Regenerate inst/extdata/example/example.rds and its GC index.
##
## The bundled example is a pre-computed BamAnalysisResult replayed by the
## Shiny app. Four chromosome pairs of different lengths each carry one event
## of a different width:
##
##   1A 160 kb  40-43 kb   balanced substitution
##   2A 220 kb 110-115 kb  addition
##   3A 180 kb  90-94 kb   deletion
##   4A 260 kb 200-207 kb  complex
##
## The low-signal gaps segment as non-passing `other` rows, so the Calls table
## also exposes neutral calls and the `needs_review` status. The four events
## pass calibration and cover all four introgression classes.

args <- commandArgs(trailingOnly = TRUE)
value_after <- function(args, flag) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) {
        return(NULL)
    }
    args[[i + 1L]]
}

gcrank_build <- value_after(args, "--gcrank-build")
if (is.null(gcrank_build) || !nzchar(gcrank_build)) {
    gcrank_build <- Sys.which("gcrank-build")
    if (!nzchar(gcrank_build)) {
        local_builder <- file.path("test_data", "gcidx", "gcrank-build")
        if (file.exists(local_builder)) {
            gcrank_build <- local_builder
        }
    }
}
if (!nzchar(gcrank_build) || !file.exists(gcrank_build)) {
    stop(
        "gcrank-build not found; pass --gcrank-build or build it from src/gcrank",
        call. = FALSE
    )
}

suppressPackageStartupMessages(library(MaMaMIA))
root <- normalizePath(file.path(getwd()), winslash = "/")
example_dir <- file.path(root, "inst", "extdata", "example")

chromosome_lengths <- c(160000L, 220000L, 180000L, 260000L)
read_length <- 50L
windows_per_chromosome <- 1000L
sequence <- strrep("A", read_length)
quality <- strrep("I", read_length)
chromosomes <- 1:4

## A 20 bp block pattern gives the GC spline enough variation in every
## component. Homologs use the same landscape so the correction is not driven
## by a donor/recipient GC difference.
variable_sequence <- function(length_bp, offset = 0L) {
    blocks <- length_bp / 20L
    gc <- ((seq_len(blocks) + offset) %% 17L) + 4L
    paste0(
        vapply(gc, function(x) {
            paste0(strrep("G", x), strrep("A", 20L - x))
        }, character(1)),
        collapse = ""
    )
}

directory <- tempfile("mamamia-example-")
dir.create(directory)
on.exit(unlink(directory, recursive = TRUE), add = TRUE)
contigs <- c(paste0("d", chromosomes), paste0("r", chromosomes))
contig_lengths <- rep(chromosome_lengths, 2L)
fasta <- file.path(directory, "example.fa")
writeLines(
    unlist(lapply(seq_along(contigs), function(i) {
        c(paste0(">", contigs[[i]]), variable_sequence(contig_lengths[[i]]))
    })),
    fasta
)

## Rebuild the bundled GC index for the four homolog pairs.
old_epoch <- Sys.getenv("SOURCE_DATE_EPOCH", unset = NA_character_)
Sys.setenv(SOURCE_DATE_EPOCH = "1700000000")
status <- system2(
    gcrank_build,
    c(
        "--fasta", fasta, "--out", example_dir, "--name", "demo",
        "--threads", "1", "--overwrite"
    ),
    stdout = FALSE, stderr = FALSE
)
if (is.na(old_epoch)) {
    Sys.unsetenv("SOURCE_DATE_EPOCH")
} else {
    Sys.setenv(SOURCE_DATE_EPOCH = old_epoch)
}
stopifnot(identical(status, 0L))

base <- lapply(chromosome_lengths, function(length_bp) {
    seq(1L, length_bp - read_length + 1L, by = 200L)
})
in_event <- function(positions, from, to) {
    positions >= from & positions <= to
}

## The donor keeps a base layer everywhere. Adding a denser layer raises the
## donor side; removing the recipient layer removes the recipient side.
events <- list(
    list(
        chromosome = 1L, from = 40000L, to = 43000L,
        donor_extra = 200L, recipient = "absent"
    ),
    list(
        chromosome = 2L, from = 110000L, to = 115000L,
        donor_extra = 20L, recipient = "base"
    ),
    list(
        chromosome = 3L, from = 90000L, to = 94000L,
        donor_extra = NULL, recipient = "absent"
    ),
    list(
        chromosome = 4L, from = 200000L, to = 207000L,
        donor_extra = 20L, recipient = "absent"
    )
)

make_reads <- function(contig, positions) {
    sprintf(
        "%s_%s_%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\tRG:Z:%s",
        "example-sample", contig, positions, contig, positions, read_length,
        sequence, quality, "example-sample"
    )
}

sam_lines <- c(
    "@HD\tVN:1.6\tSO:coordinate",
    sprintf("@SQ\tSN:%s\tLN:%d", contigs, contig_lengths),
    sprintf("@RG\tID:%s\tSM:%s", "example-sample", "example-sample")
)
for (event in events) {
    index <- event$chromosome
    donor <- base[[index]]
    recipient <- base[[index]]
    if (!is.null(event$donor_extra)) {
        donor <- sort(unique(c(
            donor,
            seq(
                event$from + if (event$donor_extra == 200L) 100L else 0L,
                event$to, by = event$donor_extra
            )
        )))
    }
    if (identical(event$recipient, "absent")) {
        recipient <- recipient[!in_event(recipient, event$from, event$to)]
    }
    sam_lines <- c(
        sam_lines,
        make_reads(paste0("d", index), donor),
        make_reads(paste0("r", index), recipient)
    )
}

sam <- file.path(directory, "example.sam")
bam <- file.path(directory, "example.bam")
writeLines(sam_lines, sam)
status <- system2(
    "samtools", c("view", "--no-PG", "-b", "-o", bam, sam),
    stdout = FALSE, stderr = FALSE
)
stopifnot(identical(status, 0L))
ensureBamIndex(bam)

config_path <- file.path(example_dir, "assemblies.yml")
config <- readMaMaMIAConfig(config_path)
pair <- resolveSpeciesAssembly(
    config, donor = "secale_like", recipient = "triticum_like"
)

preflight <- preflightBams(bam, pair$config, assembly = pair$key)
result <- analyzeBamInput(
    bam, pair$config, assembly = pair$key,
    windows_per_chromosome = windows_per_chromosome,
    calibration = list(
        min_mb = 0.001,
        length_bins_mb = c(0.001, 0.01, Inf),
        floors = c(0.3, 0.3),
        fdr = 0.2
    ),
    cores = 1L
)

passing <- result$calls[result$calls$pass %in% TRUE, , drop = FALSE]
expected <- c("substitution", "addition", "deletion", "complex")
missing <- setdiff(expected, unique(as.character(passing$class)))
if (length(missing) > 0L) {
    stop(
        "the example does not produce all introgression classes: ",
        paste(missing, collapse = ", "),
        call. = FALSE
    )
}
event_groups <- vapply(events, function(event) {
    as.character(event$chromosome)
}, character(1))
if (!all(event_groups %in% unique(as.character(passing$group)))) {
    stop(
        "the example events are not distributed across chromosomes",
        call. = FALSE
    )
}
data <- list(
    results = list("1000" = result),
    preflight = preflight,
    species_selection = list(
        recipient = "triticum_like",
        recipient_component = "A",
        donor = "secale_like",
        donor_component = "R"
    )
)
saveRDS(data, file.path(example_dir, "example.rds"), version = 2)
cat(
    "wrote", file.path(example_dir, "example.rds"),
    "with", nrow(passing), "calibrated call(s):",
    paste(sort(unique(as.character(passing$class))), collapse = ", "), "\n"
)
