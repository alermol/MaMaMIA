#!/usr/bin/env Rscript
## Regenerate inst/extdata/example/example.rds and its GC index.
##
## The bundled example is a pre-computed BamAnalysisResult replayed by the
## Shiny app. Seven chromosome pairs of different lengths demonstrate every
## plotted event category:
##
##   1A 160 kb  40-43 kb   balanced substitution (segment)
##   2A 220 kb 110-115 kb  addition (segment)
##   3A 180 kb  90-94 kb   deletion (segment)
##   4A 260 kb 200-207 kb  complex (segment)
##   5A 200 kb  full       added donor chromosome
##   6A 240 kb  full       deleted recipient chromosome
##   7A 180 kb  40-80 kb   other (weak segment)
##   8A 200 kb  full       balanced substitution
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

chromosome_lengths <- c(
    160000L, 220000L, 180000L, 260000L, 200000L, 240000L, 180000L,
    200000L
)
read_length <- 50L
windows_per_chromosome <- 1000L
sequence <- strrep("A", read_length)
quality <- strrep("I", read_length)
chromosomes <- 1:8

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
    ),
    list(
        chromosome = 5L, from = 1L, to = 200000L,
        donor_extra = 20L, recipient = "base"
    ),
    list(
        chromosome = 6L, from = 1L, to = 240000L,
        donor_extra = NULL, recipient = "absent"
    ),
    list(
        chromosome = 7L, from = 40000L, to = 80000L,
        donor_extra = 220L, recipient = "base"
    ),
    list(
        chromosome = 8L, from = 1L, to = 200000L,
        donor_extra = 200L, recipient = "absent"
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
        floors = c(0.2, 0.2),
        fdr = 0.5
    ),
    mask = list(physical_fallback_mb = 0),
    cores = 1L
)

## The weak `other` segment has no recipient depletion, so assignment leaves
## it without a host. Pin it to 7A for the demo so the Other legend colour is
## visible; its class, copy state and needs_review status are unchanged. The
## segment addition on 2A is pinned the same way: a donor-only addition has no
## recipient deletion to anchor a host, but the demo should show a segment
## addition as well as the whole-chromosome addition on 5R.
pin_to_group <- function(row, chr_id, chr_name, group) {
    result$calls$recipient_component[[row]] <<- "A"
    result$calls$recipient_chr_id[[row]] <<- chr_id
    result$calls$recipient_chr_name[[row]] <<- chr_name
    result$calls$recipient_group[[row]] <<- as.character(group)
    result$calls$recipient_start_bp[[row]] <<- result$calls$start_bp[[row]]
    result$calls$recipient_end_bp[[row]] <<- result$calls$end_bp[[row]]
    result$calls$recipient_assignment[[row]] <<- "ambiguous"
}
other <- which(
    as.character(result$calls$class) == "other" &
        result$calls$pass %in% TRUE
)
if (length(other) == 1L) {
    pin_to_group(other, "r7", "7A", 7L)
}
segment_addition <- which(
    as.character(result$calls$class) == "addition" &
        as.character(result$calls$group) == "2" &
        result$calls$pass %in% TRUE
)
if (length(segment_addition) == 1L) {
    pin_to_group(segment_addition, "r2", "2A", 2L)
}
full_substitution <- which(
    as.character(result$calls$class) == "substitution" &
        as.character(result$calls$group) == "8" &
        result$calls$pass %in% TRUE
)
if (length(full_substitution) == 1L) {
    pin_to_group(full_substitution, "r8", "8A", 8L)
}

passing <- result$calls[result$calls$pass %in% TRUE, , drop = FALSE]
expected <- c("substitution", "addition", "deletion", "complex", "other")
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
