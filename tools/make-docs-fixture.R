#!/usr/bin/env Rscript

## Generate the small BAM/GC fixture used by the documentation examples.
##
## Usage:
##   Rscript tools/make-docs-fixture.R [--force]
##       [--out inst/extdata/bam-demo]
##       [--samtools /path/to/samtools]
##       [--gcrank-build /path/to/gcrank-build]
##
## Unlike tests/testthat/fixtures/tiny.gcrank, the two contigs have variable
## GC content, so the GC model can be fitted and the complete pipeline
## (BAM -> metrics -> scoring -> segmentation -> BEDPE) runs. The
## sequences, read positions and manifest timestamp are fixed.

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
    out_dir <- file.path("inst", "extdata", "bam-demo")
}
samtools <- value_after(args, "--samtools")
if (is.null(samtools) || !nzchar(samtools)) {
    samtools <- Sys.which("samtools")
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

if (!nzchar(samtools) || !file.exists(samtools)) {
    stop("samtools not found; pass --samtools", call. = FALSE)
}
if (!nzchar(gcrank_build) || !file.exists(gcrank_build)) {
    stop(
        "gcrank-build not found; pass --gcrank-build or build it from src/gcrank",
        call. = FALSE
    )
}

if (dir.exists(out_dir) && !force) {
    if (file.exists(file.path(out_dir, "demo.bam"))) {
        stop(
            "fixture already exists in ", out_dir,
            "; pass --force to regenerate it",
            call. = FALSE
        )
    }
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run <- function(command, args) {
    status <- system2(command, args)
    if (!identical(status, 0L)) {
        stop(
            "command failed (", status, "): ",
            paste(c(command, args), collapse = " "),
            call. = FALSE
        )
    }
}

contig_length <- 20000L
block_length <- 20L
contig_names <- c("d1", "r1")

## 20 bp blocks with 4..20 G/C bases: GC varies between 20% and 100%, which
## gives the smoother enough signal within every component.
variable_sequence <- function(offset) {
    blocks <- contig_length / block_length
    gc <- ((seq_len(blocks) + offset) %% 17L) + 4L
    paste0(
        vapply(gc, function(x) {
            paste0(strrep("G", x), strrep("A", block_length - x))
        }, character(1)),
        collapse = ""
    )
}

write_fasta <- function(path, sequences) {
    lines <- unlist(lapply(names(sequences), function(name) {
        sequence <- sequences[[name]]
        starts <- seq(1L, nchar(sequence), by = 60L)
        chunks <- substring(
            sequence, starts, pmin(starts + 59L, nchar(sequence))
        )
        c(paste0(">", name), chunks)
    }))
    writeLines(lines, path)
}

fasta_path <- file.path(out_dir, "demo.fa")
write_fasta(
    fasta_path,
    list(
        d1 = variable_sequence(0L),
        r1 = variable_sequence(8L)
    )
)

read_length <- 50L
sequence <- strrep("A", read_length)
quality <- strrep("I", read_length)

## The donor has twice the coverage in the first half of the contig, which
## yields one clear donor-excess segment for the BEDPE example.
donor_positions <- c(
    seq(1L, contig_length / 2L - read_length + 1L, by = 4L),
    seq(contig_length / 2L + 1L, contig_length - read_length + 1L, by = 8L)
)
recipient_positions <- seq(1L, contig_length - read_length + 1L, by = 8L)

make_read <- function(contig, position) {
    sprintf(
        "r_%s_%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\tRG:Z:demo",
        contig, position, contig, position, read_length, sequence, quality
    )
}

sam <- c(
    "@HD\tVN:1.6\tSO:coordinate",
    sprintf("@SQ\tSN:%s\tLN:%d", contig_names, contig_length),
    "@RG\tID:demo\tSM:demo-sample\tPL:ILLUMINA",
    make_read("d1", donor_positions),
    make_read("r1", recipient_positions)
)

sam_path <- file.path(out_dir, "demo.sam")
writeLines(sam, sam_path)

bam_path <- file.path(out_dir, "demo.bam")
index_path <- paste0(bam_path, ".csi")
run(samtools, c("view", "--no-PG", "-b", "-o", bam_path, sam_path))
run(samtools, c("index", "-c", "-o", index_path, bam_path))
unlink(sam_path)

index_dir <- file.path(out_dir, "demo.gcrank")
old_epoch <- Sys.getenv("SOURCE_DATE_EPOCH", unset = NA_character_)
Sys.setenv(SOURCE_DATE_EPOCH = "1700000000")
if (force && dir.exists(index_dir)) {
    unlink(index_dir, recursive = TRUE, force = TRUE)
}
build_args <- c(
    "--fasta", fasta_path,
    "--out", out_dir,
    "--name", "demo",
    "--threads", "1"
)
if (force) {
    build_args <- c(build_args, "--overwrite")
}
run(gcrank_build, build_args)
if (is.na(old_epoch)) {
    Sys.unsetenv("SOURCE_DATE_EPOCH")
} else {
    Sys.setenv(SOURCE_DATE_EPOCH = old_epoch)
}

config <- list(
    version = 2L,
    defaults = list(
        windows_per_chromosome = 100L,
        allowed_windows = c(100L, 1000L, 10000L),
        mapq = 20L,
        exclude_flags = 1796L,
        gc_mmap = TRUE,
        gc_validate_crc = FALSE
    ),
    server = list(
        max_bams = 4L,
        max_total_bytes = 20 * 1024^2
    ),
    species = list(
        triticum_like = list(
            label = "Triticum-like (A)",
            description = "Synthetic 20 kb contig with variable GC",
            default_component = "A",
            gc_indexes = list("demo.gcrank"),
            contigs = list(list(id = "r1", name = "1A", component = "A"))
        ),
        secale_like = list(
            label = "Secale-like (R)",
            description = "Synthetic 20 kb contig with variable GC",
            default_component = "R",
            gc_indexes = list("demo.gcrank"),
            contigs = list(list(id = "d1", name = "1R", component = "R"))
        )
    )
)
yaml_lines <- strsplit(yaml::as.yaml(config), "\n", fixed = TRUE)[[1L]]
yaml_lines <- sub(": yes$", ": true", yaml_lines)
yaml_lines <- sub(": no$", ": false", yaml_lines)
writeLines(yaml_lines, file.path(out_dir, "assemblies.yml"))

message("Wrote demo fixture to ", normalizePath(out_dir, winslash = "/"))
