#!/usr/bin/env Rscript

## Generate the tiny committed BAM/GC fixtures under tests/testthat/fixtures.
##
## Usage:
##   Rscript tools/make-test-fixtures.R [--force]
##       [--out tests/testthat/fixtures]
##       [--samtools /path/to/samtools]
##       [--gcrank-build /path/to/gcrank-build]
##
## The script is deliberately deterministic: fixed sequences, fixed read
## positions and SOURCE_DATE_EPOCH for the gcrank manifest.

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
    out_dir <- file.path("tests", "testthat", "fixtures")
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
    stop(
        "samtools not found; pass --samtools or activate the conda mamamia environment",
        call. = FALSE
    )
}
if (!nzchar(gcrank_build) || !file.exists(gcrank_build)) {
    stop(
        "gcrank-build not found; pass --gcrank-build or build it from src/gcrank",
        call. = FALSE
    )
}

if (dir.exists(out_dir) && !force) {
    existing <- c("tiny.fa", "tiny.bam", "tiny.bam.csi", "assemblies.yml")
    if (any(file.exists(file.path(out_dir, existing)))) {
        stop(
            "fixtures already exist in ", out_dir,
            "; pass --force to regenerate them",
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

repeat_motif <- function(motif, n) {
    reps <- ceiling(n / nchar(motif))
    substr(paste(rep(motif, reps), collapse = ""), 1L, n)
}

write_fasta <- function(path, sequences) {
    lines <- unlist(lapply(names(sequences), function(name) {
        sequence <- sequences[[name]]
        starts <- seq(1L, nchar(sequence), by = 60L)
        chunks <- substring(
            sequence,
            starts,
            pmin(starts + 59L, nchar(sequence))
        )
        c(paste0(">", name), chunks)
    }))
    writeLines(lines, path)
}

contig_length <- 20000L
d1 <- repeat_motif("GGCCTTTTTT", contig_length) # 40% GC
r1 <- repeat_motif("GGGCCCTTTT", contig_length) # 60% GC
write_fasta(file.path(out_dir, "tiny.fa"), list(d1 = d1, r1 = r1))

read_length <- 50L
sequence <- paste(rep("A", read_length), collapse = "")
quality <- paste(rep("I", read_length), collapse = "")

make_read <- function(name, flag, contig, position, mapq) {
    paste(
        name, flag, contig, position, mapq, paste0(read_length, "M"),
        "*", 0, 0, sequence, quality, "RG:Z:tiny",
        sep = "\t"
    )
}

sam <- c(
    "@HD\tVN:1.6\tSO:coordinate",
    paste0("@SQ\tSN:d1\tLN:", contig_length),
    paste0("@SQ\tSN:r1\tLN:", contig_length),
    "@RG\tID:tiny\tSM:tiny-sample\tPL:ILLUMINA",
    "@PG\tID:make-test-fixtures\tPN:make-test-fixtures\tVN:0.0.1",
    make_read("d1_normal_1", 0L, "d1", 101L, 60L),
    make_read("d1_normal_2", 0L, "d1", 201L, 60L),
    make_read("d1_normal_3", 0L, "d1", 301L, 60L),
    make_read("d1_normal_4", 0L, "d1", 401L, 60L),
    make_read("d1_normal_5", 0L, "d1", 501L, 60L),
    make_read("d1_normal_6", 0L, "d1", 601L, 60L),
    make_read("d1_low_mapq", 0L, "d1", 701L, 10L),
    make_read("d1_secondary", 256L, "d1", 801L, 60L),
    make_read("d1_supplementary", 2048L, "d1", 901L, 60L),
    make_read("r1_normal_1", 0L, "r1", 1001L, 60L),
    make_read("r1_normal_2", 0L, "r1", 1101L, 60L),
    make_read("r1_normal_3", 0L, "r1", 1201L, 60L),
    make_read("r1_normal_4", 0L, "r1", 1301L, 60L),
    make_read("r1_duplicate", 1024L, "r1", 1401L, 60L)
)

sam_path <- file.path(out_dir, "tiny.sam")
writeLines(sam, sam_path)

bam_path <- file.path(out_dir, "tiny.bam")
run(samtools, c("view", "--no-PG", "-b", "-o", bam_path, sam_path))
run(samtools, c("index", "-c", bam_path))
unlink(sam_path)

index_dir <- file.path(out_dir, "tiny.gcrank")
old_epoch <- Sys.getenv("SOURCE_DATE_EPOCH", unset = NA_character_)
Sys.setenv(SOURCE_DATE_EPOCH = "1700000000")
if (force && dir.exists(index_dir)) {
    unlink(index_dir, recursive = TRUE, force = TRUE)
}
build_args <- c(
    "--fasta", file.path(out_dir, "tiny.fa"),
    "--out", out_dir,
    "--name", "tiny",
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
    server = list(max_bams = 4L, max_total_bytes = 1024^2),
    species = list(
        tiny_recipient = list(
            label = "Tiny synthetic recipient",
            description = "Two 20 kb contigs with 40% and 60% GC",
            default_component = "R",
            gc_indexes = list("tiny.gcrank"),
            contigs = list(list(id = "r1", name = "1R", component = "R"))
        ),
        tiny_donor = list(
            label = "Tiny synthetic donor",
            description = "Two 20 kb contigs with 40% and 60% GC",
            default_component = "D",
            gc_indexes = list("tiny.gcrank"),
            contigs = list(list(id = "d1", name = "1D", component = "D"))
        )
    )
)
yaml_lines <- strsplit(yaml::as.yaml(config), "\n", fixed = TRUE)[[1L]]
yaml_lines <- sub(": yes$", ": true", yaml_lines)
yaml_lines <- sub(": no$", ": false", yaml_lines)
writeLines(yaml_lines, file.path(out_dir, "assemblies.yml"))

message("Wrote fixtures to ", normalizePath(out_dir, winslash = "/"))
