#!/usr/bin/env Rscript

## Benchmark for the high-resolution GC-fit subsampling (T7.1).
##
## The script builds a synthetic RCA with 28 target contigs (14 donor and 14
## recipient chromosomes) at two resolutions:
##
##   28,000 windows  = 1,000 windows per contig
##   280,000 windows = 10,000 windows per contig
##
## Each case runs in a fresh R process so `/proc/self/status` VmHWM is a clean
## peak-RSS measurement. Results are appended to
## `inst/benchmarks/fit-subsampling-results.csv`.
##
## Usage:
##   Rscript inst/benchmarks/benchmark-fit-subsampling.R
##   Rscript inst/benchmarks/benchmark-fit-subsampling.R --windows=28000 \
##       --cores=8 --max-fit-windows=50000 --out=bench-out
##
## The no-argument form runs both cases in child processes and writes the CSV
## next to this script.

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(name, default = NULL) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit) == 0L) {
        return(default)
    }
    sub(paste0("^--", name, "="), "", hit[[1L]])
}

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script <- normalizePath(script[[1L]], winslash = "/", mustWork = TRUE)
out_dir <- arg_value("out", dirname(script))
windows <- as.integer(arg_value("windows", NA_character_))
cores <- as.integer(arg_value("cores", "4"))
max_fit_windows <- as.numeric(arg_value("max-fit-windows", "50000"))

if (is.na(windows)) {
    ## Parent mode: run each case in a fresh process and append to one CSV.
    out <- file.path(out_dir, "fit-subsampling-results.csv")
    unlink(out)
    for (n in c(28000L, 280000L)) {
        status <- system2(
            file.path(R.home("bin"), "Rscript"),
            c(
                script,
                paste0("--windows=", n),
                paste0("--cores=", cores),
                paste0("--max-fit-windows=", max_fit_windows),
                paste0("--out=", out_dir)
            )
        )
        if (!identical(status, 0L)) {
            stop("benchmark for ", n, " windows failed", call. = FALSE)
        }
    }
    cat("\nresults written to ", out, "\n", sep = "")
    quit(status = 0L)
}

suppressPackageStartupMessages(library(MaMaMIA))

read_vm <- function(field) {
    path <- "/proc/self/status"
    if (!file.exists(path)) {
        return(NA_real_)
    }
    hit <- grep(paste0("^", field, ":"), readLines(path, warn = FALSE),
        value = TRUE
    )
    if (length(hit) == 0L) {
        return(NA_real_)
    }
    as.numeric(sub(" kB", "", sub("^.*:\\s+", "", hit[[1L]])))
}

heap_mb <- function(g) {
    g["Ncells", "max used"] * 56 / 1e6 + g["Vcells", "max used"] * 8 / 1e6
}

## 28 target contigs: 14 donor/recipient pairs with equal window counts.
make_windows <- function(chr_id, n, value) {
    start <- seq(0, by = 1000, length.out = n)
    data.frame(
        chr_id = chr_id,
        start = start,
        stop = start + 1000,
        value = value,
        stringsAsFactors = FALSE
    )
}

make_benchmark_rca <- function(total_windows) {
    n_contigs <- 28L
    n_per <- as.integer(total_windows / n_contigs)
    if (n_per * n_contigs != total_windows) {
        stop("`windows` must be a multiple of 28", call. = FALSE)
    }

    gc <- seq(0.2, 0.8, length.out = n_per)
    gc <- pmin(
        0.95,
        pmax(0.05, gc + 0.015 * sin(seq_len(n_per) * 2 * pi / 37))
    )
    coverage <- function(offset) {
        pmax(
            1,
            round(
                30 + 18 * gc +
                    8 * sin(seq_len(n_per) * 2 * pi / 53 + offset) +
                    3 * cos(seq_len(n_per) * 2 * pi / 11 + offset)
            )
        )
    }

    donor_ids <- sprintf("d%02d", 1:14)
    rec_ids <- sprintf("r%02d", 1:14)
    donor_cov <- do.call(rbind, lapply(donor_ids, make_windows,
        n = n_per, value = coverage(0)
    ))
    donor_gc <- do.call(rbind, lapply(donor_ids, make_windows,
        n = n_per, value = gc
    ))
    rec_cov <- do.call(rbind, lapply(rec_ids, make_windows,
        n = n_per, value = coverage(pi / 3)
    ))
    rec_gc <- do.call(rbind, lapply(rec_ids, make_windows,
        n = n_per, value = gc
    ))
    meta <- data.frame(
        chr_id = c(donor_ids, rec_ids),
        chr_name = c(
            sprintf("%dD", 1:14),
            sprintf("%dR", 1:14)
        ),
        subgenome = c(rep("D", 14L), rep("R", 14L)),
        stringsAsFactors = FALSE
    )
    RCA(donor_cov, donor_gc, rec_cov, rec_gc, meta)
}

set.seed(1L)
rca <- make_benchmark_rca(windows)

before <- read_vm("VmHWM")
gc(reset = TRUE)
timing <- system.time(
    corrected <- suppressWarnings(correctReadCounts(
        rca,
        cores = cores,
        verbose = FALSE,
        cache = NULL,
        max_fit_windows = max_fit_windows
    ))
)
g <- gc()
after <- read_vm("VmHWM")

row <- data.frame(
    n_contigs = 28L,
    windows_total = nrow(rca$data),
    windows_per_contig = nrow(rca$data) / 28L,
    max_fit_windows = max_fit_windows,
    n_ideal = corrected$fit_sampling$n_ideal,
    n_fitted = corrected$fit_sampling$n_fitted,
    sampled = corrected$fit_sampling$sampled,
    cores = cores,
    elapsed_s = round(timing[["elapsed"]], 3),
    peak_rss_mb = round(after / 1024, 2),
    peak_rss_delta_mb = round((after - before) / 1024, 2),
    peak_heap_mb = round(heap_mb(g), 2),
    finite_cor_gc = all(is.finite(corrected$data$cor.gc)),
    stringsAsFactors = FALSE
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- file.path(out_dir, "fit-subsampling-results.csv")
utils::write.table(
    row,
    out,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    append = file.exists(out),
    col.names = !file.exists(out)
)
print(row, row.names = FALSE)
cat("\nresults appended to ", out, "\n", sep = "")
