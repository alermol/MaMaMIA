#!/usr/bin/env Rscript

## Golden benchmark and regression oracle for the BAM-first pipeline (T3.4).
##
## The script compares the BAM-first preprocessing with the legacy BED tables
## in test_data:
##
##   1. read counts on the legacy intervals (bedtools -counts vs start-based
##      bamsignals) and an explanation of every mismatch with samtools;
##   2. GC content on the legacy intervals at the six decimals stored in the
##      legacy BED files;
##   3. per-contig read totals on the generated exact-N windows;
##   4. optionally, the full correction/segmentation pipeline against the old
##      directory-layout pipeline.
##
## Usage:
##   Rscript tools/benchmark-bam-pipeline.R [--test-data=DIR] [--cores=4]
##       [--out=bench-out-bam] [--sample=Sq54] [--windows=1000]
##       [--resolution=1000] [--bam=PATH] [--skip-full]
##
## --resolution is an alias for --windows; --bam overrides the sample BAM path
## while the legacy oracle still comes from --sample.
##
## The script never writes into test_data. All derived files go to --out.

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
    hit <- grep(paste0("^--", name, "="), args, value = TRUE)
    if (length(hit) == 0L) {
        return(default)
    }
    sub(paste0("^--", name, "="), "", hit[[1L]])
}

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

opts <- list(
    test_data = arg_value("test-data", Sys.getenv("MAMAMIA_TEST_DATA_DIR", "test_data")),
    cores = as.integer(arg_value("cores", "4")),
    out = arg_value("out", "bench-out-bam"),
    sample = arg_value("sample", "Sq54"),
    windows = as.integer(arg_value(
        "windows", arg_value("resolution", "1000")
    )),
    bam = arg_value("bam", NULL),
    skip_full = "--skip-full" %in% args
)
if (is.na(opts$cores) || opts$cores < 1L) {
    stop("--cores must be a positive integer", call. = FALSE)
}
if (is.na(opts$windows) || opts$windows < 1L) {
    stop("--windows must be a positive integer", call. = FALSE)
}

suppressPackageStartupMessages(library(MaMaMIA))

abort <- function(...) stop(..., call. = FALSE)
check <- function(env, label, ok, detail = "") {
    env$checks[[length(env$checks) + 1L]] <- data.frame(
        check = label,
        result = if (isTRUE(ok)) "PASS" else "FAIL",
        detail = as.character(detail),
        stringsAsFactors = FALSE
    )
    invisible(ok)
}
record <- function(env, step, expr) {
    timing <- system.time(value <- force(expr))
    env$timings[[length(env$timings) + 1L]] <- data.frame(
        step = step,
        elapsed_s = round(timing[["elapsed"]], 3),
        stringsAsFactors = FALSE
    )
    value
}
hash_object <- function(x) {
    path <- tempfile()
    on.exit(unlink(path), add = TRUE)
    saveRDS(x, path, version = 2)
    unname(tools::md5sum(path))
}

## Reference-span length of a CIGAR string (M/D/N/=/X consume the reference).
cigar_ref_length <- function(cigar) {
    ops <- regmatches(cigar, gregexpr("[0-9]+[MIDNSHP=X]", cigar))
    vapply(ops, function(x) {
        if (length(x) == 0L) {
            return(0L)
        }
        sizes <- as.integer(sub("[MIDNSHP=X]$", "", x))
        kinds <- sub("^[0-9]+", "", x)
        sum(sizes[kinds %in% c("M", "D", "N", "=", "X")])
    }, integer(1))
}

root <- normalizePath(opts$test_data, winslash = "/", mustWork = TRUE)
analysis_dir <- file.path(root, "mamamia_analysis", "arsenal")
assemblies_dir <- file.path(root, "assemblies")
bam_path <- if (!is.null(opts$bam)) {
    normalizePath(opts$bam, winslash = "/", mustWork = TRUE)
} else {
    file.path(root, "bam_files", opts$sample, "cs_t2t_se_cereale.bam")
}
legacy_dir <- file.path(analysis_dir, opts$sample)
legacy_prefix <- "cs_t2t_se_cereale"
required <- c(
    bam_path,
    file.path(legacy_dir, paste0(legacy_prefix, ".w", opts$windows, ".rcov.bed")),
    file.path(legacy_dir, paste0(legacy_prefix, ".w", opts$windows, ".dcov.bed")),
    file.path(assemblies_dir, "cs_t2t_se_cereale", "chromosomes.tsv")
)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
    abort("missing test data:\n", paste(missing, collapse = "\n"))
}

dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(opts$out, winslash = "/", mustWork = TRUE)
env <- new.env(parent = emptyenv())
env$checks <- env$timings <- list()
samtools <- Sys.which("samtools")

## The real assembly registry, assembled from the test data (not from YAML in
## test_data, which is read-only).
real_assembly_config <- function(root) {
    index_paths <- c(
        normalizePath(
            file.path(root, "gcidx", "triticum-aestivum", "index.gcrank"),
            winslash = "/", mustWork = TRUE
        ),
        normalizePath(
            file.path(root, "gcidx", "secale-cereale", "index.gcrank"),
            winslash = "/", mustWork = TRUE
        )
    )
    chromosomes <- utils::read.delim(
        file.path(root, "assemblies", "cs_t2t_se_cereale", "chromosomes.tsv"),
        header = FALSE,
        col.names = c("id", "name", "subgenome"),
        stringsAsFactors = FALSE
    )
    entries <- function(x) {
        lapply(seq_len(nrow(x)), function(i) {
            list(id = x$id[[i]], name = x$name[[i]], subgenome = x$subgenome[[i]])
        })
    }
    list(
        version = 1L,
        defaults = list(
            windows_per_chromosome = as.integer(opts$windows),
            allowed_windows = sort(unique(c(100L, opts$windows, 1000L, 10000L))),
            mapq = 20L,
            exclude_flags = 2308L,
            max_fit_windows = 50000L,
            gc_mmap = TRUE,
            gc_validate_crc = FALSE,
            normalization_by_side = TRUE
        ),
        server = list(max_bams = 20L, max_total_bytes = 20 * 1024^3),
        assemblies = list(cs_t2t_se_cereale = list(
            label = "Triticum aestivum CS T2T + Secale cereale",
            description = "Golden regression fixture",
            gc_indexes = as.list(index_paths),
            contigs = list(
                recipient = entries(chromosomes[chromosomes$subgenome != "R", ]),
                donor = entries(chromosomes[chromosomes$subgenome == "R", ])
            ),
            pairs = list(
                OZ284382.1 = "CP155614.1",
                OZ284383.1 = "CP155611.1",
                OZ284384.1 = "CP155608.1",
                OZ284385.1 = "CP155605.1",
                OZ284386.1 = "CP155602.1",
                OZ284387.1 = "CP155599.1",
                OZ284388.1 = "CP155597.1"
            )
        ))
    )
}

config_path <- file.path(out_dir, "assemblies-real.yml")
writeLines(yaml::as.yaml(real_assembly_config(root)), config_path)
config <- readMaMaMIAConfig(config_path)
assembly <- config$assemblies$cs_t2t_se_cereale
assembly_key <- "cs_t2t_se_cereale"

cat("MaMaMIA BAM-first golden benchmark\n")
cat("test_data : ", root, "\n", sep = "")
cat("sample    : ", opts$sample, "\n", sep = "")
cat("windows   : ", opts$windows, "\n", sep = "")
cat("cores     : ", opts$cores, "\n", sep = "")
cat("samtools  : ", if (nzchar(samtools)) samtools else "<not found>", "\n\n", sep = "")

## Stage the BAM once; the source directory is never modified.
session <- tempfile("mamamia-golden-")
dir.create(session)
staged_bam <- MaMaMIA:::stageBam(bam_path, session)
index_path <- ensureBamIndex(staged_bam)
check(env, "CSI index created for the staged BAM", file.exists(index_path), index_path)

legacy <- readMaMaMIAInput(
    analysis_dir, assemblies_dir, opts$sample, assembly_key,
    window = opts$windows
)
cat(
    sprintf(
        "legacy tables: %d recipient + %d donor windows\n",
        nrow(legacy$rcov), nrow(legacy$dcov)
    )
)

meta <- legacy$meta
annotate_meta <- function(table, role) {
    position <- match(table$chr_id, meta$chr_id)
    data.frame(
        chr_id = table$chr_id,
        start = as.numeric(table$start),
        stop = as.numeric(table$stop),
        role = role,
        chr_name = meta$chr_name[position],
        subgenome = meta$subgenome[position],
        stringsAsFactors = FALSE
    )
}
legacy_table <- rbind(
    annotate_meta(legacy$rcov, "recipient"),
    annotate_meta(legacy$dcov, "donor")
)
legacy_lengths <- c(
    stats::setNames(
        assembly$contigs$recipient$length,
        assembly$contigs$recipient$id
    ),
    stats::setNames(
        assembly$contigs$donor$length,
        assembly$contigs$donor$id
    )
)
legacy_gr <- GenomicRanges::GRanges(
    legacy_table$chr_id,
    IRanges::IRanges(legacy_table$start + 1, legacy_table$stop),
    seqlengths = legacy_lengths[unique(legacy_table$chr_id)]
)
legacy_windows <- list(
    gr = legacy_gr,
    table = legacy_table,
    n_per_chromosome = as.integer(opts$windows)
)

## 1. Counts on the legacy intervals ---------------------------------------

cat("\n[1/6] counting reads on legacy intervals ...\n")
legacy_counts <- record(env, "countReads on legacy intervals", {
    countReads(staged_bam, legacy_gr, mapq = 20L, exclude_flags = 2308L)
})
legacy_expected <- as.integer(c(legacy$rcov$cov, legacy$dcov$cov))
delta <- legacy_expected - legacy_counts
cat(
    sprintf(
        "  mismatched windows: %d / %d (%.3f%%), max delta %d, total |delta| %d\n",
        sum(delta != 0L), length(delta), 100 * mean(delta != 0L),
        max(abs(delta)), sum(abs(delta))
    )
)
check(env, "all count deltas are non-negative", all(delta >= 0L),
    sprintf("%d negative", sum(delta < 0L)))
check(env, "mismatched windows <= 2%", mean(delta != 0L) <= 0.02,
    sprintf("%.4f%%", 100 * mean(delta != 0L)))
check(env, "maximum count delta <= 20", max(abs(delta)) <= 20L,
    max(abs(delta)))
check(
    env, "total absolute delta < 0.1% of legacy counts",
    sum(abs(delta)) / sum(legacy_expected) < 0.001,
    sprintf("%.6f%%", 100 * sum(abs(delta)) / sum(legacy_expected))
)

if (nzchar(samtools)) {
    cat("  verifying each mismatch at the window with samtools ...\n")
    mismatch <- which(delta != 0L)
    outside_counts <- integer(length(mismatch))
    for (k in seq_along(mismatch)) {
        i <- mismatch[[k]]
        start <- legacy_table$start[[i]]
        stop <- legacy_table$stop[[i]]
        region <- sprintf(
            "%s:%s-%s", legacy_table$chr_id[[i]],
            format(start + 1, scientific = FALSE),
            format(stop, scientific = FALSE)
        )
        output <- suppressWarnings(system2(
            samtools,
            c(
                "view", "-q", "20", "-F", "2308",
                staged_bam, region
            ),
            stdout = TRUE, stderr = FALSE
        ))
        if (length(output) > 0L) {
            fields <- strsplit(output, "\t", fixed = TRUE)
            positions <- vapply(
                fields, function(x) as.integer(x[[4L]]), integer(1)
            )
            flags <- vapply(
                fields, function(x) as.integer(x[[2L]]), integer(1)
            )
            cigars <- vapply(fields, `[[`, character(1), 6L)
            ends <- positions + cigar_ref_length(cigars) - 1L
            ## bamsignals counts the 5' end: leftmost for forward reads,
            ## rightmost for reverse reads.
            five_prime <- ifelse(
                bitwAnd(flags, 16L) != 0L, ends, positions
            )
            outside_counts[[k]] <- sum(
                five_prime < start + 1L | five_prime > stop
            )
        }
    }
    explained <- delta[mismatch] == outside_counts
    check(
        env, "every mismatch equals reads whose 5' end is outside the window",
        all(explained),
        sprintf(
            "%d / %d explained; max unexplained %d",
            sum(explained), length(explained),
            if (any(!explained)) max(abs(delta[mismatch][!explained] -
                outside_counts[!explained])) else 0L
        )
    )
    cat(
        sprintf(
            "  5'-end assignment explains %d / %d mismatches\n",
            sum(explained), length(explained)
        )
    )
} else {
    cat("  samtools not found; skipping the boundary explanation\n")
}

## 2. GC on the legacy intervals -------------------------------------------

cat("\n[2/6] comparing GC content on legacy intervals ...\n")
legacy_gc <- record(env, "gcrank GC on legacy intervals", {
    MaMaMIA:::extractWindowGc(legacy_windows, assembly)
})
legacy_gc_expected <- c(legacy$rgc$gc, legacy$dgc$gc)
gc_delta <- abs(legacy_gc$gc - legacy_gc_expected)
cat(sprintf(
    "  max |GC delta| %.3e; identical at 6 decimals: %d / %d\n",
    max(gc_delta), sum(round(legacy_gc$gc, 6) == round(legacy_gc_expected, 6)),
    length(gc_delta)
))
check(env, "GC matches at the six-decimal legacy precision",
    max(gc_delta) < 5e-7, sprintf("max %.3e", max(gc_delta)))

## 3. Generated windows and per-contig totals ------------------------------

cat("\n[3/6] generating exact-N windows and counting reads ...\n")
windows <- record(env, "makeAnalysisWindows", {
    makeAnalysisWindows(assembly, NULL, opts$windows)
})
generated_counts <- record(env, "countReads on generated windows", {
    countReads(staged_bam, windows$gr, mapq = 20L, exclude_flags = 2308L)
})
per_contig <- data.frame(
    chr_id = windows$table$chr_id,
    count = generated_counts,
    stringsAsFactors = FALSE
)
per_contig <- stats::aggregate(count ~ chr_id, per_contig, sum)
per_contig$count <- as.numeric(per_contig$count)
legacy_contig <- stats::aggregate(count ~ chr_id, data.frame(
    chr_id = legacy_table$chr_id,
    count = legacy_expected,
    stringsAsFactors = FALSE
), sum)
legacy_contig$count <- as.numeric(legacy_contig$count)
totals <- merge(per_contig, legacy_contig, by = "chr_id", suffixes = c(".new", ".legacy"))
totals$delta <- totals$count.legacy - totals$count.new
cat(sprintf(
    "  per-contig legacy - start-based totals: min %d, max %d, total %d\n",
    min(totals$delta), max(totals$delta), sum(totals$delta)
))
if (nzchar(samtools)) {
    samtools_totals <- vapply(totals$chr_id, function(contig) {
        as.numeric(suppressWarnings(system2(
            samtools,
            c("view", "-q", "20", "-F", "2308", "-c", staged_bam, contig),
            stdout = TRUE, stderr = FALSE
        )))
    }, numeric(1))
    check(
        env, "5'-end counts equal samtools per-contig reads",
        all(totals$count.new == samtools_totals),
        sprintf(
            "max difference %s",
            max(abs(totals$count.new - samtools_totals))
        )
    )
    check(
        env, "legacy excess equals the per-window count delta",
        sum(totals$delta) == sum(delta),
        sprintf("sum(delta) = %d", sum(delta))
    )
} else {
    cat("  samtools not found; skipping per-contig verification\n")
}

## 4. Full pipeline on generated windows -----------------------------------

rca <- isa <- NULL
legacy_rca <- legacy_isa <- NULL
legacy_interval_isa <- NULL

if (!opts$skip_full) {
    cat("\n[4/6] full BAM-first pipeline (this is the slow step) ...\n")
    pipeline <- record(env, "analyzeBamInput (full pipeline)", {
        analyzeBamInput(
            staged_bam,
            config_path,
            normalization = "none",
            windows_per_chromosome = opts$windows,
            cores = opts$cores,
            cache = file.path(out_dir, "fit-cache-bam.rds")
        )
    })
    rca <- pipeline$rca[[1L]]
    isa <- pipeline$isa[[1L]]
    check(env, "full pipeline produced one RCA and one ISA",
        inherits(rca, "RCA") && inherits(isa, "ISA"))
    check(env, "full pipeline used normalization = 'none'",
        identical(pipeline$normalization$mode, "none"))

    cat("\n[5/6] legacy directory-layout pipeline on the legacy tables ...\n")
    legacy_rca <- record(env, "legacy RCA + correction", {
        corrected <- correctReadCounts(
            do.call(RCA, legacy[c("dcov", "dgc", "rcov", "rgc", "meta")]),
            cores = opts$cores,
            verbose = FALSE,
            cache = file.path(out_dir, "fit-cache-legacy.rds")
        )
        corrected
    })
    legacy_isa <- record(env, "legacy segments", {
        segments(
            legacy_rca,
            target_pairs = assembly$pairs,
            verbose = FALSE
        )
    })

    cat("\n[6/6] BAM-first core on the same legacy intervals ...\n")
    legacy_new_tables <- record(env, "buildSampleTables on legacy intervals", {
        buildSampleTables(
            legacy_counts, legacy_gc, legacy_windows, assembly, "Sq54"
        )
    })
    legacy_interval_isa <- record(env, "BAM-first core on legacy intervals", {
        corrected <- correctReadCounts(
            do.call(
                RCA,
                legacy_new_tables[c("dcov", "dgc", "rcov", "rgc", "meta")]
            ),
            cores = opts$cores,
            verbose = FALSE,
            cache = file.path(out_dir, "fit-cache-intervals.rds")
        )
        segments(
            corrected,
            target_pairs = legacy_new_tables$target_pairs,
            verbose = FALSE
        )
    })

    legacy_corrected <- legacy_rca$data[, c(
        "chr_id", "iid", "cor.gc", "ideal"
    )]
    interval_corrected <- legacy_interval_isa$data[, c(
        "chr_id", "iid", "cor.gc", "ideal"
    )]
    merged_corrected <- merge(
        legacy_corrected, interval_corrected,
        by = c("chr_id", "iid"), suffixes = c(".legacy", ".bam")
    )
    both_ideal <- merged_corrected$ideal.legacy & merged_corrected$ideal.bam
    rca_relative <- abs(
        merged_corrected$cor.gc.legacy - merged_corrected$cor.gc.bam
    ) / pmax(abs(merged_corrected$cor.gc.legacy), 1e-12)
    rca_q99 <- as.numeric(stats::quantile(
        rca_relative[both_ideal], 0.99, na.rm = TRUE
    ))
    rca_max <- max(rca_relative[both_ideal], na.rm = TRUE)
    check(
        env, "corrected coverage agrees on ideal windows",
        nrow(merged_corrected) == nrow(interval_corrected) && rca_q99 < 5e-3,
        sprintf(
            "99%% relative %.3e, max %.3e over %d ideal windows (%d paired rows)",
            rca_q99, rca_max, sum(both_ideal), nrow(merged_corrected)
        )
    )

    segment_table <- function(x) {
        data.frame(
            chr_id.don = x$segments$chr_id.don,
            chr_id.rec = x$segments$chr_id.rec,
            loc.start = x$segments$loc.start,
            loc.end = x$segments$loc.end,
            seg.median = x$segments$seg.median,
            putative = x$segments$putative_introgression,
            stringsAsFactors = FALSE
        )
    }
    legacy_segments <- segment_table(legacy_isa)
    interval_segments <- segment_table(legacy_interval_isa)
    merged_segments <- merge(
        legacy_segments, interval_segments,
        by = c("chr_id.don", "chr_id.rec", "loc.start", "loc.end"),
        suffixes = c(".legacy", ".bam"), all = TRUE
    )
    check(
        env, "segment breakpoints are identical on the same intervals",
        nrow(merged_segments) == nrow(legacy_segments) &&
            nrow(merged_segments) == nrow(interval_segments),
        sprintf(
            "%d vs %d segments; %d matched rows",
            nrow(legacy_segments), nrow(interval_segments),
            nrow(merged_segments)
        )
    )
    if (nrow(merged_segments) == nrow(legacy_segments) &&
        nrow(merged_segments) == nrow(interval_segments)) {
        median_delta <- max(abs(
            merged_segments$seg.median.legacy -
                merged_segments$seg.median.bam
        ))
        median_scale <- max(abs(merged_segments$seg.median.legacy))
        median_relative <- if (median_scale == 0) {
            0
        } else {
            median_delta / median_scale
        }
        check(
            env, "segment medians agree on the same intervals",
            median_relative < 0.01,
            sprintf(
                "max |delta| %.4f (%.4f%% of scale)",
                median_delta, 100 * median_relative
            )
        )
        check(
            env, "putative calls agree on the same intervals",
            identical(
                merged_segments$putative.legacy,
                merged_segments$putative.bam
            ),
            sprintf(
                "%d vs %d putative",
                sum(merged_segments$putative.legacy),
                sum(merged_segments$putative.bam)
            )
        )
    }
    saveRDS(
        list(
            legacy_corrected = legacy_corrected,
            interval_corrected = interval_corrected,
            bam_segments = isa$segments,
            legacy_segments = legacy_isa$segments,
            interval_segments = legacy_interval_isa$segments
        ),
        file.path(out_dir, "comparison.rds")
    )
} else {
    cat("\n[4-6/6] skipped (--skip-full)\n")
}

## Segment-level comparison between the generated windows and the legacy
## pipeline. Windows differ by up to one window width, so exact equality is not
## expected; the number of segments, the number of putative segments and the
## total introgression length should agree closely.
segment_summary <- function(isa) {
    if (is.null(isa) || nrow(isa$segments) == 0L) {
        return(data.frame(
            chr_id.don = character(0), chr_id.rec = character(0),
            n_segments = integer(0), n_putative = integer(0),
            bp = numeric(0), stringsAsFactors = FALSE
        ))
    }
    segments <- isa$segments
    annotated <- MaMaMIA:::annotate_segment_coords(segments, isa$out)
    annotated$bp <- annotated$stop.don - annotated$start.don
    pieces <- split(
        annotated,
        list(annotated$chr_id.don, annotated$chr_id.rec),
        drop = TRUE
    )
    out <- do.call(rbind, lapply(pieces, function(x) {
        data.frame(
            chr_id.don = x$chr_id.don[[1L]],
            chr_id.rec = x$chr_id.rec[[1L]],
            n_segments = nrow(x),
            n_putative = sum(x$putative_introgression),
            bp = sum(x$bp[x$putative_introgression]),
            stringsAsFactors = FALSE
        )
    }))
    rownames(out) <- NULL
    out
}

if (!is.null(isa) && !is.null(legacy_isa)) {
    generated <- segment_summary(isa)
    reference <- segment_summary(legacy_isa)
    comparison <- merge(
        generated, reference,
        by = c("chr_id.don", "chr_id.rec"),
        suffixes = c(".bam", ".legacy"),
        all = TRUE
    )
    comparison[is.na(comparison)] <- 0
    comparison$segment_delta <- comparison$n_segments.bam - comparison$n_segments.legacy
    comparison$putative_delta <- comparison$n_putative.bam - comparison$n_putative.legacy
    comparison$bp_relative <- ifelse(
        comparison$bp.legacy == 0,
        ifelse(comparison$bp.bam == 0, 0, Inf),
        abs(comparison$bp.bam - comparison$bp.legacy) / comparison$bp.legacy
    )
    print(comparison, row.names = FALSE)
    check(env, "both pipelines found the same donor-recipient pairs",
        nrow(comparison) == nrow(generated) &&
            nrow(comparison) == nrow(reference),
        sprintf("%d pairs", nrow(comparison)))
    check(env, "segment counts differ by at most one per pair",
        all(abs(comparison$segment_delta) <= 1L),
        sprintf("max %d", max(abs(comparison$segment_delta))))
    check(env, "putative segment counts differ by at most one per pair",
        all(abs(comparison$putative_delta) <= 1L),
        sprintf("max %d", max(abs(comparison$putative_delta))))
    check(env, "introgression length agrees within 2%",
        all(comparison$bp_relative <= 0.02),
        sprintf("max %.4f", max(comparison$bp_relative)))
}

## Outputs ------------------------------------------------------------------

timings <- do.call(rbind, env$timings)
checks <- do.call(rbind, env$checks)
utils::write.csv(timings, file.path(out_dir, "TIMINGS.csv"), row.names = FALSE)
utils::write.csv(checks, file.path(out_dir, "CHECKS.csv"), row.names = FALSE)

fingerprint <- c(
    package_version = as.character(utils::packageVersion("MaMaMIA")),
    r_version = as.character(getRversion()),
    htslib_version = MaMaMIA:::.htslib_version(),
    gcrank_commit = MaMaMIA:::.MaMaMIAgcrankCommit,
    sample = opts$sample,
    windows = opts$windows,
    n_windows = nrow(windows$table),
    n_legacy_windows = nrow(legacy_table),
    n_count_mismatches = sum(delta != 0L),
    max_count_delta = max(abs(delta)),
    total_abs_count_delta = sum(abs(delta)),
    max_gc_delta = max(gc_delta),
    count_hash = hash_object(legacy_counts),
    generated_count_hash = hash_object(generated_counts),
    gc_hash = hash_object(legacy_gc$gc),
    isa_hash = if (is.null(isa)) NA_character_ else hash_object(isa$segments),
    legacy_isa_hash = if (is.null(legacy_isa)) NA_character_ else hash_object(legacy_isa$segments),
    interval_isa_hash = if (is.null(legacy_interval_isa)) {
        NA_character_
    } else {
        hash_object(legacy_interval_isa$segments)
    }
)
writeLines(
    c(
        "# MaMaMIA BAM-first golden benchmark",
        "[environment]",
        sprintf("%s: %s", names(fingerprint[1:4]), fingerprint[1:4]),
        "[counts]",
        sprintf("%s: %s", names(fingerprint[5:14]), fingerprint[5:14]),
        "[fingerprints]",
        sprintf("%s: %s", names(fingerprint[15:length(fingerprint)]),
            fingerprint[15:length(fingerprint)]),
        "[checks]",
        sprintf("%s: %s", gsub("[^A-Za-z0-9]+", "_", checks$check),
            checks$result)
    ),
    file.path(out_dir, "FINGERPRINT.txt")
)

cat("\n== checks ==\n")
print(checks, row.names = FALSE)
cat("\n== timings ==\n")
print(timings, row.names = FALSE)
cat("\noutputs written to ", out_dir, "\n", sep = "")

if (any(checks$result != "PASS")) {
    quit(status = 1L)
}
cat("all golden checks passed\n")
