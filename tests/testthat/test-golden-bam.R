## T3.4: golden regression against the legacy Sq54 BED tables.
##
## The count/GC checks run whenever test_data is available. The full
## correction/segmentation comparison is slow and opt-in:
##
##   MAMAMIA_SLOW_TESTS=1 MAMAMIA_TEST_DATA=1
##
## The same checks are implemented in tools/benchmark-bam-pipeline.R.

golden_skip_slow <- function() {
    if (!identical(Sys.getenv("MAMAMIA_SLOW_TESTS"), "1")) {
        testthat::skip("set MAMAMIA_SLOW_TESTS=1 to run the full golden pipeline")
    }
}

golden_cores <- function() {
    detected <- parallel::detectCores()
    if (is.na(detected) || detected < 1L) {
        return(1L)
    }
    as.integer(min(4L, detected))
}

golden_config <- function(root, session) {
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
    config <- list(
        version = 1L,
        defaults = list(
            windows_per_chromosome = 1000L,
            allowed_windows = c(100L, 1000L, 10000L),
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
    path <- file.path(session, "assemblies-real.yml")
    writeLines(yaml::as.yaml(config), path)
    readMaMaMIAConfig(path)
}

golden_legacy_windows <- function(legacy, assembly, n = 1000L) {
    meta <- legacy$meta
    annotate <- function(table, role) {
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
    table <- rbind(
        annotate(legacy$rcov, "recipient"),
        annotate(legacy$dcov, "donor")
    )
    lengths <- c(
        stats::setNames(
            assembly$contigs$recipient$length,
            assembly$contigs$recipient$id
        ),
        stats::setNames(
            assembly$contigs$donor$length,
            assembly$contigs$donor$id
        )
    )
    gr <- GenomicRanges::GRanges(
        table$chr_id,
        IRanges::IRanges(table$start + 1, table$stop),
        seqlengths = lengths[unique(table$chr_id)]
    )
    list(gr = gr, table = table, n_per_chromosome = as.integer(n))
}

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

golden_mismatch_explained <- function(samtools, bam, table, delta) {
    mismatch <- which(delta != 0L)
    explained <- logical(length(mismatch))
    for (k in seq_along(mismatch)) {
        i <- mismatch[[k]]
        start <- table$start[[i]]
        stop <- table$stop[[i]]
        region <- sprintf(
            "%s:%s-%s", table$chr_id[[i]],
            format(start + 1, scientific = FALSE),
            format(stop, scientific = FALSE)
        )
        output <- suppressWarnings(system2(
            samtools,
            c("view", "-q", "20", "-F", "2308", bam, region),
            stdout = TRUE, stderr = FALSE
        ))
        if (length(output) == 0L) {
            explained[[k]] <- FALSE
            next
        }
        fields <- strsplit(output, "\t", fixed = TRUE)
        positions <- vapply(fields, function(x) as.integer(x[[4L]]), integer(1))
        flags <- vapply(fields, function(x) as.integer(x[[2L]]), integer(1))
        cigars <- vapply(fields, `[[`, character(1), 6L)
        ends <- positions + cigar_ref_length(cigars) - 1L
        five_prime <- ifelse(bitwAnd(flags, 16L) != 0L, ends, positions)
        outside <- sum(five_prime < start + 1L | five_prime > stop)
        explained[[k]] <- identical(as.integer(outside), as.integer(delta[[i]]))
    }
    explained
}

golden_segment_summary <- function(isa) {
    if (is.null(isa) || nrow(isa$segments) == 0L) {
        return(data.frame(
            chr_id.don = character(0), chr_id.rec = character(0),
            n_segments = integer(0), n_putative = integer(0),
            bp = numeric(0), stringsAsFactors = FALSE
        ))
    }
    annotated <- MaMaMIA:::annotate_segment_coords(isa$segments, isa$out)
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

golden_segment_table <- function(isa) {
    data.frame(
        chr_id.don = isa$segments$chr_id.don,
        chr_id.rec = isa$segments$chr_id.rec,
        loc.start = isa$segments$loc.start,
        loc.end = isa$segments$loc.end,
        seg.median = isa$segments$seg.median,
        putative = isa$segments$putative_introgression,
        stringsAsFactors = FALSE
    )
}

golden_legacy_isa <- function(data, cache) {
    rca <- do.call(
        RCA,
        data$legacy[c("dcov", "dgc", "rcov", "rgc", "meta")]
    )
    corrected <- correctReadCounts(
        rca,
        cores = golden_cores(),
        verbose = FALSE,
        cache = cache
    )
    segments(
        corrected,
        target_pairs = data$assembly$pairs,
        verbose = FALSE
    )
}

golden_state <- new.env(parent = emptyenv())

golden_data <- function() {
    if (!is.null(golden_state$data)) {
        return(golden_state$data)
    }
    root <- skip_if_no_test_data()
    skip_if_no_samtools()
    bam_path <- file.path(
        root, "bam_files", "Sq54", "cs_t2t_se_cereale.bam"
    )
    legacy_dir <- file.path(
        root, "mamamia_analysis", "arsenal", "Sq54"
    )
    required <- c(
        bam_path,
        file.path(
            legacy_dir, "cs_t2t_se_cereale.w1000.rcov.bed"
        ),
        file.path(
            legacy_dir, "cs_t2t_se_cereale.w1000.dcov.bed"
        )
    )
    if (!all(file.exists(required))) {
        testthat::skip("the Sq54 BAM or legacy BED tables are missing")
    }

    session <- tempfile("mamamia-golden-")
    dir.create(session)
    config <- golden_config(root, session)
    assembly <- config$assemblies$cs_t2t_se_cereale
    staged_bam <- MaMaMIA:::stageBam(bam_path, session)
    ensureBamIndex(staged_bam)

    legacy <- read_golden_legacy(
        root, "Sq54", "cs_t2t_se_cereale", 1000L
    )
    legacy_windows <- golden_legacy_windows(legacy, assembly)
    legacy_counts <- countReads(
        staged_bam, legacy_windows$gr,
        mapq = 20L, exclude_flags = 2308L
    )
    legacy_expected <- as.integer(c(legacy$rcov$cov, legacy$dcov$cov))
    delta <- legacy_expected - legacy_counts
    legacy_gc <- MaMaMIA:::extractWindowGc(legacy_windows, assembly)
    legacy_gc_expected <- c(legacy$rgc$gc, legacy$dgc$gc)

    golden_state$data <- list(
        root = root,
        session = session,
        config = config,
        config_path = file.path(session, "assemblies-real.yml"),
        assembly = assembly,
        staged_bam = staged_bam,
        legacy = legacy,
        legacy_windows = legacy_windows,
        legacy_counts = legacy_counts,
        legacy_expected = legacy_expected,
        delta = delta,
        legacy_gc = legacy_gc,
        legacy_gc_expected = legacy_gc_expected,
        generated_windows = makeAnalysisWindows(assembly, NULL, 1000L),
        samtools = Sys.which("samtools")
    )
    golden_state$data
}

test_that("legacy count mismatches are explained by the 5' end rule", {
    data <- golden_data()
    delta <- data$delta

    expect_true(all(delta >= 0L))
    expect_lte(mean(delta != 0L), 0.02)
    expect_lte(max(abs(delta)), 20L)
    expect_lt(sum(abs(delta)) / sum(data$legacy_expected), 0.001)

    explained <- golden_mismatch_explained(
        data$samtools, data$staged_bam, data$legacy_windows$table, delta
    )
    expect_true(all(explained))
    sides <- rep(
        c("recipient", "donor"),
        c(nrow(data$legacy$rcov), nrow(data$legacy$dcov))
    )
    by_side <- tapply(delta != 0L, sides, sum)
    totals <- tapply(rep(TRUE, length(sides)), sides, sum)
    message(sprintf(
        "Sq54: %d / %d windows differ (recipient %d/%d, donor %d/%d), max delta %d, total |delta| %d; %d / %d explained by 5' ends",
        sum(delta != 0L), length(delta),
        by_side[["recipient"]], totals[["recipient"]],
        by_side[["donor"]], totals[["donor"]],
        max(abs(delta)), sum(abs(delta)),
        sum(explained), length(explained)
    ))
})

test_that("legacy GC matches at the six-decimal BED precision", {
    data <- golden_data()
    gc_delta <- abs(data$legacy_gc$gc - data$legacy_gc_expected)

    expect_lt(max(gc_delta), 5e-7)
    expect_identical(
        sum(round(data$legacy_gc$gc, 6) == round(data$legacy_gc_expected, 6)),
        length(gc_delta)
    )
    message(sprintf(
        "Sq54: max |GC delta| %.3e over %d windows",
        max(gc_delta), length(gc_delta)
    ))
})

test_that("generated windows follow the exact-N tiling rule", {
    data <- golden_data()
    windows <- data$generated_windows
    expect_identical(windows$n_per_chromosome, 1000L)
    expect_identical(nrow(windows$table), 28000L)

    targets <- rbind(
        data$assembly$contigs$recipient,
        data$assembly$contigs$donor
    )
    for (i in seq_len(nrow(targets))) {
        idx <- which(windows$table$chr_id == targets$id[[i]])
        expect_length(idx, 1000L)
        expected_start <- floor((0:999) * targets$length[[i]] / 1000)
        expected_stop <- floor((1:1000) * targets$length[[i]] / 1000)
        expect_identical(as.numeric(windows$table$start[idx]), expected_start)
        expect_identical(as.numeric(windows$table$stop[idx]), expected_stop)
        expect_identical(as.numeric(windows$table$start[idx][[1L]]), 0)
        expect_identical(
            as.numeric(windows$table$stop[idx][[1000L]]),
            as.numeric(targets$length[[i]])
        )
    }

    ## The legacy tables use fixed-width windows with a stretched last window;
    ## document the resulting shift instead of pretending the intervals match.
    legacy_table <- data$legacy_windows$table
    position <- match(
        paste(windows$table$chr_id, ave(
            seq_len(nrow(windows$table)), windows$table$chr_id, FUN = seq_along
        )),
        paste(legacy_table$chr_id, ave(
            seq_len(nrow(legacy_table)), legacy_table$chr_id, FUN = seq_along
        ))
    )
    shift <- abs(
        as.numeric(windows$table$start) -
            as.numeric(legacy_table$start[position])
    )
    expect_lt(max(shift), 1000)
    expect_lt(sum(shift == 0), nrow(legacy_table))
    message(sprintf(
        "Sq54: exact-N windows shift by at most %d bp relative to the legacy fixed-width tiling",
        max(shift)
    ))
})

test_that("the full BAM pipeline matches the legacy segment summary", {
    golden_skip_slow()
    data <- golden_data()
    result <- analyzeBamInput(
        data$staged_bam,
        data$config_path,
        normalization = "none",
        windows_per_chromosome = 1000L,
        cores = golden_cores(),
        cache = file.path(data$session, "fit-cache-bam.rds")
    )

    expect_named(result$rca, "Sq54")
    expect_named(result$isa, "Sq54")
    expect_identical(result$normalization$mode, "none")

    table <- result$windows$table
    counts <- numeric(nrow(table))
    counts[table$role == "donor"] <- result$counts$donor[, 1L]
    counts[table$role == "recipient"] <- result$counts$recipient[, 1L]
    per_contig <- stats::aggregate(
        count ~ chr_id,
        data.frame(chr_id = table$chr_id, count = counts),
        sum
    )
    samtools_totals <- vapply(as.character(per_contig$chr_id), function(contig) {
        as.numeric(suppressWarnings(system2(
            data$samtools,
            c(
                "view", "-q", "20", "-F", "2308",
                "-c", data$staged_bam, contig
            ),
            stdout = TRUE, stderr = FALSE
        )))
    }, numeric(1))
    expect_equal(as.numeric(per_contig$count), unname(samtools_totals))

    legacy_isa <- golden_legacy_isa(
        data, file.path(data$session, "fit-cache-legacy.rds")
    )
    generated <- golden_segment_summary(result$isa[[1L]])
    reference <- golden_segment_summary(legacy_isa)
    comparison <- merge(
        generated, reference,
        by = c("chr_id.don", "chr_id.rec"),
        suffixes = c(".bam", ".legacy"),
        all = TRUE
    )
    comparison[is.na(comparison)] <- 0

    expect_identical(nrow(comparison), nrow(generated))
    expect_identical(nrow(comparison), nrow(reference))
    expect_lte(max(abs(comparison$n_segments.bam - comparison$n_segments.legacy)), 1L)
    expect_lte(max(abs(comparison$n_putative.bam - comparison$n_putative.legacy)), 1L)
    relative_bp <- ifelse(
        comparison$bp.legacy == 0,
        ifelse(comparison$bp.bam == 0, 0, Inf),
        abs(comparison$bp.bam - comparison$bp.legacy) / comparison$bp.legacy
    )
    expect_lte(max(relative_bp), 0.02)
    message(sprintf(
        "Sq54: %d pairs; segment count delta <= %d, putative delta <= %d, introgression bp relative delta %.3e",
        nrow(comparison),
        max(abs(comparison$n_segments.bam - comparison$n_segments.legacy)),
        max(abs(comparison$n_putative.bam - comparison$n_putative.legacy)),
        max(relative_bp)
    ))
})

test_that("the BAM-first core reproduces the legacy correction on the same intervals", {
    golden_skip_slow()
    data <- golden_data()
    legacy_isa <- golden_legacy_isa(
        data, file.path(data$session, "fit-cache-legacy.rds")
    )
    tables <- buildSampleTables(
        data$legacy_counts, data$legacy_gc, data$legacy_windows,
        data$assembly, "Sq54"
    )
    corrected <- correctReadCounts(
        do.call(RCA, tables[c("dcov", "dgc", "rcov", "rgc", "meta")]),
        cores = golden_cores(),
        verbose = FALSE,
        cache = file.path(data$session, "fit-cache-intervals.rds")
    )
    isa <- segments(
        corrected,
        target_pairs = tables$target_pairs,
        verbose = FALSE
    )

    legacy_corrected <- legacy_isa$data[, c("chr_id", "iid", "cor.gc", "ideal")]
    interval_corrected <- isa$data[, c("chr_id", "iid", "cor.gc", "ideal")]
    merged <- merge(
        legacy_corrected, interval_corrected,
        by = c("chr_id", "iid"), suffixes = c(".legacy", ".bam")
    )
    both_ideal <- merged$ideal.legacy & merged$ideal.bam
    relative <- abs(
        merged$cor.gc.legacy - merged$cor.gc.bam
    ) / pmax(abs(merged$cor.gc.legacy), 1e-12)
    expect_identical(nrow(merged), nrow(interval_corrected))
    expect_lt(
        as.numeric(stats::quantile(relative[both_ideal], 0.99, na.rm = TRUE)),
        5e-3
    )

    legacy_segments <- golden_segment_table(legacy_isa)
    interval_segments <- golden_segment_table(isa)
    merged_segments <- merge(
        legacy_segments, interval_segments,
        by = c("chr_id.don", "chr_id.rec", "loc.start", "loc.end"),
        suffixes = c(".legacy", ".bam"),
        all = TRUE
    )
    expect_identical(nrow(merged_segments), nrow(legacy_segments))
    expect_identical(nrow(merged_segments), nrow(interval_segments))
    expect_identical(
        merged_segments$putative.legacy,
        merged_segments$putative.bam
    )
    median_delta <- max(abs(
        merged_segments$seg.median.legacy - merged_segments$seg.median.bam
    ))
    median_scale <- max(abs(merged_segments$seg.median.legacy))
    expect_lt(if (median_scale == 0) 0 else median_delta / median_scale, 0.01)
    message(sprintf(
        "Sq54: %d segments; breakpoints identical; max median delta %.4f (%.4f%% of scale)",
        nrow(merged_segments), median_delta,
        100 * if (median_scale == 0) 0 else median_delta / median_scale
    ))
})
