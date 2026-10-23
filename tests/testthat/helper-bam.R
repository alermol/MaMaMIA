## Helpers for the tiny committed BAM/gcrank fixtures and gated real data.

tiny_fixture_dir <- function() {
    testthat::test_path("fixtures")
}

tiny_fixture <- function(...) {
    file.path(tiny_fixture_dir(), ...)
}

skip_if_no_tiny_fixtures <- function() {
    required <- c(
        "tiny.fa",
        "tiny.bam",
        "tiny.bam.csi",
        "tiny.gcrank",
        "tiny.gcrank/contigs.tsv",
        "tiny.gcrank/manifest.json",
        "assemblies.yml"
    )
    missing <- required[!file.exists(tiny_fixture(required))]
    if (length(missing) > 0L) {
        testthat::skip(paste(
            "tiny fixtures are missing:",
            paste(missing, collapse = ", ")
        ))
    }
}

tiny_config_path <- function() {
    tiny_fixture("assemblies.yml")
}

tiny_bam_path <- function() {
    tiny_fixture("tiny.bam")
}

tiny_gcrank_index <- function() {
    tiny_fixture("tiny.gcrank")
}

tiny_contigs <- function() {
    readContigsTsv(file.path(tiny_gcrank_index(), "contigs.tsv"))
}

have_samtools <- function() {
    nzchar(Sys.which("samtools"))
}

skip_if_no_samtools <- function() {
    if (!have_samtools()) {
        testthat::skip("samtools is not available")
    }
}

## Build a BAM from SAM text for header-error tests. Indexing is optional
## because bamHeaderInfo() never needs an index.
write_test_bam <- function(lines, path, index = FALSE) {
    sam <- tempfile(fileext = ".sam")
    on.exit(unlink(sam), add = TRUE)
    writeLines(lines, sam)
    status <- system2(
        "samtools",
        c("view", "--no-PG", "-b", "-o", path, sam),
        stdout = FALSE, stderr = FALSE
    )
    if (!identical(status, 0L)) {
        testthat::skip("samtools could not create the test BAM")
    }
    if (index) {
        status <- system2(
            "samtools", c("index", "-c", path),
            stdout = FALSE, stderr = FALSE
        )
        if (!identical(status, 0L)) {
            testthat::skip("samtools could not index the test BAM")
        }
    }
    path
}

## The committed tiny fixture has exactly 40% and 60% GC, so the ZINB smooth
## cannot be fitted on it. These helpers build a small assembly whose GC varies
## between 20% and 100% in 20 bp blocks, plus a densely covered BAM, entirely in
## session temp. They need samtools and the gcrank builder.

gcrank_build_binary <- function() {
    candidate <- Sys.which("gcrank-build")
    if (nzchar(candidate)) {
        return(candidate)
    }
    candidate <- file.path(
        mamamia_test_data_dir(), "gcidx", "gcrank-build"
    )
    if (file.exists(candidate)) {
        return(candidate)
    }
    ""
}

variable_gc_sequence <- function(length_bp, offset) {
    blocks <- length_bp / 20L
    gc_counts <- ((seq_len(blocks) + offset) %% 17L) + 4L
    paste0(vapply(gc_counts, function(gc) {
        paste0(strrep("G", gc), strrep("A", 20L - gc))
    }, character(1)), collapse = "")
}

variable_gc_sam <- function(sample, contig_length = 20000L,
                            read_length = 50L, step = 4L) {
    sequence <- strrep("A", read_length)
    quality <- strrep("I", read_length)
    positions <- seq(1L, contig_length - read_length + 1L, by = step)
    make_read <- function(contig, position) {
        sprintf(
            "r_%s_%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s",
            contig, position, contig, position, read_length,
            sequence, quality
        )
    }
    c(
        "@HD\tVN:1.6\tSO:coordinate",
        sprintf("@SQ\tSN:d1\tLN:%d", contig_length),
        sprintf("@SQ\tSN:r1\tLN:%d", contig_length),
        sprintf("@RG\tID:%s\tSM:%s", sample, sample),
        make_read("d1", positions),
        make_read("r1", positions)
    )
}

add_variable_gc_bam <- function(fixture, sample, filename = NULL, step = 4L) {
    if (is.null(filename)) {
        filename <- paste0(sample, ".bam")
    }
    contig_length <- if (is.null(fixture$contig_length)) {
        20000L
    } else {
        fixture$contig_length
    }
    path <- file.path(fixture$directory, filename)
    write_test_bam(
        variable_gc_sam(
            sample, contig_length = contig_length, step = step
        ),
        path,
        index = FALSE
    )
    path
}

## One coordinate-sorted SAM containing several samples, each with its own
## @RG ID and RG tag on every read. Reads are emitted in coordinate order so
## samtools can index the resulting BAM directly.
multi_sample_variable_gc_sam <- function(samples, contig_length = 20000L,
                                         read_length = 50L,
                                         steps = NULL) {
    if (is.null(steps)) {
        steps <- rep(4L, length(samples))
    }
    if (length(steps) == 1L) {
        steps <- rep(steps, length(samples))
    }
    if (length(steps) != length(samples)) {
        stop("steps must have length one or length(samples)", call. = FALSE)
    }
    sequence <- strrep("A", read_length)
    quality <- strrep("I", read_length)
    reads <- list()
    for (sample_index in seq_along(samples)) {
        sample <- samples[[sample_index]]
        positions <- seq(
            1L, contig_length - read_length + 1L, by = steps[[sample_index]]
        )
        for (contig in c("d1", "r1")) {
            reads[[length(reads) + 1L]] <- data.frame(
                contig = contig,
                position = positions,
                sample = sample,
                stringsAsFactors = FALSE
            )
        }
    }
    reads <- do.call(rbind, reads)
    reads <- reads[order(
        match(reads$contig, c("d1", "r1")), reads$position, reads$sample
    ), , drop = FALSE]
    read_lines <- sprintf(
        "r_%s_%s_%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\tRG:Z:%s",
        reads$sample, reads$contig, reads$position, reads$contig,
        reads$position, read_length, sequence, quality, reads$sample
    )
    c(
        "@HD\tVN:1.6\tSO:coordinate",
        sprintf("@SQ\tSN:d1\tLN:%d", contig_length),
        sprintf("@SQ\tSN:r1\tLN:%d", contig_length),
        sprintf("@RG\tID:%s\tSM:%s", samples, samples),
        read_lines
    )
}

add_multi_sample_gc_bam <- function(fixture, samples, filename = "multi.bam",
                                    steps = NULL) {
    contig_length <- if (is.null(fixture$contig_length)) {
        20000L
    } else {
        fixture$contig_length
    }
    path <- file.path(fixture$directory, filename)
    write_test_bam(
        multi_sample_variable_gc_sam(
            samples,
            contig_length = contig_length,
            steps = steps
        ),
        path,
        index = TRUE
    )
    path
}

variable_gc_fixture <- function(sample = "synthetic", windows = 100L,
                                contig_length = 20000L,
                                allowed_windows = c(100L, 1000L, 10000L),
                                species = FALSE) {
    skip_if_no_samtools()
    builder <- gcrank_build_binary()
    if (!nzchar(builder)) {
        testthat::skip("gcrank-build is not available")
    }

    directory <- tempfile("variable-gc-")
    dir.create(directory)
    fasta <- file.path(directory, "variable.fa")
    writeLines(
        c(
            ">d1", variable_gc_sequence(contig_length, 0L),
            ">r1", variable_gc_sequence(contig_length, 8L)
        ),
        fasta
    )
    status <- system2(
        builder,
        c(
            "--fasta", fasta, "--out", directory,
            "--name", "variable", "--threads", "1"
        ),
        stdout = FALSE, stderr = FALSE
    )
    if (!identical(status, 0L)) {
        unlink(directory, recursive = TRUE)
        testthat::skip("gcrank-build could not create the variable-GC index")
    }

    config <- list(
        version = 1L,
        defaults = list(
            windows_per_chromosome = as.integer(windows),
            allowed_windows = as.integer(allowed_windows),
            mapq = 20L,
            exclude_flags = 2308L,
            max_fit_windows = 1000L,
            gc_mmap = TRUE,
            gc_validate_crc = FALSE,
            normalization_by_side = TRUE
        ),
        server = list(max_bams = 4L, max_total_bytes = 20 * 1024^2)
    )
    index <- file.path(directory, "variable.gcrank")
    if (species) {
        config$species <- list(
            recipient = list(
                label = "Variable-GC recipient",
                default_subgenome = "R",
                gc_indexes = list(index),
                contigs = list(list(
                    id = "r1", name = "1R", subgenome = "R"
                ))
            ),
            donor = list(
                label = "Variable-GC donor",
                default_subgenome = "D",
                gc_indexes = list(index),
                contigs = list(list(
                    id = "d1", name = "1D", subgenome = "D"
                ))
            )
        )
    } else {
        config$assemblies <- list(
            variable = list(
                label = "Variable-GC synthetic assembly",
                description = "Two 20 kb contigs with variable GC",
                gc_indexes = list(index),
                contigs = list(
                    recipient = list(list(
                        id = "r1", name = "1R", subgenome = "R"
                    )),
                    donor = list(list(
                        id = "d1", name = "1D", subgenome = "D"
                    ))
                ),
                pairs = list(d1 = "r1")
            )
        )
    }
    config_path <- file.path(directory, "assemblies.yml")
    writeLines(yaml::as.yaml(config), config_path)

    fixture <- list(
        directory = directory,
        config_path = config_path,
        sample = sample,
        contig_length = as.integer(contig_length)
    )
    fixture$bam_path <- add_variable_gc_bam(fixture, sample)
    fixture
}

## A synthetic species pair with two recipient subgenomes (A and B) that
## share one donor subgenome (R). The BAM contains d1, r1 and r2, so switching
## the recipient subgenome is a valid assembly change. Used to test that the
## first analysis precomputes both combinations.
multi_subgenome_fixture <- function(samples = c("sample1", "sample2"),
                                    contig_length = 20000L) {
    skip_if_no_samtools()
    builder <- gcrank_build_binary()
    if (!nzchar(builder)) {
        testthat::skip("gcrank-build is not available")
    }

    directory <- tempfile("multi-subgenome-")
    dir.create(directory)
    fasta <- file.path(directory, "multi.fa")
    writeLines(
        c(
            ">d1", variable_gc_sequence(contig_length, 0L),
            ">r1", variable_gc_sequence(contig_length, 8L),
            ">r2", variable_gc_sequence(contig_length, 3L)
        ),
        fasta
    )
    status <- system2(
        builder,
        c(
            "--fasta", fasta, "--out", directory,
            "--name", "multi", "--threads", "1"
        ),
        stdout = FALSE, stderr = FALSE
    )
    if (!identical(status, 0L)) {
        unlink(directory, recursive = TRUE)
        testthat::skip("gcrank-build could not create the multi-subgenome index")
    }

    config <- list(
        version = 1L,
        defaults = list(
            windows_per_chromosome = 100L,
            allowed_windows = c(100L, 1000L, 10000L),
            mapq = 20L,
            exclude_flags = 2308L,
            max_fit_windows = 1000L,
            gc_mmap = TRUE,
            gc_validate_crc = FALSE,
            normalization_by_side = TRUE
        ),
        server = list(max_bams = 4L, max_total_bytes = 20 * 1024^2),
        species = list(
            recipient = list(
                label = "Multi-subgenome recipient",
                default_subgenome = "A",
                gc_indexes = list(file.path(directory, "multi.gcrank")),
                contigs = list(
                    list(id = "r1", name = "1A", subgenome = "A"),
                    list(id = "r2", name = "1B", subgenome = "B")
                )
            ),
            donor = list(
                label = "Multi-subgenome donor",
                default_subgenome = "R",
                gc_indexes = list(file.path(directory, "multi.gcrank")),
                contigs = list(
                    list(id = "d1", name = "1R", subgenome = "R")
                )
            )
        )
    )
    config_path <- file.path(directory, "assemblies.yml")
    writeLines(yaml::as.yaml(config), config_path)

    read_length <- 50L
    positions <- seq(1L, contig_length - read_length + 1L, by = 4L)
    reads <- list()
    for (sample in samples) {
        for (contig in c("d1", "r1", "r2")) {
            reads[[length(reads) + 1L]] <- data.frame(
                contig = contig,
                position = positions,
                sample = sample,
                stringsAsFactors = FALSE
            )
        }
    }
    reads <- do.call(rbind, reads)
    reads <- reads[order(
        match(reads$contig, c("d1", "r1", "r2")),
        reads$position, reads$sample
    ), , drop = FALSE]
    sam <- c(
        "@HD\tVN:1.6\tSO:coordinate",
        sprintf("@SQ\tSN:%s\tLN:%d", c("d1", "r1", "r2"), contig_length),
        sprintf("@RG\tID:%s\tSM:%s", samples, samples),
        sprintf(
            "r_%s_%s_%d\t0\t%s\t%d\t60\t%dM\t*\t0\t0\t%s\t%s\tRG:Z:%s",
            reads$sample, reads$contig, reads$position, reads$contig,
            reads$position, read_length,
            strrep("A", read_length), strrep("I", read_length), reads$sample
        )
    )
    bam_path <- write_test_bam(
        sam, file.path(directory, "multi.bam"), index = TRUE
    )

    list(
        directory = directory,
        config_path = config_path,
        bam_path = bam_path,
        samples = samples,
        contig_length = as.integer(contig_length)
    )
}

## Read the legacy directory-layout BED tables used only as the golden oracle.
## The package no longer ships a directory-layout reader (T7.5); this helper
## keeps the gated BAM-vs-legacy regression self-contained.
read_golden_legacy <- function(root, sample = "Sq54",
                               reference = "cs_t2t_se_cereale",
                               window = 1000L) {
    read_bed <- function(path, value_col) {
        table <- utils::read.delim(
            path, header = FALSE, comment.char = "#",
            colClasses = "character", stringsAsFactors = FALSE,
            check.names = FALSE
        )
        table <- table[, 1:4, drop = FALSE]
        names(table) <- c("chr_id", "start", "stop", value_col)
        table$chr_id <- trimws(table$chr_id)
        for (col in c("start", "stop", value_col)) {
            table[[col]] <- as.numeric(trimws(table[[col]]))
        }
        table
    }

    analysis <- file.path(root, "mamamia_analysis", "arsenal", sample)
    assembly_dir <- file.path(root, "assemblies", reference)
    dcov <- read_bed(
        file.path(analysis, sprintf("%s.w%d.dcov.bed", reference, window)),
        "cov"
    )
    rcov <- read_bed(
        file.path(analysis, sprintf("%s.w%d.rcov.bed", reference, window)),
        "cov"
    )
    dgc <- read_bed(
        file.path(assembly_dir, sprintf("donor.w%d.gc.bed", window)), "gc"
    )
    rgc <- read_bed(
        file.path(assembly_dir, sprintf("recipient.w%d.gc.bed", window)), "gc"
    )
    meta <- utils::read.delim(
        file.path(assembly_dir, "chromosomes.tsv"),
        header = FALSE,
        col.names = c("chr_id", "chr_name", "subgenome"),
        colClasses = "character",
        stringsAsFactors = FALSE
    )
    meta <- meta[meta$chr_id %in% c(dcov$chr_id, rcov$chr_id), , drop = FALSE]
    rownames(meta) <- NULL

    list(dcov = dcov, dgc = dgc, rcov = rcov, rgc = rgc, meta = meta)
}

mamamia_test_data_dir <- function() {
    configured <- Sys.getenv("MAMAMIA_TEST_DATA_DIR", unset = "")
    if (nzchar(configured)) {
        return(normalizePath(configured, winslash = "/", mustWork = FALSE))
    }
    normalizePath(
        testthat::test_path("..", "..", "test_data"),
        winslash = "/",
        mustWork = FALSE
    )
}

skip_if_no_test_data <- function() {
    path <- mamamia_test_data_dir()
    have_environment <- identical(Sys.getenv("MAMAMIA_TEST_DATA"), "1")
    if (!dir.exists(path)) {
        testthat::skip(paste(
            "test_data not found at", path,
            "(set MAMAMIA_TEST_DATA_DIR or MAMAMIA_TEST_DATA=1 on a data host)"
        ))
    }
    if (!have_environment && !dir.exists(file.path(path, "bam_files", "Sq54"))) {
        testthat::skip("real BAM fixtures are not available")
    }
    invisible(path)
}
