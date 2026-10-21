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
    path <- file.path(fixture$directory, filename)
    write_test_bam(variable_gc_sam(sample, step = step), path, index = FALSE)
    path
}

variable_gc_fixture <- function(sample = "synthetic", windows = 100L) {
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
            ">d1", variable_gc_sequence(20000L, 0L),
            ">r1", variable_gc_sequence(20000L, 8L)
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
            allowed_windows = c(100L, 1000L, 10000L),
            mapq = 20L,
            exclude_flags = 2308L,
            max_fit_windows = 1000L,
            gc_mmap = TRUE,
            gc_validate_crc = FALSE,
            normalization_by_side = TRUE
        ),
        server = list(max_bams = 4L, max_total_bytes = 20 * 1024^2),
        assemblies = list(
            variable = list(
                label = "Variable-GC synthetic assembly",
                description = "Two 20 kb contigs with variable GC",
                gc_indexes = list(file.path(directory, "variable.gcrank")),
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
    )
    config_path <- file.path(directory, "assemblies.yml")
    writeLines(yaml::as.yaml(config), config_path)

    fixture <- list(
        directory = directory,
        config_path = config_path,
        sample = sample
    )
    fixture$bam_path <- add_variable_gc_bam(fixture, sample)
    fixture
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
