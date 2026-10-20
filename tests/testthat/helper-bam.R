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
