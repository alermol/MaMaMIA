#!/usr/bin/env Rscript
## T7.4 real-data release gate.
##
## Runs the BAM-first pipeline end to end against the real Sq54/Sq81/Sq84
## BAMs and the wheat + rye gcrank indexes:
##
##   1. Sq54 alone at 1,000 windows per contig.
##   2. Sq54 + Sq81 + Sq84 with median normalization.
##   3. The same cohort with Sq54 as the reference sample.
##   4. Sq54 at 10,000 windows per contig with max_fit_windows = 50000.
##
## It also resolves the role-free species registry in both directions to prove
## that donor/recipient are assigned at run time and that swapping works.
##
## Usage (from the repository root, after R CMD INSTALL):
##
##   MAMAMIA_TEST_DATA=1 Rscript tests/integration/run-bam-e2e.R
##   Rscript tests/integration/run-bam-e2e.R \
##       --test-data=/srv/mamamia/test_data --cores=8
##
## The script exits 0 only when every check passes. It prints the process peak
## RSS (VmHWM) at the end; fork-based workers are not included, so watch the
## container memory as well during the high-resolution run.

suppressPackageStartupMessages(library(MaMaMIA))

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
    prefix <- paste0("--", name, "=")
    hit <- args[startsWith(args, prefix)]
    if (length(hit) == 0L) {
        return(default)
    }
    substring(hit[[1L]], nchar(prefix) + 1L)
}

has_flag <- function(name) paste0("--", name) %in% args

test_data <- arg_value(
    "test-data",
    Sys.getenv("MAMAMIA_TEST_DATA_DIR", unset = "test_data")
)
if (!dir.exists(test_data)) {
    message("SKIP: test_data not found at ", test_data)
    message("Set MAMAMIA_TEST_DATA_DIR or pass --test-data=<path>.")
    quit(status = 0L)
}
test_data <- normalizePath(test_data, winslash = "/", mustWork = TRUE)

cores <- as.integer(arg_value("cores", "8"))
if (is.na(cores) || cores < 1L) {
    stop("--cores must be a positive integer", call. = FALSE)
}
detected <- parallel::detectCores()
if (!is.na(detected) && detected >= 1L) {
    cores <- min(cores, as.integer(detected))
}

resolution <- as.integer(arg_value("resolution", "1000"))
high_resolution <- as.integer(arg_value("high-resolution", "10000"))
max_fit_windows <- as.numeric(arg_value("max-fit-windows", "50000"))
skip_high_res <- has_flag("skip-high-res")

if (is.na(resolution) || resolution < 1L ||
    is.na(high_resolution) || high_resolution < 1L) {
    stop("resolutions must be positive integers", call. = FALSE)
}
if (is.na(max_fit_windows) || max_fit_windows < 1) {
    stop("--max-fit-windows must be a positive number", call. = FALSE)
}

session <- tempfile("mamamia-e2e-")
dir.create(session, recursive = TRUE)
on.exit(unlink(session, recursive = TRUE, force = TRUE), add = TRUE)

## Checks --------------------------------------------------------------------

checks <- list()
failures <- 0L

check <- function(scenario, label, ok, detail = "") {
    ok <- isTRUE(ok)
    checks[[length(checks) + 1L]] <<- data.frame(
        scenario = scenario,
        check = label,
        result = if (ok) "PASS" else "FAIL",
        detail = as.character(detail),
        stringsAsFactors = FALSE
    )
    cat(sprintf(
        "  %-4s %s%s\n",
        if (ok) "PASS" else "FAIL", label,
        if (nzchar(detail)) paste0(" (", detail, ")") else ""
    ))
    if (!ok) {
        failures <<- failures + 1L
    }
    invisible(ok)
}

peak_rss <- function() {
    status <- tryCatch(
        readLines("/proc/self/status", warn = FALSE),
        error = function(e) character()
    )
    hit <- grep("^VmHWM:", status, value = TRUE)
    if (length(hit) == 0L) {
        return(NA_real_)
    }
    as.numeric(sub("\\s*kB$", "", sub("^VmHWM:\\s*", "", hit[[1L]])))
}

rss_mb <- function() peak_rss() / 1024

format_detail <- function(detail) {
    if (is.null(detail) || length(detail) == 0L) {
        return("")
    }
    fields <- vapply(names(detail), function(name) {
        value <- detail[[name]]
        text <- if (length(value) == 1L) {
            as.character(value)
        } else {
            paste0(length(value), " values")
        }
        paste0(name, "=", text)
    }, character(1))
    paste(fields, collapse = " ")
}

## Configuration -------------------------------------------------------------

read_chromosomes <- function(root) {
    utils::read.delim(
        file.path(
            root, "assemblies", "cs_t2t_se_cereale", "chromosomes.tsv"
        ),
        header = FALSE,
        col.names = c("id", "name", "subgenome"),
        stringsAsFactors = FALSE
    )
}

contig_entries <- function(chromosomes) {
    lapply(seq_len(nrow(chromosomes)), function(i) {
        list(
            id = chromosomes$id[[i]],
            name = chromosomes$name[[i]],
            subgenome = chromosomes$subgenome[[i]]
        )
    })
}

gc_indexes <- function(root) {
    c(
        normalizePath(
            file.path(root, "gcidx", "triticum-aestivum", "index.gcrank"),
            winslash = "/", mustWork = TRUE
        ),
        normalizePath(
            file.path(root, "gcidx", "secale-cereale", "index.gcrank"),
            winslash = "/", mustWork = TRUE
        )
    )
}

defaults <- function() {
    list(
        windows_per_chromosome = resolution,
        allowed_windows = c(100L, 1000L, 10000L),
        mapq = 20L,
        exclude_flags = 2308L,
        max_fit_windows = max_fit_windows,
        gc_mmap = TRUE,
        gc_validate_crc = FALSE,
        normalization_by_side = TRUE
    )
}

server <- function() {
    list(
        max_bams = 20L,
        max_total_bytes = 20 * 1024^3
    )
}

## The full 28-contig wheat + rye assembly, used for the analysis scenarios.
real_config <- function(root) {
    chromosomes <- read_chromosomes(root)
    indexes <- gc_indexes(root)
    pairs <- list(
        OZ284382.1 = "CP155614.1",
        OZ284383.1 = "CP155611.1",
        OZ284384.1 = "CP155608.1",
        OZ284385.1 = "CP155605.1",
        OZ284386.1 = "CP155602.1",
        OZ284387.1 = "CP155599.1",
        OZ284388.1 = "CP155597.1"
    )
    config <- list(
        version = 1L,
        defaults = defaults(),
        server = server(),
        assemblies = list(cs_t2t_se_cereale = list(
            label = "Triticum aestivum CS T2T + Secale cereale",
            description = "T7.4 real-data integration fixture",
            gc_indexes = as.list(indexes),
            contigs = list(
                recipient = contig_entries(
                    chromosomes[chromosomes$subgenome != "R", , drop = FALSE]
                ),
                donor = contig_entries(
                    chromosomes[chromosomes$subgenome == "R", , drop = FALSE]
                )
            ),
            pairs = pairs
        ))
    )
    path <- file.path(session, "assemblies-real.yml")
    writeLines(yaml::as.yaml(config), path)
    readMaMaMIAConfig(path)
}

## The species-only registry: no predefined assemblies, roles assigned at run
## time, both directions resolvable.
species_config <- function(root) {
    chromosomes <- read_chromosomes(root)
    indexes <- gc_indexes(root)
    config <- list(
        version = 1L,
        defaults = defaults(),
        server = server(),
        species = list(
            triticum_aestivum = list(
                label = "Triticum aestivum",
                description = "Wheat A-genome contigs of the CS T2T assembly",
                gc_indexes = list(indexes[[1L]]),
                contigs = contig_entries(
                    chromosomes[chromosomes$subgenome != "R", , drop = FALSE]
                )
            ),
            secale_cereale = list(
                label = "Secale cereale",
                description = "Rye R-genome contigs",
                gc_indexes = list(indexes[[2L]]),
                contigs = contig_entries(
                    chromosomes[chromosomes$subgenome == "R", , drop = FALSE]
                )
            )
        )
    )
    path <- file.path(session, "assemblies-species.yml")
    writeLines(yaml::as.yaml(config), path)
    readMaMaMIAConfig(path)
}

config_path <- file.path(session, "assemblies-real.yml")
real_config(test_data)

## Species registry ----------------------------------------------------------

cat("\n== 0. Role-free species registry ==\n")
species <- species_config(test_data)
forward <- resolveSpeciesAssembly(
    species,
    donor = "secale_cereale", recipient = "triticum_aestivum",
    donor_subgenome = "R", recipient_subgenome = "A"
)
reverse <- resolveSpeciesAssembly(
    species,
    donor = "triticum_aestivum", recipient = "secale_cereale",
    donor_subgenome = "A", recipient_subgenome = "R"
)
check(
    "species", "forward resolution pairs wheat and rye",
    length(forward$assembly$pairs) == 7L,
    sprintf("%d pairs", length(forward$assembly$pairs))
)
check(
    "species", "roles are run-time assignments",
    identical(forward$donor_species, "secale_cereale") &&
        identical(forward$recipient_species, "triticum_aestivum") &&
        identical(reverse$donor_species, "triticum_aestivum") &&
        identical(reverse$recipient_species, "secale_cereale")
)
check(
    "species", "swapping donor and recipient reverses every pair",
    setequal(
        paste(names(forward$assembly$pairs), unname(forward$assembly$pairs)),
        paste(unname(reverse$assembly$pairs), names(reverse$assembly$pairs))
    ) &&
        !identical(forward$key, reverse$key),
    sprintf("%s <-> %s", forward$key, reverse$key)
)

## BAM staging ----------------------------------------------------------------

bam_names <- c("Sq54", "Sq81", "Sq84")
bam_paths <- vapply(bam_names, function(sample) {
    path <- file.path(
        test_data, "bam_files", sample, "cs_t2t_se_cereale.bam"
    )
    if (!file.exists(path)) {
        stop("Missing real-data BAM: ", path, call. = FALSE)
    }
    MaMaMIA:::stageBam(path, session)
}, character(1))
names(bam_paths) <- bam_names

run_analysis <- function(label, bams, normalization,
                         reference_sample = NULL,
                         windows = resolution,
                         max_fit = max_fit_windows) {
    cat("\n==", label, "==\n")
    warnings <- character()
    started <- proc.time()
    result <- withCallingHandlers(
        analyzeBamInput(
            bams,
            config_path,
            windows_per_chromosome = windows,
            reference_sample = reference_sample,
            normalization = normalization,
            max_fit_windows = max_fit,
            cores = cores,
            cache = file.path(
                session,
                paste0("fit-cache-", gsub("[^A-Za-z0-9]+", "-", label), ".rds")
            ),
            progress = function(stage, detail = list()) {
                cat(sprintf(
                    "  [%s]%s\n", stage,
                    if (length(detail) == 0L) {
                        ""
                    } else {
                        paste0(" ", format_detail(detail))
                    }
                ))
            }
        ),
        warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    elapsed <- unname((proc.time() - started)[["elapsed"]])
    cat(sprintf(
        "  elapsed %.1f s, peak RSS %.0f MB\n", elapsed, rss_mb()
    ))
    list(
        result = result,
        elapsed = elapsed,
        warnings = warnings
    )
}

check_warnings <- function(scenario, warnings) {
    expected <- grepl("^Total analysis windows", warnings)
    check(
        scenario, "only documented performance warnings",
        all(expected),
        if (length(warnings) == 0L) {
            "none"
        } else {
            sprintf("%d warning(s)", length(warnings))
        }
    )
}

## 1. Single BAM -------------------------------------------------------------

single <- run_analysis(
    "1. Sq54 single BAM, 1,000 windows", bam_paths[["Sq54"]], "none"
)
r <- single$result
check("single", "result class", inherits(r, "BamAnalysisResult"))
check("single", "sample name", identical(r$samples$sample, "Sq54"))
check(
    "single", "28,000 windows",
    nrow(r$windows$table) == 28000L,
    sprintf("%d", nrow(r$windows$table))
)
check(
    "single", "window counts by side",
    identical(dim(r$counts$donor), c(7000L, 1L)) &&
        identical(dim(r$counts$recipient), c(21000L, 1L)),
    sprintf(
        "%dx%d donor, %dx%d recipient",
        nrow(r$counts$donor), ncol(r$counts$donor),
        nrow(r$counts$recipient), ncol(r$counts$recipient)
    )
)
check(
    "single", "GC is complete and finite",
    nrow(r$gc) == 28000L && all(is.finite(r$gc$gc)) &&
        all(r$gc$gc >= 0 & r$gc$gc <= 1),
    sprintf("%d windows", nrow(r$gc))
)
check(
    "single", "single BAM is not normalized",
    identical(r$normalization$mode, "none") &&
        is.null(r$normalization$reference) &&
        !("dose" %in% names(r$rca$Sq54$data)) &&
        identical(r$isa$Sq54$meta$value_col, "cor.gc")
)
check(
    "single", "segmentation produced introgressions",
    nrow(r$isa$Sq54$segments) > 0L &&
        sum(r$isa$Sq54$segments$putative_introgression) > 0L,
    sprintf(
        "%d segments, %d putative",
        nrow(r$isa$Sq54$segments),
        sum(r$isa$Sq54$segments$putative_introgression)
    )
)
check(
    "single", "provenance records the frozen semantics",
    identical(r$provenance$normalization_mode, "none") &&
        identical(r$provenance$count_semantics, "bamsignals_start_based") &&
        identical(r$provenance$gc_source, "gcrank") &&
        isTRUE(all.equal(
            as.numeric(r$provenance$max_fit_windows), max_fit_windows
        )) &&
        identical(as.integer(r$provenance$n_windows), 28000L) &&
        identical(as.integer(r$provenance$n_samples), 1L)
)
single_bedpe <- file.path(session, "single-selected.bedpe")
suppressMessages(writeCohortToBedpe(r, single_bedpe))
single_table <- utils::read.delim(
    single_bedpe, header = TRUE, sep = "\t", check.names = FALSE
)
check(
    "single", "BEDPE export carries the sample column",
    "sample" %in% names(single_table) &&
        nrow(single_table) > 0L &&
        all(single_table$sample == "Sq54"),
    sprintf("%d rows, %d columns", nrow(single_table), ncol(single_table))
)
check_warnings("single", single$warnings)
single <- NULL
r <- NULL
gc()

## 2. Median cohort ----------------------------------------------------------

cohort_bams <- bam_paths[c("Sq54", "Sq81", "Sq84")]
median_run <- run_analysis(
    "2. Sq54 + Sq81 + Sq84, median normalization",
    cohort_bams, "median"
)
r <- median_run$result
check(
    "median", "three samples in input order",
    identical(r$samples$sample, c("Sq54", "Sq81", "Sq84")),
    paste(r$samples$sample, collapse = ", ")
)
check(
    "median", "cohort count matrices",
    identical(dim(r$counts$donor), c(7000L, 3L)) &&
        identical(dim(r$counts$recipient), c(21000L, 3L)),
    sprintf(
        "%dx%d donor, %dx%d recipient",
        nrow(r$counts$donor), ncol(r$counts$donor),
        nrow(r$counts$recipient), ncol(r$counts$recipient)
    )
)
check(
    "median", "normalization metadata",
    identical(r$normalization$mode, "median") &&
        is.null(r$normalization$reference) &&
        isTRUE(r$normalization$by_side) &&
        identical(r$normalization$pseudocount, 0.5)
)
dose_finite <- vapply(
    names(r$isa),
    function(sample) all(is.finite(r$rca[[sample]]$data$dose)),
    logical(1)
)
check(
    "median", "every sample has a finite dose",
    all(dose_finite),
    paste(names(dose_finite)[!dose_finite], collapse = ", ")
)
value_cols <- vapply(
    names(r$isa),
    function(sample) r$isa[[sample]]$meta$value_col,
    character(1)
)
check(
    "median", "CBS consumes the dose column",
    all(value_cols == "dose"),
    paste(unique(value_cols), collapse = ", ")
)
segment_finite <- vapply(
    names(r$isa),
    function(sample) all(is.finite(r$isa[[sample]]$segments$seg.median)),
    logical(1)
)
check(
    "median", "segment medians are finite",
    all(segment_finite) &&
        sum(r$samples$n_putative_introgressions) > 0L,
    sprintf(
        "%d segments, %d putative",
        sum(r$samples$n_segments),
        sum(r$samples$n_putative_introgressions)
    )
)
median_bedpe <- file.path(session, "median-cohort.bedpe")
suppressMessages(writeCohortToBedpe(r, median_bedpe))
median_table <- utils::read.delim(
    median_bedpe, header = TRUE, sep = "\t", check.names = FALSE
)
check(
    "median", "combined BEDPE identifies multiple samples",
    "sample" %in% names(median_table) &&
        all(unique(median_table$sample) %in% names(r$isa)) &&
        length(unique(median_table$sample)) >= 2L,
    sprintf(
        "%d rows from: %s",
        nrow(median_table),
        paste(sort(unique(median_table$sample)), collapse = ", ")
    )
)
check_warnings("median", median_run$warnings)
median_run <- NULL
r <- NULL
gc()

## 3. Reference cohort -------------------------------------------------------

reference_run <- run_analysis(
    "3. Sq54 as the reference sample",
    cohort_bams, "reference", reference_sample = "Sq54"
)
r <- reference_run$result
check(
    "reference", "reference mode and sample",
    identical(r$normalization$mode, "reference") &&
        identical(r$normalization$reference, "Sq54")
)
check(
    "reference", "reference flags",
    identical(r$samples$reference, c(TRUE, FALSE, FALSE)),
    paste(r$samples$reference, collapse = ", ")
)
check(
    "reference", "reference sample dose is exactly one",
    all(r$rca$Sq54$data$dose == 1)
)
non_reference <- setdiff(names(r$isa), "Sq54")
other_finite <- vapply(
    non_reference,
    function(sample) all(is.finite(r$rca[[sample]]$data$dose)),
    logical(1)
)
check(
    "reference", "non-reference doses are finite",
    all(other_finite),
    paste(non_reference[!other_finite], collapse = ", ")
)
check(
    "reference", "segment medians are finite",
    all(vapply(
        names(r$isa),
        function(sample) all(is.finite(r$isa[[sample]]$segments$seg.median)),
        logical(1)
    )),
    sprintf("%d segments", sum(r$samples$n_segments))
)
check_warnings("reference", reference_run$warnings)
reference_run <- NULL
r <- NULL
gc()

## 4. High-resolution single BAM --------------------------------------------

if (skip_high_res) {
    message("\n--skip-high-res: high-resolution scenario skipped")
} else {
    high <- run_analysis(
        sprintf(
            "4. Sq54 at %s windows, max_fit_windows = %s",
            format(high_resolution, scientific = FALSE),
            format(max_fit_windows, scientific = FALSE)
        ),
        bam_paths[["Sq54"]], "none",
        windows = high_resolution, max_fit = max_fit_windows
    )
    r <- high$result
    expected_windows <- 28L * high_resolution
    check(
        "high-res", "exact high-resolution window count",
        nrow(r$windows$table) == expected_windows,
        sprintf("%d", nrow(r$windows$table))
    )
    check(
        "high-res", "high-resolution count matrices",
        ncol(r$counts$donor) == 1L && ncol(r$counts$recipient) == 1L &&
            nrow(r$counts$donor) == 7L * high_resolution &&
            nrow(r$counts$recipient) == 21L * high_resolution,
        sprintf(
            "%d donor, %d recipient",
            nrow(r$counts$donor), nrow(r$counts$recipient)
        )
    )
    sampling <- r$rca$Sq54$fit_sampling
    check(
        "high-res", "fit windows are subsampled deterministically",
        isTRUE(sampling$sampled) &&
            sampling$n_ideal > max_fit_windows &&
            sampling$n_fitted <= max_fit_windows &&
            !is.na(sampling$seed),
        sprintf(
            "%d ideal, %d fitted, seed %s",
            sampling$n_ideal, sampling$n_fitted, sampling$seed
        )
    )
    check(
        "high-res", "fit cap is recorded in provenance",
        isTRUE(all.equal(
            as.numeric(sampling$max_fit_windows), max_fit_windows
        )) &&
            isTRUE(all.equal(
                as.numeric(r$provenance$max_fit_windows), max_fit_windows
            )) &&
            identical(
                as.integer(r$provenance$n_windows), as.integer(expected_windows)
            )
    )
    check(
        "high-res", "segmentation still succeeds",
        nrow(r$isa$Sq54$segments) > 0L &&
            all(is.finite(r$isa$Sq54$segments$seg.median)),
        sprintf("%d segments", nrow(r$isa$Sq54$segments))
    )
    check_warnings("high-res", high$warnings)
    high <- NULL
    r <- NULL
    gc()
}

## Summary -------------------------------------------------------------------

table <- do.call(rbind, checks)
cat("\n== checks ==\n")
print(table[, c("scenario", "check", "result", "detail")], row.names = FALSE)
cat(sprintf(
    "\n%d checks, %d failure(s); peak RSS %.0f MB\n",
    nrow(table), failures, rss_mb()
))
if (failures > 0L) {
    quit(status = 1L)
}
cat("REAL-DATA INTEGRATION GATE PASSED\n")
