#!/usr/bin/env Rscript

## Validate a production MaMaMIA species/assembly registry before it is
## mounted into the ShinyProxy container.
##
##   Rscript tools/validate-config.R /srv/mamamia/config/assemblies.yml
##
## This checks everything the application checks at startup (YAML shape,
## species groups and defaults, gcrank format version and contig tables,
## target-contig coverage, pairs, windows) and adds deployment checks: index
## payload files and the optional server-side BAM root. Exit status 0 means the
## registry is usable.

usage <- "Usage: Rscript tools/validate-config.R <assemblies.yml>"

args <- commandArgs(trailingOnly = TRUE)
args <- args[nzchar(args)]
if (length(args) != 1L) {
    message(usage)
    quit(status = 2L)
}
config_path <- args[[1L]]
options(width = 200L)

if (!requireNamespace("MaMaMIA", quietly = TRUE)) {
    if (requireNamespace("pkgload", quietly = TRUE)) {
        pkgload::load_all(".", quiet = TRUE)
    } else {
        message(
            "MaMaMIA is not installed. Run `uvr sync` and install the package, ",
            "or run this script from the source tree with pkgload available."
        )
        quit(status = 2L)
    }
}

fail <- function(fmt, ...) {
    message("ERROR: ", sprintf(fmt, ...))
    quit(status = 1L)
}

config <- tryCatch(
    MaMaMIA::readMaMaMIAConfig(config_path),
    error = function(e) e
)
if (inherits(config, "condition")) {
    fail("%s", conditionMessage(config))
}

cat("Registry:", normalizePath(config_path, mustWork = TRUE), "\n")
config_dir <- dirname(normalizePath(config_path, mustWork = TRUE))
cat("Config version:", config$version, "\n")
cat(
    "Defaults: windows/contig =", config$defaults$windows_per_chromosome,
    "| allowed =", paste(config$defaults$allowed_windows, collapse = ", "),
    "| mapq =", config$defaults$mapq,
    "| max_fit_windows =", config$defaults$max_fit_windows, "\n"
)
cat("Server limits: max_bams =", config$server$max_bams,
    "| max_total_bytes =", format(config$server$max_total_bytes, scientific = FALSE),
    "\n\n")

summary_table <- MaMaMIA::listMaMaMIAAssemblies(config)
if (nrow(summary_table) > 0L) {
    print(summary_table, row.names = FALSE)
    cat("\n")
}

species_table <- MaMaMIA::listMaMaMIASpecies(config)
if (nrow(species_table) > 0L) {
    cat("Species:\n")
    print(species_table, row.names = FALSE)
    cat("\n")
}

check_index <- function(index_path, context) {
    required <- c(
        "contigs.tsv", "manifest.json", "genome.bv", "genome.rb"
    )
    missing <- required[!file.exists(file.path(index_path, required))]
    if (length(missing) > 0L) {
        fail(
            "%s: gc index %s is missing %s",
            context, index_path, paste(missing, collapse = ", ")
        )
    }
    size <- sum(file.info(file.path(index_path, required))$size)
    cat(
        "  gc index:", index_path,
        sprintf("(%.1f MiB)", size / 1024^2), "\n"
    )
}

for (key in names(config$species)) {
    species <- config$species[[key]]
    cat("Species:", key, "-", species$label, "\n")
    cat(
        "  contigs:", nrow(species$contigs),
        "| subgenomes:",
        paste(unique(species$contigs$subgenome), collapse = ", "),
        "| default:",
        if (is.null(species$default_subgenome)) {
            "none"
        } else {
            species$default_subgenome
        },
        "\n"
    )
    for (index_path in species$gc_indexes) {
        check_index(index_path, sprintf("species '%s'", key))
    }
    if (anyDuplicated(species$contigs$id) != 0L) {
        fail("species '%s' configures a contig more than once", key)
    }
}

for (key in names(config$assemblies)) {
    assembly <- config$assemblies[[key]]
    donor <- assembly$contigs$donor
    recipient <- assembly$contigs$recipient
    cat("Assembly:", key, "-", assembly$label, "\n")
    cat(
        "  contigs:", nrow(donor), "donor /", nrow(recipient), "recipient",
        "| pairs:", length(assembly$pairs), "\n"
    )
    cat(
        "  target length:",
        format(sum(donor$length) + sum(recipient$length), scientific = FALSE),
        "bp\n"
    )

    for (index_path in assembly$gc_indexes) {
        check_index(index_path, sprintf("assembly '%s'", key))
    }

    duplicate <- anyDuplicated(c(donor$id, recipient$id))
    if (duplicate != 0L) {
        fail("assembly '%s' uses a contig id on both sides", key)
    }
    if (any(donor$subgenome %in% recipient$subgenome)) {
        fail("assembly '%s' shares a subgenome label between sides", key)
    }
}

bam_root <- config$server$bam_root
if (is.null(bam_root)) {
    cat("BAM root: not configured (uploads only)\n")
} else if (!dir.exists(bam_root)) {
    fail("server.bam_root does not exist or is not a directory: %s", bam_root)
} else {
    bams <- list.files(
        bam_root, pattern = "[.]bam$", recursive = TRUE, full.names = TRUE
    )
    cat("BAM root:", bam_root, "|", length(bams), "BAM file(s)\n")
    if (length(bams) == 0L) {
        warning(
            "server.bam_root contains no BAM files; the picker will be empty",
            call. = FALSE
        )
    }
    if (file.access(bam_root, mode = 2L) == 0L) {
        warning(
            "server.bam_root is writable; mount it read-only in production",
            call. = FALSE
        )
    }
}

cat("\nRegistry is valid.\n")
