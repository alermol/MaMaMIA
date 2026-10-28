#!/usr/bin/env Rscript

## Validate a production MaMaMIA species/assembly registry before it is
## mounted into the ShinyProxy container.
##
##   Rscript tools/validate-config.R /srv/mamamia/config/assemblies.yml
##
## This checks everything the application checks at startup (YAML shape,
## species groups and defaults, gcrank format version and contig tables,
## target-contig coverage, pairs, windows) and adds deployment checks for the
## index payload files. Exit status 0 means the registry is usable.

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
cat("Config version:", config$version, "\n")
cat(
    "Defaults: windows/contig =", config$defaults$windows_per_chromosome,
    "| allowed =", paste(config$defaults$allowed_windows, collapse = ", "),
    "| mapq =", config$defaults$mapq,
    "| metric =", config$defaults$metric,
    "| score =", config$defaults$score, "\n"
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
        "| components:",
        paste(unique(species$contigs$component), collapse = ", "),
        "| default:",
        if (is.null(species$default_component)) {
            "none"
        } else {
            species$default_component
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

## The registry is species-only. Check that every ordered species pair the
## application can select resolves to an assembly with at least one paired
## chromosome group, and report the contigs that would be ignored.
species_keys <- names(config$species)
checked <- 0L
for (recipient_key in species_keys) {
    for (donor_key in species_keys) {
        if (identical(recipient_key, donor_key)) {
            next
        }
        pair <- tryCatch(
            MaMaMIA::resolveSpeciesAssembly(
                config, donor = donor_key, recipient = recipient_key
            ),
            error = function(e) e
        )
        if (inherits(pair, "condition")) {
            cat(
                "Pair", recipient_key, "<-", donor_key,
                ": not usable (", conditionMessage(pair), ")\n"
            )
            next
        }
        assembly <- pair$assembly
        checked <- checked + 1L
        cat(
            "Pair", recipient_key, "<-", donor_key, ":",
            nrow(assembly$contigs$donor), "donor /",
            nrow(assembly$contigs$recipient), "recipient contigs,",
            length(assembly$pairs), "groups,",
            nrow(pair$unpaired), "unpaired\n"
        )
        if (length(assembly$pairs) == 0L) {
            fail(
                "pair '%s' + '%s' has no paired chromosome group",
                recipient_key, donor_key
            )
        }
    }
}
if (checked == 0L) {
    fail("no species pair resolves to an assembly")
}

cat("\nRegistry is valid.\n")
