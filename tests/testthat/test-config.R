## T0.3: YAML registry loading and validation.

write_index_fixture <- function(dir,
                                ids,
                                lengths,
                                format_version = "0.1.0",
                                total_bases = NULL,
                                num_contigs = NULL,
                                offsets = NULL) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    if (is.null(offsets)) {
        offsets <- cumsum(lengths) - lengths
    }
    utils::write.table(
        data.frame(ids, lengths, offsets),
        file.path(dir, "contigs.tsv"),
        sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    if (is.null(total_bases)) {
        total_bases <- sum(lengths)
    }
    if (is.null(num_contigs)) {
        num_contigs <- length(ids)
    }
    manifest <- list(
        builder_version = "0.1.0",
        created = "2026-01-01T00:00:00Z",
        files = list(),
        format_version = format_version,
        index = list(
            basic_block_bits = 512L, bit_order = "msb0",
            endianness = "little", remove_ambiguous_supported = TRUE,
            superblock_bits = 4096L, total_gc = 0, total_valid = 0
        ),
        semantics = list(base_encoding = "test"),
        source = list(
            file_bytes = total_bases, num_contigs = num_contigs,
            path = "tiny.fa", sha256 = "test", total_bases = total_bases
        )
    )
    writeLines(
        jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE),
        file.path(dir, "manifest.json")
    )
    dir
}

contig_entry <- function(id, name, subgenome) {
    list(id = id, name = name, subgenome = subgenome)
}

assembly_fixture <- function(gc_indexes,
                             donor,
                             recipient,
                             pairs = NULL,
                             label = "Test assembly",
                             description = "Synthetic test assembly") {
    out <- list(
        label = label,
        description = description,
        gc_indexes = as.list(gc_indexes),
        contigs = list(recipient = recipient, donor = donor)
    )
    if (!is.null(pairs)) {
        out$pairs <- as.list(pairs)
    }
    out
}

write_config_fixture <- function(dir,
                                 assemblies,
                                 defaults = list(),
                                 server = list(),
                                 species = NULL,
                                 version = 1L) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    config <- list()
    if (!is.null(version)) {
        config$version <- version
    }
    if (length(defaults) > 0L) {
        config$defaults <- defaults
    }
    if (length(server) > 0L) {
        config$server <- server
    }
    if (!is.null(species)) {
        config$species <- species
    }
    if (!is.null(assemblies)) {
        config$assemblies <- assemblies
    }
    path <- file.path(dir, "assemblies.yml")
    yaml::write_yaml(config, path)
    path
}

tiny_assembly <- function(indexes = "idx",
                           donor = list(contig_entry("d1", "1D", "D")),
                           recipient = list(contig_entry("r1", "1R", "R")),
                           pairs = c(d1 = "r1"),
                           label = "Tiny assembly") {
    assembly_fixture(
        gc_indexes = indexes,
        donor = donor,
        recipient = recipient,
        pairs = pairs,
        label = label
    )
}

species_fixture <- function(gc_indexes, contigs, label = "Test species",
                            description = "Synthetic species",
                            default_subgenome = NULL) {
    out <- list(
        label = label,
        description = description,
        gc_indexes = as.list(gc_indexes),
        contigs = contigs
    )
    if (!is.null(default_subgenome)) {
        out$default_subgenome <- default_subgenome
    }
    out
}

test_that("readMaMaMIAConfig() normalizes a valid registry", {
    root <- tempfile("mamamia-config-valid-")
    write_index_fixture(file.path(root, "idx"), c("d1", "r1"), c(20000, 20000))
    path <- write_config_fixture(
        root,
        assemblies = list(tiny = tiny_assembly()),
        defaults = list(windows_per_chromosome = 100L),
        server = list(max_bams = 3, max_total_bytes = 1024)
    )

    config <- readMaMaMIAConfig(path)

    expect_type(config, "list")
    expect_named(config, c(
        "version", "defaults", "species", "assemblies", "server"
    ))
    expect_identical(config$version, 1L)
    expect_identical(config$defaults$windows_per_chromosome, 100L)
    expect_identical(config$defaults$allowed_windows, c(100L, 1000L, 10000L))
    expect_identical(config$defaults$mapq, 20L)
    expect_identical(config$defaults$exclude_flags, 2308L)
    expect_identical(config$defaults$max_fit_windows, 50000L)
    expect_true(config$defaults$gc_mmap)
    expect_false(config$defaults$gc_validate_crc)
    expect_true(config$defaults$normalization_by_side)
    expect_identical(config$server$max_bams, 3L)
    expect_identical(config$server$max_total_bytes, 1024)
    expect_false("bam_root" %in% names(config$server))

    assembly <- config$assemblies$tiny
    expect_named(assembly, c(
        "key", "label", "description", "gc_indexes", "contigs", "pairs"
    ))
    expect_identical(assembly$key, "tiny")
    expect_identical(assembly$label, "Tiny assembly")
    expect_identical(assembly$description, "Synthetic test assembly")
    expect_identical(
        assembly$gc_indexes,
        normalizePath(file.path(root, "idx"), winslash = "/", mustWork = TRUE)
    )
    expect_identical(assembly$pairs, c(d1 = "r1"))
    expect_named(assembly$contigs, c("recipient", "donor"))
    expect_named(
        assembly$contigs$donor,
        c("id", "name", "subgenome", "length", "index_path")
    )
    expect_identical(assembly$contigs$donor$id, "d1")
    expect_identical(assembly$contigs$donor$length, 20000)
    expect_identical(
        assembly$contigs$donor$index_path,
        normalizePath(file.path(root, "idx"), winslash = "/", mustWork = TRUE)
    )
    expect_true(!is.null(attr(config, "config_path")))
})

test_that("readMaMaMIAConfig() resolves relative index paths", {
    root <- tempfile("mamamia-config-relative-")
    write_index_fixture(file.path(root, "indexes", "a"), "d1", 30000)
    write_index_fixture(file.path(root, "indexes", "b"), "r1", 30000)
    path <- write_config_fixture(
        root,
        assemblies = list(tiny = tiny_assembly(indexes = c(
            "indexes/a", "indexes/b"
        )))
    )

    config <- readMaMaMIAConfig(path)

    expect_identical(
        config$assemblies$tiny$gc_indexes,
        normalizePath(
            file.path(root, "indexes", c("a", "b")),
            winslash = "/", mustWork = TRUE
        )
    )
})

test_that("listMaMaMIAAssemblies() returns one row per assembly", {
    root <- tempfile("mamamia-config-list-")
    write_index_fixture(file.path(root, "idx"), c("d1", "r1"), c(20000, 20000))
    path <- write_config_fixture(root, assemblies = list(
        tiny = tiny_assembly(),
        tiny2 = tiny_assembly(label = "Second")
    ))

    config <- readMaMaMIAConfig(path)
    table <- listMaMaMIAAssemblies(config)

    expect_s3_class(table, "data.frame")
    expect_named(table, c(
        "key", "label", "description", "n_donor", "n_recipient",
        "n_pairs", "n_gc_indexes"
    ))
    expect_identical(table$key, c("tiny", "tiny2"))
    expect_identical(table$label, c("Tiny assembly", "Second"))
    expect_identical(table$n_donor, c(1L, 1L))
    expect_identical(table$n_recipient, c(1L, 1L))
    expect_identical(table$n_pairs, c(1L, 1L))
    expect_identical(table$n_gc_indexes, c(1L, 1L))
    expect_identical(listMaMaMIAAssemblies(path)$key, c("tiny", "tiny2"))
})

test_that("config defaults are validated", {
    root <- tempfile("mamamia-config-defaults-")
    write_index_fixture(file.path(root, "idx"), c("d1", "r1"), c(20000, 20000))

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "mapq"),
            assemblies = list(tiny = tiny_assembly()),
            defaults = list(mapq = 256)
        )),
        "defaults\\.mapq"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "mapq-negative"),
            assemblies = list(tiny = tiny_assembly()),
            defaults = list(mapq = -1)
        )),
        "defaults\\.mapq"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "flags"),
            assemblies = list(tiny = tiny_assembly()),
            defaults = list(exclude_flags = 4096)
        )),
        "defaults\\.exclude_flags"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "windows"),
            assemblies = list(tiny = tiny_assembly()),
            defaults = list(windows_per_chromosome = 500)
        )),
        "windows_per_chromosome"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "allowed-duplicate"),
            assemblies = list(tiny = tiny_assembly()),
            defaults = list(allowed_windows = c(100, 100, 1000))
        )),
        "must not contain duplicates"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "allowed-zero"),
            assemblies = list(tiny = tiny_assembly()),
            defaults = list(allowed_windows = c(0, 1000))
        )),
        "positive integers"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "fit"),
            assemblies = list(tiny = tiny_assembly()),
            defaults = list(max_fit_windows = 0)
        )),
        "max_fit_windows"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "logical"),
            assemblies = list(tiny = tiny_assembly()),
            defaults = list(gc_mmap = "yes")
        )),
        "gc_mmap"
    )
})

test_that("assembly structure is validated", {
    root <- tempfile("mamamia-config-structure-")
    index <- write_index_fixture(file.path(root, "idx"), c("d1", "r1"), c(20000, 20000))

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-indexes"),
            assemblies = list(tiny = local({
                x <- tiny_assembly(indexes = index)
                x$gc_indexes <- NULL
                x
            }))
        )),
        "gc_indexes must not be empty"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "duplicate-index"),
            assemblies = list(tiny = tiny_assembly(indexes = c(index, index)))
        )),
        "lists the same index twice"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-label"),
            assemblies = list(tiny = local({
                x <- tiny_assembly(indexes = index)
                x$label <- NULL
                x
            }))
        )),
        "assemblies\\.tiny\\.label"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "empty-label"),
            assemblies = list(tiny = local({
                x <- tiny_assembly(indexes = index)
                x$label <- "  "
                x
            }))
        )),
        "assemblies\\.tiny\\.label"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-contigs"),
            assemblies = list(tiny = local({
                x <- tiny_assembly(indexes = index)
                x$contigs <- NULL
                x
            }))
        )),
        "contigs must contain"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-recipient"),
            assemblies = list(tiny = local({
                x <- tiny_assembly(indexes = index)
                x$contigs$recipient <- NULL
                x
            }))
        )),
        "at least one recipient"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-donor"),
            assemblies = list(tiny = local({
                x <- tiny_assembly(indexes = index)
                x$contigs$donor <- NULL
                x
            }))
        )),
        "at least one donor"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "unknown-role"),
            assemblies = list(tiny = local({
                x <- tiny_assembly(indexes = index)
                x$contigs$other <- list(contig_entry("x1", "1X", "X"))
                x
            }))
        )),
        "unknown contigs role"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-id"),
            assemblies = list(tiny = tiny_assembly(
                indexes = index,
                donor = list(list(name = "1D", subgenome = "D"))
            ))
        )),
        "contigs\\.donor\\[\\[1\\]\\]\\.id"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-subgenome"),
            assemblies = list(tiny = tiny_assembly(
                indexes = index,
                donor = list(list(id = "d1", name = "1D"))
            ))
        )),
        "contigs\\.donor\\[\\[1\\]\\]\\.subgenome"
    )

    duplicate_names <- list(version = 1L, assemblies = list(
        tiny = tiny_assembly(indexes = index),
        tiny = tiny_assembly(indexes = index)
    ))
    expect_error(
        validateMaMaMIAConfig(duplicate_names),
        "Assembly keys must be unique"
    )
    expect_error(
        validateMaMaMIAConfig(list(version = 1L, assemblies = list(tiny = "oops"))),
        "must be a mapping"
    )
})

test_that("contigs.tsv content is validated", {
    root <- tempfile("mamamia-config-contigs-tsv-")

    bad_offsets <- file.path(root, "bad-offsets")
    write_index_fixture(
        bad_offsets, c("d1", "r1"), c(20000, 20000),
        offsets = c(0, 1)
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "bad-offsets-config"),
            assemblies = list(tiny = tiny_assembly(indexes = bad_offsets))
        )),
        "not the cumulative lengths"
    )

    duplicate_ids <- file.path(root, "duplicate-ids")
    write_index_fixture(
        duplicate_ids, c("d1", "d1"), c(20000, 20000)
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "duplicate-ids-config"),
            assemblies = list(tiny = tiny_assembly(indexes = duplicate_ids))
        )),
        "Duplicate contig"
    )
})

test_that("large numeric limits are parsed without integer overflow", {
    root <- tempfile("mamamia-config-large-")
    write_index_fixture(file.path(root, "idx"), c("d1", "r1"), c(20000, 20000))
    path <- write_config_fixture(
        root,
        assemblies = list(tiny = tiny_assembly()),
        server = list(max_total_bytes = 1024)
    )
    lines <- readLines(path, warn = FALSE)
    lines <- sub(
        "max_total_bytes:.*",
        "max_total_bytes: 5368709120",
        lines
    )
    writeLines(lines, path)

    config <- readMaMaMIAConfig(path)

    expect_identical(config$server$max_total_bytes, 5368709120)
})

test_that("config version and top-level fields are required", {
    root <- tempfile("mamamia-config-top-")
    write_index_fixture(file.path(root, "idx"), c("d1", "r1"), c(20000, 20000))

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-version"),
            assemblies = list(tiny = tiny_assembly()),
            version = NULL
        )),
        "`version` is required"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "bad-version"),
            assemblies = list(tiny = tiny_assembly()),
            version = 2L
        )),
        "Unsupported MaMaMIA config version"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-assemblies"),
            assemblies = list()
        )),
        "at least one `species` or one `assemblies`"
    )
})

test_that("missing and malformed GC indexes are rejected", {
    root <- tempfile("mamamia-config-index-")
    write_index_fixture(file.path(root, "idx"), c("d1", "r1"), c(20000, 20000))

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "missing"),
            assemblies = list(tiny = tiny_assembly(indexes = "nope"))
        )),
        "does not exist"
    )

    missing_contigs <- file.path(root, "idx-no-contigs")
    write_index_fixture(missing_contigs, "d1", 20000)
    unlink(file.path(missing_contigs, "contigs.tsv"))
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-contigs"),
            assemblies = list(tiny = tiny_assembly(indexes = missing_contigs))
        )),
        "has no contigs.tsv"
    )

    missing_manifest <- file.path(root, "idx-no-manifest")
    write_index_fixture(missing_manifest, "d1", 20000)
    unlink(file.path(missing_manifest, "manifest.json"))
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "no-manifest"),
            assemblies = list(tiny = tiny_assembly(indexes = missing_manifest))
        )),
        "has no manifest.json"
    )

    bad_format <- file.path(root, "idx-bad-format")
    write_index_fixture(bad_format, "d1", 20000, format_version = "9.9.9")
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "bad-format"),
            assemblies = list(tiny = tiny_assembly(indexes = bad_format))
        )),
        "unsupported gcrank format_version"
    )

    bad_total <- file.path(root, "idx-bad-total")
    write_index_fixture(bad_total, "d1", 20000, total_bases = 1)
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "bad-total"),
            assemblies = list(tiny = tiny_assembly(indexes = bad_total))
        )),
        "source.total_bases"
    )

    bad_contigs <- file.path(root, "idx-bad-contig-count")
    write_index_fixture(bad_contigs, "d1", 20000, num_contigs = 2)
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "bad-contig-count"),
            assemblies = list(tiny = tiny_assembly(indexes = bad_contigs))
        )),
        "source.num_contigs"
    )
})

test_that("contigs must be present exactly once and roles must be disjoint", {
    root <- tempfile("mamamia-config-contigs-")
    index_a <- file.path(root, "a")
    index_b <- file.path(root, "b")
    only_d1 <- file.path(root, "only-d1")
    write_index_fixture(index_a, c("d1", "r1"), c(20000, 20000))
    write_index_fixture(index_b, c("d1", "r2"), c(20000, 20000))

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "duplicate-index"),
            assemblies = list(tiny = tiny_assembly(indexes = c(index_a, index_b)))
        )),
        "appears in more than one gc_index"
    )

    write_index_fixture(only_d1, "d1", 20000)
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "missing-contig"),
            assemblies = list(tiny = tiny_assembly(indexes = only_d1))
        )),
        "not found in any configured gc_index"
    )

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "shared-subgenome"),
            assemblies = list(tiny = tiny_assembly(
                indexes = index_a,
                recipient = list(contig_entry("r1", "1R", "D"))
            ))
        )),
        "subgenome labels must be disjoint"
    )

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "shared-id"),
            assemblies = list(tiny = tiny_assembly(
                indexes = index_a,
                donor = list(contig_entry("r1", "1D", "D"))
            ))
        )),
        "configured in both roles"
    )
})

test_that("contig lengths are checked against resolution and 32-bit limits", {
    root <- tempfile("mamamia-config-lengths-")
    short_index <- file.path(root, "short")
    long_index <- file.path(root, "long")

    write_index_fixture(short_index, c("d1", "r1"), c(5000, 20000))
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "short-config"),
            assemblies = list(tiny = tiny_assembly(indexes = short_index))
        )),
        "shorter than the largest allowed window"
    )

    write_index_fixture(
        long_index,
        c("d1", "r1"),
        c(2^31, 20000)
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "long-config"),
            assemblies = list(tiny = tiny_assembly(indexes = long_index))
        )),
        "32-bit positions"
    )
})

test_that("pairs are optional but must be one-to-one and known", {
    root <- tempfile("mamamia-config-pairs-")
    index <- file.path(root, "idx")
    write_index_fixture(
        index,
        c("d1", "d2", "r1"),
        c(20000, 20000, 20000)
    )
    donor <- list(
        contig_entry("d1", "1D", "D"),
        contig_entry("d2", "2D", "D")
    )

    no_pairs <- readMaMaMIAConfig(write_config_fixture(
        file.path(root, "none"),
        assemblies = list(tiny = assembly_fixture(
            index, donor = donor,
            recipient = list(contig_entry("r1", "1R", "R")),
            pairs = NULL
        ))
    ))
    expect_identical(
        no_pairs$assemblies$tiny$pairs,
        stats::setNames(character(0), character(0))
    )

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "unknown-donor"),
            assemblies = list(tiny = assembly_fixture(
                index, donor = donor,
                recipient = list(contig_entry("r1", "1R", "R")),
                pairs = c(x = "r1")
            ))
        )),
        "unknown donor"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "unknown-recipient"),
            assemblies = list(tiny = assembly_fixture(
                index, donor = donor,
                recipient = list(contig_entry("r1", "1R", "R")),
                pairs = c(d1 = "x")
            ))
        )),
        "unknown recipient"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "many-to-one"),
            assemblies = list(tiny = assembly_fixture(
                index, donor = donor,
                recipient = list(contig_entry("r1", "1R", "R")),
                pairs = c(d1 = "r1", d2 = "r1")
            ))
        )),
        "not one-to-one"
    )
})

test_that("server limits are validated", {
    root <- tempfile("mamamia-config-server-")
    write_index_fixture(file.path(root, "idx"), c("d1", "r1"), c(20000, 20000))

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "bams"),
            assemblies = list(tiny = tiny_assembly()),
            server = list(max_bams = 0)
        )),
        "server\\.max_bams"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "bytes"),
            assemblies = list(tiny = tiny_assembly()),
            server = list(max_total_bytes = 0)
        )),
        "server\\.max_total_bytes"
    )
    ## Server-side BAM directories were removed; BAMs are uploaded by users.
    ## An old registry key is ignored rather than rejected.
    write_index_fixture(
        file.path(root, "legacy", "idx"), c("d1", "r1"), c(20000, 20000)
    )
    legacy <- readMaMaMIAConfig(write_config_fixture(
        file.path(root, "legacy"),
        assemblies = list(tiny = tiny_assembly()),
        server = list(bam_root = "/srv/mamamia/bams")
    ))
    expect_false("bam_root" %in% names(legacy$server))
})

test_that("duplicate YAML keys fail during parsing", {
    root <- tempfile("mamamia-config-duplicate-")
    dir.create(root, recursive = TRUE)
    path <- file.path(root, "assemblies.yml")
    writeLines(c(
        "version: 1",
        "version: 2",
        "assemblies:"
    ), path)

    expect_error(readMaMaMIAConfig(path), "Duplicate map key")
})

test_that("explicit chromosome groups override the name", {
    root <- tempfile("mamamia-config-species-group-")
    index <- file.path(root, "idx")
    write_index_fixture(index, c("r1", "d1"), rep(20000, 2))
    recipient <- species_fixture(
        index,
        list(list(id = "r1", name = "chrOne", subgenome = "A", group = "1")),
        label = "Recipient"
    )
    donor <- species_fixture(
        index,
        list(list(id = "d1", name = "scaffoldA", subgenome = "R", group = "1")),
        label = "Donor"
    )
    path <- write_config_fixture(
        root,
        assemblies = NULL,
        species = list(recipient = recipient, donor = donor)
    )
    config <- readMaMaMIAConfig(path)
    expect_identical(
        config$species$donor$contigs$group, "1"
    )
    resolved <- resolveSpeciesAssembly(config, "donor", "recipient")
    expect_identical(resolved$assembly$pairs, c(d1 = "r1"))
})

test_that("a missing config file fails clearly", {
    expect_error(
        readMaMaMIAConfig(file.path(tempfile(), "missing.yml")),
        "does not exist"
    )
})

test_that("a species-only registry resolves a donor/recipient pair", {
    root <- tempfile("mamamia-config-species-")
    write_index_fixture(
        file.path(root, "idx"),
        c("w1", "w2", "b1", "r1", "r2"),
        rep(20000, 5)
    )
    wheat <- species_fixture(
        "idx",
        list(
            contig_entry("w1", "1A", "A"),
            contig_entry("w2", "2A", "A"),
            contig_entry("b1", "1B", "B")
        ),
        label = "Wheat",
        default_subgenome = "A"
    )
    rye <- species_fixture(
        "idx",
        list(
            contig_entry("r1", "1R", "R"),
            contig_entry("r2", "2R", "R")
        ),
        label = "Rye",
        default_subgenome = "R"
    )
    path <- write_config_fixture(
        root, assemblies = NULL, species = list(wheat = wheat, rye = rye)
    )

    config <- readMaMaMIAConfig(path)
    expect_identical(names(config$assemblies), character(0))
    expect_named(
        config$species$wheat,
        c(
            "key", "label", "description", "default_subgenome",
            "gc_indexes", "contigs"
        )
    )
    expect_named(
        config$species$wheat$contigs,
        c("id", "name", "subgenome", "group", "length", "index_path")
    )
    expect_identical(config$species$wheat$contigs$group, c("1", "2", "1"))

    table <- listMaMaMIASpecies(config)
    expect_named(table, c(
        "key", "label", "description", "subgenomes",
        "default_subgenome", "n_contigs", "n_gc_indexes"
    ))
    expect_identical(table$key, c("wheat", "rye"))
    expect_identical(table$subgenomes, c("A, B", "R"))
    expect_identical(table$default_subgenome, c("A", "R"))
    expect_identical(table$n_contigs, c(3L, 2L))

    resolved <- resolveSpeciesAssembly(config, "rye", "wheat")
    expect_identical(resolved$key, "wheat__A__rye__R")
    expect_identical(resolved$assembly$pairs, c(r1 = "w1", r2 = "w2"))
    expect_identical(resolved$assembly$contigs$donor$id, c("r1", "r2"))
    expect_identical(resolved$assembly$contigs$recipient$id, c("w1", "w2"))
    expect_identical(nrow(resolved$unpaired), 0L)
    expect_identical(
        names(resolved$config$assemblies), "wheat__A__rye__R"
    )
})

test_that("species assemblies can be presets or use explicit pairs", {
    root <- tempfile("mamamia-config-species-assembly-")
    write_index_fixture(
        file.path(root, "idx"),
        c("w1", "w2", "r1", "r2"),
        rep(20000, 4)
    )
    wheat <- species_fixture(
        "idx",
        list(
            contig_entry("w1", "1A", "A"),
            contig_entry("w2", "2A", "A")
        ),
        label = "Wheat",
        default_subgenome = "A"
    )
    rye <- species_fixture(
        "idx",
        list(
            contig_entry("r1", "1R", "R"),
            contig_entry("r2", "2R", "R")
        ),
        label = "Rye",
        default_subgenome = "R"
    )
    path <- write_config_fixture(
        root,
        species = list(wheat = wheat, rye = rye),
        assemblies = list(
            derived = list(
                label = "Derived wheat + rye",
                recipient = "wheat",
                donor = "rye"
            ),
            explicit = list(
                label = "Explicit wheat + rye",
                recipient = "wheat",
                donor = "rye",
                pairs = list(r1 = "w2")
            )
        )
    )

    config <- readMaMaMIAConfig(path)
    expect_identical(
        config$assemblies$derived$pairs, c(r1 = "w1", r2 = "w2")
    )
    expect_identical(
        config$assemblies$explicit$pairs, c(r1 = "w2")
    )
    expect_identical(config$assemblies$explicit$contigs$donor$id, "r1")
    expect_identical(
        listMaMaMIAAssemblies(config)$n_pairs, c(2L, 1L)
    )
})

test_that("polyploid species need a subgenome choice or a default", {
    root <- tempfile("mamamia-config-species-subgenome-")
    index <- file.path(root, "idx")
    write_index_fixture(
        index,
        c("w1", "w2", "r1"),
        rep(20000, 3)
    )
    wheat <- species_fixture(
        index,
        list(
            contig_entry("w1", "1A", "A"),
            contig_entry("w2", "1B", "B")
        ),
        label = "Wheat"
    )
    rye <- species_fixture(
        index,
        list(contig_entry("r1", "1R", "R")),
        label = "Rye"
    )
    path <- write_config_fixture(
        root, assemblies = NULL, species = list(wheat = wheat, rye = rye)
    )
    config <- readMaMaMIAConfig(path)

    expect_error(
        resolveSpeciesAssembly(config, "rye", "wheat"),
        "several subgenomes"
    )
    a <- resolveSpeciesAssembly(
        config, "rye", "wheat", recipient_subgenome = "A"
    )
    b <- resolveSpeciesAssembly(
        config, "rye", "wheat", recipient_subgenome = "B"
    )
    expect_identical(a$assembly$pairs, c(r1 = "w1"))
    expect_identical(b$assembly$pairs, c(r1 = "w2"))
    expect_identical(
        a$assembly$contigs$recipient$subgenome, "A"
    )

    wheat$default_subgenome <- "B"
    path <- write_config_fixture(
        file.path(root, "default"),
        assemblies = NULL,
        species = list(wheat = wheat, rye = rye)
    )
    config <- readMaMaMIAConfig(path)
    expect_identical(
        resolveSpeciesAssembly(config, "rye", "wheat")$assembly$contigs$recipient$id,
        "w2"
    )
})

test_that("unpaired species contigs are reported and never analysed", {
    root <- tempfile("mamamia-config-species-unpaired-")
    write_index_fixture(
        file.path(root, "idx"),
        c("w1", "r1", "r2", "r3"),
        rep(20000, 4)
    )
    wheat <- species_fixture(
        "idx", list(contig_entry("w1", "1A", "A")), label = "Wheat"
    )
    rye <- species_fixture(
        "idx",
        list(
            contig_entry("r1", "1R", "R"),
            contig_entry("r2", "2R", "R"),
            contig_entry("r3", "unplaced", "R")
        ),
        label = "Rye"
    )
    path <- write_config_fixture(
        root, assemblies = NULL, species = list(wheat = wheat, rye = rye)
    )
    resolved <- resolveSpeciesAssembly(readMaMaMIAConfig(path), "rye", "wheat")

    expect_identical(resolved$assembly$pairs, c(r1 = "w1"))
    expect_identical(resolved$assembly$contigs$donor$id, "r1")
    expect_true(all(c("r2", "r3") %in% resolved$unpaired$contig))
    expect_true(all(resolved$unpaired$role == "donor"))
})

test_that("incompatible species pairs fail with a clear message", {
    root <- tempfile("mamamia-config-species-invalid-")
    write_index_fixture(
        file.path(root, "idx"),
        c("w1", "d1"),
        rep(20000, 2)
    )
    wheat <- species_fixture(
        "idx", list(contig_entry("w1", "1A", "A")), label = "Wheat"
    )
    same_label <- species_fixture(
        "idx", list(contig_entry("d1", "1A", "A")), label = "Other"
    )
    path <- write_config_fixture(
        root,
        assemblies = NULL,
        species = list(wheat = wheat, other = same_label)
    )
    config <- readMaMaMIAConfig(path)
    expect_error(
        resolveSpeciesAssembly(config, "other", "wheat"),
        "subgenome labels must be disjoint"
    )
    expect_error(
        resolveSpeciesAssembly(config, "wheat", "wheat"),
        "two different species"
    )
    expect_error(
        resolveSpeciesAssembly(config, "missing", "wheat"),
        "Unknown donor species"
    )
})

test_that("species configs reject structural mistakes", {
    root <- tempfile("mamamia-config-species-errors-")
    index <- file.path(root, "idx")
    write_index_fixture(
        index,
        c("w1", "r1"),
        rep(20000, 2)
    )
    wheat <- species_fixture(
        index, list(contig_entry("w1", "1A", "A")), label = "Wheat"
    )

    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "one-species"),
            assemblies = NULL, species = list(wheat = wheat)
        )),
        "at least two"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "unknown-reference"),
            species = list(wheat = wheat),
            assemblies = list(
                broken = list(
                    label = "Broken", recipient = "wheat", donor = "rye"
                )
            )
        )),
        "unknown donor species"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "default"),
            assemblies = NULL,
            species = list(
                wheat = species_fixture(
                    index, list(contig_entry("w1", "1A", "A")),
                    default_subgenome = "B"
                ),
                rye = wheat
            )
        )),
        "default_subgenome"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "duplicate-group"),
            assemblies = NULL,
            species = list(
                wheat = species_fixture(
                    index,
                    list(
                        contig_entry("w1", "1A", "A"),
                        contig_entry("r1", "1A", "A")
                    )
                ),
                rye = wheat
            )
        )),
        "more than one contig with chromosome group"
    )
    expect_error(
        readMaMaMIAConfig(write_config_fixture(
            file.path(root, "missing-contig"),
            assemblies = NULL,
            species = list(
                wheat = species_fixture(
                    index, list(contig_entry("nope", "1A", "A"))
                ),
                rye = wheat
            )
        )),
        "not found in any configured gc_index"
    )
})
