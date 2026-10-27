## YAML assembly registry: loading, path resolution and validation.
##
## The registry is mounted read-only in production and is never edited by the
## application. Validation happens once, before any BAM work starts, so that a
## broken index or an inconsistent contig map fails with a message that names
## the assembly, contig and file.
##
## Version 2 is ploidy-aware: a species declares its subgenomes per contig and
## an assembly selects one or more recipient and donor subgenomes. Contigs are
## paired by chromosome group; the scoring layer aggregates homoeologs within a
## group instead of requiring a one-to-one contig map.

.MaMaMIASupportedConfigVersions <- 2L
.MaMaMIASupportedGcrankFormatVersions <- "0.1.0"

.config_error <- function(fmt, ...) {
    message <- sprintf(fmt, ...)
    condition <- structure(
        class = c("MaMaMIAConfigError", "MaMaMIAError", "error", "condition"),
        list(message = message, call = NULL)
    )
    stop(condition)
}

.config_defaults <- function() {
    list(
        windows_per_chromosome = 1000L,
        allowed_windows = c(100L, 1000L, 10000L),
        mapq = 20L,
        exclude_flags = 1796L,
        gc_mmap = TRUE,
        gc_validate_crc = FALSE,
        metric = "breadth",
        score = "paired_log_ratio",
        donor_aggregation = "max",
        recipient_aggregation = "mean",
        mask = list(
            panel = NULL,
            blacklist = NULL,
            physical_fallback_mb = 1,
            recurrence_threshold = 0.4
        ),
        calibration = list(
            length_bins_mb = c(5, 10, 25, 50, 100, Inf),
            fdr = 0.01,
            single_sample_floors = c(2.0, 1.5, 1.2, 1.0, 0.5)
        )
    )
}

.config_server_defaults <- function() {
    list(
        max_bams = 20L,
        max_total_bytes = 5 * 1024^3
    )
}

.config_string <- function(x, what) {
    if (length(x) != 1L || is.na(x) || is.list(x) || !is.atomic(x)) {
        .config_error("`%s` must be a single non-empty value", what)
    }
    value <- trimws(as.character(x))
    if (!nzchar(value)) {
        .config_error("`%s` must be a single non-empty value", what)
    }
    value
}

.config_integer <- function(x, what, min = 0, max = .Machine$integer.max) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
        !is.finite(x) || x != floor(x)) {
        .config_error("`%s` must be a single integer", what)
    }
    if (x < min || x > max) {
        .config_error(
            "`%s` must be between %s and %s (got %s)",
            what, format(min, scientific = FALSE),
            format(max, scientific = FALSE), format(x, scientific = FALSE)
        )
    }
    as.integer(x)
}

.config_number <- function(x, what, min = 0, max = Inf) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
        !is.finite(x) || x < min || x > max) {
        .config_error("`%s` must be a single number in [%s, %s]", what,
            format(min, scientific = FALSE),
            format(max, scientific = FALSE)
        )
    }
    as.numeric(x)
}

.config_logical <- function(x, what) {
    if (length(x) != 1L || is.na(x) || !is.logical(x)) {
        .config_error("`%s` must be true or false", what)
    }
    x
}

.config_path <- function(path, base_dir) {
    path <- path.expand(path)
    absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
    if (!absolute) {
        path <- file.path(base_dir, path)
    }
    normalizePath(path, winslash = "/", mustWork = FALSE)
}

.config_path_or_null <- function(path, what, base_dir) {
    if (is.null(path)) {
        return(NULL)
    }
    value <- .config_string(path, what)
    .config_path(value, base_dir)
}

.read_contigs_tsv <- function(path) {
    if (!file.exists(path)) {
        .config_error("Missing contigs.tsv: %s", path)
    }
    contigs <- tryCatch(
        utils::read.delim(
            path,
            header = FALSE,
            col.names = c("id", "length", "offset"),
            colClasses = c("character", "numeric", "numeric"),
            quote = "",
            comment.char = "",
            stringsAsFactors = FALSE
        ),
        error = function(e) {
            .config_error("Could not read %s: %s", path, conditionMessage(e))
        }
    )
    if (nrow(contigs) == 0L) {
        .config_error("No contigs found in %s", path)
    }
    if (anyNA(contigs$id) || any(!nzchar(contigs$id))) {
        .config_error("Empty contig name in %s", path)
    }
    if (anyNA(contigs$length) || any(contigs$length <= 0) ||
        any(contigs$length != floor(contigs$length))) {
        .config_error("Contig lengths in %s must be positive integers", path)
    }
    if (anyNA(contigs$offset) || any(contigs$offset < 0) ||
        any(contigs$offset != floor(contigs$offset))) {
        .config_error("Contig offsets in %s must be non-negative integers", path)
    }
    if (anyDuplicated(contigs$id)) {
        duplicate <- unique(contigs$id[duplicated(contigs$id)])[[1L]]
        .config_error("Duplicate contig '%s' in %s", duplicate, path)
    }
    expected_offsets <- cumsum(contigs$length) - contigs$length
    if (any(contigs$offset != expected_offsets)) {
        .config_error(
            "Contig offsets in %s are not the cumulative lengths of the preceding contigs",
            path
        )
    }
    contigs
}

.read_gcrank_index <- function(index_path, assembly_key) {
    if (!dir.exists(index_path)) {
        .config_error(
            "Assembly '%s': gc_index '%s' does not exist or is not a directory",
            assembly_key, index_path
        )
    }
    contigs_path <- file.path(index_path, "contigs.tsv")
    manifest_path <- file.path(index_path, "manifest.json")
    if (!file.exists(contigs_path)) {
        .config_error(
            "Assembly '%s': gc_index '%s' has no contigs.tsv (%s)",
            assembly_key, index_path, contigs_path
        )
    }
    if (!file.exists(manifest_path)) {
        .config_error(
            "Assembly '%s': gc_index '%s' has no manifest.json (%s)",
            assembly_key, index_path, manifest_path
        )
    }

    manifest <- tryCatch(
        jsonlite::fromJSON(manifest_path, simplifyVector = TRUE),
        error = function(e) {
            .config_error(
                "Assembly '%s': could not read %s: %s",
                assembly_key, manifest_path, conditionMessage(e)
            )
        }
    )
    format_version <- manifest$format_version
    if (length(format_version) != 1L || is.na(format_version) ||
        !is.character(format_version)) {
        .config_error(
            "Assembly '%s': %s has no string format_version",
            assembly_key, manifest_path
        )
    }
    if (!format_version %in% .MaMaMIASupportedGcrankFormatVersions) {
        .config_error(
            "Assembly '%s': gc_index '%s' uses unsupported gcrank format_version '%s' (supported: %s)",
            assembly_key, index_path, format_version,
            paste(.MaMaMIASupportedGcrankFormatVersions, collapse = ", ")
        )
    }

    contigs <- .read_contigs_tsv(contigs_path)
    total_bases <- manifest$source$total_bases
    if (length(total_bases) != 1L || is.na(total_bases) ||
        !is.numeric(total_bases) || !is.finite(total_bases)) {
        .config_error(
            "Assembly '%s': %s has no numeric source.total_bases",
            assembly_key, manifest_path
        )
    }
    expected_bases <- sum(contigs$length)
    if (total_bases != expected_bases) {
        .config_error(
            "Assembly '%s': source.total_bases in %s is %s but contigs.tsv sums to %s",
            assembly_key, manifest_path,
            format(total_bases, scientific = FALSE),
            format(expected_bases, scientific = FALSE)
        )
    }
    num_contigs <- manifest$source$num_contigs
    if (length(num_contigs) == 1L && is.numeric(num_contigs) &&
        !is.na(num_contigs) && num_contigs != nrow(contigs)) {
        .config_error(
            "Assembly '%s': source.num_contigs in %s is %s but contigs.tsv has %s rows",
            assembly_key, manifest_path,
            format(num_contigs, scientific = FALSE), nrow(contigs)
        )
    }

    contigs$index_path <- index_path
    contigs
}

.collect_index_contigs <- function(gc_indexes, assembly_key) {
    pieces <- lapply(
        gc_indexes,
        .read_gcrank_index,
        assembly_key = assembly_key
    )
    contigs <- do.call(rbind, pieces)
    rownames(contigs) <- NULL
    duplicated_ids <- unique(contigs$id[duplicated(contigs$id)])
    if (length(duplicated_ids) > 0L) {
        id <- duplicated_ids[[1L]]
        paths <- unique(contigs$index_path[contigs$id == id])
        .config_error(
            "Assembly '%s': contig '%s' appears in more than one gc_index (%s)",
            assembly_key, id, paste(paths, collapse = ", ")
        )
    }
    contigs
}

.is_contig_mapping <- function(x) {
    is.list(x) && !is.null(names(x)) &&
        any(names(x) %in% c("id", "name", "subgenome"))
}

## First number in a contig display name, used as the chromosome group when the
## registry does not declare one (for example "1A" and "1R" both give "1").
#' @noRd
chr_number <- function(chr_name) {
    has_number <- grepl("[0-9]", chr_name)
    number <- sub("^[^0-9]*([0-9]+).*$", "\\1", chr_name)
    number[!has_number] <- NA_character_
    number
}

.config_copy_number <- function(x, what) {
    if (is.null(x)) {
        return(2)
    }
    .config_integer(x, what, min = 0L, max = 1000L)
}

.species_contig_entries <- function(x, species_key) {
    if (is.null(x) || length(x) == 0L) {
        return(data.frame(
            id = character(0), name = character(0),
            subgenome = character(0), group = character(0),
            copy_number = numeric(0), stringsAsFactors = FALSE
        ))
    }
    entries <- if (is.data.frame(x)) {
        lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
    } else if (.is_contig_mapping(x)) {
        list(x)
    } else if (is.list(x)) {
        x
    } else {
        .config_error(
            "Species '%s': contigs must be a list of mappings", species_key
        )
    }

    rows <- lapply(seq_along(entries), function(i) {
        entry <- entries[[i]]
        where <- sprintf("Species '%s': contigs[[%d]]", species_key, i)
        if (!is.list(entry)) {
            .config_error("%s must be a mapping with id, name and subgenome", where)
        }
        id <- .config_string(entry$id, paste0(where, ".id"))
        name <- if (is.null(entry$name)) id else {
            .config_string(entry$name, paste0(where, ".name"))
        }
        subgenome <- .config_string(
            entry$subgenome, paste0(where, ".subgenome")
        )
        group <- if (is.null(entry$group)) {
            NA_character_
        } else {
            .config_string(entry$group, paste0(where, ".group"))
        }
        data.frame(
            id = id, name = name, subgenome = subgenome, group = group,
            copy_number = .config_copy_number(
                entry$copy_number, paste0(where, ".copy_number")
            ),
            stringsAsFactors = FALSE
        )
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL

    if (anyDuplicated(out$id)) {
        duplicate <- unique(out$id[duplicated(out$id)])[[1L]]
        .config_error(
            "Species '%s': contig '%s' is configured more than once",
            species_key, duplicate
        )
    }

    derived <- chr_number(out$name)
    missing_group <- is.na(out$group) | !nzchar(out$group)
    out$group[missing_group] <- derived[missing_group]

    for (subgenome in unique(out$subgenome)) {
        idx <- out$subgenome == subgenome & !is.na(out$group)
        if (anyDuplicated(out$group[idx])) {
            group <- unique(out$group[idx][duplicated(out$group[idx])])[[1L]]
            .config_error(
                "Species '%s': subgenome '%s' has more than one contig with chromosome group '%s'",
                species_key, subgenome, group
            )
        }
    }
    out
}

.normalize_species <- function(species, key, config_dir) {
    if (!is.list(species)) {
        .config_error("Species '%s' must be a mapping", key)
    }
    label <- .config_string(species$label, sprintf("species.%s.label", key))
    description <- if (is.null(species$description)) "" else {
        x <- species$description
        if (length(x) != 1L || is.na(x) || is.list(x) || !is.atomic(x)) {
            .config_error("`species.%s.description` must be a string", key)
        }
        trimws(as.character(x))
    }

    gc_indexes <- .config_gc_indexes(
        species$gc_indexes, sprintf("Species '%s'", key), config_dir
    )
    contigs <- .species_contig_entries(species$contigs, key)
    if (nrow(contigs) == 0L) {
        .config_error("Species '%s' must have at least one contig", key)
    }

    index_contigs <- .collect_index_contigs(gc_indexes, key)
    contigs <- .attach_contig_lengths(
        contigs, index_contigs, sprintf("Species '%s'", key), "configured"
    )
    too_long <- contigs$length >= 2^31
    if (any(too_long)) {
        contig <- contigs$id[too_long][[1L]]
        .config_error(
            "Species '%s': contig '%s' is %s bp; 32-bit positions require lengths below 2^31",
            key, contig,
            format(contigs$length[too_long][[1L]], scientific = FALSE)
        )
    }

    available <- unique(contigs$subgenome)
    default_subgenome <- if (is.null(species$default_subgenome)) NULL else {
        value <- .config_string(
            species$default_subgenome,
            sprintf("species.%s.default_subgenome", key)
        )
        if (!value %in% available) {
            .config_error(
                "Species '%s': default_subgenome '%s' is not one of: %s",
                key, value, paste(available, collapse = ", ")
            )
        }
        value
    }

    list(
        key = key,
        label = label,
        description = description,
        default_subgenome = default_subgenome,
        gc_indexes = unname(gc_indexes),
        contigs = contigs[, c(
            "id", "name", "subgenome", "group", "copy_number", "length",
            "index_path"
        )]
    )
}

.config_gc_indexes <- function(gc_indexes, what, config_dir) {
    if (is.null(gc_indexes) || length(gc_indexes) == 0L) {
        .config_error("%s: gc_indexes must not be empty", what)
    }
    if (!is.character(gc_indexes)) {
        gc_indexes <- unlist(gc_indexes, use.names = FALSE)
    }
    if (!is.character(gc_indexes) || anyNA(gc_indexes) ||
        any(!nzchar(trimws(gc_indexes)))) {
        .config_error(
            "%s: gc_indexes must be a list of non-empty paths", what
        )
    }
    gc_indexes <- vapply(gc_indexes, .config_path, character(1),
        base_dir = config_dir
    )
    if (anyDuplicated(gc_indexes)) {
        .config_error("%s: gc_indexes lists the same index twice", what)
    }
    unname(gc_indexes)
}

.attach_contig_lengths <- function(contigs, index_contigs, what, role) {
    position <- match(contigs$id, index_contigs$id)
    if (anyNA(position)) {
        missing_ids <- contigs$id[is.na(position)]
        .config_error(
            "%s: %s contig(s) not found in any configured gc_index: %s",
            what, role, paste(missing_ids, collapse = ", ")
        )
    }
    contigs$length <- index_contigs$length[position]
    contigs$index_path <- index_contigs$index_path[position]
    contigs
}

.check_contig_lengths <- function(contigs, what, allowed_windows) {
    max_resolution <- max(allowed_windows)
    too_short <- contigs$length < max_resolution
    if (any(too_short)) {
        contig <- contigs$id[too_short][[1L]]
        length_bp <- contigs$length[too_short][[1L]]
        .config_error(
            "%s: contig '%s' is %s bp, shorter than the largest allowed window (%s bp)",
            what, contig, format(length_bp, scientific = FALSE),
            format(max_resolution, scientific = FALSE)
        )
    }
    too_long <- contigs$length >= 2^31
    if (any(too_long)) {
        contig <- contigs$id[too_long][[1L]]
        .config_error(
            "%s: contig '%s' is %s bp; 32-bit positions require lengths below 2^31",
            what, contig,
            format(contigs$length[too_long][[1L]], scientific = FALSE)
        )
    }
    invisible(contigs)
}

.species_resolve_subgenomes <- function(species, requested, role, context) {
    available <- unique(as.character(species$contigs$subgenome))
    if (is.null(requested) || length(requested) == 0L) {
        if (!is.null(species$default_subgenome)) {
            return(species$default_subgenome)
        }
        return(available)
    }
    if (!is.character(requested)) {
        requested <- unlist(requested, use.names = FALSE)
    }
    requested <- trimws(as.character(requested))
    if (length(requested) == 0L || anyNA(requested) || any(!nzchar(requested))) {
        .config_error("%s: `%s_subgenomes` must be a non-empty character vector",
            context, role)
    }
    unknown <- setdiff(requested, available)
    if (length(unknown) > 0L) {
        .config_error(
            "%s: species '%s' has no %s subgenome(s) %s; available: %s",
            context, species$key, role,
            paste(sprintf("'%s'", unknown), collapse = ", "),
            paste(available, collapse = ", ")
        )
    }
    unique(requested)
}

.species_derived_pairs <- function(donor, recipient) {
    donor_group <- donor$group
    recipient_group <- recipient$group
    pairs <- character(0)
    unpaired <- data.frame(
        role = character(0), contig = character(0), reason = character(0),
        stringsAsFactors = FALSE
    )

    for (i in seq_len(nrow(donor))) {
        group <- donor_group[[i]]
        if (is.na(group) || !nzchar(group)) {
            unpaired <- rbind(unpaired, data.frame(
                role = "donor", contig = donor$id[[i]],
                reason = "no chromosome group", stringsAsFactors = FALSE
            ))
            next
        }
        matched <- which(!is.na(recipient_group) & recipient_group == group)
        if (length(matched) == 0L) {
            unpaired <- rbind(unpaired, data.frame(
                role = "donor", contig = donor$id[[i]],
                reason = sprintf("no recipient contig in group %s", group),
                stringsAsFactors = FALSE
            ))
            next
        }
        ## Group-level scoring uses every recipient homoeolog; this map names a
        ## representative so one-to-one consumers still have a pairing.
        pairs[[donor$id[[i]]]] <- recipient$id[[matched[[1L]]]]
    }

    used <- unname(pairs)
    unused <- recipient$id[!recipient$id %in% used]
    for (id in unused) {
        group <- recipient$group[match(id, recipient$id)]
        unpaired <- rbind(unpaired, data.frame(
            role = "recipient", contig = id,
            reason = if (is.na(group) || !nzchar(group)) {
                "no chromosome group"
            } else {
                sprintf("no donor contig in group %s", group)
            },
            stringsAsFactors = FALSE
        ))
    }

    list(pairs = pairs, unpaired = unpaired)
}

.sample_composition <- function(contigs, subgenomes, requested, context) {
    if (is.null(requested)) {
        values <- vapply(subgenomes, function(s) {
            x <- contigs$copy_number[contigs$subgenome == s]
            if (length(x) == 0L) 2 else as.numeric(x[[1L]])
        }, numeric(1))
        return(stats::setNames(values, subgenomes))
    }
    values <- unlist(requested, use.names = TRUE)
    if (is.null(values) || is.null(names(values)) || anyNA(names(values)) ||
        any(!nzchar(names(values)))) {
        .config_error(
            "%s: `default_sample_composition` must be a named numeric mapping",
            context
        )
    }
    unknown <- setdiff(names(values), subgenomes)
    if (length(unknown) > 0L) {
        .config_error(
            "%s: `default_sample_composition` names unknown recipient subgenome(s): %s",
            context, paste(unknown, collapse = ", ")
        )
    }
    out <- vapply(seq_along(values), function(i) {
        .config_integer(
            values[[i]],
            sprintf("%s: default_sample_composition.%s", context, names(values)[[i]]),
            min = 0L, max = 1000L
        )
    }, integer(1))
    names(out) <- names(values)
    missing <- setdiff(subgenomes, names(out))
    if (length(missing) > 0L) {
        defaults <- vapply(missing, function(s) {
            x <- contigs$copy_number[contigs$subgenome == s]
            if (length(x) == 0L) 2 else as.numeric(x[[1L]])
        }, numeric(1))
        out <- c(out, stats::setNames(defaults, missing))
    }
    out[subgenomes]
}

.build_species_assembly <- function(recipient_species, donor_species, key = NULL,
                                    label = NULL, description = NULL,
                                    recipient_subgenomes = NULL,
                                    donor_subgenomes = NULL,
                                    recipient_subgenome = NULL,
                                    donor_subgenome = NULL,
                                    default_sample_composition = NULL,
                                    allowed_windows) {
    if (identical(recipient_species$key, donor_species$key)) {
        .config_error(
            "A species pair needs two different species; both sides are '%s'",
            recipient_species$key
        )
    }
    context <- sprintf(
        "Species pair '%s' (recipient) + '%s' (donor)",
        recipient_species$key, donor_species$key
    )
    if (is.null(recipient_subgenomes) && !is.null(recipient_subgenome)) {
        recipient_subgenomes <- recipient_subgenome
    }
    if (is.null(donor_subgenomes) && !is.null(donor_subgenome)) {
        donor_subgenomes <- donor_subgenome
    }
    recipient_subgenomes <- .species_resolve_subgenomes(
        recipient_species, recipient_subgenomes, "recipient", context
    )
    donor_subgenomes <- .species_resolve_subgenomes(
        donor_species, donor_subgenomes, "donor", context
    )

    recipient_all <- recipient_species$contigs[
        recipient_species$contigs$subgenome %in% recipient_subgenomes, ,
        drop = FALSE
    ]
    donor_all <- donor_species$contigs[
        donor_species$contigs$subgenome %in% donor_subgenomes, ,
        drop = FALSE
    ]
    if (nrow(recipient_all) == 0L || nrow(donor_all) == 0L) {
        .config_error(
            "%s: no contigs in the selected subgenomes (%s/%s)",
            context, paste(recipient_subgenomes, collapse = ","),
            paste(donor_subgenomes, collapse = ",")
        )
    }

    derived <- .species_derived_pairs(donor_all, recipient_all)
    pairs <- derived$pairs
    unpaired <- derived$unpaired
    if (length(pairs) == 0L) {
        .config_error(
            "%s: no contig pairs could be derived by chromosome group; add `group` values to the species contigs",
            context
        )
    }

    donor <- donor_all
    recipient <- recipient_all
    if (anyDuplicated(c(recipient$id, donor$id))) {
        duplicate <- unique(c(recipient$id, donor$id)[
            duplicated(c(recipient$id, donor$id))
        ])[[1L]]
        .config_error(
            "%s: contig '%s' is present on both sides",
            context, duplicate
        )
    }

    shared <- intersect(
        unique(recipient$subgenome), unique(donor$subgenome)
    )
    if (length(shared) > 0L) {
        .config_error(
            "%s: donor and recipient subgenome labels must be disjoint; shared: %s",
            context, paste(shared, collapse = ", ")
        )
    }
    .check_contig_lengths(
        rbind(recipient, donor), context, allowed_windows
    )

    if (is.null(key)) {
        key <- paste(
            recipient_species$key, paste(recipient_subgenomes, collapse = "+"),
            donor_species$key, paste(donor_subgenomes, collapse = "+"),
            sep = "__"
        )
    } else {
        key <- .config_string(key, "assembly key")
    }
    if (is.null(label)) {
        label <- sprintf("%s + %s", recipient_species$label, donor_species$label)
    } else {
        label <- .config_string(label, sprintf("assemblies.%s.label", key))
    }
    if (is.null(description)) {
        description <- sprintf(
            paste0(
                "Species pair: %s [%s] (recipient) + %s [%s] (donor); ",
                "contigs paired by chromosome group"
            ),
            recipient_species$label,
            paste(recipient_subgenomes, collapse = ","),
            donor_species$label,
            paste(donor_subgenomes, collapse = ",")
        )
    } else {
        description <- .config_string(
            description, sprintf("assemblies.%s.description", key)
        )
    }

    composition <- .sample_composition(
        recipient, recipient_subgenomes, default_sample_composition, context
    )
    columns <- c(
        "id", "name", "subgenome", "group", "copy_number", "length",
        "index_path"
    )
    assembly <- list(
        key = key,
        label = label,
        description = description,
        gc_indexes = unique(c(recipient$index_path, donor$index_path)),
        contigs = list(
            recipient = recipient[, columns, drop = FALSE],
            donor = donor[, columns, drop = FALSE]
        ),
        pairs = pairs,
        recipient_subgenomes = recipient_subgenomes,
        donor_subgenomes = donor_subgenomes,
        sample_composition = composition
    )
    list(
        key = key,
        assembly = assembly,
        unpaired = unpaired,
        recipient_subgenomes = recipient_subgenomes,
        donor_subgenomes = donor_subgenomes,
        recipient_subgenome = recipient_subgenomes[[1L]],
        donor_subgenome = donor_subgenomes[[1L]]
    )
}

.validate_metric <- function(x, what) {
    value <- .config_string(x, what)
    allowed <- c("breadth", "count", "combined")
    if (!value %in% allowed) {
        .config_error("`%s` must be one of: %s", what, paste(allowed, collapse = ", "))
    }
    value
}

.validate_defaults <- function(defaults, config_dir = getwd()) {
    if (is.null(defaults)) {
        defaults <- list()
    }
    if (!is.list(defaults)) {
        .config_error("`defaults` must be a mapping")
    }
    supported <- c(
        "windows_per_chromosome", "allowed_windows", "mapq", "exclude_flags",
        "gc_mmap", "gc_validate_crc", "metric", "score", "donor_aggregation",
        "recipient_aggregation", "mask", "calibration"
    )
    unknown <- setdiff(names(defaults), supported)
    if (length(unknown) > 0L) {
        .config_error(
            "Unknown `defaults` key(s): %s. Supported keys are: %s. See MIGRATION.md for the version 1 key changes.",
            paste(unknown, collapse = ", "),
            paste(supported, collapse = ", ")
        )
    }
    out <- .config_defaults()

    if (!is.null(defaults$allowed_windows)) {
        allowed <- defaults$allowed_windows
        if (length(allowed) == 0L || !is.numeric(allowed) ||
            anyNA(allowed) || any(!is.finite(allowed)) ||
            any(allowed != floor(allowed)) || any(allowed <= 0)) {
            .config_error("`defaults.allowed_windows` must be positive integers")
        }
        if (anyDuplicated(allowed)) {
            .config_error("`defaults.allowed_windows` must not contain duplicates")
        }
        if (any(allowed > .Machine$integer.max)) {
            .config_error("`defaults.allowed_windows` exceeds the 32-bit window limit")
        }
        out$allowed_windows <- sort(as.integer(allowed))
    }
    if (!is.null(defaults$windows_per_chromosome)) {
        out$windows_per_chromosome <- .config_integer(
            defaults$windows_per_chromosome,
            "defaults.windows_per_chromosome",
            min = 1L
        )
    }
    if (!out$windows_per_chromosome %in% out$allowed_windows) {
        .config_error(
            "`defaults.windows_per_chromosome` (%s) must be one of `defaults.allowed_windows` (%s)",
            out$windows_per_chromosome,
            paste(out$allowed_windows, collapse = ", ")
        )
    }
    if (!is.null(defaults$mapq)) {
        out$mapq <- .config_integer(defaults$mapq, "defaults.mapq", 0L, 255L)
    }
    if (!is.null(defaults$exclude_flags)) {
        out$exclude_flags <- .config_integer(
            defaults$exclude_flags, "defaults.exclude_flags", 0L, 4095L
        )
    }
    if (!is.null(defaults$gc_mmap)) {
        out$gc_mmap <- .config_logical(defaults$gc_mmap, "defaults.gc_mmap")
    }
    if (!is.null(defaults$gc_validate_crc)) {
        out$gc_validate_crc <- .config_logical(
            defaults$gc_validate_crc, "defaults.gc_validate_crc"
        )
    }
    if (!is.null(defaults$metric)) {
        out$metric <- .validate_metric(defaults$metric, "defaults.metric")
    }
    if (!is.null(defaults$score)) {
        out$score <- .config_string(defaults$score, "defaults.score")
        if (!out$score %in% c("paired_log_ratio", "donor_only")) {
            .config_error(
                "`defaults.score` must be one of: paired_log_ratio, donor_only"
            )
        }
    }
    if (!is.null(defaults$donor_aggregation)) {
        out$donor_aggregation <- .config_string(
            defaults$donor_aggregation, "defaults.donor_aggregation"
        )
        if (!out$donor_aggregation %in% c("mean", "max", "sum")) {
            .config_error(
                "`defaults.donor_aggregation` must be one of: mean, max, sum"
            )
        }
    }
    if (!is.null(defaults$recipient_aggregation)) {
        out$recipient_aggregation <- .config_string(
            defaults$recipient_aggregation, "defaults.recipient_aggregation"
        )
        if (!identical(out$recipient_aggregation, "mean")) {
            .config_error("`defaults.recipient_aggregation` must be 'mean'")
        }
    }

    if (!is.null(defaults$mask)) {
        if (!is.list(defaults$mask)) {
            .config_error("`defaults.mask` must be a mapping")
        }
        mask <- defaults$mask
        unknown <- setdiff(
            names(mask),
            c("panel", "blacklist", "physical_fallback_mb",
              "recurrence_threshold")
        )
        if (length(unknown) > 0L) {
            .config_error("`defaults.mask` has unknown key(s): %s",
                paste(unknown, collapse = ", "))
        }
        out$mask$panel <- .config_path_or_null(
            mask$panel, "defaults.mask.panel", config_dir
        )
        out$mask$blacklist <- .config_path_or_null(
            mask$blacklist, "defaults.mask.blacklist", config_dir
        )
        if (!is.null(mask$physical_fallback_mb)) {
            out$mask$physical_fallback_mb <- .config_number(
                mask$physical_fallback_mb, "defaults.mask.physical_fallback_mb",
                min = 0
            )
        }
        if (!is.null(mask$recurrence_threshold)) {
            out$mask$recurrence_threshold <- .config_number(
                mask$recurrence_threshold,
                "defaults.mask.recurrence_threshold", min = 0, max = 1
            )
        }
    }
    if (!is.null(defaults$calibration)) {
        if (!is.list(defaults$calibration)) {
            .config_error("`defaults.calibration` must be a mapping")
        }
        calibration <- defaults$calibration
        unknown <- setdiff(
            names(calibration),
            c("length_bins_mb", "fdr", "single_sample_floors")
        )
        if (length(unknown) > 0L) {
            .config_error("`defaults.calibration` has unknown key(s): %s",
                paste(unknown, collapse = ", "))
        }
        if (!is.null(calibration$length_bins_mb)) {
            breaks <- as.numeric(unlist(calibration$length_bins_mb))
            if (length(breaks) < 2L || anyNA(breaks) || any(!is.finite(breaks[-length(breaks)])) ||
                any(diff(breaks) <= 0)) {
                .config_error(
                    "`defaults.calibration.length_bins_mb` must be increasing"
                )
            }
            out$calibration$length_bins_mb <- breaks
        }
        if (!is.null(calibration$fdr)) {
            out$calibration$fdr <- .config_number(
                calibration$fdr, "defaults.calibration.fdr",
                min = .Machine$double.eps, max = 1 - .Machine$double.eps
            )
        }
        if (!is.null(calibration$single_sample_floors)) {
            floors <- as.numeric(unlist(calibration$single_sample_floors))
            if (anyNA(floors) || any(!is.finite(floors)) || any(floors < 0)) {
                .config_error(
                    "`defaults.calibration.single_sample_floors` must be non-negative numbers"
                )
            }
            out$calibration$single_sample_floors <- floors
        }
    }
    out
}

.validate_server <- function(server) {
    if (is.null(server)) {
        server <- list()
    }
    if (!is.list(server)) {
        .config_error("`server` must be a mapping")
    }
    out <- .config_server_defaults()
    if (!is.null(server$max_bams)) {
        out$max_bams <- .config_integer(server$max_bams, "server.max_bams", min = 1L)
    }
    if (!is.null(server$max_total_bytes)) {
        out$max_total_bytes <- .config_number(
            server$max_total_bytes, "server.max_total_bytes", min = 1
        )
    }
    out
}

#' Read and validate a MaMaMIA assembly registry
#'
#' Reads the YAML configuration, resolves relative paths against the directory
#' containing the YAML file, and validates the species registry, assemblies, GC
#' indexes, contigs and pairs before returning a normalized configuration list.
#'
#' Version 2 is required. A version 1 file raises a migration error that points
#' to `MIGRATION.md`.
#'
#' @param path Path to the YAML configuration file.
#'
#' @return A normalized `MaMaMIAConfig` list with `version`, `defaults`,
#'   `species`, `assemblies` and `server` elements.
#'
#' @examples
#' demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
#' config <- readMaMaMIAConfig(file.path(demo, "assemblies.yml"))
#' names(config$species)
#' config$defaults$windows_per_chromosome
#'
#' @export
readMaMaMIAConfig <- function(path) {
    if (length(path) != 1L || is.na(path) || !is.character(path) ||
        !nzchar(path)) {
        .config_error("`path` must be a single non-empty file path")
    }
    if (!file.exists(path)) {
        .config_error("MaMaMIA config file does not exist: %s", path)
    }
    path <- normalizePath(path, winslash = "/", mustWork = TRUE)
    config <- tryCatch(
        yaml::read_yaml(
            path,
            ## Parse YAML integers as doubles: max_total_bytes (5 GiB by
            ## default) exceeds R's 32-bit integer range and would otherwise
            ## become NA with a coercion warning.
            handlers = list(int = function(value) as.numeric(value))
        ),
        error = function(e) {
            .config_error(
                "Could not read MaMaMIA config '%s': %s",
                path, conditionMessage(e)
            )
        }
    )
    if (!is.list(config)) {
        .config_error("MaMaMIA config '%s' must contain a YAML mapping", path)
    }
    attr(config, "config_path") <- path
    validateMaMaMIAConfig(config)
}

#' Validate a MaMaMIA configuration
#'
#' Takes a configuration list, either read from YAML or assembled in memory,
#' and returns a normalized configuration. Relative GC index paths are resolved
#' against the directory recorded by `readMaMaMIAConfig()`, or the current
#' working directory for a plain list.
#'
#' @param config A configuration list with `version`, optional `defaults`,
#'   optional `server`, optional `assemblies`, and optional `species`. At least
#'   one assembly or two species are required.
#'
#' @return The validated, normalized configuration, invisibly.
#'
#' @export
validateMaMaMIAConfig <- function(config) {
    if (!is.list(config)) {
        .config_error("`config` must be a list")
    }
    config_path <- attr(config, "config_path", exact = TRUE)
    config_dir <- if (is.null(config_path)) {
        getwd()
    } else {
        dirname(config_path)
    }

    if (is.null(config$version)) {
        .config_error("`version` is required in the MaMaMIA config")
    }
    version <- suppressWarnings(as.integer(config$version))
    if (is.na(version)) {
        .config_error("`version` must be an integer")
    }
    if (version == 1L) {
        .config_error(
            paste0(
                "MaMaMIA config version 1 is no longer supported. ",
                "Change `version: 1` to `version: 2` and see MIGRATION.md ",
                "for the schema changes (multi-subgenome ploidy, metric, score, ",
                "mask and calibration defaults)."
            )
        )
    }
    if (!version %in% .MaMaMIASupportedConfigVersions) {
        .config_error(
            "Unsupported MaMaMIA config version %s (supported: %s)",
            version, paste(.MaMaMIASupportedConfigVersions, collapse = ", ")
        )
    }

    defaults <- .validate_defaults(config$defaults, config_dir)
    server <- .validate_server(config$server)

    species <- config$species
    if (is.null(species)) {
        species <- list()
    }
    if (!is.list(species)) {
        .config_error("`species` must be a named mapping")
    }
    if (length(species) > 0L &&
        (is.null(names(species)) || anyNA(names(species)) ||
            any(!nzchar(names(species))))) {
        .config_error("`species` must be a named mapping")
    }
    if (anyDuplicated(names(species))) {
        .config_error("Species keys must be unique")
    }
    normalized_species <- lapply(seq_along(species), function(i) {
        .normalize_species(
            species[[i]], names(species)[[i]], config_dir
        )
    })
    names(normalized_species) <- if (length(normalized_species) == 0L) {
        character(0)
    } else {
        names(species)
    }

    if (!is.null(config$assemblies) && length(config$assemblies) > 0L) {
        .config_error(
            paste0(
                "`assemblies` is not part of the registry schema: declare ",
                "`species` only and pair them at run time with ",
                "resolveSpeciesAssembly(). See MIGRATION.md."
            )
        )
    }
    if (length(normalized_species) < 2L) {
        .config_error(
            "The config must define at least two `species` entries to form a donor/recipient pair"
        )
    }

    out <- list(
        version = version,
        defaults = defaults,
        species = normalized_species,
        assemblies = list(),
        server = server
    )
    attr(out, "config_path") <- config_path
    invisible(out)
}

#' List the assemblies resolved in a MaMaMIA configuration
#'
#' Returns one row per assembly. A freshly read registry has none: assemblies
#' are created by [resolveSpeciesAssembly()] and stored in the returned
#' `config`, so pass that resolved configuration to list them.
#'
#' @param config A validated configuration list, usually the `config` element
#'   returned by [resolveSpeciesAssembly()], or a path to a YAML file.
#'
#' @return A data frame with `key`, `label`, `description`, `n_donor`,
#'   `n_recipient`, `n_pairs` and `n_gc_indexes`.
#'
#' @export
listMaMaMIAAssemblies <- function(config) {
    if (is.character(config)) {
        config <- readMaMaMIAConfig(config)
    }
    if (!is.list(config) || is.null(config$assemblies)) {
        .config_error("`config` must be a validated MaMaMIA configuration")
    }
    assemblies <- config$assemblies
    keys <- names(assemblies)
    out <- data.frame(
        key = keys,
        label = vapply(assemblies, function(x) x$label, character(1)),
        description = vapply(assemblies, function(x) {
            if (is.null(x$description)) "" else as.character(x$description)
        }, character(1)),
        n_donor = vapply(assemblies, function(x) {
            nrow(x$contigs$donor)
        }, integer(1)),
        n_recipient = vapply(assemblies, function(x) {
            nrow(x$contigs$recipient)
        }, integer(1)),
        n_pairs = vapply(assemblies, function(x) {
            length(x$pairs)
        }, integer(1)),
        n_gc_indexes = vapply(assemblies, function(x) {
            length(x$gc_indexes)
        }, integer(1)),
        recipient_subgenomes = vapply(assemblies, function(x) {
            paste(x$recipient_subgenomes, collapse = ",")
        }, character(1)),
        donor_subgenomes = vapply(assemblies, function(x) {
            paste(x$donor_subgenomes, collapse = ",")
        }, character(1)),
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
}

#' List the species in a MaMaMIA configuration
#'
#' Returns one row per configured species for the application's donor and
#' recipient selectors.
#'
#' @param config A validated configuration list or a path to a YAML file.
#'
#' @return A data frame with `key`, `label`, `description`, `subgenomes`,
#'   `default_subgenome`, `n_contigs` and `n_gc_indexes`.
#'
#' @export
listMaMaMIASpecies <- function(config) {
    if (is.character(config)) {
        config <- readMaMaMIAConfig(config)
    }
    if (!is.list(config) || is.null(config$species)) {
        .config_error("`config` must be a validated MaMaMIA configuration")
    }
    species <- config$species
    keys <- names(species)
    out <- data.frame(
        key = keys,
        label = vapply(species, function(x) x$label, character(1)),
        description = vapply(species, function(x) {
            if (is.null(x$description)) "" else as.character(x$description)
        }, character(1)),
        subgenomes = vapply(species, function(x) {
            paste(unique(x$contigs$subgenome), collapse = ", ")
        }, character(1)),
        default_subgenome = vapply(species, function(x) {
            if (is.null(x$default_subgenome)) "" else x$default_subgenome
        }, character(1)),
        n_contigs = vapply(species, function(x) {
            nrow(x$contigs)
        }, integer(1)),
        n_gc_indexes = vapply(species, function(x) {
            length(x$gc_indexes)
        }, integer(1)),
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL
    out
}

#' Build an assembly from a recipient and a donor species
#'
#' Selects one or more subgenomes per species and pairs their contigs by
#' chromosome group. The group is taken from the optional `group` field of a
#' species contig, or from the first number in its `name` (for example `1A` and
#' `1R` are both group 1). All contigs in the selected subgenomes are included;
#' the scoring layer compares donor and recipient homoeologs within a group.
#'
#' `recipient_subgenomes` and `donor_subgenomes` are optional when a species has
#' a single subgenome or a `default_subgenome`; otherwise every subgenome is
#' selected. The singular `recipient_subgenome` / `donor_subgenome` arguments
#' are accepted as one-element shortcuts.
#'
#' @param config A validated configuration list or a path to a YAML file.
#' @param donor Donor species key.
#' @param recipient Recipient species key.
#' @param donor_subgenomes Optional donor subgenome labels.
#' @param recipient_subgenomes Optional recipient subgenome labels.
#' @param donor_subgenome Optional single donor subgenome (shortcut).
#' @param recipient_subgenome Optional single recipient subgenome (shortcut).
#' @param key Optional assembly key. Defaults to a key built from the species
#'   and subgenome labels.
#' @param label Optional assembly label.
#' @param description Optional assembly description.
#' @param default_sample_composition Optional named recipient copy-number
#'   vector for the derived assembly.
#'
#' @return A list with `key`, `assembly`, `config` (the input configuration
#'   with the derived assembly registered), `unpaired` (contigs excluded from
#'   the pair, with a reason), `donor_species`, `recipient_species`,
#'   `donor_subgenomes` and `recipient_subgenomes`.
#'
#' @export
resolveSpeciesAssembly <- function(config, donor, recipient,
                                   donor_subgenomes = NULL,
                                   recipient_subgenomes = NULL,
                                   donor_subgenome = NULL,
                                   recipient_subgenome = NULL,
                                   key = NULL,
                                   label = NULL, description = NULL,
                                   default_sample_composition = NULL) {
    if (is.character(config)) {
        config <- readMaMaMIAConfig(config)
    }
    if (!is.list(config) || is.null(config$species)) {
        .config_error("`config` must be a validated MaMaMIA configuration")
    }
    if (length(config$species) == 0L) {
        .config_error("`config` has no `species` registry")
    }
    donor <- .config_string(donor, "donor")
    recipient <- .config_string(recipient, "recipient")
    if (!donor %in% names(config$species)) {
        .config_error(
            "Unknown donor species '%s'; configured species: %s",
            donor, paste(names(config$species), collapse = ", ")
        )
    }
    if (!recipient %in% names(config$species)) {
        .config_error(
            "Unknown recipient species '%s'; configured species: %s",
            recipient, paste(names(config$species), collapse = ", ")
        )
    }

    built <- .build_species_assembly(
        recipient_species = config$species[[recipient]],
        donor_species = config$species[[donor]],
        key = key,
        label = label,
        description = description,
        recipient_subgenomes = recipient_subgenomes,
        donor_subgenomes = donor_subgenomes,
        recipient_subgenome = recipient_subgenome,
        donor_subgenome = donor_subgenome,
        default_sample_composition = default_sample_composition,
        allowed_windows = config$defaults$allowed_windows
    )

    if (is.null(config$assemblies)) {
        config$assemblies <- list()
    }
    existing <- config$assemblies[[built$key]]
    if (!is.null(existing) && !identical(existing, built$assembly)) {
        .config_error(
            "Assembly key '%s' already exists with a different definition; pass `key` to choose another key",
            built$key
        )
    }
    config$assemblies[[built$key]] <- built$assembly

    list(
        key = built$key,
        assembly = built$assembly,
        config = config,
        unpaired = built$unpaired,
        donor_species = donor,
        recipient_species = recipient,
        donor_subgenomes = built$donor_subgenomes,
        recipient_subgenomes = built$recipient_subgenomes,
        donor_subgenome = built$donor_subgenome,
        recipient_subgenome = built$recipient_subgenome
    )
}
