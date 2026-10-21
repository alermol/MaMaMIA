## contigs.tsv parsing and exact-N window generation (T2.1, T2.2).
##
## Window tables are 0-based half-open. GenomicRanges objects are 1-based
## inclusive. Every target contig is tiled into exactly n_per_chromosome
## windows using start_i = floor((i - 1) * L / n), stop_i = floor(i * L / n).

.window_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIAWindowError", fmt, ...)
}

.window_integer <- function(x, what, minimum = 1L) {
    if (length(x) != 1L || is.na(x) || !is.numeric(x) ||
        !is.finite(x) || x != floor(x) || x < minimum ||
        x > .Machine$integer.max) {
        .window_error("`%s` must be a single integer >= %d", what, minimum)
    }
    as.integer(x)
}

#' Read a gcrank `contigs.tsv` file
#'
#' Reads the three-column `name`, `length`, `offset` table shipped with a
#' gcrank index. Contig names must be unique, lengths positive integers and
#' offsets the cumulative lengths of the preceding contigs (the index
#' concatenates contigs without gaps).
#'
#' @param path Path to `contigs.tsv`.
#'
#' @return A data frame with `name`, `length` and `offset`.
#' @export
readContigsTsv <- function(path) {
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(path)) {
        .window_error("`path` must be a single non-empty file path")
    }
    path <- path.expand(path)
    if (!file.exists(path)) {
        .window_error("Contigs file does not exist: %s", path)
    }

    lines <- tryCatch(
        readLines(path, warn = FALSE),
        error = function(e) {
            .window_error("Could not read %s: %s", path, conditionMessage(e))
        }
    )
    lines <- sub("\r$", "", lines)
    lines <- lines[nzchar(trimws(lines))]
    if (length(lines) == 0L) {
        .window_error("No contigs found in %s", path)
    }

    fields <- strsplit(lines, "\t", fixed = TRUE)
    widths <- lengths(fields)
    if (any(widths != 3L)) {
        row <- which(widths != 3L)[[1L]]
        .window_error(
            "%s line %d has %d tab-separated columns; expected exactly 3 (name, length, offset)",
            path, row, widths[[row]]
        )
    }

    name <- vapply(fields, `[[`, character(1), 1L)
    length_text <- vapply(fields, `[[`, character(1), 2L)
    offset_text <- vapply(fields, `[[`, character(1), 3L)

    if (anyNA(name) || any(!nzchar(trimws(name)))) {
        .window_error("Empty contig name in %s", path)
    }
    length <- suppressWarnings(as.numeric(length_text))
    offset <- suppressWarnings(as.numeric(offset_text))
    if (anyNA(length) || any(!is.finite(length))) {
        row <- which(is.na(length) | !is.finite(length))[[1L]]
        .window_error(
            "Contig length on line %d of %s is not numeric: '%s'",
            row, path, length_text[[row]]
        )
    }
    if (anyNA(offset) || any(!is.finite(offset))) {
        row <- which(is.na(offset) | !is.finite(offset))[[1L]]
        .window_error(
            "Contig offset on line %d of %s is not numeric: '%s'",
            row, path, offset_text[[row]]
        )
    }
    if (any(length <= 0) || any(length != floor(length))) {
        .window_error("Contig lengths in %s must be positive integers", path)
    }
    if (any(offset < 0) || any(offset != floor(offset))) {
        .window_error("Contig offsets in %s must be non-negative integers", path)
    }
    if (anyDuplicated(name)) {
        duplicate <- unique(name[duplicated(name)])[[1L]]
        .window_error("Duplicate contig '%s' in %s", duplicate, path)
    }

    expected_offsets <- c(0, utils::head(cumsum(length), -1L))
    if (any(offset != expected_offsets)) {
        row <- which(offset != expected_offsets)[[1L]]
        .window_error(
            "Contig offsets in %s must be the cumulative lengths of the preceding contigs (first mismatch on line %d)",
            path, row
        )
    }

    data.frame(
        name = name,
        length = as.numeric(length),
        offset = as.numeric(offset),
        stringsAsFactors = FALSE
    )
}

.normalize_contigs_table <- function(contigs, source_label) {
    if (is.character(contigs)) {
        if (length(contigs) != 1L) {
            .window_error("`contigs` must be one contigs.tsv path or a data frame")
        }
        return(readContigsTsv(contigs))
    }
    if (!is.data.frame(contigs)) {
        .window_error(
            "`contigs` must be a contigs.tsv path or a data frame with name and length"
        )
    }
    if (!"length" %in% names(contigs)) {
        .window_error("%s must have a `length` column", source_label)
    }
    if ("name" %in% names(contigs)) {
        name <- contigs$name
    } else if ("id" %in% names(contigs)) {
        name <- contigs$id
    } else {
        .window_error("%s must have a `name` or `id` column", source_label)
    }
    name <- as.character(name)
    length <- suppressWarnings(as.numeric(contigs$length))
    if (anyNA(name) || any(!nzchar(name))) {
        .window_error("%s contains an empty contig name", source_label)
    }
    if (anyNA(length) || any(!is.finite(length)) || any(length <= 0)) {
        .window_error("%s contains a missing or non-positive contig length", source_label)
    }
    if (anyDuplicated(name)) {
        duplicate <- unique(name[duplicated(name)])[[1L]]
        .window_error("%s contains duplicate contig '%s'", source_label, duplicate)
    }
    data.frame(
        name = name,
        length = length,
        stringsAsFactors = FALSE
    )
}

#' Flatten and validate assembly target contigs
#'
#' Returns one row per configured target in assembly order: recipient contigs
#' first, then donor contigs. When a contigs table is supplied (one or more
#' gcrank indexes concatenated), every target must be present exactly once with
#' the same length as in the configuration.
#'
#' @param assembly An `AssemblyConfig` list.
#' @param contigs Optional contigs table or path.
#'
#' @return A data frame with `chr_id`, `chr_name`, `subgenome`, `length`,
#'   `role` and `index_path`.
#' @noRd
.assembly_target_contigs <- function(assembly, contigs = NULL) {
    if (!is.list(assembly) || !is.list(assembly$contigs)) {
        .window_error("`assembly` must be an AssemblyConfig with contigs")
    }
    key <- if (is.null(assembly$key)) "<unknown>" else assembly$key
    pieces <- lapply(c("recipient", "donor"), function(role) {
        x <- assembly$contigs[[role]]
        if (!is.data.frame(x) || nrow(x) == 0L) {
            .window_error(
                "Assembly '%s' must configure at least one %s contig",
                key, role
            )
        }
        required <- c("id", "name", "subgenome", "length")
        if (!all(required %in% names(x))) {
            .window_error(
                "Assembly '%s' %s contigs must have columns: %s",
                key, role, paste(required, collapse = ", ")
            )
        }
        data.frame(
            chr_id = as.character(x$id),
            chr_name = as.character(x$name),
            subgenome = as.character(x$subgenome),
            length = as.numeric(x$length),
            role = role,
            index_path = if ("index_path" %in% names(x)) {
                as.character(x$index_path)
            } else {
                rep(NA_character_, nrow(x))
            },
            stringsAsFactors = FALSE
        )
    })
    targets <- do.call(rbind, pieces)
    rownames(targets) <- NULL
    if (anyDuplicated(targets$chr_id)) {
        duplicate <- unique(targets$chr_id[duplicated(targets$chr_id)])[[1L]]
        .window_error(
            "Assembly '%s' configures contig '%s' in more than one role",
            key, duplicate
        )
    }
    if (anyNA(targets$length) || any(targets$length <= 0)) {
        .window_error(
            "Assembly '%s' has a contig with a missing or non-positive length",
            key
        )
    }

    if (!is.null(contigs)) {
        source_label <- if (is.character(contigs) && length(contigs) == 1L) {
            sprintf("contigs.tsv '%s'", contigs)
        } else {
            "the supplied contigs table"
        }
        table <- .normalize_contigs_table(contigs, source_label)
        position <- match(targets$chr_id, table$name)
        missing <- targets$chr_id[is.na(position)]
        if (length(missing) > 0L) {
            .window_error(
                "Assembly '%s': target contig(s) not found in %s: %s",
                key, source_label, paste(missing, collapse = ", ")
            )
        }
        observed <- table$length[position]
        mismatch <- which(observed != targets$length)
        if (length(mismatch) > 0L) {
            first <- mismatch[[1L]]
            .window_error(
                "Assembly '%s': contig '%s' is %s bp in %s but %s bp in the assembly config",
                key, targets$chr_id[[first]],
                format(observed[[first]], scientific = FALSE),
                source_label,
                format(targets$length[[first]], scientific = FALSE)
            )
        }
    }
    targets
}

#' Build exact-N analysis windows
#'
#' Tiles every target contig into exactly `n_per_chromosome` contiguous
#' half-open windows. Windows are returned in assembly target order (recipient
#' contigs first, then donor), never in lexical chromosome order.
#'
#' @param assembly An `AssemblyConfig` list from a validated configuration.
#' @param contigs A `contigs.tsv` path or a data frame with `name` and `length`
#'   (for example [readContigsTsv()], optionally concatenated across indexes).
#' @param n_per_chromosome Number of windows per target contig.
#'
#' @return A list with `gr` (1-based inclusive `GRanges`), `table` (0-based
#'   half-open data frame with `chr_id`, `start`, `stop`, `role`, `chr_name`,
#'   `subgenome`) and `n_per_chromosome`.
#'
#' @examples
#' demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
#' config <- readMaMaMIAConfig(file.path(demo, "assemblies.yml"))
#' pair <- resolveSpeciesAssembly(
#'     config, donor = "secale_like", recipient = "triticum_like"
#' )
#' contigs <- readContigsTsv(file.path(demo, "demo.gcrank", "contigs.tsv"))
#' windows <- makeAnalysisWindows(pair$assembly, contigs, n_per_chromosome = 100L)
#' nrow(windows$table)
#' table(windows$table$role)
#'
#' @export
makeAnalysisWindows <- function(assembly, contigs, n_per_chromosome = 1000L) {
    n <- .window_integer(n_per_chromosome, "n_per_chromosome")
    targets <- .assembly_target_contigs(assembly, contigs)
    key <- if (is.null(assembly$key)) "<unknown>" else assembly$key

    rows <- vector("list", nrow(targets))
    for (i in seq_len(nrow(targets))) {
        length_bp <- targets$length[[i]]
        chr_id <- targets$chr_id[[i]]
        if (length_bp < n) {
            .window_error(
                "Assembly '%s': contig '%s' is %s bp, shorter than the requested %d windows",
                key, chr_id,
                format(length_bp, scientific = FALSE), n
            )
        }
        index <- seq_len(n)
        start <- floor((index - 1) * length_bp / n)
        stop <- floor(index * length_bp / n)
        if (any(stop <= start) || start[[1L]] != 0 ||
            stop[[n]] != length_bp) {
            .window_error(
                "Assembly '%s': could not tile contig '%s' (%s bp) into %d non-empty windows",
                key, chr_id, format(length_bp, scientific = FALSE), n
            )
        }
        rows[[i]] <- data.frame(
            chr_id = rep(chr_id, n),
            start = as.integer(start),
            stop = as.integer(stop),
            role = rep(targets$role[[i]], n),
            chr_name = rep(targets$chr_name[[i]], n),
            subgenome = rep(targets$subgenome[[i]], n),
            stringsAsFactors = FALSE
        )
    }
    win_table <- do.call(rbind, rows)
    rownames(win_table) <- NULL

    for (i in seq_len(nrow(targets))) {
        idx <- which(win_table$chr_id == targets$chr_id[[i]])
        if (length(idx) != n || idx[[1L]] != (i - 1L) * n + 1L) {
            .window_error(
                "Assembly '%s': contig '%s' does not have exactly %d contiguous windows",
                key, targets$chr_id[[i]], n
            )
        }
        start <- win_table$start[idx]
        stop <- win_table$stop[idx]
        if (start[[1L]] != 0 ||
            stop[[n]] != targets$length[[i]] ||
            any(start[-1L] != stop[-n]) ||
            any(stop <= start)) {
            .window_error(
                "Assembly '%s': windows for contig '%s' are not a contiguous half-open tiling",
                key, targets$chr_id[[i]]
            )
        }
    }

    ids <- unique(win_table$chr_id)
    lengths <- targets$length[match(ids, targets$chr_id)]
    gr <- GenomicRanges::GRanges(
        seqnames = win_table$chr_id,
        ranges = IRanges::IRanges(
            start = win_table$start + 1L,
            end = win_table$stop
        ),
        seqlengths = stats::setNames(as.numeric(lengths), ids)
    )

    list(
        gr = gr,
        table = win_table,
        n_per_chromosome = n
    )
}
