## Per-sample tables for the RCA constructor (T3.2).
##
## `buildSampleTables()` is the single place where the generated windows, the
## GC table and the read counts are joined into the `SampleTables` contract.
## Windows are 0-based half-open; counts and GC values must be aligned to
## `windows$table` element by element. The function validates the resulting
## tables through the same `RCA()` constructor used by the legacy workflow, so
## a malformed sample fails before the expensive GC-bias fit.

.sample_tables_error <- function(fmt, ...) {
    .mamamia_error("MaMaMIASampleError", fmt, ...)
}

.sample_string <- function(x, what) {
    if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
        .sample_tables_error("`%s` must be a single non-empty string", what)
    }
    x
}

## Normalize the AnalysisWindows input to a plain table plus the window count.
.sample_window_table <- function(windows) {
    if (is.list(windows) && !is.data.frame(windows) &&
        is.data.frame(windows$table)) {
        table <- windows$table
        declared <- windows$n_per_chromosome
    } else if (is.data.frame(windows)) {
        table <- windows
        declared <- NULL
    } else {
        .mamamia_error(
            "MaMaMIAWindowError",
            "`windows` must be an AnalysisWindows list or a data frame"
        )
    }

    required <- c("chr_id", "start", "stop", "role", "chr_name", "subgenome")
    missing <- setdiff(required, names(table))
    if (length(missing) > 0L) {
        .mamamia_error(
            "MaMaMIAWindowError",
            "The window table must have columns: %s",
            paste(required, collapse = ", ")
        )
    }
    if (nrow(table) == 0L) {
        .mamamia_error("MaMaMIAWindowError", "The window table is empty")
    }

    table <- as.data.frame(table, stringsAsFactors = FALSE)
    for (column in c("chr_id", "role", "chr_name", "subgenome")) {
        if (is.factor(table[[column]])) {
            table[[column]] <- as.character(table[[column]])
        }
        table[[column]] <- as.character(table[[column]])
        if (anyNA(table[[column]]) || any(!nzchar(table[[column]]))) {
            .mamamia_error(
                "MaMaMIAWindowError",
                "The window table has an empty `%s` value", column
            )
        }
    }

    for (column in c("start", "stop")) {
        value <- suppressWarnings(as.numeric(table[[column]]))
        if (anyNA(value) || any(!is.finite(value))) {
            .mamamia_error(
                "MaMaMIAWindowError",
                "The window table `%s` column must be numeric", column
            )
        }
        if (any(value != floor(value))) {
            .mamamia_error(
                "MaMaMIAWindowError",
                "The window table `%s` column must contain whole numbers",
                column
            )
        }
        table[[column]] <- value
    }
    if (any(table$start < 0) || any(table$stop <= table$start)) {
        bad <- which(table$start < 0 | table$stop <= table$start)[[1L]]
        .mamamia_error(
            "MaMaMIAWindowError",
            "Invalid window %s:%s-%s: start must be non-negative and smaller than stop",
            table$chr_id[[bad]], table$start[[bad]], table$stop[[bad]]
        )
    }
    if (!all(table$role %in% c("donor", "recipient"))) {
        bad <- which(!table$role %in% c("donor", "recipient"))[[1L]]
        .mamamia_error(
            "MaMaMIAWindowError",
            "Unknown window role '%s' for contig '%s'; expected donor or recipient",
            table$role[[bad]], table$chr_id[[bad]]
        )
    }

    ## A contig belongs to exactly one role and has one stable display name and
    ## subgenome label. Overlapping windows are left to RCA(), which validates
    ## the tiling with the same rules as the legacy workflow.
    per_contig <- split(
        seq_len(nrow(table)),
        factor(table$chr_id, levels = unique(table$chr_id))
    )
    for (idx in per_contig) {
        first <- idx[[1L]]
        if (length(unique(table$role[idx])) != 1L ||
            length(unique(table$chr_name[idx])) != 1L ||
            length(unique(table$subgenome[idx])) != 1L) {
            .mamamia_error(
                "MaMaMIAWindowError",
                "Contig '%s' must have one role, name and subgenome across its windows",
                table$chr_id[[first]]
            )
        }
    }

    window_counts <- lengths(per_contig)
    if (length(unique(window_counts)) != 1L) {
        .mamamia_error(
            "MaMaMIAWindowError",
            "All contigs in the window table must have the same number of windows"
        )
    }
    n_per <- as.integer(window_counts[[1L]])
    if (!is.null(declared)) {
        if (length(declared) != 1L || is.na(declared) ||
            !is.numeric(declared) || declared != floor(declared) ||
            as.integer(declared) != n_per) {
            .mamamia_error(
                "MaMaMIAWindowError",
                "`windows$n_per_chromosome` does not match the window table"
            )
        }
    }

    list(table = table, n_per_chromosome = n_per)
}

## Match a per-window value vector to the window table.
.sample_values <- function(values, table, what) {
    n <- nrow(table)
    if (is.data.frame(values)) {
        required <- c("chr_id", "start", "stop", what)
        missing <- setdiff(required, names(values))
        if (length(missing) > 0L) {
            .sample_tables_error(
                "`%s` must have columns: %s", what,
                paste(required, collapse = ", ")
            )
        }
        if (nrow(values) != n) {
            .sample_tables_error(
                "`%s` has %d rows but the window table has %d",
                what, nrow(values), n
            )
        }
        aligned <- identical(as.character(values$chr_id), as.character(table$chr_id)) &&
            isTRUE(all.equal(as.numeric(values$start), as.numeric(table$start))) &&
            isTRUE(all.equal(as.numeric(values$stop), as.numeric(table$stop)))
        if (!aligned) {
            .sample_tables_error(
                "`%s` is not aligned to the window table (chr_id/start/stop differ)",
                what
            )
        }
        values <- values[[what]]
    }

    if (is.list(values) && !is.data.frame(values)) {
        out <- rep(NA_real_, n)
        for (role in c("donor", "recipient")) {
            idx <- which(table$role == role)
            part <- values[[role]]
            if (is.matrix(part) || is.data.frame(part)) {
                if (ncol(part) != 1L) {
                    .sample_tables_error(
                        "`%s$%s` must be a vector or a one-column matrix",
                        what, role
                    )
                }
                part <- part[, 1L]
            }
            if (!is.numeric(part) || length(part) != length(idx)) {
                .sample_tables_error(
                    "`%s$%s` must have %d values (one per %s window)",
                    what, role, length(idx), role
                )
            }
            out[idx] <- as.numeric(part)
        }
        if (anyNA(out)) {
            .sample_tables_error(
                "`%s` must provide donor and recipient values", what
            )
        }
        values <- out
    }

    if (!is.numeric(values) || length(values) != n) {
        .sample_tables_error(
            "`%s` must be a numeric vector of length %d (one value per window)",
            what, n
        )
    }
    values <- as.numeric(values)
    if (anyNA(values) || any(!is.finite(values))) {
        .sample_tables_error("`%s` must not contain missing values", what)
    }
    values
}

.sample_meta <- function(table) {
    first <- !duplicated(table$chr_id)
    meta <- data.frame(
        chr_id = table$chr_id[first],
        chr_name = table$chr_name[first],
        subgenome = table$subgenome[first],
        stringsAsFactors = FALSE
    )
    rownames(meta) <- NULL
    meta
}

## Use the explicit assembly pairs when present; otherwise fall back to the
## chromosome-number rule used by the legacy workflow.
.sample_target_pairs <- function(assembly, meta, don_ids, rec_ids, sample) {
    pairs <- assembly$pairs
    if (!is.null(pairs) && length(pairs) > 0L) {
        flat <- unlist(pairs, use.names = TRUE)
        if (is.null(names(flat)) || anyNA(names(flat)) ||
            any(!nzchar(names(flat)))) {
            .sample_tables_error(
                "Sample '%s': assembly pairs must be a named donor-to-recipient mapping",
                sample
            )
        }
        values <- as.character(flat)
        names(values) <- names(flat)
        if (anyNA(values) || any(!nzchar(values))) {
            .sample_tables_error(
                "Sample '%s': assembly pairs must not contain empty values",
                sample
            )
        }
        if (anyDuplicated(names(values))) {
            .sample_tables_error(
                "Sample '%s': assembly pairs has a duplicate donor", sample
            )
        }
        if (anyDuplicated(unname(values))) {
            .sample_tables_error(
                "Sample '%s': assembly pairs is not one-to-one", sample
            )
        }
        unknown_donor <- setdiff(names(values), don_ids)
        if (length(unknown_donor) > 0L) {
            .sample_tables_error(
                "Sample '%s': assembly pairs refers to unknown donor contig(s): %s",
                sample, paste(unknown_donor, collapse = ", ")
            )
        }
        unknown_recipient <- setdiff(unname(values), rec_ids)
        if (length(unknown_recipient) > 0L) {
            .sample_tables_error(
                "Sample '%s': assembly pairs refers to unknown recipient contig(s): %s",
                sample, paste(unknown_recipient, collapse = ", ")
            )
        }
        return(values)
    }

    derived <- tryCatch(
        defaultTargetPairs(meta, don_ids, rec_ids),
        error = function(e) e
    )
    if (inherits(derived, "error") || length(derived) == 0L) {
        .sample_tables_error(
            "Sample '%s': no explicit assembly pairs and no donor/recipient pair could be derived by chromosome number%s",
            sample,
            if (inherits(derived, "error")) {
                paste0(": ", conditionMessage(derived))
            } else {
                ""
            }
        )
    }
    derived
}

#' Build the per-sample tables consumed by [RCA()]
#'
#' Joins generated windows, GC content and read counts into the `SampleTables`
#' contract. Donor and recipient windows are split by `windows$table$role`;
#' each coverage table has exactly four columns (`chr_id`, `start`, `stop`,
#' `cov`) and each GC table exactly four columns (`chr_id`, `start`, `stop`,
#' `gc`). The result is validated by running [RCA()], which must return an
#' uncorrected object.
#'
#' Donor-recipient pairs come from `assembly$pairs` when configured. Otherwise
#' [defaultTargetPairs()] derives them from the chromosome numbers of the
#' display names.
#'
#' @param counts Read counts aligned to `windows$table`, either a numeric
#'   vector (one value per window) or a list with `donor` and `recipient`
#'   vectors in window order.
#' @param gc GC content aligned to `windows$table`, either a data frame with
#'   columns `chr_id`, `start`, `stop`, `gc` (as returned by the internal
#'   `extractWindowGc()`) or a numeric vector (one value per window).
#' @param windows An `AnalysisWindows` list from [makeAnalysisWindows()], or a
#'   window data frame with columns `chr_id`, `start`, `stop`, `role`,
#'   `chr_name` and `subgenome`.
#' @param assembly An `AssemblyConfig` list.
#' @param sample Sample name.
#'
#' @return A `SampleTables` list with `dcov`, `dgc`, `rcov`, `rgc`, `meta`,
#'   `don_chr_ids`, `rec_chr_ids`, `sample`, `reference`, `window` and
#'   `target_pairs`.
#' @export
buildSampleTables <- function(counts, gc, windows, assembly, sample) {
    sample <- .sample_string(sample, "sample")
    normalized <- .sample_window_table(windows)
    table <- normalized$table
    n_per <- normalized$n_per_chromosome

    if (!is.list(assembly) || !is.list(assembly$contigs)) {
        .sample_tables_error(
            "Sample '%s': `assembly` must be an AssemblyConfig with contigs",
            sample
        )
    }
    targets <- .assembly_target_contigs(assembly)
    position <- match(table$chr_id, targets$chr_id)
    if (anyNA(position)) {
        unknown <- unique(table$chr_id[is.na(position)])[[1L]]
        .sample_tables_error(
            "Sample '%s': window contig '%s' is not configured in assembly '%s'",
            sample, unknown,
            if (is.null(assembly$key)) "<unknown>" else assembly$key
        )
    }
    if (any(table$role != targets$role[position])) {
        bad <- which(table$role != targets$role[position])[[1L]]
        .sample_tables_error(
            "Sample '%s': contig '%s' is a %s target but appears as %s windows",
            sample, table$chr_id[[bad]], targets$role[position[[bad]]],
            table$role[[bad]]
        )
    }

    values <- .sample_values(counts, table, "cov")
    if (any(values < 0)) {
        .sample_tables_error(
            "Sample '%s': read counts must be non-negative", sample
        )
    }
    if (any(abs(values - round(values)) > 1e-6)) {
        .sample_tables_error(
            "Sample '%s': read counts must be whole numbers", sample
        )
    }
    gc_values <- .sample_values(gc, table, "gc")
    if (any(gc_values < 0 | gc_values > 1)) {
        .sample_tables_error(
            "Sample '%s': GC content must be inside [0, 1]", sample
        )
    }

    donor <- table$role == "donor"
    recipient <- table$role == "recipient"
    if (!any(donor)) {
        .sample_tables_error(
            "Sample '%s': the window table has no donor windows", sample
        )
    }
    if (!any(recipient)) {
        .sample_tables_error(
            "Sample '%s': the window table has no recipient windows", sample
        )
    }

    make_table <- function(idx, value, column) {
        out <- data.frame(
            chr_id = table$chr_id[idx],
            start = as.integer(table$start[idx]),
            stop = as.integer(table$stop[idx]),
            value = as.numeric(value[idx]),
            stringsAsFactors = FALSE
        )
        names(out) <- c("chr_id", "start", "stop", column)
        rownames(out) <- NULL
        out
    }

    dcov <- make_table(donor, values, "cov")
    dgc <- make_table(donor, gc_values, "gc")
    rcov <- make_table(recipient, values, "cov")
    rgc <- make_table(recipient, gc_values, "gc")
    meta <- .sample_meta(table)

    don_chr_ids <- unique(dcov$chr_id)
    rec_chr_ids <- unique(rcov$chr_id)
    shared <- intersect(don_chr_ids, rec_chr_ids)
    if (length(shared) > 0L) {
        .sample_tables_error(
            "Sample '%s': contig(s) appear as both donor and recipient: %s",
            sample, paste(shared, collapse = ", ")
        )
    }

    target_pairs <- .sample_target_pairs(
        assembly, meta, don_chr_ids, rec_chr_ids, sample
    )

    tables <- list(
        dcov = dcov,
        dgc = dgc,
        rcov = rcov,
        rgc = rgc,
        meta = meta
    )
    rca <- tryCatch(
        do.call(RCA, tables),
        error = function(e) {
            .sample_tables_error(
                "Sample '%s': could not build an RCA object: %s",
                sample, conditionMessage(e)
            )
        }
    )
    if (!identical(rca$corrected, FALSE)) {
        .sample_tables_error(
            "Sample '%s': buildSampleTables() must produce an uncorrected RCA object",
            sample
        )
    }

    structure(
        c(
            tables,
            list(
                don_chr_ids = don_chr_ids,
                rec_chr_ids = rec_chr_ids,
                sample = sample,
                reference = if (is.null(assembly$key)) {
                    NA_character_
                } else {
                    as.character(assembly$key)
                },
                window = as.integer(n_per),
                target_pairs = target_pairs
            )
        ),
        class = c("SampleTables", "list")
    )
}
