## Match a BAM header against the YAML assembly registry (T1.3).
##
## A match is exact only when every configured donor and recipient target is
## present with the expected length. Extra unplaced contigs in the BAM header
## are ignored. When several assemblies match exactly, the most specific one
## (largest target set) wins; a remaining tie is ambiguous and stops preflight.

.assembly_targets_for_match <- function(assembly, config_key = NULL) {
    key <- if (is.null(config_key)) {
        if (is.null(assembly$key)) "<unknown>" else assembly$key
    } else {
        config_key
    }
    if (!is.list(assembly) || !is.list(assembly$contigs)) {
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "Assembly '%s' has no contig definitions", key
        )
    }
    pieces <- lapply(c("recipient", "donor"), function(role) {
        contigs <- assembly$contigs[[role]]
        if (!is.data.frame(contigs) || nrow(contigs) == 0L) {
            .mamamia_error(
                "MaMaMIAAssemblyMatchError",
                "Assembly '%s' must configure at least one %s contig",
                key, role
            )
        }
        required <- c("id", "name", "subgenome", "length")
        if (!all(required %in% names(contigs))) {
            .mamamia_error(
                "MaMaMIAAssemblyMatchError",
                "Assembly '%s' %s contigs must have columns: %s",
                key, role, paste(required, collapse = ", ")
            )
        }
        data.frame(
            chr_id = as.character(contigs$id),
            chr_name = as.character(contigs$name),
            subgenome = as.character(contigs$subgenome),
            length = as.numeric(contigs$length),
            role = role,
            stringsAsFactors = FALSE
        )
    })
    targets <- do.call(rbind, pieces)
    rownames(targets) <- NULL
    if (anyDuplicated(targets$chr_id)) {
        duplicate <- unique(targets$chr_id[duplicated(targets$chr_id)])[[1L]]
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "Assembly '%s' configures contig '%s' more than once",
            key, duplicate
        )
    }
    if (anyNA(targets$length) || any(targets$length <= 0)) {
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "Assembly '%s' has a contig with a missing or non-positive length",
            key
        )
    }
    targets
}

.match_candidate <- function(header_targets, assembly, config_key,
                             header_path) {
    targets <- .assembly_targets_for_match(assembly, config_key)
    position <- match(targets$chr_id, header_targets$name)
    present <- !is.na(position)
    observed <- rep(NA_real_, nrow(targets))
    observed[present] <- header_targets$length[position[present]]
    mismatch <- present & observed != targets$length

    missing <- targets$chr_id[!present]
    length_mismatches <- data.frame(
        chr_id = targets$chr_id[mismatch],
        role = targets$role[mismatch],
        observed = observed[mismatch],
        expected = targets$length[mismatch],
        stringsAsFactors = FALSE
    )
    n_targets <- nrow(targets)
    n_matched <- sum(present & !mismatch)
    exact <- length(missing) == 0L && nrow(length_mismatches) == 0L

    detail <- character(0)
    if (length(missing) > 0L) {
        detail <- c(
            detail,
            paste0("missing: ", paste(missing, collapse = ", "))
        )
    }
    if (nrow(length_mismatches) > 0L) {
        detail <- c(
            detail,
            paste0(
                "length mismatches: ",
                paste(
                    sprintf(
                        "%s (observed %s, expected %s)",
                        length_mismatches$chr_id,
                        format(length_mismatches$observed, scientific = FALSE),
                        format(length_mismatches$expected, scientific = FALSE)
                    ),
                    collapse = ", "
                )
            )
        )
    }
    summary <- if (exact) {
        sprintf(
            "BAM '%s' matches assembly '%s' exactly (%d/%d targets)",
            header_path, config_key, n_matched, n_targets
        )
    } else {
        sprintf(
            "BAM '%s' matches %d/%d targets of assembly '%s'%s",
            header_path, n_matched, n_targets, config_key,
            if (length(detail) == 0L) "" else paste0("; ", paste(detail, collapse = "; "))
        )
    }

    list(
        assembly = config_key,
        assembly_config = assembly,
        header_path = header_path,
        exact = exact,
        missing = missing,
        length_mismatches = length_mismatches,
        score = n_matched / n_targets,
        n_targets = n_targets,
        n_matched = n_matched,
        summary = summary
    )
}

.candidate_table <- function(candidates) {
    out <- data.frame(
        key = vapply(candidates, `[[`, character(1), "assembly"),
        label = vapply(candidates, function(x) {
            label <- x$assembly_config$label
            if (is.null(label)) "" else as.character(label)[[1L]]
        }, character(1)),
        n_targets = vapply(candidates, `[[`, integer(1), "n_targets"),
        n_matched = vapply(candidates, `[[`, integer(1), "n_matched"),
        n_missing = vapply(candidates, function(x) {
            as.integer(length(x$missing))
        }, integer(1)),
        n_length_mismatches = vapply(candidates, function(x) {
            as.integer(nrow(x$length_mismatches))
        }, integer(1)),
        score = vapply(candidates, `[[`, numeric(1), "score"),
        exact = vapply(candidates, `[[`, logical(1), "exact"),
        stringsAsFactors = FALSE
    )
    out[order(-out$score, -out$n_targets, out$key), , drop = FALSE]
}

.public_match <- function(candidate, candidates_table) {
    list(
        assembly = candidate$assembly,
        assembly_config = candidate$assembly_config,
        header_path = candidate$header_path,
        candidates = candidates_table,
        matched = candidate$exact,
        exact = candidate$exact,
        missing = candidate$missing,
        length_mismatches = candidate$length_mismatches,
        score = candidate$score,
        summary = candidate$summary
    )
}

.check_header_match <- function(candidate, header_path) {
    if (nrow(candidate$length_mismatches) > 0L) {
        first <- candidate$length_mismatches[1L, , drop = FALSE]
        .mamamia_error(
            "MaMaMIAHeaderMismatchError",
            "BAM '%s' target '%s' has length %s; assembly '%s' expects %s",
            header_path, first$chr_id,
            format(first$observed, scientific = FALSE),
            candidate$assembly,
            format(first$expected, scientific = FALSE)
        )
    }
    if (length(candidate$missing) > 0L) {
        donor <- candidate$assembly_config$contigs$donor
        recipient <- candidate$assembly_config$contigs$recipient
        header_names <- candidate$header_target_names
        donor_complete <- all(donor$id %in% header_names)
        recipient_complete <- all(recipient$id %in% header_names)
        if (xor(donor_complete, recipient_complete)) {
            .mamamia_error(
                "MaMaMIAAssemblyMatchError",
                "BAM '%s' covers only one side of assembly '%s' (missing: %s)",
                header_path, candidate$assembly,
                paste(candidate$missing, collapse = ", ")
            )
        }
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "BAM '%s' is missing target(s) required by assembly '%s': %s",
            header_path, candidate$assembly,
            paste(candidate$missing, collapse = ", ")
        )
    }
    if (!candidate$exact) {
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "BAM '%s' does not match assembly '%s': %s",
            header_path, candidate$assembly, candidate$summary
        )
    }
    invisible(candidate)
}

#' Match a BAM header to a configured assembly
#'
#' Compares the target names and lengths in a BAM header with every assembly
#' in a validated MaMaMIA configuration. A candidate is exact only when all
#' configured donor and recipient targets are present with the expected
#' lengths; extra header contigs are ignored.
#'
#' If `assembly` is supplied, that assembly is validated instead of being
#' auto-selected. Without it, exactly one exact match is selected; when several
#' assemblies match exactly, the one with the largest target set (the most
#' specific) wins, and a remaining tie raises an ambiguity error.
#'
#' @param header A `BamHeaderInfo` list from [bamHeaderInfo()], a path to a BAM,
#'   or a list with a `targets` data frame.
#' @param config A validated configuration list, or a path to a YAML config.
#' @param assembly Optional assembly key to validate instead of auto-selecting.
#'
#' @return A list with `assembly` (the selected key), `assembly_config`,
#'   `candidates` (one row per evaluated assembly), `matched`, `exact`,
#'   `missing`, `length_mismatches`, `score` and a human-readable `summary`.
#'
#' @examples
#' demo <- system.file("extdata", "bam-demo", package = "MaMaMIA")
#' config <- readMaMaMIAConfig(file.path(demo, "assemblies.yml"))
#' pair <- resolveSpeciesAssembly(
#'     config, donor = "secale_like", recipient = "triticum_like"
#' )
#' header <- bamHeaderInfo(file.path(demo, "demo.bam"))
#' match <- matchAssembly(header, pair$config)
#' match$assembly
#'
#' @export
matchAssembly <- function(header, config, assembly = NULL) {
    if (is.character(header)) {
        header <- bamHeaderInfo(header)
    }
    if (!is.list(header) || is.null(header$targets) ||
        !is.data.frame(header$targets) ||
        !all(c("name", "length") %in% names(header$targets))) {
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "`header` must be a BamHeaderInfo object with a targets data frame"
        )
    }
    if (is.character(config)) {
        config <- readMaMaMIAConfig(config)
    }
    if (!is.list(config) || !is.list(config$assemblies) ||
        length(config$assemblies) == 0L) {
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "`config` must be a validated MaMaMIA configuration with at least one assembly"
        )
    }
    header_path <- if (is.null(header$path)) "<unknown BAM>" else header$path
    header_targets <- data.frame(
        name = as.character(header$targets$name),
        length = as.numeric(header$targets$length),
        stringsAsFactors = FALSE
    )
    if (anyNA(header_targets$length) ||
        any(!is.finite(header_targets$length)) ||
        any(header_targets$length < 0)) {
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "BAM '%s' has a missing or invalid target length",
            header_path
        )
    }
    if (anyDuplicated(header_targets$name)) {
        duplicate <- unique(header_targets$name[
            duplicated(header_targets$name)
        ])[[1L]]
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "BAM '%s' has duplicate target sequence name '%s'",
            header_path, duplicate
        )
    }

    candidates <- lapply(seq_along(config$assemblies), function(i) {
        candidate <- .match_candidate(
            header_targets,
            config$assemblies[[i]],
            names(config$assemblies)[[i]],
            header_path
        )
        candidate$header_target_names <- header_targets$name
        candidate
    })
    names(candidates) <- names(config$assemblies)
    candidates_table <- .candidate_table(candidates)

    if (!is.null(assembly)) {
        if (!is.character(assembly) || length(assembly) != 1L ||
            is.na(assembly) || !nzchar(assembly)) {
            .mamamia_error(
                "MaMaMIAAssemblyMatchError",
                "`assembly` must be a single configured assembly key"
            )
        }
        if (!assembly %in% names(candidates)) {
            .mamamia_error(
                "MaMaMIAAssemblyMatchError",
                "Unknown assembly '%s'; configured assemblies: %s",
                assembly, paste(names(candidates), collapse = ", ")
            )
        }
        candidate <- candidates[[assembly]]
        .check_header_match(candidate, header_path)
        return(.public_match(candidate, candidates_table))
    }

    exact_keys <- names(candidates)[vapply(
        candidates, `[[`, logical(1), "exact"
    )]
    if (length(exact_keys) == 1L) {
        return(.public_match(candidates[[exact_keys]], candidates_table))
    }
    if (length(exact_keys) > 1L) {
        sizes <- vapply(exact_keys, function(key) {
            candidates[[key]]$n_targets
        }, integer(1))
        largest <- exact_keys[sizes == max(sizes)]
        if (length(largest) == 1L) {
            return(.public_match(candidates[[largest]], candidates_table))
        }
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "BAM '%s' matches more than one assembly: %s",
            header_path, paste(sort(largest), collapse = ", "),
            data = list(candidates = candidates_table)
        )
    }

    top <- candidates_table[1L, , drop = FALSE]
    .mamamia_error(
        "MaMaMIAAssemblyMatchError",
        "BAM '%s' does not match any configured assembly; closest match: '%s' (%d/%d targets, %.1f%%%s)",
        header_path, top$key, top$n_matched, top$n_targets,
        100 * top$score,
        if (top$n_length_mismatches > 0L) {
            sprintf(", %d length mismatch(es)", top$n_length_mismatches)
        } else {
            ""
        },
        data = list(candidates = candidates_table)
    )
}

#' Match several BAM headers to one assembly
#'
#' Matches each header with [matchAssembly()] and requires all inputs to select
#' the same assembly key. This is the multi-BAM rule used by preflight.
#'
#' @param headers A list of `BamHeaderInfo` objects or BAM paths.
#' @param config A validated configuration list or YAML path.
#' @param assembly Optional assembly key to validate for every input.
#'
#' @return A list of match results.
#' @noRd
matchAssemblyCohort <- function(headers, config, assembly = NULL) {
    if (is.character(headers)) {
        headers <- as.list(headers)
    }
    if (is.list(headers) && !is.null(headers$targets)) {
        headers <- list(headers)
    }
    if (!is.list(headers) || length(headers) == 0L) {
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "`headers` must contain at least one BAM header or path"
        )
    }
    matches <- lapply(headers, matchAssembly, config = config,
        assembly = assembly
    )
    keys <- vapply(matches, `[[`, character(1), "assembly")
    if (length(unique(keys)) > 1L) {
        detail <- vapply(seq_along(matches), function(i) {
            sprintf("'%s' -> '%s'", matches[[i]]$header_path, keys[[i]])
        }, character(1))
        .mamamia_error(
            "MaMaMIAAssemblyMatchError",
            "BAM inputs match different assemblies: %s",
            paste(detail, collapse = "; ")
        )
    }
    matches
}
