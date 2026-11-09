## One runnable analysis configuration.
##
## A registry describes species, not assemblies. Preflight and the pipeline
## need a resolved assembly, so callers otherwise repeat the same dance:
## read the registry, choose a species pair, resolve the assembly, then carry
## both the key and the updated config into every call. An analysis context
## captures that result once. The pipeline still accepts a plain resolved
## configuration for backwards compatibility; when it receives a context it
## also knows which assembly the caller selected.

#' Is this a resolved analysis context?
#'
#' @param x Any object.
#' @return `TRUE` when `x` came from [openAnalysisConfig()].
isMaMaMIAAnalysisContext <- function(x) {
    inherits(x, "MaMaMIAAnalysisContext")
}

#' Configuration carried by an analysis context
#'
#' @param context A context from [openAnalysisConfig()], or a plain resolved
#'   configuration (returned unchanged).
#' @return The resolved configuration list.
analysisContextConfig <- function(context) {
    if (isMaMaMIAAnalysisContext(context)) {
        return(context$config)
    }
    context
}

#' Assembly selected by an analysis context
#'
#' @param context A context from [openAnalysisConfig()].
#' @return The assembly key, or `NULL` for a plain configuration without a
#'   context.
analysisContextAssemblyKey <- function(context) {
    if (isMaMaMIAAnalysisContext(context)) {
        return(context$assembly_key)
    }
    NULL
}

.context_error <- function(fmt, ...) {
    .config_error(fmt, ...)
}

.context_registry <- function(config) {
    if (is.character(config)) {
        config <- readMaMaMIAConfig(config)
    }
    if (!is.list(config) || is.null(config$species)) {
        .context_error("`config` must be a validated MaMaMIA configuration")
    }
    config
}

.context_assembly <- function(config, key) {
    if (length(config$assemblies) == 0L) {
        .context_error(
            "`config` has no resolved assembly; pass `donor` and `recipient`"
        )
    }
    assembly <- config$assemblies[[key]]
    if (is.null(assembly)) {
        .context_error(
            "Unknown assembly '%s'; configured assemblies: %s",
            key, paste(names(config$assemblies), collapse = ", ")
        )
    }
    assembly
}

#' Open a runnable analysis context
#'
#' Resolves a registry and a species pair once, returning the configuration
#' with the derived assembly registered plus the assembly key, labels and
#' components that preflight and [analyzeBamInput()] need. It accepts a YAML
#' path or an already-validated configuration list. A resolved configuration
#' with exactly one assembly, or an explicit `assembly` key, may omit the
#' species pair.
#'
#' @param config Path to a YAML registry or a validated configuration list.
#' @param donor Optional donor species key.
#' @param recipient Optional recipient species key.
#' @param assembly Optional assembly key. When `donor` and `recipient` are
#'   also given it is used as the derived assembly key.
#' @param donor_components Optional donor component labels.
#' @param recipient_components Optional recipient component labels.
#' @param donor_component Optional single donor component (shortcut).
#' @param recipient_component Optional single recipient component (shortcut).
#' @param key Optional assembly key (alias for `assembly` when deriving).
#' @param label Optional assembly label.
#' @param description Optional assembly description.
#' @param default_sample_composition Optional named recipient copy-number
#'   vector for the derived assembly.
#'
#' @return A `MaMaMIAAnalysisContext` list with `config`, `assembly_key`,
#'   `assembly`, `donor_species`, `recipient_species`,
#'   `donor_components`, `recipient_components` and `unpaired`.
openAnalysisConfig <- function(config, donor = NULL, recipient = NULL,
                               assembly = NULL,
                               donor_components = NULL,
                               recipient_components = NULL,
                               donor_component = NULL,
                               recipient_component = NULL,
                               key = NULL,
                               label = NULL, description = NULL,
                               default_sample_composition = NULL) {
    config <- .context_registry(config)
    if (!is.null(key) && !is.null(assembly) && !identical(key, assembly)) {
        .context_error("`key` and `assembly` name different assemblies")
    }
    if (is.null(key)) {
        key <- assembly
    }

    has_donor <- !is.null(donor)
    has_recipient <- !is.null(recipient)
    if (has_donor != has_recipient) {
        .context_error(
            "Pass both `donor` and `recipient`, or neither"
        )
    }

    if (has_donor) {
        resolved <- resolveSpeciesAssembly(
            config,
            donor = donor,
            recipient = recipient,
            donor_components = donor_components,
            recipient_components = recipient_components,
            donor_component = donor_component,
            recipient_component = recipient_component,
            key = key,
            label = label,
            description = description,
            default_sample_composition = default_sample_composition
        )
        context <- list(
            config = resolved$config,
            assembly_key = resolved$key,
            assembly = resolved$assembly,
            donor_species = resolved$donor_species,
            recipient_species = resolved$recipient_species,
            donor_components = resolved$donor_components,
            recipient_components = resolved$recipient_components,
            unpaired = resolved$unpaired
        )
    } else {
        if (is.null(key)) {
            keys <- names(config$assemblies)
            if (length(keys) == 1L) {
                key <- keys[[1L]]
            } else if (length(keys) == 0L) {
                .context_error(
                    "`config` has no resolved assembly; pass `donor` and `recipient`"
                )
            } else {
                .context_error(
                    "`config` has several resolved assemblies; pass `assembly` or a species pair"
                )
            }
        }
        assembly_config <- .context_assembly(config, key)
        context <- list(
            config = config,
            assembly_key = key,
            assembly = assembly_config,
            donor_species = NULL,
            recipient_species = NULL,
            donor_components = assembly_config$donor_components,
            recipient_components = assembly_config$recipient_components,
            unpaired = data.frame(
                role = character(0), contig = character(0),
                reason = character(0), stringsAsFactors = FALSE
            )
        )
    }

    structure(
        context,
        class = c("MaMaMIAAnalysisContext", "list")
    )
}

#' @export
print.MaMaMIAAnalysisContext <- function(x, ...) {
    cat(sprintf(
        "<MaMaMIAAnalysisContext> assembly '%s' (%s)\n",
        x$assembly_key, x$assembly$label
    ))
    cat(sprintf(
        "  recipient: %s [%s]\n",
        x$recipient_species %||% "<unknown>",
        paste(x$recipient_components, collapse = ",")
    ))
    cat(sprintf(
        "  donor: %s [%s]\n",
        x$donor_species %||% "<unknown>",
        paste(x$donor_components, collapse = ",")
    ))
    invisible(x)
}

## Internal normalizers: the pipeline accepts either a context or a plain
## configuration, and both preflight and the pipeline need the config alone.
.analysis_context <- function(config) {
    if (isMaMaMIAAnalysisContext(config)) config else NULL
}

.analysis_config <- function(config) {
    analysisContextConfig(config)
}
