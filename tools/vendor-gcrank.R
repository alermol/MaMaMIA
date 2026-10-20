#!/usr/bin/env Rscript

## Vendor the pinned gcrank checkout into src/gcrank.
##
## Usage:
##   Rscript tools/vendor-gcrank.R --source /path/to/gcrank [--commit <sha>]
##                                 [--dest src/gcrank] [--force]
##
## Only the library headers and sources (plus LICENSE and CHANGELOG.md) are
## copied. The command line client, tests, examples and build artefacts of the
## upstream checkout are deliberately excluded.

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
    "Usage: Rscript tools/vendor-gcrank.R --source <dir>",
    "[--commit <sha>] [--dest src/gcrank] [--force]"
)

value_after <- function(args, flag) {
    i <- match(flag, args)
    if (is.na(i)) {
        return(NULL)
    }
    if (i == length(args)) {
        stop(flag, " requires a value", call. = FALSE)
    }
    args[[i + 1L]]
}

has_flag <- function(args, flag) flag %in% args

source_dir <- value_after(args, "--source")
commit <- value_after(args, "--commit")
dest_dir <- value_after(args, "--dest")
force <- has_flag(args, "--force")

if (is.null(source_dir) || !nzchar(source_dir)) {
    stop(usage, call. = FALSE)
}
if (is.null(dest_dir) || !nzchar(dest_dir)) {
    dest_dir <- "src/gcrank"
}

if (!dir.exists(source_dir)) {
    stop("gcrank source directory does not exist: ", source_dir, call. = FALSE)
}
source_abs <- normalizePath(source_dir, mustWork = TRUE)
dest_abs <- normalizePath(dest_dir, mustWork = FALSE)

if (identical(source_abs, dest_abs)) {
    stop(
        "--source and --dest must differ; vendor from a separate checkout",
        call. = FALSE
    )
}

if (dir.exists(dest_dir) && !force) {
    stop(
        "destination already exists: ", dest_dir,
        "; pass --force to replace it",
        call. = FALSE
    )
}

## Prefer the source checkout's own metadata when --commit is omitted.
if (is.null(commit) || !nzchar(commit)) {
    commit <- tryCatch(
        {
            out <- suppressWarnings(system2(
                "git",
                c("-C", shQuote(source_abs), "rev-parse", "HEAD"),
                stdout = TRUE,
                stderr = FALSE
            ))
            if (length(out) == 1L && grepl("^[0-9a-f]{40}$", out)) {
                out
            } else {
                "unknown"
            }
        },
        error = function(e) "unknown"
    )
}

copy_files <- function(from, to, pattern) {
    if (!dir.exists(from)) {
        stop("missing upstream directory: ", from, call. = FALSE)
    }
    files <- list.files(
        from,
        pattern = pattern,
        recursive = TRUE,
        full.names = TRUE,
        all.files = FALSE
    )
    files <- files[file.info(files)$isdir %in% FALSE]
    if (length(files) == 0L) {
        stop("no files matched in ", from, call. = FALSE)
    }
    prefix <- paste0(from, .Platform$file.sep)
    for (file in files) {
        relative <- sub(paste0("^", prefix), "", file, fixed = FALSE)
        target <- file.path(to, relative)
        dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
        if (!file.copy(file, target, overwrite = TRUE)) {
            stop("could not copy ", file, " to ", target, call. = FALSE)
        }
    }
    invisible(length(files))
}

if (dir.exists(dest_dir)) {
    unlink(dest_dir, recursive = TRUE, force = TRUE)
}
dir.create(file.path(dest_dir, "include"), recursive = TRUE)
dir.create(file.path(dest_dir, "src"), recursive = TRUE)

n_headers <- copy_files(
    file.path(source_abs, "include"),
    file.path(dest_dir, "include"),
    "\\.(h|hpp)$"
)
n_sources <- copy_files(
    file.path(source_abs, "src"),
    file.path(dest_dir, "src"),
    "\\.(cpp|cc|cxx|h|hpp)$"
)

for (file in c("LICENSE", "CHANGELOG.md")) {
    from <- file.path(source_abs, file)
    if (!file.exists(from)) {
        stop("missing upstream file: ", from, call. = FALSE)
    }
    if (!file.copy(from, file.path(dest_dir, file), overwrite = TRUE)) {
        stop("could not copy ", from, call. = FALSE)
    }
}

epoch <- Sys.getenv("SOURCE_DATE_EPOCH", unset = "")
stamp <- if (nzchar(epoch)) {
    format(as.POSIXct(as.numeric(epoch), origin = "1970-01-01", tz = "UTC"),
        "%Y-%m-%dT%H:%M:%SZ",
        tz = "UTC"
    )
} else {
    format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}

notes <- c(
    "upstream = https://github.com/alermol/gcrank.git",
    paste0("commit = ", commit),
    paste0("vendored_at = ", stamp),
    "vendored_paths = include/, src/, LICENSE, CHANGELOG.md",
    "excluded_paths = cli/, tests/, examples/, cmake/, data-raw/, CMakeLists.txt",
    paste0("headers = ", n_headers),
    paste0("sources = ", n_sources)
)
writeLines(notes, "tools/vendor-gcrank.txt", useBytes = TRUE)

message(
    "Vendored gcrank ", commit, " into ", dest_dir, " (",
    n_headers, " headers, ", n_sources, " sources)"
)
