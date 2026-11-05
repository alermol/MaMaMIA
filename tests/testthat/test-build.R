## Build and dependency contract for the web service.

app_root <- function() {
    candidates <- c(testthat::test_path("..", ".."), getwd())
    for (candidate in candidates) {
        if (file.exists(file.path(candidate, "VERSION")) &&
            file.exists(file.path(candidate, "R", "load.R"))) {
            return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
        }
    }
    testthat::skip("application root not found")
}

read_uvr_sections <- function(root) {
    lines <- readLines(file.path(root, "uvr.toml"), warn = FALSE)
    section <- NA_character_
    out <- list(`dependencies` = character(0), `dev-dependencies` = character(0))
    for (line in lines) {
        trimmed <- trimws(sub("#.*$", "", line))
        if (!nzchar(trimmed)) {
            next
        }
        header <- regmatches(trimmed, regexec("^\\[(.+)\\]$", trimmed))[[1L]]
        if (length(header) == 2L) {
            section <- header[[2L]]
            next
        }
        if (!is.na(section) && section %in% names(out)) {
            match <- regmatches(
                trimmed,
                regexec("^([A-Za-z0-9._]+)\\s*=", trimmed)
            )[[1L]]
            if (length(match) == 2L) {
                out[[section]] <- c(out[[section]], match[[2L]])
            }
        }
    }
    lapply(out, function(x) sort(unique(x)))
}

test_that("the recorded pipeline version matches VERSION", {
    root <- app_root()
    expect_identical(
        .MaMaMIAPipelineVersion,
        trimws(readLines(file.path(root, "VERSION"), n = 1L))
    )
})

test_that("uvr.toml declares the Bioconductor and dashboard packages", {
    root <- app_root()
    dependencies <- read_uvr_sections(root)
    lines <- readLines(file.path(root, "uvr.toml"), warn = FALSE)
    bioc <- c("DNAcopy", "GenomicRanges", "IRanges", "Rhtslib")
    for (package in bioc) {
        line <- lines[grepl(
            paste0("^\\s*", package, "\\s*="),
            lines
        )]
        expect_length(line, 1L)
        expect_match(line, "bioc\\s*=\\s*true")
    }
    expect_true(all(c(
        "shiny", "bslib", "DT", "htmltools", "htmlwidgets",
        "shinyWidgets", "gridlayout"
    ) %in% dependencies$dependencies))
})

test_that("the native bridge is compiled as C++23", {
    root <- app_root()
    makevars <- readLines(file.path(root, "src", "Makevars.in"), warn = FALSE)
    expect_true(any(grepl("^\\s*CXX_STD\\s*=\\s*CXX23\\s*$", makevars)))
    expect_true(any(grepl("\\bOBJECTS\\b", makevars)))
    expect_true(any(grepl("gcrank/src/", makevars)))
})

test_that("the shared library loads and reports htslib", {
    version <- htslib_version_cpp()
    expect_type(version, "character")
    expect_length(version, 1L)
    expect_match(version, "^[0-9]+\\.[0-9]+")
})

test_that("gcrank is vendored with attribution and a pinned commit", {
    root <- app_root()
    vendored <- c(
        file.path(root, "src", "gcrank", "include", "gcrank", "gc_rank.h"),
        file.path(root, "src", "gcrank", "src", "gc_rank.cpp"),
        file.path(root, "src", "gcrank", "LICENSE"),
        file.path(root, "src", "gcrank", "CHANGELOG.md"),
        file.path(root, "COPYRIGHTS"),
        file.path(root, "tools", "vendor-gcrank.txt")
    )
    expect_true(all(file.exists(vendored)))

    note <- readLines(file.path(root, "tools", "vendor-gcrank.txt"), warn = FALSE)
    expect_true(any(grepl("^commit = [0-9a-f]{40}$", note)))
    expect_true(any(grepl("^upstream = https://github.com/alermol/gcrank.git$", note)))
    expect_false(dir.exists(file.path(root, "src", "gcrank", "tests")))
    expect_false(dir.exists(file.path(root, "src", "gcrank", "cli")))
})
