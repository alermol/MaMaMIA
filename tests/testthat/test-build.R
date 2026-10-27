## T0.2: build system and dependency contract.

package_root <- function() {
    candidates <- c(
        testthat::test_path("..", ".."),
        getwd()
    )
    for (candidate in candidates) {
        if (file.exists(file.path(candidate, "DESCRIPTION"))) {
            return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
        }
    }
    testthat::skip("package root with DESCRIPTION not found")
}

read_dcf_field <- function(root, field) {
    description <- read.dcf(file.path(root, "DESCRIPTION"))
    if (!field %in% colnames(description) || is.na(description[1L, field])) {
        return(character(0))
    }
    items <- strsplit(description[1L, field], ",", fixed = TRUE)[[1L]]
    items <- trimws(items)
    items <- sub("\\s*\\(.*$", "", items)
    sort(unique(items[nzchar(items)]))
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

test_that("uvr.toml mirrors the DESCRIPTION dependency sets", {
    root <- package_root()
    dependencies <- read_uvr_sections(root)
    runtime <- sort(unique(c(
        read_dcf_field(root, "Imports"),
        read_dcf_field(root, "LinkingTo")
    )))
    development <- read_dcf_field(root, "Suggests")

    expect_setequal(dependencies$dependencies, runtime)
    expect_setequal(dependencies$`dev-dependencies`, development)
})

test_that("Bioconductor packages are declared as such in uvr.toml", {
    root <- package_root()
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
})

test_that("the package is compiled as C++23", {
    root <- package_root()
    expect_match(
        read.dcf(file.path(root, "DESCRIPTION"))[1L, "SystemRequirements"],
        "C\\+\\+23"
    )
    makevars <- readLines(file.path(root, "src", "Makevars"), warn = FALSE)
    expect_true(any(grepl("^\\s*CXX_STD\\s*=\\s*CXX23\\s*$", makevars)))
    expect_true(any(grepl("\\bOBJECTS\\b", makevars)))
    expect_true(any(grepl("gcrank/src/", makevars)))
})

test_that("the recorded pipeline version matches DESCRIPTION", {
    root <- package_root()
    expect_identical(
        MaMaMIA:::.MaMaMIAPipelineVersion,
        unname(read.dcf(file.path(root, "DESCRIPTION"))[1L, "Version"])
    )
})

test_that("the shared library loads and reports htslib", {
    version <- MaMaMIA:::htslib_version_cpp()
    expect_type(version, "character")
    expect_length(version, 1L)
    expect_match(version, "^[0-9]+\\.[0-9]+")
})

test_that("gcrank is vendored with attribution and a pinned commit", {
    root <- package_root()
    vendored <- c(
        file.path(root, "src", "gcrank", "include", "gcrank", "gc_rank.h"),
        file.path(root, "src", "gcrank", "src", "gc_rank.cpp"),
        file.path(root, "src", "gcrank", "LICENSE"),
        file.path(root, "src", "gcrank", "CHANGELOG.md"),
        file.path(root, "inst", "COPYRIGHTS"),
        file.path(root, "tools", "vendor-gcrank.txt")
    )
    expect_true(all(file.exists(vendored)))

    note <- readLines(file.path(root, "tools", "vendor-gcrank.txt"), warn = FALSE)
    expect_true(any(grepl("^commit = [0-9a-f]{40}$", note)))
    expect_true(any(grepl("^upstream = https://github.com/alermol/gcrank.git$", note)))
    expect_false(dir.exists(file.path(root, "src", "gcrank", "tests")))
    expect_false(dir.exists(file.path(root, "src", "gcrank", "cli")))
})
