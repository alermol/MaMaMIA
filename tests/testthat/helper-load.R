## Make the source application available when a test file is run directly
## through testthat::test_dir() or testthat::test_file().

## The dictionaries and their tests are UTF-8. Under LC_ALL=C the C locale
## cannot compare or print Cyrillic reliably, so widen LC_CTYPE while keeping
## numeric/time collation unchanged. Fall back silently when no UTF-8 locale
## is installed.
if (!grepl("UTF-8|utf8", Sys.getlocale("LC_CTYPE"), ignore.case = TRUE)) {
    for (candidate in c("C.UTF-8", "C.utf8", "en_US.UTF-8", "en_US.utf8")) {
        if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", candidate)))) {
            break
        }
    }
}

if (!exists("mamamia_load", mode = "function")) {
    root <- Sys.getenv("MAMAMIA_ROOT", unset = getwd())
    root <- normalizePath(root, winslash = "/", mustWork = FALSE)
    while (!file.exists(file.path(root, "R", "load.R")) &&
        !identical(dirname(root), root)) {
        root <- dirname(root)
    }
    source(file.path(root, "R", "load.R"), local = FALSE)
    mamamia_load(root)
}
