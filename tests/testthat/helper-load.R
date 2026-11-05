## Make the source application available when a test file is run directly
## through testthat::test_dir() or testthat::test_file().

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
