windows <- data.frame(chr_id = "c1", iid = 1:4, cov = 10:13, gc = c(0.40, 0.41, 0.42, 0.43), subgenome = "A")

expect_key_change <- function(key, windows = get("windows", parent.frame()), cores = 4L, call = MaMaMIA:::coverage_model_call(cores)) {
    expect_false(identical(key, MaMaMIA:::fit_cache_key(windows, cores, call)))
}

test_that("the fit cache key follows data, cores, versions and model", {
    call <- MaMaMIA:::coverage_model_call(4L)
    key <- MaMaMIA:::fit_cache_key(windows, 4L, call)
    expect_type(key, "character"); expect_length(key, 1L)
    expect_identical(key, MaMaMIA:::fit_cache_key(windows, 4L, call))

    expect_key_change(key, cores = 8L)
    for (mutated in list(transform(windows, cov = cov + 1), transform(windows, gc = gc + 1e-3), transform(windows, subgenome = "B"), windows[-1L, ])) {
        expect_key_change(key, mutated)
    }

    wider <- as.list(call); wider[[2L]] <- quote(cov ~ s(gc, k = 15) + subgenome)
    no_zi <- as.list(call); no_zi$ziformula <- quote(~ 1)
    relabelled <- as.list(call); relabelled$data <- quote(anything)
    expect_key_change(key, call = as.call(wider))
    expect_key_change(key, call = as.call(no_zi))
    expect_identical(key, MaMaMIA:::fit_cache_key(windows, 4L, as.call(relabelled)))
})

test_that("cached fits are reused only when the key and format match", {
    path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
    expect_null(MaMaMIA:::read_cached_fit(NULL, "key")); expect_null(MaMaMIA:::read_cached_fit(path, "key"))
    fake <- structure(list(fixef = 1), class = "glmmTMB")
    MaMaMIA:::write_cached_fit(path, "key", fake)
    expect_true(file.exists(path)); expect_identical(MaMaMIA:::read_cached_fit(path, "key"), fake)
    expect_identical(readRDS(path)$format, MaMaMIA:::CACHE_LAYOUT_VERSION)
    expect_null(MaMaMIA:::read_cached_fit(path, "other"))
    writeLines("not an rds file", path); expect_null(MaMaMIA:::read_cached_fit(path, "key"))
    saveRDS(list(format = MaMaMIA:::CACHE_LAYOUT_VERSION + 1L, key = "key", fit = fake), path)
    expect_null(MaMaMIA:::read_cached_fit(path, "key"))
    MaMaMIA:::write_cached_fit(NULL, "key", fake)
    expect_false(file.exists(paste0(path, ".tmp")))
})
