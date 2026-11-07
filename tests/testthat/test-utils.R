## Shared utility helpers.

test_that("hash_object hashes values in memory and is deterministic", {
    first <- hash_object(list(a = 1, b = "x"))
    second <- hash_object(list(a = 1, b = "x"))
    expect_identical(first, second)
    expect_match(first, "^[0-9a-f]{32}$")
    expect_false(identical(first, hash_object(list(a = 2, b = "x"))))
})

test_that("the default BAM budget leaves one available core free", {
    expect_identical(mamamia_default_cores(8), 7L)
    expect_identical(mamamia_default_cores(2), 1L)
    expect_identical(mamamia_default_cores(1), 1L)
    expect_identical(mamamia_default_cores(NA_real_), 1L)
    expect_identical(mamamia_default_cores(1000), 256L)
})

test_that("direct BAM counting uses the default budget", {
    withr::local_options(
        MaMaMIA.metrics_workers = NULL,
        MaMaMIA.htslib_threads = NULL
    )
    expect_identical(.bam_htslib_threads(), mamamia_default_cores())
})
