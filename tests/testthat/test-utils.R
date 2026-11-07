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

test_that("natural sorting orders digit runs as numbers", {
    expect_identical(
        mamamia_natural_sort(c("sample10", "sample2", "sample1")),
        c("sample1", "sample2", "sample10")
    )
    expect_identical(
        mamamia_natural_sort(c("Sq100", "Sq9", "Sq80b", "Sq80a", "Sq80")),
        c("Sq9", "Sq80", "Sq80a", "Sq80b", "Sq100")
    )
    expect_identical(
        mamamia_natural_sort(c("b2", "A10", "a2")),
        c("a2", "A10", "b2")
    )
})

test_that("natural sorting keeps missing values last and ties stable", {
    expect_identical(
        mamamia_natural_sort(c("s2", NA, "s1")),
        c("s1", "s2", NA)
    )
    expect_identical(
        mamamia_natural_sort(c("s01", "s1")),
        c("s01", "s1")
    )
})

test_that("bam_sample_choices returns samples in natural order", {
    result <- list(samples = data.frame(
        sample = c("sample10", "sample2", "sample1"),
        stringsAsFactors = FALSE
    ))
    expect_identical(
        bam_sample_choices(result),
        c("sample1", "sample2", "sample10")
    )
    expect_identical(bam_sample_choices(NULL), character())
})

test_that("the bundle worker budget respects the container memory limit", {
    gib <- 1024^3
    expect_identical(mamamia_memory_worker_cap(16 * gib), 6L)
    expect_identical(mamamia_memory_worker_cap(64 * gib), 30L)
    expect_identical(mamamia_memory_worker_cap(4 * gib), 1L)
    expect_identical(mamamia_memory_worker_cap(2 * gib), 1L)
    expect_true(is.infinite(mamamia_memory_worker_cap(Inf)))
})
