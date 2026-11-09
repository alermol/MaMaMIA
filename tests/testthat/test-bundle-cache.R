## One-entry result-bundle cache and its request key.

test_that("bundle_cache_serve() copies a repeated request without rebuilding", {
    cache <- bundle_cache_state()
    on.exit(bundle_cache_clear(cache), add = TRUE)
    first <- tempfile(fileext = ".zip")
    second <- tempfile(fileext = ".zip")
    builds <- 0L
    build <- function(path) {
        builds <<- builds + 1L
        writeLines("zip-content", path)
    }
    key <- list(analysis_id = "run-1", sample = "S")

    bundle_cache_serve(cache, key, first, build)
    expect_identical(builds, 1L)
    expect_true(file.exists(first))

    bundle_cache_serve(cache, key, second, build)
    expect_identical(builds, 1L)
    expect_identical(readLines(first), readLines(second))
})

test_that("a different key rebuilds and replaces the cached archive", {
    cache <- bundle_cache_state()
    on.exit(bundle_cache_clear(cache), add = TRUE)
    builds <- 0L
    build <- function(path) {
        builds <<- builds + 1L
        writeLines(paste0("zip-", builds), path)
    }
    one <- tempfile(fileext = ".zip")
    two <- tempfile(fileext = ".zip")
    bundle_cache_serve(cache, list(key = 1), one, build)
    cached_one <- file.path(cache$path)
    bundle_cache_serve(cache, list(key = 2), two, build)
    expect_identical(builds, 2L)
    expect_false(file.exists(cached_one))
    expect_true(file.exists(cache$path))
})

test_that("a failed build is not cached and is retried", {
    cache <- bundle_cache_state()
    on.exit(bundle_cache_clear(cache), add = TRUE)
    attempts <- 0L
    build <- function(path) {
        attempts <<- attempts + 1L
        if (attempts == 1L) {
            stop("build failed", call. = FALSE)
        }
        writeLines("zip-content", path)
    }
    out <- tempfile(fileext = ".zip")
    expect_error(
        bundle_cache_serve(cache, list(key = 1), out, build),
        "build failed"
    )
    expect_null(bundle_cache_get(cache, list(key = 1)))
    bundle_cache_serve(cache, list(key = 1), out, build)
    expect_identical(attempts, 2L)
    expect_true(file.exists(out))
})

test_that("empty archives are rejected and cleared", {
    cache <- bundle_cache_state()
    empty <- tempfile(fileext = ".zip")
    file.create(empty)
    expect_error(
        bundle_cache_store(cache, list(key = 1), empty),
        "empty"
    )
    good <- tempfile(fileext = ".zip")
    writeLines("zip-content", good)
    bundle_cache_store(cache, list(key = 1), good)
    stored <- cache$path
    bundle_cache_clear(cache)
    expect_null(cache$key)
    expect_null(cache$path)
    expect_false(file.exists(stored))
})

test_that("bundle_request_key() is structured and changes with every input", {
    result <- bam_example_result(1000)
    calls <- data.frame(.call_id = c(2L, 1L), stringsAsFactors = FALSE)
    request <- bundle_request(
        result, "example-sample", "calls",
        calls = calls, groups = c("2", "1"),
        hide_deletions = FALSE, hide_unplaced = TRUE
    )
    key <- bundle_request_key(request)
    expect_identical(key$call_ids, c(1L, 2L))
    expect_identical(key$groups, c("1", "2"))
    expect_true(key$hide_unplaced)
    expect_false(key$hide_deletions)
    expect_true(nzchar(key$analysis_id))
    expect_true(nzchar(key$bundle_version))

    changed <- request
    changed$hide_unplaced <- FALSE
    expect_false(identical(key, bundle_request_key(changed)))

    changed <- request
    changed$calls$.call_id <- c(1L, 3L)
    expect_false(identical(key, bundle_request_key(changed)))

    changed <- request
    changed$sample <- "other"
    expect_false(identical(key, bundle_request_key(changed)))
})

test_that("bundle_analysis_id() prefers the pipeline id", {
    result <- bam_example_result(1000)
    expect_identical(
        bundle_analysis_id(result),
        result$provenance$analysis_id
    )
    old <- result
    old$provenance$analysis_id <- NULL
    old$provenance$started_at <- "2026-01-01T00:00:00+0000"
    expect_match(bundle_analysis_id(old), "2026-01-01")
})
