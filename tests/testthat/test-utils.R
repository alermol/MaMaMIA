## Shared utility helpers.

test_that("hash_object hashes values in memory and is deterministic", {
    first <- hash_object(list(a = 1, b = "x"))
    second <- hash_object(list(a = 1, b = "x"))
    expect_identical(first, second)
    expect_match(first, "^[0-9a-f]{32}$")
    expect_false(identical(first, hash_object(list(a = 2, b = "x"))))
})
