test_that("defaultTargetPairs() matches chromosome numbers within subgenomes", {
    meta <- make_meta(
        c("d1", "d2", "d3", "r1", "r2", "r3"),
        c("1At", "2At", "3G", "1A", "2A", "1B"),
        c("At", "At", "G", "A", "A", "B")
    )
    pairs <- defaultTargetPairs(
        meta, c("d1", "d2", "d3"), c("r1", "r2", "r3"),
        don_subgenomes = "At", rec_subgenomes = "A"
    )
    expect_equal(pairs, c(d1 = "r1", d2 = "r2"))

    ## A second recipient subgenome makes every number ambiguous
    expect_error(
        defaultTargetPairs(meta, "d1", c("r1", "r2", "r3")),
        "No donor->recipient pairs"
    )
    expect_error(
        defaultTargetPairs(meta, "d1", "r1", don_subgenomes = "Z"),
        "No donor or recipient chromosomes match"
    )
})
