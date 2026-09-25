test_that("subset() restricts an RCA object by subgenome", {
    x <- small_rca()

    # Only the donor side is filtered; recipients are left unchanged.
    s <- subset(x, don_subgenomes = "At")
    expect_setequal(unique(s$data$chr_id), c("d1", "r1", "r2"))
    expect_setequal(s$meta$don_chr_ids, "d1")
    expect_setequal(s$meta$rec_chr_ids, c("r1", "r2"))

    both <- subset(x, don_subgenomes = "G", rec_subgenomes = "B")
    expect_setequal(unique(both$data$chr_id), c("d2", "r2"))
})

test_that("subset() restricts an RCA object by chromosome ID", {
    x <- small_rca()
    s <- subset(x, don_chromlist = "d2", rec_chromlist = "r2")

    expect_setequal(unique(s$data$chr_id), c("d2", "r2"))
    expect_setequal(s$meta$don_chr_ids, "d2")
    expect_setequal(s$meta$rec_chr_ids, "r2")
})

test_that("subset() validates its arguments", {
    x <- small_rca()

    expect_error(subset(x), "at least one")
    expect_error(subset(x, don_subgenomes = "nope"), "Unknown donor subgenome")
    expect_error(subset(x, rec_subgenomes = "nope"), "Unknown recipient subgenome")
    expect_error(subset(x, don_chromlist = "nope"), "Unknown donor chromosome")
    expect_error(subset(x, rec_chromlist = "nope"), "Unknown recipient chromosome")
})

test_that("reverseWindows() flips value vectors without touching the window grid", {
    x <- small_rca()
    # GC is a ramp per chromosome, so a reversal is observable.
    x$data$cor.gc <- x$data$gc

    reversed <- reverseWindows(x, chr_ids = "d1")

    expect_equal(reversed$data$chr_id, x$data$chr_id)
    expect_equal(reversed$data$start, x$data$start)
    expect_equal(reversed$data$stop, x$data$stop)
    expect_equal(reversed$data$iid, x$data$iid)
    expect_equal(reversed$meta, x$meta)

    d1 <- x$data$chr_id == "d1"
    expect_equal(reversed$data$gc[d1], rev(x$data$gc[d1]))
    expect_equal(reversed$data$cor.gc[d1], rev(x$data$cor.gc[d1]))

    # Chromosomes that were not selected are untouched.
    d2 <- x$data$chr_id == "d2"
    expect_equal(reversed$data$gc[d2], x$data$gc[d2])
    expect_equal(reversed$data$cor.gc[d2], x$data$cor.gc[d2])
})

test_that("reverseWindows() reverses every selected chromosome", {
    x <- small_rca()
    x$data$cor.gc <- x$data$gc

    reversed <- reverseWindows(x, chr_ids = c("d1", "r2"))
    for (id in c("d1", "r2")) {
        idx <- x$data$chr_id == id
        expect_equal(reversed$data$gc[idx], rev(x$data$gc[idx]))
    }
    untouched <- x$data$chr_id %in% c("d2", "r1")
    expect_equal(reversed$data$gc[untouched], x$data$gc[untouched])
})

test_that("reverseWindows() validates its arguments", {
    x <- small_rca()

    expect_error(reverseWindows(x, chr_ids = "nope"), "not present")
    expect_warning(
        reverseWindows(x, chr_ids = character(0)),
        "No chromosomes requested"
    )
    expect_error(reverseWindows(list(), chr_ids = "d1"))
})
