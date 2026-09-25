test_that("correctReadCounts() fits, annotates and summarizes GC-corrected coverage", {
    skip_on_cran()
    res <- triticum_fit()
    expect_s3_class(res, "RCA"); expect_true(res$corrected); expect_s3_class(res$fit, "glmmTMB")
    expect_true(all(c("valid", "ideal", "cor.gc") %in% names(res$data)))
    expect_type(res$data$valid, "logical"); expect_type(res$data$ideal, "logical")
    expect_true(is.numeric(res$data$cor.gc)); expect_true(all(is.finite(res$data$cor.gc))); expect_true(any(res$data$ideal))
    expect_named(res$outliers, c("gc_lower_bound", "gc_upper_bound", "cov_upper_bound"))
    expect_true(res$outliers$gc_lower_bound < res$outliers$gc_upper_bound)
    expect_true(summary(res)$corrected); expect_output(print(summary(res)), "GC-bias corrected: Yes")
})

test_that("correctReadCounts() validates input", {
    expect_error(correctReadCounts(list()), "not RCA")
})

test_that("a cached fit is reused and invalidated by a data change", {
    skip_on_cran()
    x <- do.call(RCA, triticum_one_pair())
    path <- tempfile(fileext = ".rds"); on.exit(unlink(path), add = TRUE)
    cold <- correctReadCounts(x, cores = 1L, verbose = FALSE, cache = path)
    warm <- correctReadCounts(x, cores = 1L, verbose = FALSE, cache = path)
    expect_true(file.exists(path))
    expect_identical(warm$data[c("cor.gc", "ideal")], cold$data[c("cor.gc", "ideal")])
    expect_identical(glmmTMB::fixef(warm$fit), glmmTMB::fixef(cold$fit))
    expect_true(correctReadCounts(reverseWindows(x, unique(x$data$chr_id)[1L]), cores = 1L, verbose = FALSE, cache = path)$corrected)
})

test_that("the fast ZINB prediction reproduces predict() and falls back safely", {
    skip_on_cran()
    res <- triticum_fit()
    gc_ref <- stats::median(res$data$gc[res$data$ideal], na.rm = TRUE)
    fast <- MaMaMIA:::predict_zinb_response(res$fit, res$data$gc, res$data$subgenome, gc_ref)
    nd <- data.frame(gc = res$data$gc, subgenome = res$data$subgenome)
    ref <- list(
        actual = as.numeric(stats::predict(res$fit, newdata = nd, type = "response")),
        ref = as.numeric(stats::predict(res$fit, newdata = transform(nd, gc = gc_ref), type = "response"))
    )
    expect_equal(fast$actual, ref$actual, tolerance = 1e-8)
    expect_equal(fast$ref, ref$ref, tolerance = 1e-8)
    expect_length(fast$actual, nrow(res$data)); expect_length(fast$ref, nrow(res$data))
    expect_true(all(unlist(tapply(fast$ref, res$data$subgenome, function(z) diff(range(z)))) < 1e-8))

    broken <- res$fit; broken$fit$parfull <- NULL
    fallback <- MaMaMIA:::predict_zinb_response(broken, res$data$gc, res$data$subgenome, gc_ref)
    expect_equal(fast$actual, fallback$actual, tolerance = 1e-8)
    expect_equal(fast$ref, fallback$ref, tolerance = 1e-8)
})

test_that("correctReadCounts() output feeds segmentation", {
    skip_on_cran()
    res <- triticum_fit()
    isa <- segments(res, target_pairs = stats::setNames(res$meta$rec_chr_ids, res$meta$don_chr_ids), verbose = FALSE)
    expect_s3_class(isa, "ISA")
    expect_true(all(isa$segments$chr_id.don %in% res$meta$don_chr_ids))
    expect_true(all(isa$segments$chr_id.rec %in% res$meta$rec_chr_ids))
    expect_true(all(is.finite(isa$segments$seg.median)))
})
