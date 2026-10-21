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

test_that("fit subsampling respects the cap and the per-stratum minimum", {
    skip_on_cran()
    x <- do.call(RCA, triticum_one_pair())
    sampled <- correctReadCounts(
        x,
        cores = 1L,
        verbose = FALSE,
        max_fit_windows = 500L
    )
    sampling <- sampled$fit_sampling

    expect_true(sampling$sampled)
    expect_identical(sampling$max_fit_windows, 500)
    expect_identical(sampling$n_fitted, 500L)
    expect_gt(sampling$n_ideal, sampling$n_fitted)
    expect_identical(sum(sampling$strata$n_fitted), 500L)
    expect_true(all(
        sampling$strata$n_fitted >=
            pmin(sampling$strata$n_ideal, 10L)
    ))
    expect_true(all(sampling$strata$n_fitted <= sampling$strata$n_ideal))
    expect_identical(
        names(sampling$strata),
        c("subgenome", "gc_decile", "n_ideal", "n_fitted")
    )
    expect_identical(anyDuplicated(paste(
        sampling$strata$subgenome, sampling$strata$gc_decile
    )), 0L)
    expect_true(all(is.finite(sampled$data$cor.gc)))
    expect_output(print(summary(sampled)), "Fit windows: 500 of")
})

test_that("fit subsampling is reproducible and leaves the RNG alone", {
    skip_on_cran()
    x <- do.call(RCA, triticum_one_pair())
    set.seed(42L)
    before <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    first <- correctReadCounts(
        x, cores = 1L, verbose = FALSE, max_fit_windows = 500L
    )
    expect_identical(
        get(".Random.seed", envir = globalenv(), inherits = FALSE),
        before
    )
    second <- correctReadCounts(
        x, cores = 1L, verbose = FALSE, max_fit_windows = 500L
    )
    expect_identical(second$fit_sampling, first$fit_sampling)
    expect_equal(second$data$cor.gc, first$data$cor.gc, tolerance = 1e-10)
    expect_equal(
        glmmTMB::fixef(second$fit), glmmTMB::fixef(first$fit),
        tolerance = 1e-8
    )
})

test_that("the fit cache follows max_fit_windows", {
    skip_on_cran()
    x <- do.call(RCA, triticum_one_pair())
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path), add = TRUE)

    first <- correctReadCounts(
        x, cores = 1L, verbose = FALSE, cache = path,
        max_fit_windows = 500L
    )
    key <- readRDS(path)$key
    expect_identical(first$fit_sampling$n_fitted, 500L)

    warm <- correctReadCounts(
        x, cores = 1L, verbose = FALSE, cache = path,
        max_fit_windows = 500L
    )
    expect_identical(readRDS(path)$key, key)
    expect_equal(warm$data$cor.gc, first$data$cor.gc, tolerance = 1e-10)

    larger <- correctReadCounts(
        x, cores = 1L, verbose = FALSE, cache = path,
        max_fit_windows = 700L
    )
    expect_false(identical(readRDS(path)$key, key))
    expect_identical(larger$fit_sampling$n_fitted, 700L)
})

test_that("a subsampled fit stays close to the full fit", {
    skip_on_cran()
    full <- triticum_fit()
    sampled <- correctReadCounts(
        do.call(RCA, triticum_one_pair()),
        cores = 1L,
        verbose = FALSE,
        max_fit_windows = 600L
    )
    common <- full$data$valid & sampled$data$valid
    expect_gt(sum(common), 1000L)
    expect_true(all(is.finite(sampled$data$cor.gc[common])))
    ## The deterministic subsample seed is part of the fit-cache key, so a
    ## package version bump selects a different 600-window subset. Keep this
    ## sanity bound generous (the median relative error below is the tighter
    ## check); observed values are around 0.88-0.99 depending on the subset.
    expect_gt(
        stats::cor(full$data$cor.gc[common], sampled$data$cor.gc[common]),
        0.85
    )
    relative <- abs(
        full$data$cor.gc[common] - sampled$data$cor.gc[common]
    ) / pmax(abs(full$data$cor.gc[common]), 1e-8)
    expect_lt(stats::median(relative), 0.05)
})

test_that("a subsampled fit feeds segmentation", {
    skip_on_cran()
    sampled <- correctReadCounts(
        do.call(RCA, triticum_one_pair()),
        cores = 1L,
        verbose = FALSE,
        max_fit_windows = 600L
    )
    isa <- segments(
        sampled,
        target_pairs = stats::setNames(
            sampled$meta$rec_chr_ids,
            sampled$meta$don_chr_ids
        ),
        verbose = FALSE
    )
    expect_s3_class(isa, "ISA")
    expect_gt(nrow(isa$segments), 0L)
    expect_true(all(is.finite(isa$segments$seg.median)))
})

test_that("fit subsampling validates its cap and reports too-small caps", {
    skip_on_cran()
    x <- do.call(RCA, triticum_one_pair())
    expect_error(
        correctReadCounts(x, cores = 1L, verbose = FALSE, max_fit_windows = 0),
        "positive integer or Inf"
    )
    expect_error(
        correctReadCounts(x, cores = 1L, verbose = FALSE, max_fit_windows = 1.5),
        "positive integer or Inf"
    )
    expect_error(
        correctReadCounts(x, cores = 1L, verbose = FALSE, max_fit_windows = NA),
        "positive integer or Inf"
    )
    expect_error(
        correctReadCounts(x, cores = 1L, verbose = FALSE, max_fit_windows = 1L),
        "at least"
    )
})
