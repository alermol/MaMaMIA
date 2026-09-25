#!/usr/bin/env Rscript

abort <- function(...) stop(..., call. = FALSE)
arg <- function(name, default = NA_character_) {
    hit <- grep(paste0("^--", name, "="), commandArgs(TRUE), value = TRUE)
    if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else default
}
opts <- list(
    cores = as.integer(arg("cores", "4")),
    reps = as.integer(arg("repeat", "1")),
    out = arg("out", "bench-out"),
    compare = arg("compare")
)
if (is.na(opts$cores) || opts$cores < 1L) abort("--cores must be a positive integer")
if (is.na(opts$reps) || opts$reps < 1L) abort("--repeat must be a positive integer")

read_vm <- function(field) {
    f <- "/proc/self/status"
    if (!file.exists(f)) return(NA_real_)
    hit <- grep(paste0("^", field, ":"), readLines(f, warn = FALSE), value = TRUE)
    if (length(hit)) as.numeric(sub(" kB", "", sub("^.*:\\s+", "", hit[[1L]]))) else NA_real_
}
heap_mb <- function(g) g["Ncells", "max used"] * 56 / 1e6 + g["Vcells", "max used"] * 8 / 1e6
hash_object <- function(x) {
    f <- tempfile(); on.exit(unlink(f), add = TRUE)
    saveRDS(x, f, version = 2); unname(tools::md5sum(f))
}
record <- function(env, section, step, expr, size = TRUE) {
    h0 <- read_vm("VmHWM"); gc(reset = TRUE)
    tm <- system.time(value <- force(expr)); g <- gc(); h1 <- read_vm("VmHWM")
    env$rows[[length(env$rows) + 1L]] <- data.frame(
        section, step,
        elapsed_s = round(tm[["elapsed"]], 3),
        peak_heap_MB = round(heap_mb(g), 2),
        peak_rss_MB = round(if (anyNA(c(h0, h1))) NA_real_ else (h1 - h0) / 1024, 2),
        vmhwm_MB = round(h1 / 1024, 2),
        obj_MB = round(if (size) as.numeric(utils::object.size(value)) / 1e6 else NA_real_, 3)
    )
    invisible(value)
}
check <- function(env, label, ok, detail = "") {
    env$checks[[length(env$checks) + 1L]] <- data.frame(
        check = label, result = if (isTRUE(ok)) "PASS" else "FAIL", detail = detail
    )
    invisible(ok)
}
print_section <- function(env, section) {
    x <- do.call(rbind, env$rows)
    x <- x[x$section == section, setdiff(names(x), "section"), drop = FALSE]
    if (nrow(x)) print(x, row.names = FALSE)
}

run_pipeline <- function(env, cores, section = "pipeline") {
    record(env, section, "data(triticum)", utils::data("triticum", package = "MaMaMIA", envir = environment()), FALSE)
    rca <- record(env, section, "RCA()", RCA(don_cov, don_gc, rec_cov, rec_gc, meta))
    corrected <- record(env, section, sprintf("correctReadCounts(cores=%d)", cores), correctReadCounts(rca, cores = cores, verbose = FALSE))
    subset_cd <- record(env, section, "subset(At / A)", subset(corrected, don_subgenomes = "At", rec_subgenomes = "A"))
    record(env, section, "reverseWindows(RCA)", reverseWindows(subset_cd, chr_ids = "OY997261.1"))
    pairs <- stats::setNames(subset_cd$meta$rec_chr_ids, subset_cd$meta$don_chr_ids)
    isa <- record(env, section, "segments()", segments(subset_cd, target_pairs = pairs, verbose = FALSE))
    record(env, section, "summary(ISA)", summary(isa, unit = "Mb", digits = 2))
    bedpe <- tempfile(fileext = ".bedpe")
    record(env, section, "writeToBedpe()", writeToBedpe(isa, bedpe), FALSE)
    record(env, section, "plot(ISA, mirror)", plot(isa, plot.type = "mirror"), FALSE)
    record(env, section, "plot(ISA, diff)", plot(isa, plot.type = "diff"), FALSE)
    list(rca = rca, corrected = corrected, subset_cd = subset_cd, isa = isa, bedpe = bedpe)
}

correction_internals <- function(env, cores) {
    rca <- record(env, "correction", "data + RCA()", {
        utils::data("triticum", package = "MaMaMIA", envir = environment())
        RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
    }, FALSE)
    dat <- rca$data
    dat <- record(env, "correction", "  filter + quantiles", {
        dat$valid <- dat$cov >= 0 & dat$gc > 0
        gc_q <- stats::quantile(dat$gc[dat$valid], c(0.01, 0.99))
        cov_q <- stats::quantile(dat$cov[dat$valid], 0.99)
        dat$ideal <- dat$valid & dat$gc >= gc_q[1] & dat$gc <= gc_q[2] & dat$cov <= cov_q
        dat
    })
    ideal <- dat[dat$ideal, ]
    fit <- record(env, "correction", "  glmmTMB fit", glmmTMB::glmmTMB(
        cov ~ s(gc, k = 10) + subgenome, ziformula = ~ s(gc, k = 10) + subgenome,
        family = glmmTMB::nbinom2(), data = ideal, REML = TRUE,
        control = glmmTMB::glmmTMBControl(parallel = list(n = cores))
    ))
    gc_ref <- stats::median(dat$gc[dat$ideal], na.rm = TRUE)
    pred <- record(env, "correction", "  predict (fast path)", MaMaMIA:::predict_zinb_response(fit, dat$gc, dat$subgenome, gc_ref))
    dat <- record(env, "correction", "  cor.gc + post-filter", {
        dat$cor.gc <- dat$cov * pred$ref / (pred$actual + 1e-8)
        dat$ideal <- dat$ideal & dat$cor.gc < stats::quantile(dat$cor.gc, 0.99, na.rm = TRUE)
        dat
    })
    list(fit = fit, dat = dat, gc_ref = gc_ref, pred = pred)
}

predict_reference <- function(env, fit, dat, gc_ref) {
    nd <- data.frame(gc = dat$gc, subgenome = dat$subgenome)
    record(env, "predict", "  glmmTMB::predict() reference", {
        list(
            actual = as.numeric(stats::predict(fit, newdata = nd, type = "response")),
            ref = as.numeric(stats::predict(fit, newdata = transform(nd, gc = gc_ref), type = "response"))
        )
    })
}

fingerprint <- function(res, cores) {
    nums <- res$isa$segments[vapply(res$isa$segments, is.numeric, logical(1L))]
    c(
        cores = cores,
        n_ideal = sum(res$corrected$data$ideal),
        n_segments = nrow(res$isa$segments),
        n_putative = sum(res$isa$segments$putative_introgression),
        r_cor_gc = hash_object(signif(res$corrected$data$cor.gc, 10L)),
        r_fixef = hash_object(signif(unlist(glmmTMB::fixef(res$corrected$fit)), 10L)),
        r_segments = hash_object(lapply(nums, signif, 10L)),
        h_cor_gc = hash_object(res$corrected$data$cor.gc),
        h_ideal = hash_object(res$corrected$data$ideal),
        h_outliers = hash_object(res$corrected$outliers),
        h_fixef = hash_object(unlist(glmmTMB::fixef(res$corrected$fit))),
        h_segments = hash_object(res$isa$segments),
        h_bedpe = hash_object(readLines(res$bedpe))
    )
}
read_section <- function(path, section) {
    x <- readLines(path, warn = FALSE); i <- match(paste0("[", section, "]"), x)
    if (is.na(i) || i == length(x)) return(character())
    x <- x[(i + 1L):length(x)]; j <- which(startsWith(x, "[")); if (length(j)) x <- x[seq_len(j[1L] - 1L)]
    x <- x[grepl(":", x)]; stats::setNames(sub("^[^:]*:\\s*", "", x), sub(":.*$", "", x))
}
write_outputs <- function(env, fp, opts) {
    dir.create(opts$out, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(do.call(rbind, env$rows), file.path(opts$out, "RESULTS.csv"), row.names = FALSE)
    meta <- c(package_version = as.character(utils::packageVersion("MaMaMIA")), R_version = as.character(getRversion()), platform = R.version$platform, glmmTMB_version = as.character(utils::packageVersion("glmmTMB")), DNAcopy_version = as.character(utils::packageVersion("DNAcopy")), reps = opts$reps, seed = 1L)
    checks <- do.call(rbind, env$checks)
    writeLines(c("# MaMaMIA pipeline fingerprint", "[env]", sprintf("%s: %s", names(meta), meta), "[outputs]", sprintf("%s: %s", names(fp), fp), "[checks]", sprintf("%s: %s", gsub(":", "_", checks$check), checks$result)), file.path(opts$out, "FINGERPRINT.txt"))
}

main <- function() {
    suppressPackageStartupMessages(library(MaMaMIA))
    env <- new.env(parent = emptyenv()); env$rows <- env$checks <- list(); set.seed(1L)
    cat("== pipeline ==\n"); res <- run_pipeline(env, opts$cores); print_section(env, "pipeline")
    cat("\n== correctReadCounts() internals ==\n"); internals <- correction_internals(env, opts$cores); print_section(env, "correction")
    cat("\n== prediction reference ==\n"); ref <- predict_reference(env, internals$fit, internals$dat, internals$gc_ref); print_section(env, "predict")

    fast <- internals$pred
    max_rel <- max(abs(fast$actual - ref$actual) / pmax(abs(ref$actual), 1e-300))
    check(env, "fast prediction matches glmmTMB::predict()", max_rel < 1e-6, sprintf("max rel %.3e", max_rel))
    rows <- do.call(rbind, env$rows)
    whole <- rows$elapsed_s[rows$step == sprintf("correctReadCounts(cores=%d)", opts$cores)]
    parts <- sum(rows$elapsed_s[rows$section == "correction" & grepl("^  ", rows$step)])
    check(env, "piecewise internals account for the whole call", abs(parts - whole) / whole < 0.15, sprintf("whole %.2f s vs pieces %.2f s", whole, parts))
    repeat_fit <- suppressMessages(correctReadCounts(res$rca, cores = opts$cores, verbose = FALSE))
    check(env, "correctReadCounts() is repeatable", identical(hash_object(repeat_fit$data$cor.gc), hash_object(res$corrected$data$cor.gc)))

    pairs <- stats::setNames(res$subset_cd$meta$rec_chr_ids, res$subset_cd$meta$don_chr_ids)
    set.seed(42L); s1 <- segments(res$subset_cd, target_pairs = pairs, verbose = FALSE)
    set.seed(99L); s2 <- segments(res$subset_cd, target_pairs = pairs, verbose = FALSE)
    check(env, "segments() is independent of caller RNG", identical(s1$segments, s2$segments), sprintf("%d segments", nrow(s1$segments)))
    set.seed(123L); baseline <- stats::runif(3); set.seed(123L); invisible(segments(res$subset_cd, target_pairs = pairs, verbose = FALSE))
    check(env, "segments() leaves RNG stream untouched", identical(baseline, stats::runif(3)))
    f1 <- tempfile(fileext = ".bedpe"); f2 <- tempfile(fileext = ".bedpe")
    suppressMessages(writeToBedpe(res$isa, f1)); suppressMessages(writeToBedpe(res$isa, f2))
    check(env, "writeToBedpe() is byte-stable", identical(readLines(f1), readLines(f2)), sprintf("%d data rows", length(readLines(f1)) - 1L))

    fps <- list(fingerprint(res, opts$cores))
    for (i in seq_len(opts$reps - 1L)) {
        set.seed(1L); cat(sprintf("repeat %d of %d ...\n", i + 1L, opts$reps))
        fps[[i + 1L]] <- fingerprint(run_pipeline(env, opts$cores, paste0("repeat", i)), opts$cores)
    }
    fp <- fps[[1L]]
    check(env, sprintf("pipeline output identical across %d run(s)", opts$reps), all(vapply(fps[-1L], identical, logical(1L), fp)))

    if (!is.na(opts$compare)) {
        ref_fp <- read_section(opts$compare, "outputs"); if (!length(ref_fp)) abort("no [outputs] section in ", opts$compare)
        diff <- data.frame(key = union(names(fp), names(ref_fp)), stringsAsFactors = FALSE)
        diff$previous <- ref_fp[diff$key]; diff$current <- fp[diff$key]; diff$match <- diff$previous == diff$current
        info <- diff$key %in% c("h_cor_gc", "h_fixef", "h_segments")
        if (any(!diff$match)) print(diff[!diff$match, ], row.names = FALSE)
        check(env, "outputs unchanged versus reference fingerprint", !any(!diff$match & !info), sprintf("%d gating keys differ", sum(!diff$match & !info)))
    }

    print(do.call(rbind, env$checks), row.names = FALSE)
    write_outputs(env, fp, opts)
    if (any(vapply(env$checks, function(x) x$result != "PASS", logical(1L)))) quit(status = 1L)
    cat("\nall consistency checks passed\n")
}
main()
