## Shared synthetic tables for segmentation, calibration and pipeline tests.

synthetic_score_table <- function(n = 200L, groups = c("1", "2"),
                                  sample = "S1", mode = "single",
                                  event_group = "1", event_from = 50L,
                                  event_to = 90L, event_effect = 1.5,
                                  noise = 0.05, seed = 1L) {
    set.seed(seed)
    rows <- list()
    for (group in groups) {
        iid <- seq.int(0L, n - 1L)
        score <- stats::rnorm(n, 0, noise)
        if (!is.null(event_group) && identical(group, event_group)) {
            score[iid >= event_from & iid <= event_to] <-
                score[iid >= event_from & iid <= event_to] + event_effect
        }
        rows[[length(rows) + 1L]] <- data.frame(
            sample = sample, group = group, iid = iid,
            start = iid * 1e6, stop = (iid + 1L) * 1e6,
            mode = mode, metric = "breadth", score_type = "paired_log_ratio",
            score = score,
            donor_mean = 0.1 * exp(score / 2),
            recipient_mean = rep(0.1, n),
            donor_component = "S", n_donor = 1L, n_recipient = 1L,
            donor_copies = NA_real_, donor_fraction = NA_real_,
            chr_name = paste0(group, "S"), length = n * 1e6,
            n_per_chromosome = n,
            gc_donor = 0.45, gc_recipient = 0.45,
            log_donor = log(0.1) + score,
            log_recipient = rep(log(0.1), n),
            donor_residual = score,
            recipient_residual = stats::rnorm(n, 0, noise),
            lambda_copy = 0.1,
            stringsAsFactors = FALSE
        )
    }
    do.call(rbind, rows)
}

synthetic_window_metrics <- function(n = 200L, groups = c("1", "2"),
                                     samples = "S1",
                                     recipient_components = c("A"),
                                     donor_components = "S",
                                     donor_lambda = 0.1,
                                     recipient_lambda = 0.1,
                                     gc_effect = 0,
                                     seed = 1L) {
    set.seed(seed)
    rows <- list()
    for (sample in samples) {
        for (group in groups) {
            iid <- seq.int(0L, n - 1L)
            gc <- 0.4 + 0.1 * sin(iid / 20)
            for (component in donor_components) {
                lambda <- donor_lambda * exp(
                    gc_effect * ((gc - 0.5) / 0.1)^2
                ) * exp(stats::rnorm(n, 0, 0.02))
                rows[[length(rows) + 1L]] <- data.frame(
                    sample = sample, chr_id = paste0("d", group, component),
                    chr_name = paste0(group, component), group = group,
                    role = "donor", component = component, copy_number = 2,
                    iid = iid, start = iid * 1e6, stop = (iid + 1) * 1e6,
                    length = 1e6, count = round(lambda * 1e6 / 100),
                    covered = pmin(1e6, round((1 - exp(-lambda)) * 1e6)),
                    breadth = 1 - exp(-lambda), gc = gc,
                    lambda_breadth = lambda, lambda_count = lambda,
                    read_length = 100, stringsAsFactors = FALSE
                )
            }
            for (component in recipient_components) {
                lambda <- recipient_lambda * exp(
                    gc_effect * ((gc - 0.5) / 0.1)^2
                ) * exp(stats::rnorm(n, 0, 0.02))
                rows[[length(rows) + 1L]] <- data.frame(
                    sample = sample, chr_id = paste0("r", group, component),
                    chr_name = paste0(group, component), group = group,
                    role = "recipient", component = component,
                    copy_number = 2, iid = iid, start = iid * 1e6,
                    stop = (iid + 1) * 1e6, length = 1e6,
                    count = round(lambda * 1e6 / 100),
                    covered = pmin(1e6, round((1 - exp(-lambda)) * 1e6)),
                    breadth = 1 - exp(-lambda), gc = gc,
                    lambda_breadth = lambda, lambda_count = lambda,
                    read_length = 100, stringsAsFactors = FALSE
                )
            }
        }
    }
    do.call(rbind, rows)
}
