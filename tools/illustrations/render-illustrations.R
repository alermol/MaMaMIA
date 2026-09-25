# Render optional vignette figures into vignettes/figures/.
# Run from package root; use --install for optional plotting dependencies.

pkgs <- c("ggplot2", "dplyr", "ggtext", "gggenes", "ggdist", "ggraph", "tidygraph", "circlize", "patchwork", "ragg")
args <- commandArgs(TRUE)
if ("--install" %in% args) install.packages(pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)], repos = "https://packagemanager.posit.co/cran/latest")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing illustration packages: ", paste(missing, collapse = ", "), call. = FALSE)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

OUT <- "vignettes/figures"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PAL <- c(don = "#d6604d", rec = "#4393c3", ink = "#333333", soft = "#9e9e9e", bg = "#f0f0f0")
DATA_CACHE <- "/tmp/mamamia-illus.rds"; FIT_CACHE <- "/tmp/mamamia-fit.rds"

th <- function(base = 11) theme_minimal(base_size = base) + theme(
    panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey93"),
    plot.title = ggtext::element_markdown(face = "bold", size = base + 2),
    plot.subtitle = ggtext::element_markdown(colour = "grey30", size = base - 1, lineheight = 1.05),
    plot.caption = ggtext::element_markdown(colour = "grey40", hjust = 0, size = base - 2),
    legend.position = "bottom"
)
save_fig <- function(p, file, w, h) {
    f <- file.path(OUT, file)
    ggsave(f, p, width = w, height = h, dpi = 200, bg = "white", device = ragg::agg_png)
    cat("wrote", f, "\n")
}
seg_lwd <- function(x) ifelse(x, 2.4, 0.4)
pair_filter <- function(df, i) df[df$chr_id.don == pairs$chr_id.don[i] & df$chr_id.rec == pairs$chr_id.rec[i], , drop = FALSE]

build_data <- function() {
    suppressPackageStartupMessages(library(MaMaMIA))
    data("triticum", package = "MaMaMIA", envir = environment())
    cd <- correctReadCounts(RCA(don_cov, don_gc, rec_cov, rec_gc, meta), cores = 4L, verbose = FALSE, cache = FIT_CACHE)
    sc <- subset(cd, don_subgenomes = "At", rec_subgenomes = "A")
    isa <- segments(sc, target_pairs = stats::setNames(sc$meta$rec_chr_ids, sc$meta$don_chr_ids), verbose = FALSE)
    saveRDS(list(cd = cd, sc = sc, isa = isa), DATA_CACHE)
}
if (!file.exists(DATA_CACHE) || "--rebuild-data" %in% args) build_data()
x <- readRDS(DATA_CACHE); cd <- x$cd; sc <- x$sc; isa <- x$isa
pairs <- isa$meta$target_pairs
chrs <- unique(sc$data[, c("chr_id", "chr_name")]); chrs$len <- as.numeric(tapply(sc$data$stop, sc$data$chr_id, max)[chrs$chr_id])
don_names <- chrs$chr_name[match(pairs$chr_id.don, chrs$chr_id)]
rec_names <- chrs$chr_name[match(pairs$chr_id.rec, chrs$chr_id)]

pair_out <- function(i) pair_filter(isa$out, i)
pair_seg <- function(i) pair_filter(isa$segments, i)
segment_coords <- function(i) {
    s <- pair_seg(i); o <- pair_out(i)
    if (!nrow(s)) return(s)
    transform(s,
        x0 = o$start.don[match(loc.start, o$iid)] / 1e6,
        x1 = o$stop.don[match(loc.end, o$iid)] / 1e6,
        value = s[[isa$meta$metrics]], lwd = seg_lwd(putative_introgression)
    )
}
profile_plot <- function(pts, seg, title, subtitle, xlab = "Window along the chromosome") {
    ggplot(pts, aes(.data$x, .data$diff)) +
        geom_hline(yintercept = 0, colour = "grey55", linetype = "dashed") +
        geom_point(size = 0.25, colour = "grey35", alpha = 0.7) +
        geom_segment(data = seg, aes(x = .data$x0, xend = .data$x1, y = .data$value, yend = .data$value, colour = .data$putative_introgression, linewidth = .data$lwd)) +
        scale_linewidth_identity() +
        scale_colour_manual(NULL, values = c("FALSE" = "grey60", "TRUE" = PAL[["don"]]), labels = c("CBS segment", "putative introgression")) +
        labs(title = title, subtitle = subtitle, x = xlab, y = "Donor - recipient coverage") + th()
}

fig_mechanism <- function(i = 1L) {
    o <- pair_out(i); s <- pair_seg(i); call <- s[s$putative_introgression, , drop = FALSE]
    xr <- range(o$iid); lo <- min(call$loc.start); hi <- max(call$loc.end)
    arrow_df <- data.frame(xmin = c(xr[1], lo), xmax = c(xr[2], hi), y = 1, forward = TRUE, type = c("recipient background", "introgressed segment"))
    schema <- ggplot(arrow_df, aes(xmin = .data$xmin, xmax = .data$xmax, y = .data$y, forward = .data$forward, fill = .data$type)) +
        geom_gene_arrow(arrowhead_height = unit(7, "mm"), arrowhead_width = unit(7, "mm")) +
        scale_fill_manual(NULL, values = c("recipient background" = "grey85", "introgressed segment" = PAL[["don"]])) +
        annotate("text", x = xr[1] + diff(xr) * 0.30, y = 1, label = "recipient background", size = 3, colour = "grey25") +
        annotate("text", x = mean(c(lo, hi)), y = 1.62, label = "introgressed segment", size = 3.2, colour = PAL[["don"]]) +
        annotate("segment", x = mean(c(lo, hi)), xend = mean(c(lo, hi)), y = 1.50, yend = 1.22, colour = PAL[["don"]], arrow = arrow(length = unit(0.10, "cm"))) +
        coord_cartesian(xlim = xr, ylim = c(0.6, 1.9), expand = FALSE) + theme_genes() +
        theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank(), plot.margin = margin(0, 12, 0, 12))
    cov <- data.frame(x = o$iid, donor = o$cor.gc.don, recipient = -o$cor.gc.rec)
    ymax <- max(abs(unlist(cov[-1]))) * 1.08
    prof <- ggplot(cov) +
        geom_rect(data = data.frame(xmin = lo, xmax = hi), aes(xmin = .data$xmin, xmax = .data$xmax), ymin = -ymax, ymax = ymax, fill = PAL[["don"]], alpha = 0.10, inherit.aes = FALSE) +
        geom_area(aes(.data$x, .data$donor), fill = PAL[["don"]], alpha = 0.55) +
        geom_area(aes(.data$x, .data$recipient), fill = PAL[["rec"]], alpha = 0.55) +
        geom_hline(yintercept = 0, colour = "grey45", linewidth = 0.3) +
        annotate("text", x = xr[1] + diff(xr) * 0.02, y = c(ymax, -ymax) * 0.78, label = c("donor coverage", "recipient coverage"), hjust = 0, size = 3.2, colour = c(PAL[["don"]], PAL[["rec"]])) +
        annotate("text", x = mean(c(lo, hi)), y = ymax * 0.88, label = "donor above recipient:\nthis interval is called", size = 2.9, colour = PAL[["ink"]], lineheight = 1.1) +
        coord_cartesian(xlim = xr, ylim = c(-ymax, ymax), expand = FALSE) +
        scale_y_continuous(labels = function(v) round(abs(v))) +
        labs(x = "Window along the chromosome", y = "GC-corrected coverage\n(donor above, recipient below)") + th()
    (schema / prof) + plot_layout(heights = c(1, 4.2)) +
        plot_annotation(title = paste0("An introgression shows up as donor excess: *", don_names[i], "* against *", rec_names[i], "*")) &
        theme(plot.title = ggtext::element_markdown(face = "bold", size = 13), plot.margin = margin(6, 12, 6, 8))
}

fig_workflow <- function() {
    labels <- c(
        "Sample\nrecipient chromosome carrying\na donor-derived segment",
        "Reads mapped competitively\nto a concatenated donor + recipient reference",
        "Windowed tables\nread counts and GC content, one table per genome",
        "GC-bias correction\nzero-inflated negative binomial model, giving cor.gc",
        "Difference profile\ndonor minus recipient, on aligned windows",
        "Circular binary segmentation\nchange points along the profile",
        "Putative introgressions\nsegments whose difference passes the threshold"
    )
    nodes <- data.frame(name = seq_along(labels), label = labels, kind = c("input", rep("step", 5), "result"), x = 1, y = rev(seq_along(labels)))
    g <- tidygraph::tbl_graph(nodes = nodes, edges = data.frame(from = 1:6, to = 2:7), node_key = "name")
    ggraph::ggraph(ggraph::create_layout(g, layout = "manual", x = nodes$x, y = nodes$y)) +
        ggraph::geom_edge_link(arrow = arrow(length = unit(2.6, "mm"), type = "closed"), end_cap = ggraph::circle(2.5, "mm"), colour = "grey55", edge_width = 0.5) +
        ggraph::geom_node_label(aes(label = .data$label, fill = .data$kind), colour = "grey15", size = 2.9, lineheight = 1.05, label.r = unit(2, "mm"), label.padding = unit(2.4, "mm")) +
        scale_fill_manual(NULL, values = c(input = "#dbe9f6", step = "grey96", result = "#fdecea")) +
        coord_cartesian(xlim = c(0.4, 1.6), ylim = c(0.3, 7.7), expand = FALSE) +
        labs(title = "What MaMaMIA does to the sequencing data", caption = "Each step is a package function; the last output is the call to inspect.") +
        theme_void(base_size = 11) + theme(plot.title = ggtext::element_markdown(face = "bold", size = 13, hjust = 0.5), plot.caption = ggtext::element_markdown(colour = "grey40", hjust = 0.5), legend.position = "none")
}

fig_circos <- function() {
    regions <- data.frame(chr = c(don_names, rec_names), start = 0, end = chrs$len[match(c(pairs$chr_id.don, pairs$chr_id.rec), chrs$chr_id)])
    diff_all <- do.call(rbind, lapply(seq_along(don_names), function(i) {
        o <- pair_out(i)
        rbind(data.frame(chr = don_names[i], start = o$start.don, end = o$stop.don, value = o$diff), data.frame(chr = rec_names[i], start = o$start.rec, end = o$stop.rec, value = o$diff))
    }))
    lim <- stats::quantile(abs(diff_all$value), 0.99, na.rm = TRUE); diff_all$value <- pmax(pmin(diff_all$value, lim), -lim)
    f <- file.path(OUT, "03-circos.png"); agg_png(f, width = 2000, height = 1800, res = 200, background = "white"); on.exit(invisible(dev.off()), add = TRUE)
    par(mar = c(1, 1, 1, 1)); circos.par(start.degree = 180, gap.degree = 2, points.overflow.warning = FALSE, track.margin = c(0.004, 0.004), cell.padding = c(0, 0, 0, 0))
    circos.genomicInitialize(regions, plotType = NULL)
    circos.track(ylim = c(0, 1), bg.border = NA, track.height = mm_h(14), panel.fun = function(x, y) {
        s <- CELL_META$sector.index
        circos.rect(CELL_META$xlim[1], 0.15, CELL_META$xlim[2], 0.85, col = if (s %in% don_names) "#737373" else "#c9c9c9", border = "grey30", lwd = 0.4)
        circos.text(mean(CELL_META$xlim), 0.5, s, facing = "bending.inside", niceFacing = TRUE, cex = 0.7, col = "white")
    })
    circos.genomicTrack(diff_all, ylim = c(-lim, lim), track.height = mm_h(26), bg.border = "grey85", panel.fun = function(region, value, ...) {
        circos.lines(CELL_META$xlim, c(0, 0), col = "grey55", lwd = 0.5)
        circos.genomicRect(region, value, ytop.column = 1, ybottom = 0, col = colorRamp2(c(-lim, 0, lim), c("#2166ac", "#f7f7f7", "#b2182b"))(value[[1]]), border = NA)
    })
    called <- isa$segments[isa$segments$putative_introgression, , drop = FALSE]
    for (i in seq_len(nrow(called))) {
        s <- called[i, ]; o <- pair_filter(isa$out, match(s$chr_id.don, pairs$chr_id.don))
        circos.link(chrs$chr_name[match(s$chr_id.don, chrs$chr_id)], c(o$start.don[match(s$loc.start, o$iid)], o$stop.don[match(s$loc.end, o$iid)]), chrs$chr_name[match(s$chr_id.rec, chrs$chr_id)], c(o$start.rec[match(s$loc.start, o$iid)], o$stop.rec[match(s$loc.end, o$iid)]), col = "#cb181d33", border = "#cb181d")
    }
    circos.clear()
    grid::grid.text("At/A chromosomes: coverage difference and the segments called", y = 0.975, gp = grid::gpar(fontsize = 13, fontface = "bold"))
    grid::grid.text("outer ring: donor (dark) and recipient (light) chromosomes | inner ring: donor - recipient coverage | ribbons: putative introgressions", y = 0.945, gp = grid::gpar(fontsize = 8, col = "grey35"))
    cat("wrote", f, "\n")
}

fig_gc_bias <- function() {
    d <- cd$data[cd$data$ideal, ]; d$subgenome <- as.character(d$subgenome)
    bin_one <- function(dd) { b <- cut(dd$gc, seq(min(dd$gc), max(dd$gc), length.out = 13), include.lowest = TRUE); data.frame(gc = tapply(dd$gc, b, mean), value = tapply(dd$cov, b, mean), se = tapply(dd$cov, b, function(v) stats::sd(v) / sqrt(length(v))), n = tapply(dd$cov, b, length), subgenome = dd$subgenome[1]) }
    pts <- do.call(rbind, lapply(sort(unique(d$subgenome)), function(s) bin_one(d[d$subgenome == s, ])))
    grid <- expand.grid(gc = seq(min(d$gc), max(d$gc), length.out = 200), subgenome = sort(unique(d$subgenome)))
    grid$value <- MaMaMIA:::eta_zinb_response(cd$fit, grid)
    ggplot(pts, aes(.data$gc, .data$value)) +
        geom_errorbar(aes(ymin = .data$value - .data$se, ymax = .data$value + .data$se), width = 0.0012, colour = "grey45", linewidth = 0.35) +
        geom_point(aes(size = .data$n), colour = "grey15", fill = "white", shape = 21) +
        geom_line(data = grid, colour = PAL[["don"]], linewidth = 0.9) +
        scale_size_continuous("windows per bin", range = c(0.8, 3.4), breaks = c(50, 500, 1500)) +
        facet_wrap(~subgenome, scales = "free_y", nrow = 1) +
        labs(title = "The GC response the model removes", subtitle = "Binned ideal-window coverage versus the fitted ZINB mean; all subgenomes peak near 45-46% GC.", caption = "correctReadCounts() divides each window by this fitted response, putting equal true coverage on one scale.", x = "GC content of the window", y = "Coverage") + th() +
        theme(legend.position = "right", strip.text = element_text(face = "bold"))
}

fig_profile <- function(i = 1L) {
    o <- pair_out(i); s <- segment_coords(i)
    profile_plot(data.frame(x = o$start.don / 1e6, diff = o$diff), s, paste0("Segmentation of one pair: *", don_names[i], "* / ", rec_names[i]), paste0("Bars are the ", isa$meta$metrics, " of each CBS segment; the red one passes the threshold"), "Donor position (Mb)")
}

fig_orientation <- function(i = 7L) {
    collect <- function(obj, label) {
        o <- pair_filter(obj$out, i); s <- pair_filter(obj$segments, i)
        s$value <- s[[obj$meta$metrics]]; s$x0 <- s$loc.start; s$x1 <- s$loc.end; s$lwd <- seg_lwd(s$putative_introgression)
        list(pts = data.frame(x = o$iid, diff = o$diff, orientation = label), seg = transform(s, orientation = label))
    }
    a <- collect(isa, "as assembled"); b <- collect(MaMaMIA::reverseWindows(isa, chr_ids = pairs$chr_id.don[i]), "donor chromosome reversed")
    pts <- rbind(a$pts, b$pts); seg <- rbind(a$seg, b$seg); lv <- c("as assembled", "donor chromosome reversed")
    pts$orientation <- factor(pts$orientation, lv); seg$orientation <- factor(seg$orientation, lv)
    profile_plot(pts, seg, paste0("Assembly orientation moves the signal: *", don_names[i], "* / ", rec_names[i]), sprintf("The donor coverage vector is reversed before segmentation; calls change from %d/%d to %d/%d segments.", sum(a$seg$putative_introgression), nrow(a$seg), sum(b$seg$putative_introgression), nrow(b$seg))) + facet_wrap(~orientation, ncol = 1)
}

fig_arrowmap <- function() {
    rows <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
        s <- segment_coords(i); s <- s[s$putative_introgression, , drop = FALSE]
        one <- function(nm, len) rbind(data.frame(chr = nm, start = 0, end = len, forward = TRUE, type = "chromosome"), if (nrow(s)) data.frame(chr = nm, start = s$x0 * 1e6, end = s$x1 * 1e6, forward = TRUE, type = "putative introgression"))
        rbind(one(don_names[i], chrs$len[match(pairs$chr_id.don[i], chrs$chr_id)]), one(rec_names[i], chrs$len[match(pairs$chr_id.rec[i], chrs$chr_id)]))
    }))
    rows$chr <- factor(rows$chr, unique(rows$chr)); rows$type <- factor(rows$type, c("chromosome", "putative introgression"))
    ggplot(rows, aes(xmin = .data$start / 1e6, xmax = .data$end / 1e6, y = .data$chr, fill = .data$type, forward = .data$forward)) +
        geom_gene_arrow(arrowhead_height = unit(2.2, "mm"), arrowhead_width = unit(2.2, "mm")) +
        scale_fill_manual(NULL, values = c("chromosome" = "grey80", "putative introgression" = PAL[["don"]])) + theme_genes() +
        labs(title = "Introgression map: each chromosome as a gene-arrow track", subtitle = "The same 14 chromosomes as in the other views; red arrows are the segments passing the threshold", x = "Position (Mb)", y = NULL) +
        theme(plot.title = ggtext::element_markdown(face = "bold", size = 12), plot.subtitle = ggtext::element_markdown(colour = "grey30", size = 9), legend.position = "bottom")
}

fig_distributions <- function() {
    diffs <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) data.frame(pair = paste0(don_names[i], " / ", rec_names[i]), diff = pair_out(i)$diff)))
    diffs$pair <- factor(diffs$pair, names(sort(tapply(diffs$diff, diffs$pair, stats::median))))
    ggplot(diffs, aes(.data$diff, .data$pair)) +
        ggdist::stat_halfeye(.width = c(0.5, 0.95), fill = PAL[["rec"]], colour = "grey30", slab_alpha = 0.75, point_colour = "grey15", point_size = 1.1, normalize = "groups") +
        geom_vline(xintercept = 0, colour = "grey40", linetype = "dashed") +
        labs(title = "Where the donor-recipient difference actually sits", subtitle = "Distribution of the per-window difference for each pair, ordered by median; the point is the median with 50% and 95% intervals", x = "Donor - recipient coverage", y = NULL) + th()
}

save_fig(fig_mechanism(), "01-mechanism.png", 9, 6)
save_fig(fig_workflow(), "02-workflow.png", 8.5, 7)
fig_circos()
save_fig(fig_gc_bias(), "04-gc-bias.png", 9, 6.5)
save_fig(fig_profile(), "05-pair-profile.png", 9, 4)
save_fig(fig_orientation(), "06-orientation.png", 9, 6)
save_fig(fig_arrowmap(), "07-arrowmap.png", 9, 7)
save_fig(fig_distributions(), "08-distributions.png", 9, 5)
