\name{plot.ISA}
\alias{plot.ISA}
\title{Plot method for ISA objects}
\usage{
\method{plot}{ISA}(x, plot.type = c("mirror", "diff"), all_segments = FALSE, ...)
}
\arguments{
\item{x}{An object of class \code{ISA}, as returned by \code{segments()}.}

\item{plot.type}{Character string selecting the plot style. One of
\code{"mirror"} (donor coverage above zero and recipient coverage below) or
\code{"diff"} (coverage difference with horizontal segment summary metrics).}

\item{all_segments}{Logical; if \code{FALSE} (default), only segments marked
\code{putative_introgression} are drawn on difference plots. If
\code{TRUE}, all segments are shown.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{plot()}.}
}
\value{
A \code{ggplot2::ggplot()} object.
}
\description{
Draws pairwise donor-recipient profiles from an \code{ISA} object, either as a
mirrored plot of the selected per-window values (normalized dosage when the
ISA was segmented on \code{dose}, otherwise GC-corrected coverage) or as a
difference profile with segment summaries (mean or median, per
\code{x$meta$metrics}) overlaid.
}
\examples{
data(triticum)
rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
\donttest{
## Note: correctReadCounts() can take a while on the full dataset
rca <- correctReadCounts(rca, cores = 4L, verbose = FALSE)
}

## For demonstration purposes, assume that GC bias correction was performed, 
## by copying the 'cov' column to 'cor.gc'.
rca$data$valid <- TRUE
gc_q <- quantile(rca$data$gc[rca$data$valid], c(0.01, 0.99))
cov_q <- quantile(rca$data$cov[rca$data$valid], 0.99)
rca$data$ideal <- rca$data$valid & 
                    (rca$data$gc >= gc_q[1]) & 
                    (rca$data$gc <= gc_q[2]) & 
                    (rca$data$cov <= cov_q)
rca$corrected <- TRUE
rca$outliers <- list(
        gc_lower_bound = gc_q[1],
        gc_upper_bound = gc_q[2],
        cov_upper_bound = cov_q
    )
rca$data$cor.gc <- rca$data$cov

isa <- segments(
  rca,
  target_pairs = c("OY997261.1" = "CP155614.1"),
  verbose = FALSE
)
plot(isa, plot.type = "mirror")
plot(isa, plot.type = "diff")
plot(isa, plot.type = "diff", all_segments = TRUE)
}
