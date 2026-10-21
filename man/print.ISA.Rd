\name{print.ISA}
\alias{print.ISA}
\title{Print method for ISA objects}
\usage{
\method{print}{ISA}(x, ...)
}
\arguments{
\item{x}{An object of class \code{ISA}, as returned by \code{segments()}.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{print()}.}
}
\value{
\code{x}, invisibly.
}
\description{
Prints a one-line overview of an \code{ISA} object (pair count, segment
count, putative introgression count, active metric, and the value column used
for segmentation).
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
print(isa)
}
