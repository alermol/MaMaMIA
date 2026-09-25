\name{segments}
\alias{segments}
\title{Segment pairwise donor-recipient coverage differences}
\usage{
segments(
  RCA,
  target_pairs,
  min_segment_median = 0L,
  min_segment_mean = 0L,
  threshold_type = c("median", "mean"),
  alpha = 0.001,
  min_width = 2,
  undo_SD = 1,
  seed = 1L,
  verbose = TRUE
)
}
\arguments{
\item{RCA}{An object of class \code{RCA} that has been processed with
\code{correctReadCounts()} (\code{corrected} must be \code{TRUE}).}

\item{target_pairs}{Named vector mapping donor chromosome IDs (\code{names})
to recipient chromosome IDs (\code{values}). Only these pairs are segmented.
Must be non-empty; there is no default all-pairs expansion.}

\item{min_segment_median}{Minimum segment median difference used to flag
putative introgressions when \code{threshold_type = "median"}
(default \code{0L}).}

\item{min_segment_mean}{Minimum segment mean difference used to flag
putative introgressions when \code{threshold_type = "mean"}
(default \code{0L}).}

\item{threshold_type}{Which segment summary statistic to use for flagging
putative introgressions: \code{"median"} (default) or \code{"mean"}.}

\item{alpha}{Significance level for the CBS change-point test passed to
\code{DNAcopy::segment()} (default \code{0.001}). Must satisfy
\code{0 < alpha <= 1}. Smaller values yield fewer segments.}

\item{min_width}{Minimum number of windows retained after segmentation
(default \code{2}). Segments shorter than this are dropped. Must be
\code{>= 2}. Note that \code{DNAcopy::segment()} itself is always called
with \code{min.width = 2}; this argument is applied as a post-filter.}

\item{undo_SD}{Number of standard deviations between segment means required
to keep a split under \code{undo.splits = "sdundo"} (default \code{1}).
Must satisfy \code{0 < undo_SD <= 10}. Larger values prune more splits.}

\item{seed}{Integer seed used for the permutation-based p-values computed by
\code{DNAcopy::segment()} (default \code{1L}). CBS is stochastic, so fixing the
seed makes repeated calls on the same object reproducible. The caller's RNG
state is restored afterwards.}

\item{verbose}{Logical; if \code{TRUE} (default), print a progress message.}
}
\value{
An object of class \code{ISA}: a list with components
\describe{
  \item{\code{data}}{\code{RCA} window data restricted to chromosomes
    involved in \code{target_pairs}.}
  \item{\code{out}}{Per-window pairwise data (donor and recipient columns,
    plus \code{diff = cor.gc.don - cor.gc.rec}), restricted to ideal windows.}
  \item{\code{segments}}{Segment summary table for all pairs, with a logical
    \code{putative_introgression} column
    (\code{seg.median} or \code{seg.mean} at or above the chosen threshold),
    annotated with donor and recipient chromosome metadata. Segments shorter
    than \code{min_width} windows are excluded.}
  \item{\code{param}}{Resolved segmentation settings: \code{alpha},
    \code{undo.SD}, and \code{seed}.}
  \item{\code{meta}}{Chromosome metadata restricted to chromosomes in
    \code{target_pairs}, plus \code{don_chr_ids}, \code{rec_chr_ids},
    \code{target_pairs} (a data frame with columns \code{chr_id.don} and
    \code{chr_id.rec}), \code{metrics} (\code{"seg.median"} or
    \code{"seg.mean"}), and \code{metrics_threshold}.}
}
}
\description{
For selected donor-recipient chromosome pairs, computes the difference in
GC-corrected coverage across aligned ideal windows and runs circular binary
segmentation to identify putative introgressed segments. Returns an object
of class \code{ISA} (Introgression Segment Array).
}
\details{
\code{target_pairs} must be supplied explicitly. There is no default
all-pairs expansion of donor and recipient chromosomes.

Pairs with fewer than 2 ideal windows (after filtering) cannot be segmented
by CBS; they are skipped with a warning and contribute no rows to
\code{segments}. Their (possibly empty) pairwise window tables are still
included in \code{out}.
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
isa
}
