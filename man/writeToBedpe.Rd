\name{writeToBedpe}
\alias{writeToBedpe}
\title{Write putative introgressions to a BEDPE file}
\usage{
writeToBedpe(ISA, file, sample = NULL)
}
\arguments{
\item{ISA}{An object of class \code{ISA}, as returned by \code{segments()}.}

\item{file}{Path to the output file. Written as a tab-separated table with a
header row.}

\item{sample}{Optional sample name. When supplied, it is written to an extra
\code{sample} column so that single-sample and cohort exports can be merged
without losing the sample identity. \code{NULL} preserves the legacy column
layout byte for byte.}
}
\value{
Called for its side effect of writing \code{file}. Invisibly returns
\code{NULL}. Prints a message with the output path when finished.
}
\description{
Exports putative introgressed segments from an \code{ISA} object to a BEDPE
(Browser Extensible Data Paired-End) file. Each row is one donor-recipient
segment pair flagged as a putative introgression.
}
\details{
Only segments with \code{putative_introgression = TRUE} are written.
Coordinates are taken from the donor and recipient windows that span each
segment (\code{start}/\code{stop} of the first and last ideal windows).
If no putative introgressions are present, a warning is issued and a
header-only file is written. With \code{sample} supplied, the \code{sample} column is
appended last; the frozen multi-sample format (D11, resolved question 5).

Output columns (in order):
\itemize{
  \item \code{chrom1}, \code{start1}, \code{end1}: donor chromosome ID and
    interval
  \item \code{chrom2}, \code{start2}, \code{end2}: recipient chromosome ID
    and interval
  \item \code{name}: segment label (\code{PIS1}, \code{PIS2}, \ldots)
  \item \code{score}: active segment metric from \code{ISA$meta$metrics}
    (\code{seg.median} or \code{seg.mean})
  \item \code{strand1}, \code{strand2}: always \code{"."}
  \item \code{chrom_name1}, \code{chrom_name2}: human-readable chromosome
    names
  \item \code{subgen_name1}, \code{subgen_name2}: subgenome labels
  \item \code{sample}: sample name, only when \code{sample} is supplied
}
}
\seealso{
\code{\link[=writeCohortToBedpe]{writeCohortToBedpe()}}
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
\donttest{
writeToBedpe(isa, file = 'introgressions.bedpe')
}
}
