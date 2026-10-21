\name{summary.RCA}
\alias{summary.RCA}
\alias{print.summary.RCA}
\title{Summary method for RCA objects}
\usage{
\method{summary}{RCA}(object, ...)

\method{print}{summary.RCA}(x, ...)
}
\arguments{
\item{object}{An object of class \code{RCA}.}

\item{x}{An object of class \code{summary.RCA}.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generics.}
}
\value{
An object of class \code{summary.RCA}: a list with
\code{n_chromosomes}, \code{n_subgenomes},
\code{chromosomes_per_subgenome} (named integer vector of chromosome counts
per subgenome label), \code{corrected}, \code{normalization} (the
\code{mode}, \code{reference}, \code{by_side} and \code{pseudocount} list added
by \code{normalizeReadCounts()}, or \code{NULL}), \code{outliers} (the GC
and coverage bounds from \code{correctReadCounts()}, or \code{NULL} if
correction has not been applied), and \code{fit_sampling} (the fit-window
subsampling metadata from \code{correctReadCounts()}, or \code{NULL} if
correction has not been applied). Printing the object shows a
human-readable summary, including outlier thresholds, the fit-window count and
the normalization mode when available.
}
\description{
Builds a structured summary of chromosome and subgenome counts for an
\code{RCA} object, whether GC-bias correction has been applied, and whether
the object has been normalized to dosage.
}
\examples{
data(triticum)
rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
summary(rca)
s <- summary(rca)
s$n_chromosomes
s$outliers
}
