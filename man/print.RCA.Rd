\name{print.RCA}
\alias{print.RCA}
\title{Print method for RCA objects}
\usage{
\method{print}{RCA}(x, ...)
}
\arguments{
\item{x}{An object of class \code{RCA}.}

\item{...}{Additional arguments (ignored); present for compatibility with
the generic \code{print()}.}
}
\value{
\code{x}, invisibly.
}
\description{
Prints a one-line overview of an \code{RCA} object: total chromosome count,
donor and recipient counts, window count, whether GC-bias correction has been
applied, and the normalization mode when \code{normalizeReadCounts()} has
run.
}
\examples{
data(triticum)
rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
print(rca)
}
