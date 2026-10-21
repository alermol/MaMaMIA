\name{reverseWindows.RCA}
\alias{reverseWindows.RCA}
\alias{reverseWindows}
\alias{reverseWindows.ISA}
\title{Reverse the effective orientation of selected chromosomes}
\usage{
\method{reverseWindows}{RCA}(x, chr_ids, ...)

reverseWindows(x, chr_ids, ...)

\method{reverseWindows}{ISA}(x, chr_ids, ...)
}
\arguments{
\item{x}{An object of class \code{RCA} or \code{ISA}.}

\item{chr_ids}{Character vector of chromosome IDs in the object's meta table
(\code{chr_id} column) whose orientation should be reversed. The short display
name (\code{chr_name}), used mainly for plotting labels, must not be used here.}

\item{...}{Further arguments passed to methods.}
}
\value{
The modified object, of the same class as \code{x}.
}
\description{
When comparing two assemblies, homologous chromosomes can be assembled in
opposite orientations. As a result, the per-window coverage and GC vectors of
such chromosomes run in opposite directions along the genome.

\code{reverseWindows()} flips the effective orientation of the requested
chromosomes so that windows can be compared position by position. It keeps
the window indexes (\code{iid}, \code{start}, \code{stop}, and the chromosome identifier)
unchanged, and instead reverses the value vectors (coverage, GC-corrected
coverage, normalized dosage when present, and the per-window validity flags)
within each selected chromosome. The (modified) object is returned.

For an \code{ISA} object, segmentation is \strong{recomputed} for every
donor-recipient pair that contains a reversed chromosome, because these
pairs' difference profiles have changed orientation. The difference is
recomputed from the column recorded in \code{ISA$meta$value_col} (\code{dose}
when normalization was applied, otherwise \code{cor.gc}). Segments of
unaffected pairs are preserved unchanged.
}
\examples{
\dontrun{
rev_rca <- reverseWindows(rca, chr_ids = c("OY997261.1", "CP155614.1"))
rev_isa <- reverseWindows(isa, chr_ids = "OY997261.1")
}
}
