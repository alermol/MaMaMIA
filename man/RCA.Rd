\name{RCA}
\alias{RCA}
\title{Create an RCA (Read Count Array) object}
\usage{
RCA(dcov, dgc, rcov, rgc, meta, presorted = TRUE)
}
\arguments{
\item{dcov}{A data frame, tibble, or matrix of donor genome coverage. Must
contain exactly 4 columns (column names are ignored): chromosome ID, start
coordinate of the interval (0-based), end coordinate of the interval
(exclusive), and number of reads mapped in the interval.}

\item{dgc}{A data frame, tibble, or matrix of donor genome GC-content. Must
contain exactly 4 columns (column names are ignored): chromosome ID, start
coordinate (0-based), end coordinate (exclusive), and GC-content of the
interval (values between 0 and 1). Window coordinates must match
\code{dcov} exactly (same chromosomes and start/stop pairs).}

\item{rcov}{A data frame, tibble, or matrix of recipient genome coverage, in
the same format as \code{dcov}.}

\item{rgc}{A data frame, tibble, or matrix of recipient genome GC-content, in
the same format as \code{dgc}. Window coordinates must match \code{rcov}.}

\item{meta}{A data frame, tibble, or matrix of chromosome metadata with
exactly 3 columns (column names are ignored): \code{chr_id} (chromosome IDs
matching those in coverage and GC data for both donor and recipient),
\code{chr_name} (human-readable chromosome names), and \code{subgenome}
(subgenome labels). Must not contain missing values or duplicated
\code{chr_id} values. Donor and recipient chromosome IDs must be disjoint,
and together they must equal the set of \code{chr_id} values in \code{meta}.
Donor and recipient subgenome labels must also be disjoint.}

\item{presorted}{Logical; if \code{TRUE} (default), coverage and GC tables
are used in the order provided. If \code{FALSE}, each table is sorted by
chromosome ID, start, and stop before validation and joining.}
}
\value{
An object of class \code{RCA}: a list with components
\describe{
  \item{\code{data}}{Combined coverage, GC-content, and metadata per window,
    including per-chromosome window index \code{iid} (1, 2, \ldots, \eqn{n}
    along each chromosome). Downstream pairwise comparison matches donor and
    recipient windows by \code{iid}, not by shared genomic coordinates. After
    \code{normalizeReadCounts()}, the table also has a \code{dose} column.}
  \item{\code{meta}}{List with chromosome metadata (\code{meta}) and
    donor/recipient chromosome ID vectors (\code{don_chr_ids},
    \code{rec_chr_ids}).}
  \item{\code{corrected}}{Logical flag; \code{FALSE} until
    \code{correctReadCounts()} is run.}
  \item{\code{normalization}}{Optional list (\code{mode}, \code{reference},
    \code{by_side}, \code{pseudocount}) added by
    \code{normalizeReadCounts()}.}
  \item{\code{fit_sampling}}{Optional list added by
    \code{correctReadCounts()}: \code{sampled}, \code{max_fit_windows},
    \code{n_ideal}, \code{n_fitted}, \code{seed} and a \code{strata} data
    frame describing the subgenome x GC-decile strata used by the fit.}
}
}
\description{
Constructs an object of class \code{RCA} that combines donor and recipient
windowed coverage and GC-content data with chromosome metadata. The resulting
object is used as input for GC-bias correction and subsequent CBS-based
introgression segmentation.
}
\details{
Chromosomes must be partitioned into the \strong{same number of bins}
regardless of physical chromosome length. Equal bin counts are required so
that donor and recipient profiles can be compared position-wise (by window
index) under circular binary segmentation. Bins therefore need not span equal
genomic distances, and donor and recipient windows need not share the same
start/stop coordinates—only the same number of ordered windows per
chromosome.

Within each input table, windows on a chromosome must form contiguous
half-open intervals (no gaps or overlaps) in row order. Within each genome,
coverage and GC tables must use identical window coordinates (\code{dcov}
with \code{dgc}, and \code{rcov} with \code{rgc}); joins that drop rows
cause an error. \code{RCA()} also requires disjoint chromosome IDs and
disjoint subgenome labels between donor and recipient, and stops if any
chromosome has a different window count from the others.
}
\examples{
data(triticum)
rca <- RCA(
  dcov = don_cov,
  dgc = don_gc,
  rcov = rec_cov,
  rgc = rec_gc,
  meta = meta
)
rca
}
