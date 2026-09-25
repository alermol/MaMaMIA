\name{exploreIntrogressions}
\alias{exploreIntrogressions}
\title{Explore introgression calls in an interactive dashboard}
\usage{
exploreIntrogressions(
  ISA,
  chr_ids = NULL,
  host = "127.0.0.1",
  port = NULL,
  launch.browser = interactive()
)
}
\arguments{
\item{ISA}{An object of class \code{ISA}, as returned by \code{segments()}.}

\item{chr_ids}{Optional chromosome IDs to show initially. Defaults to all
chromosomes present in \code{ISA}.}

\item{host}{Host passed to \code{shiny::runApp()}.}

\item{port}{Port passed to \code{shiny::runApp()}. Use \code{NULL} to let
Shiny choose.}

\item{launch.browser}{Logical passed to \code{shiny::runApp()}.}
}
\value{
The return value of \code{shiny::runApp()}. Called for its side effect of
launching the dashboard.
}
\description{
Opens a small Shiny application for interactively inspecting an \code{ISA}
object returned by \code{segments()}. The app uses a BioCircos genomic circos
plot, searchable chromosome and chromosome-pair filters, a sortable \code{DT}
segment table, and exports the currently selected table rows with
\code{writeToBedpe()}.
}
\details{
Dashboard packages are optional MaMaMIA dependencies. Install them with
\code{install.packages(c("shiny", "bslib", "shinyWidgets", "DT", "BioCircos"))}
before calling this function.
}
\examples{
\dontrun{
data(triticum)
rca <- RCA(don_cov, don_gc, rec_cov, rec_gc, meta)
rca <- correctReadCounts(rca, cores = 4L, verbose = FALSE)
isa <- segments(rca, target_pairs = c("OY997261.1" = "CP155614.1"), verbose = FALSE)
exploreIntrogressions(isa)
}
}
