\name{ReplotGSEA}
\alias{ReplotGSEA}
\title{
Replot output from the desktop tool javaGSEA.
}
\description{
Creates GSEA plots from the output of the desktop tool javaGSEA.
}
\usage{
ReplotGSEA(path, gene.set, class.name)
}
\arguments{
  \item{path}{path to the output folder of javaGSEA}
  \item{gene.set}{name of the gene set to be plotted. Approximate matching is
  used to determine the right gene set.}
  \item{class.name}{name of the variable to which gene expression has been
  matched (i.e., the description of the variable that is specified in the .cls
  file or that is used to pre-rank the data; example: "drug treatment").}
}
\details{
\code{\link{ReplotGSEA}} will recreate a GSEA plot from the output of the
desktop tool javaGSEA.
}
\author{
Thomas Kuilman (t.kuilman@nki.nl)
}
\examples{
## Generate heatmaps:

\dontrun{ReplotGSEA(path = path, gene.set = "V$AP1_Q2",
class.name = "Drug treatment")}