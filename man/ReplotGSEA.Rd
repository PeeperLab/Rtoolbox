\name{replotGSEA}
\alias{replotGSEA}
\title{
Replot output from the desktop tool javaGSEA.
}
\description{
Creates GSEA plots from the output of the desktop tool javaGSEA.
}
\usage{
replotGSEA(path, gene.set, class.name, metric.range, enrichment.score.range)
}
\arguments{
  \item{path}{path to the output folder of javaGSEA}
  \item{gene.set}{name of the gene set to be plotted. Approximate matching is
  used to determine the right gene set.}
  \item{class.name}{name of the variable to which gene expression has been
  matched (i.e., the description of the variable that is specified in the .cls
  file or that is used to pre-rank the data; example: "drug treatment").}
  \item{metric.range}{optional; the range of the ranking metric that is plotted.
  Defaults to [min(ranking metric), max(ranking metric)]}
  \item{enrichment.score.range}{optional; the range of the enrichment score that is
  plotted. Defaults to [min(enrichment score), max(enrichment score)]}
}
\details{
\code{\link{replotGSEA}} will recreate a GSEA plot from the output of the
desktop tool javaGSEA.
}
\author{
Thomas Kuilman (t.kuilman@nki.nl)
}
\examples{
## Recreate GSEA plot:

\dontrun{replotGSEA(path = path, gene.set = "V$AP1_Q2",
                    class.name = "Drug treatment", metric.range = c(-1, 1))}
}
