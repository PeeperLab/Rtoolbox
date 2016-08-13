\name{OverviewPlot}
\alias{OverviewPlot}
\alias{c}
\alias{colorRampPalette}
\title{
Create heatmaps of segmentation values from a DNAcopy object.
}
\description{
Creates color-based plots of the segmentation values from a provided DNAcopy object.
}
\usage{
OverviewPlot(DNAcopy.object, samples, range.CNA = c(-2, 2), 
              color.palette = colorRampPalette(c("blue", "white", "red"))(49))
}
\arguments{
  \item{DNAcopy.object}{DNAcopy object with the segmentation values that need to be plotted.}
  \item{samples}{character vector with the names of the samples that need to be plotted. Names need to be as given in the DNAcopy object. This is an optional argument; all samples will be plotted by default.}
  \item{range.CNA}{the range of the segmentation values that will be plotted. Segmentation values outside of range.CNA will be capped to the minimum and maximum range.CNA values. This is an optional argument; the default range is \code{\link{c}}(-2, 2).}
  \item{color.palette}{the color palette that will be used for heatmaps. This is an optional argument; the default color palette is \code{\link{colorRampPalette}}(c("blue", "white", "red"))(49).}
}
\details{
\code{\link{OverviewPlot}} will print heatmaps of the segmentation values in a DNAcopy object. The samples argument can be used to limit the amount of plotted samples. The range.CNA vector determines where the segmentation values are capped and the range of segmentation values that is plotted. Finally, the color.palette argument can be used to create heatmaps with custom colors. All arguments except DNA.copy.object are optional.
}
\author{
Thomas Kuilman (t.kuilman@nki.nl)
}
\examples{
## Generate heatmaps:

\dontrun{OverviewPlot(segment.CNA.object) ## Plot using default settings.}
\dontrun{OverviewPlot(segment.CNA.object, samples = unique(segment.CNA.object$output$ID)[1:3]) ## Plot only first three samples in segment.CNA.object.}
\dontrun{OverviewPlot(segment.CNA.object, range.CNA = c(-1, 1)) ## Plot with values from -1 to 1 (outside of this range values are capped).}
\dontrun{OverviewPlot(segment.CNA.object, color.palette = colorRampPalette(c("blue", "white", "white", "red"))(49)) ## Plot with extra white in the middle of the color palette.}
}
