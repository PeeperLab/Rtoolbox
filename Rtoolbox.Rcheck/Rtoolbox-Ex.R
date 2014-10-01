pkgname <- "Rtoolbox"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('Rtoolbox')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("OverviewPlot")
### * OverviewPlot

flush(stderr()); flush(stdout())

### Name: OverviewPlot
### Title: Create heatmaps of segmentation values from a DNAcopy object.
### Aliases: OverviewPlot c colorRampPalette

### ** Examples

## Generate heatmaps with \code{OverviewPlot}.

## Not run: OverviewPlot(segment.CNA.object) ## Plot using default settings.
## Not run: OverviewPlot(segment.CNA.object, samples = unique(segment.CNA.object$output$ID)[1:3]) ## Plot only first three samples in segment.CNA.object.
## Not run: OverviewPlot(segment.CNA.object, range.CNA = c(-1, 1)) ## Plot with values from -1 to 1 (outside of this range values are capped).
## Not run: OverviewPlot(segment.CNA.object, color.palette = colorRampPalette(c("blue", "white", "red"))(49)) ## Plot without the extra white in the middle of the palette.



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
