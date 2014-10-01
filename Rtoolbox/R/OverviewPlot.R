OverviewPlot <- function(DNAcopy.object, samples, range.CNA = c(-2, 2), color.palette = colorRampPalette(c("blue", "white", "red"))(49)) {
	
	## Use all samples by default
	if(missing(samples)) {
		samples <- unique(DNAcopy.object$output$ID)
	}

	## Select samples
	DNAcopy.object$output <- DNAcopy.object$output[DNAcopy.object$output$ID %in% samples, ]
	
	## Cap range
	DNAcopy.object$output$seg.mean[DNAcopy.object$output$seg.mean < range.CNA[1]] <- range.CNA[1]
	DNAcopy.object$output$seg.mean[DNAcopy.object$output$seg.mean > range.CNA[2]] <- range.CNA[2]
	
	## Reshape data.frame according to sample name
	# reshape2::dcast -> to wider data.frame; right-hand of tilde needs to be a 'measured variable' (i.e., needs to go into columns)
	# reshape2::melt -> to narrower data.frame
	# Names are changed by dcast according to level order -> change order levels (!)
	order.samples <- unique(DNAcopy.object$output$ID)
	DNAcopy.object$output$ID <- factor(DNAcopy.object$output$ID, levels = order.samples)
	DNAcopy.object.cast <- reshape2::dcast(data = DNAcopy.object$output, formula = ID + chrom + loc.start + loc.end + num.mark + seg.mean ~ ID
		, value.var = "num.mark")
	DNAcopy.object.cast[is.na(DNAcopy.object.cast)] <- 0
	
	## Calculate number of samples
	sample.number <- ncol(DNAcopy.object.cast) - 6
	
	## Collapse data
	DNAcopy.object.cast.aggregate <- aggregate(DNAcopy.object.cast[, c("chrom", "num.mark")], by = list(DNAcopy.object.cast$chrom), FUN = sum)
	
	## Calculate scaling factors
	range.factor <- range.CNA[2] - range.CNA[1]
	
	## Create overviewPlot
	par(mfrow = c(1, 2 + sample.number), mar = c(7, 0, 1, 0))
	barplot(matrix(DNAcopy.object.cast.aggregate[, "num.mark"], dimnames = list(NULL, "Chr")), col = c("black", "white"), border = NA
		, yaxt = "n", xlim = c(0, 0.1), width = 0.1, las = 2)
	for(i in 1:sample.number) {
		barplot(as.matrix(DNAcopy.object.cast[, 6 + i, drop = FALSE])
			, col = color.palette[1 + (DNAcopy.object.cast$seg.mean - range.CNA[1]) / range.factor * 48], border = NA
			, yaxt = "n", xlim = c(0, 0.1), width = 0.1, las = 2)
	}
	barplot(matrix(rep(1, length(color.palette)), ncol = 1, dimnames = list(NULL, paste(range.CNA, collapse = " to ")))
		, col = color.palette, border = NA, yaxt = "n", beside = FALSE, xlim = c(0, 0.1), width = 0.1, las = 2)
	
}
