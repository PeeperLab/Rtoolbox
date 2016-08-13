## Script by Thomas Kuilman
## path argument: path to output folder of analysis (e.g. PATH/my_analysis.GseaPreranked.1470948568349)
## gene.set argument: name of the gene set (e.g. V$AP1_Q2)
## class.name: the name of the class / variable to which genes have been correlated (e.g. drug-treatment)

replotGSEA <- function(path, gene.set, class.name) {
	
	if(missing(path)) {
		stop("Path argument is required")
	}
	if (!file.exists(path)) {
		stop("The path folder could not be found. Please change the path")
	}
	if(missing(gene.set)) {
		stop("Gene set argument is required")
	}

	## Load .rnk data
	path.rnk <- list.files(path = file.path(path, "edb"),
	                       pattern = ".rnk$", full.names = TRUE)
	gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
	colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
	
	## Load .edb data
	path.edb <- list.files(path = file.path(path, "edb"),
	                       pattern = ".edb$", full.names = TRUE)
	gsea.edb <- read.delim(file = path.edb,
	                       header = FALSE, stringsAsFactors = FALSE)
	gsea.edb <- unlist(gsea.edb)
	gsea.metric <- gsea.edb[grep("METRIC=", gsea.edb)]
	gsea.metric <- unlist(strsplit(gsea.metric, " "))
	gsea.metric <- gsea.metric[grep("METRIC=", gsea.metric)]
	gsea.metric <- gsub("METRIC=", "", gsea.metric)
	gsea.edb <- gsea.edb[grep("<DTG", gsea.edb)]
	
	# Select the right gene set
	gsea.edb <- gsea.edb[grep(gene.set, gsea.edb, fixed = TRUE)]
	if(length(gsea.edb) == 0) {
		stop("The gene set name was not found, please provide a correct name")
	}
	
	# Get template name
	gsea.edb <- gsub(".*TEMPLATE=(.*)", "\\1", gsea.edb)
	gsea.edb <- unlist(strsplit(gsea.edb, " "))
	gsea.template <- gsea.edb[1]
	
	# Get gene set name
	gsea.gene.set <- gsea.edb[2]
	gsea.gene.set <- gsub("GENESET=gene_sets.gmt#", "", gsea.gene.set)
	
	# Get enrichment score
	gsea.enrichment.score <- gsea.edb[3]
	gsea.enrichment.score <- gsub("ES=", "", gsea.enrichment.score)
	
	# Get gene set name
	gsea.normalized.enrichment.score <- gsea.edb[4]
	gsea.normalized.enrichment.score <- gsub("NES=", "",
	                                         gsea.normalized.enrichment.score)
	
	# Get nominal p-value
	gsea.p.value <- gsea.edb[5]
	gsea.p.value <- gsub("NP=", "", gsea.p.value)
	gsea.p.value <- as.numeric(gsea.p.value)
	
	# Get FDR
	gsea.fdr <- gsea.edb[6]
	gsea.fdr <- gsub("FDR=", "", gsea.fdr)
	gsea.fdr <- as.numeric(gsea.fdr)
	
	# Get hit indices
	gsea.edb <- gsea.edb[grep("HIT_INDICES=", gsea.edb):length(gsea.edb)]
	gsea.hit.indices <- gsea.edb[seq_len(grep("ES_PROFILE=", gsea.edb) - 1)]
	gsea.hit.indices <- gsub("HIT_INDICES=", "", gsea.hit.indices)
	gsea.hit.indices <- as.integer(gsea.hit.indices)
	
	# Get ES profile
	gsea.edb <- gsea.edb[grep("ES_PROFILE=", gsea.edb):length(gsea.edb)]
	gsea.es.profile <- gsea.edb[seq_len(grep("RANK_AT_ES=", gsea.edb) - 1)]
	gsea.es.profile <- gsub("ES_PROFILE=", "", gsea.es.profile)
	gsea.es.profile <- as.numeric(gsea.es.profile)
	
	
	## Create GSEA plot
	# Save default for resetting
	def.par <- par(no.readonly = TRUE)
	
	# Create a new device of appropriate size
	dev.new(width = 4, height = 5)
	
	# Create a division of the device
	gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))
	layout.show(gsea.layout)
	
	# Create plots
	par(mar = c(0, 5, 2, 2))
	plot(c(1, gsea.hit.indices), c(0, gsea.es.profile), type = "l", col = "red",
	     lwd = 1.5, xaxt = "n", xaxs = "i", xlab = "",
	     ylab = "Enrichment score (ES)",
	     main = list(gsea.gene.set, font = 1, cex = 1),
	     panel.first = abline(h = seq(round(min(gsea.es.profile), digits = 1),
	                                  max(gsea.es.profile), 0.1),
	                          col = "gray95", lty = 2))
	plot.coordinates <- par("usr")
	if(gsea.enrichment.score < 0) {
		text(length(gsea.rnk$metric) * 0.01, plot.coordinates[3] * 0.98,
		     paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
		           gsea.enrichment.score, "\nNormalized ES:",
		           gsea.normalized.enrichment.score), adj = c(0, 0))
	} else {
		text(length(gsea.rnk$metric) * 0.99, plot.coordinates[4] - ((plot.coordinates[4] - plot.coordinates[3]) * 0.03),
		     paste("Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
		           gsea.enrichment.score, "\nNormalized ES:",
		           gsea.normalized.enrichment.score, "\n"), adj = c(1, 1))
	}
	
	par(mar = c(0, 5, 0, 2))
	plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
	     ylab = "", xlim = c(1, length(gsea.rnk$metric)))
	abline(v = gsea.hit.indices, lwd = 0.75)
	
	par(mar = c(0, 5, 0, 2))
	plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
	     ylab = "", xlim = c(1, length(gsea.rnk$metric)))
	rank.colors <- ceiling((gsea.rnk$metric + 1) * 128)
	rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
	abline(v = seq_len(nrow(gsea.rnk)), lwd = 0.1, col = rank.colors)
	text(length(gsea.rnk$metric) / 2, 0,
	     labels = ifelse(!missing(class.name), class.name, gsea.template))
	text(length(gsea.rnk$metric) * 0.01, 0, "Positive", adj = c(0, NA))
	text(length(gsea.rnk$metric) * 0.99, 0, "Negative", adj = c(1, NA))
	
	par(mar = c(5, 5, 0, 2))
	plot(gsea.rnk$metric, col = "lightgrey", type = "h", lwd = 0.1, xaxs = "i",
	     xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
	     ylim = c(-1, 1), yaxs = "i",
	     ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric),
	     panel.first = abline(h = seq(-0.5, 0.5, 0.5), col = "gray95", lty = 2))
	
	# Reset to default
	par(def.par)

}