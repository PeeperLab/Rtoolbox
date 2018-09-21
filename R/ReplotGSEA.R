## Script by Thomas Kuilman
## path argument: path to output folder of analysis (e.g. PATH/my_analysis.GseaPreranked.1470948568349)
## gene.set argument: name of the gene set (e.g. V$AP1_Q2).
## It is used in a grep command, so multiple matching is possible.
## Also, R regular expressions can be handled, e.g. "IL2[0-9]$"
## Leading "V$" from gene set names are stripped to allow using the grep command.
## In case of multiple grep matches a warning is given and the first option is plotted.
## class.name: the name of the class / variable to which genes have been correlated (e.g. drug-treatment)
## metric.range: the range of the metric; defaults to [min(DEFINED RANGE), max(DEFINED RANGE)]
## enrichment.score.range: the range of the enrichment score; defaults to [min(ENRICHMENT SCORE), max(ENRICHMENT SCORE)]

replotGSEA <- function(path, gene.set, class.name, metric.range,
                       enrichment.score.range) {
  
  if (missing(path)) {
    stop("Path argument is required")
  }
  if (!file.exists(path)) {
    stop("The path folder could not be found. Please change the path")
  }
  if (missing(gene.set)) {
    stop("Gene set argument is required")
  }

  ## Load .rnk data
  path.rnk <- list.files(path = file.path(path, "edb"),
                         pattern = ".rnk$", full.names = TRUE)
  gsea.rnk <- read.delim(file = path.rnk, header = FALSE)
  colnames(gsea.rnk) <- c("hgnc.symbol", "metric")
  if (missing(metric.range)) {
    metric.range <- c(min(gsea.rnk$metric), max(gsea.rnk$metric))
  }  
  
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
  if (length(gsea.edb) == 0) {
    stop(paste("The gene set name was not found, please provide",
               "a correct name"))
  }
  if (length(grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)) > 1) {
    warning(paste("More than 1 gene set matched the gene.set",
                  "argument; the first match is plotted"))
  }
  gsea.edb <- gsea.edb[grep(paste0(gsub(".\\$(.*$)", "\\1", gene.set), " "), gsea.edb)[1]]
  
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
  
  # Set enrichment score range
  if (missing(enrichment.score.range)) {
    enrichment.score.range <- c(min(gsea.es.profile), max(gsea.es.profile))
  }
  
  
  ## Create GSEA plot
  # Save default for resetting
  def.par <- par(no.readonly = TRUE)
  
  # Create a new device of appropriate size
  dev.new(width = 3, height = 3)
  
  # Create a division of the device
  gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))
  layout.show(gsea.layout)
  
  # Create plots
  par(mar = c(0, 5, 2, 2))
  plot(c(1, gsea.hit.indices, length(gsea.rnk$metric)),
       c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
       xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
       ylim = enrichment.score.range,
       main = list(gsea.gene.set, font = 1, cex = 1),
       panel.first = {
          abline(h = seq(round(enrichment.score.range[1], digits = 1),
                         enrichment.score.range[2], 0.1),
                 col = "gray95", lty = 2)
         abline(h = 0, col = "gray50", lty = 2)
       })
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
  rank.colors <- gsea.rnk$metric - metric.range[1]
  rank.colors <- rank.colors / (metric.range[2] - metric.range[1])
  rank.colors <- ceiling(rank.colors * 255 + 1)
  tryCatch({
    rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
  }, error = function(e) {
    stop("Please use the metric.range argument to provide a metric range that",
         "includes all metric values")
  })
  # Use rle to prevent too many objects
  rank.colors <- rle(rank.colors)
  barplot(matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n", xlim = c(1, length(gsea.rnk$metric)))
  box()
  text(length(gsea.rnk$metric) / 2, 0.7,
       labels = ifelse(!missing(class.name), class.name, gsea.template))
  text(length(gsea.rnk$metric) * 0.01, 0.7, "Positive", adj = c(0, NA))
  text(length(gsea.rnk$metric) * 0.99, 0.7, "Negative", adj = c(1, NA))
  
  par(mar = c(5, 5, 0, 2))
  rank.metric <- rle(round(gsea.rnk$metric, digits = 2))
  plot(gsea.rnk$metric, type = "n", xaxs = "i",
	     xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
	     ylim = metric.range, yaxs = "i",
	     ylab = if(gsea.metric == "None") {"Ranking metric"} else {gsea.metric},
	     panel.first = abline(h = seq(metric.range[1] / 2,
	                                  metric.range[2] - metric.range[1] / 4,
	                                  metric.range[2] / 2), col = "gray95", lty = 2))

  barplot(rank.metric$values, col = "lightgrey", lwd = 0.1, xaxs = "i",
       xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk$metric)),
       ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,
       ylab = ifelse(gsea.metric == "None", "Ranking metric", gsea.metric), space = 0, add = TRUE)
  box()
  
  # Reset to default
  par(def.par)

}