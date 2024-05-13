# Load and install missing packages ---------------------------------------
packages <- c("ggplot2", "httr", "jsonlite", "RCurl")
package.check <- lapply(
    packages,
    FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    }
)

# Bubble plot of gene contributions across cell types - Figure 5a -------------

MCCClust <- read.delim("COVID/results_COVIDPBMCStatus/Status_clusterResults.tsv", row.names = 1)[,"MCC",drop=F]
cellTypes <- rownames(MCCClust)[order(MCCClust$MCC, decreasing = T)]
MCCClust <- MCCClust[cellTypes,,drop=F]

phenoCOVID <- read.delim("COVID/data/Phenodata.tsv")
cellTypesTop <- rownames(MCCClust[MCCClust[,1] > 0.1,,drop=F])
cellTypesTop <- sort(cellTypesTop)
clusterContributions <- list()
for (cluster in cellTypesTop) {
    clustContrib <- read.delim(paste0("COVID/results_COVIDPBMCStatus/gene_contributions/geneContributions_cluster_", cluster, ".tsv"),
                               row.names = 1, check.names = F)
    clusterContributions[[cluster]] <- clustContrib
}

genesCOVID <- lapply(clusterContributions, rowMeans)
genesCOVID <- do.call(cbind, genesCOVID)
rownames(genesCOVID) <- rownames(clustContrib)
genesCOVID <- genesCOVID[,cellTypesTop]

genesCOVIDSign <- genesCOVID
genesCOVIDSign[genesCOVIDSign > 0] <- "Positive"
genesCOVIDSign[genesCOVIDSign < 0] <- "Negative"

genesCOVIDRank <- apply(-genesCOVID, 2, rank) # Negative to rank inversely
rownames(genesCOVIDRank) = rownames(genesCOVID)

selectedGenes <- c()

for (cluster in cellTypesTop) {
    topGenes <- rownames(genesCOVID)[order(rank(-abs(genesCOVID[,cluster])), decreasing = F)[1:3]]
    selectedGenes <- unique(c(selectedGenes, topGenes))
}

genesCOVIDRank[genesCOVIDRank > (nrow(genesCOVIDRank)/2)] <- genesCOVIDRank[genesCOVIDRank > (nrow(genesCOVIDRank)/2)] - nrow(genesCOVIDRank) - 1
genesCOVIDRank <- abs(genesCOVIDRank)
genesCOVIDRank[genesCOVIDRank > 100] <- NA

n = length(selectedGenes)
m <- length(cellTypesTop)

# Convert the matrices to data frames, and add a column for the feature names
df1 <- data.frame(selectedGenes=rep(selectedGenes, m), cellTypesTop=rep(cellTypesTop, each=n),
                  Sign=as.vector(genesCOVIDSign[selectedGenes,]), Rank=-as.vector(as.matrix(genesCOVIDRank[selectedGenes,])))

feature_order <- factor(selectedGenes, levels = selectedGenes)
clusters_order <- factor(cellTypesTop, levels = cellTypesTop)

# Now we can use ggplot2 to create the bubble plot
p <- ggplot(df1, aes(x=cellTypesTop, y=selectedGenes, size=Rank, color=Sign)) +
    geom_point(alpha=0.8) +
    scale_color_manual(values = c("#006837", "#d73027")) +
    scale_x_discrete(limits = clusters_order) +
    scale_y_discrete(limits = rev(feature_order)) +
    theme_bw() +
    annotate("rect", xmin=0,xmax=5.5,ymin=-Inf,ymax=Inf, 
             fill="#8DD3C7", alpha=0.2, color=NA) +
    annotate("rect",xmin=5.5,xmax=6.5,ymin=-Inf,ymax=Inf,
             fill="#FFFFB3", alpha=0.2, color=NA) +
    annotate("rect",xmin=6.5,xmax=7.5,ymin=-Inf,ymax=Inf,
             fill="#BEBADA", alpha=0.2, color=NA) +
    annotate("rect",xmin=7.5,xmax=12.5,ymin=-Inf,ymax=Inf,
             fill="#FB8072", alpha=0.2, color=NA) +
    annotate("rect",xmin=12.5,xmax=13.5,ymin=-Inf,ymax=Inf,
             fill="#80B1D3", alpha=0.2, color=NA) +
    annotate("rect",xmin=13.5,xmax=24.5,ymin=-Inf,ymax=Inf,
             fill="#FDB462", alpha=0.2, color=NA) +
    annotate("rect",xmin=24.5,xmax=34.5,ymin=-Inf,ymax=Inf,
             fill="#B3DE69", alpha=0.2, color=NA) +
    annotate("rect",xmin=34.5,xmax=Inf,ymin=-Inf,ymax=Inf,
             fill="#FCCDE5", alpha=0.2, color=NA) +
    theme(text = element_text(size = 30), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank(), axis.title.y = element_blank(),
          legend.position="top")

ggsave("figures/figure5a.pdf", p, width=1200, height=1200, units="px", scale=5)


# Enrichment analysis - Figure 5b -----------------------------------------

genesEnrichment <- c()

for (cluster in cellTypesTop) {
    topGenes <- rownames(genesCOVID)[order(rank(-abs(genesCOVID[,cluster])), decreasing = F)[1:10]]
    genesEnrichment <- unique(c(genesEnrichment, topGenes))
}


source("./COVID/GC4libR.R")
enrichmentResults <- launchAnalysis(organism = "Homo sapiens",
                                    inputType = "genes",
                                    inputQuery = genesEnrichment,
                                    annotationsDBs = c("WikiPathways"),
                                    inputCoannotation = "no",
                                    inputName1 = "enrichmentAnalysis",
                                    inputEmail = "")

enrichmentResults <- enrichmentResults$stats_tables$`enrichmentAnalysis-WikiPathways`
enrichmentResults$pval_adj <- as.numeric(enrichmentResults$pval_adj)
enrichmentResults$pval <- as.numeric(enrichmentResults$pval)
enrichmentResults$relative_enrichment <- as.numeric(enrichmentResults$relative_enrichment)
enrichmentResults <- enrichmentResults[enrichmentResults$pval_adj < 0.05,]
enrichmentResults$logp <- -log10(enrichmentResults$pval)
enrichmentResults <- enrichmentResults[order(enrichmentResults$logp, decreasing = T),]


# Reduce long descriptions
enrichmentResults$description <- gsub(" pathway", "", enrichmentResults$description)
enrichmentResults$description <- sapply(enrichmentResults$description, function(x) {
    if (nchar(x) > 45) {return(paste0(substr(x, start=1, stop=45), "..."))} else{return(x)}
})

feature_order <- factor(enrichmentResults$description, levels = enrichmentResults$description)


p <- ggplot(enrichmentResults, aes(x=logp, y=description, color=relative_enrichment, size = logp)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = "Relative\nEnrichment",
                           guide=guide_colorbar(reverse=TRUE)) +
    scale_y_discrete(limits = rev(feature_order)) +
    ylab(NULL) + 
    scale_size(range=c(3, 15)) +
    theme(text = element_text(size = 30),
          legend.position="top")

ggsave("figures/figure5b.pdf", p, width=1000, height=1200, units="px", scale=5)