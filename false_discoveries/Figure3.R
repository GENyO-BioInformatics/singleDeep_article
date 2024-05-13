# Load and install missing packages ---------------------------------------

packages <- c("reshape2", "ggplot2", "Hmisc")
package.check <- lapply(
    packages,
    FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    }
)


# Prepare data ------------------------------------------------------------

pseudobulkSLE <- read.delim("SLE/data/Science/pseudobulk_Science_raw_sum.tsv", check.names = F, row.names = 1)
MCCClust <- read.delim("SLE/results_SLE/Status_clusterResults.tsv", row.names = 1)[,"MCC",drop=F]
cellTypes <- rownames(MCCClust)[order(MCCClust$MCC, decreasing = T)]

topExprGenes <- list()
for (cellType in cellTypes) {
    samplesCellType <- grep(paste0("_", cellType), rownames(pseudobulkSLE), value = T)
    exprCellType <- colSums(pseudobulkSLE[samplesCellType,])
    cutValues = cut2(exprCellType, g=20, onlycuts = T)
    topGenes <- names(exprCellType)[exprCellType >= cutValues[20]]
    topExprGenes[[cellType]] <- topGenes
}

# Get contributions by strategy and cell type
clusterContributionsStrategies <- list()
strategies <- c("Local + scale", "Local", "Global + scale", "Global")
folderResults <- list("Local + scale" = "SLE/results_SLE",
                      "Local" = "false_discoveries/results_NoScaled_local", 
                      "Global + scale" = "false_discoveries/results_Scaled_global",
                      "Global" = "false_discoveries/results_NoScaled_global")

for (strategy in strategies) {
    clusterContributions <- list()
    for (cluster in cellTypes) {
        
        clustContrib <- read.delim(paste0(folderResults[[strategy]], "/gene_contributions/geneContributions_cluster_", cluster, ".tsv"),
                                   row.names = 1, check.names = F)
        clusterContributions[[cluster]] <- clustContrib
    }
    clusterContributionsStrategies[[strategy]] <- clusterContributions
}


# Figure 3a ---------------------------------------------------------------

# Mean contribution of top genes
meanContribTopGenesStrategies <- list()
for (strategy in strategies) {
    for (cellType in cellTypes) {
        contributionsCellType <- sort(abs(rowMeans(clusterContributionsStrategies[[strategy]][[cellType]])), decreasing = T)[1:100]
        meanContribTopGenesStrategies[[strategy]] <- c(meanContribTopGenesStrategies[[strategy]],
                                                       sum(names(contributionsCellType) %in% topExprGenes[[cellType]]))
    }
}

datFolds <- stack(as.data.frame(meanContribTopGenesStrategies, check.names=F))

dat <- stack(as.data.frame(t(sapply(meanContribTopGenesStrategies, mean))))

# Calculate mean and sd for each cell type
for (strategy in strategies) {
    values <- datFolds[datFolds$ind == strategy, "values"]
    sdValues <- sd(values)
    dat[dat$ind == strategy, "sd"] <- sdValues
}

p <- ggplot(dat, aes(x=ind, y=values, fill=values)) +
    geom_bar( stat = "identity") +
    geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_continuous(low="#0AA363", high="#F05A3E") +
    xlab("Method") +
    ylab("Top expressed genes (%)") +
    theme_classic() +
    theme(legend.position = "none", axis.text.x=element_text(angle=45, vjust=0.95, hjust = 0.95))

ggsave("figures/Figure3a.pdf", p, scale=1.2, width = 1.75, height = 2.6)


# Figure 3b ---------------------------------------------------------------

# Percentage of top 5% expressed genes into top 100 contributing genes
topGenesContributionsStrategies <- list()
for (strategy in strategies) {
    for (cellType in cellTypes) {
        contributionsCellType <- abs(rowMeans(clusterContributionsStrategies[[strategy]][[cellType]]))
        topGenes <- topExprGenes[[cellType]][topExprGenes[[cellType]] %in% names(contributionsCellType)]
        topGenesContributionsStrategies[[strategy]] <- c(topGenesContributionsStrategies[[strategy]],
                                                         mean(contributionsCellType[topGenes]))
    }
}


datFolds <- stack(as.data.frame(topGenesContributionsStrategies, check.names=F))

dat <- stack(as.data.frame(t(sapply(topGenesContributionsStrategies, mean))))

# Calculate mean and sd for each cell type
for (strategy in strategies) {
    values <- datFolds[datFolds$ind == strategy, "values"]
    sdValues <- sd(values)
    dat[dat$ind == strategy, "sd"] <- sdValues
}

p <- ggplot(dat, aes(x=ind, y=values, fill=values)) +
    geom_bar( stat = "identity") +
    geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_continuous(low="#0AA363", high="#F05A3E") +
    xlab("Method") +
    ylab("Contribution") +
    theme_classic() +
    theme(legend.position = "none", axis.text.x=element_text(angle=45, vjust=0.95, hjust = 0.95))

ggsave("figures/Figure3b.pdf", p, scale=1.2, width = 1.75, height = 2.6)



# Figure 3c ---------------------------------------------------------------

ISG_genes <- c("IFI44L", "XAF1", "MX1", "OAS1")
ranksISG <- matrix(NA, nrow=4, ncol=4, dimnames = list(ISG_genes, strategies))

for (gene in ISG_genes) {
    meanRankGene <- list()
    for (strategy in strategies) {
        for (cellType in cellTypes) {
            contributionsRank <- rank(-rowMeans(clusterContributionsStrategies[[strategy]][[cellType]]))
            meanRankGene[[strategy]] <- c(meanRankGene[[strategy]], contributionsRank[gene])
        }
        ranksISG[gene, strategy] <- median(meanRankGene[[strategy]])
    }
}

dat <- stack(as.data.frame(t(ranksISG), check.names=F))
dat$Method <- factor(rep(strategies, 4), levels=strategies)
dat$rank <- dat$values
dat$values <- 1/dat$values

p <- ggplot(dat, aes(x=ind, y=values, fill=Method)) +
    geom_bar(position = "dodge", stat = "identity", color = "black") +
    geom_text(aes(label = rank), position = position_dodge(width = 0.9), vjust = -0.5, size=2) +
    xlab("Gene") +
    ylab("Inverse of median contributing position") +
    theme_classic() +
    scale_fill_jco() +
    labs(fill="Method") +
    theme(legend.position = "top", legend.key.size = unit(0.02, "npc"),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))

ggsave("figures/Figure3c.pdf", p, scale=1.2, width = 3.5, height = 2.4)
