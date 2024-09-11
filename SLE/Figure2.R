
# Load and install missing packages ---------------------------------------

packages <- c("metrica", "ggplot2", "ggsci", "Hmisc", "pheatmap", "paletteer")
package.check <- lapply(
    packages,
    FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    }
)


# Internal validation - Figure 2a -----------------------------------------

resultsTable <- read.delim("SLE/results_SLE/Status_testResults.tsv", row.names = 1)
resultsInternal <- list()
for (MLModel in c("LR", "SVM", "RF", "KNN", "NB", "LDA", "FNN", "DT")) {
    resultsInternal[[MLModel]] <- read.delim(paste0("SLE/results_SLE/",
                                                      MLModel, "_Whole/testResults.tsv"),
                                               row.names = 1)
}

resultsInternal <- do.call(cbind, resultsInternal)
perfMergedInternal <- data.frame(cbind(resultsTable, resultsInternal))
perfMergedInternal["MCC",] <- (unlist(perfMergedInternal["MCC",]) + 1) / 2
colnames(perfMergedInternal) <- c("singleDeep", "LR", "SVM", "RF", "KNN", "NB", "LDA", "FNN", "DT")


dataBarplot <- stack(perfMergedInternal)
dataBarplot$metric <- factor(rep(c("Accuracy", "Precision", "Recall", "F1", "normMCC"), 9))

p <- ggplot(dataBarplot, aes(x=factor(metric, c("normMCC", "Accuracy", "Precision", "Recall", "F1")), y=values, fill=ind)) +
    geom_bar(position = "dodge", stat = "identity", color = "black") +
    xlab("Metric") +
    ylab("Value") +
    theme_classic() +
    scale_fill_jco() +
    ylim(0,1) +
    labs(fill="Method") +
    theme(legend.position = "top", legend.key.size = unit(0.02, "npc"),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))

ggsave("figures/figure2a.pdf", p, scale=1.2, width = 3.5, height = 2.4)


# External validation - Table 1 -----------------------------------------

validation_predictions <- read.delim("SLE/results_SLE/validation_predictions.tsv", row.names = 1)
validation_real <- read.delim("SLE/data/SLE_pediatrics/Phenodata.tsv")
validation_real$StatusInt = ifelse(validation_real$Status == "Healthy", 0, 1)
validation_predictions$real <- validation_real[rownames(validation_predictions), "StatusInt"]

metrics <- c("accuracy", "precision", "recall", "fscore", "mcc")
resultsTable <- metrics_summary(obs = validation_predictions$real,
                                pred = validation_predictions$label_predicted, 
                                type="classification", pos_level = "1",
                                metrics_list = metrics)
rownames(resultsTable) <- resultsTable[,1]
resultsTable <- resultsTable[,-1, drop=F]
resultsTable["normMCC",] <- (resultsTable["mcc",] + 1) / 2

resultsExternal <- list()
for (MLModel in c("LR", "SVM", "RF", "KNN", "NB", "LDA", "FNN", "DT")) {
    validation_predictions <- read.delim(paste0("SLE/results_SLE/", MLModel, "_Whole/pediatrics_prediction.tsv"), row.names = 1)
    validation_predictions$real <- validation_real[rownames(validation_predictions), "StatusInt"]
    resultsModel <- metrics_summary(obs = validation_predictions$real,
                                    pred = validation_predictions$label_predicted, 
                                    type="classification", pos_level = "1",
                                    metrics_list = metrics)
    rownames(resultsModel) <- resultsModel[,1]
    resultsModel <- resultsModel[,-1, drop=F]
    resultsModel["normMCC",] <- (resultsModel["mcc",] + 1) / 2
    resultsExternal[[MLModel]] <- resultsModel
}
resultsExternal[["CloudPred"]] <- read.delim("./Cloudpred/results_validation_Cloudpred.txt", row.names = 1)
resultsExternal[["ProtoCell4P"]] <- read.delim("./ProtoCell4P/results_validation_Protocell4P.txt", row.names = 1)
resultsExternal[["ScRAT"]] <- read.delim("./scRAT/results_validation_scRAT.txt", row.names = 1)
resultsExternal <- do.call(cbind, resultsExternal)
perfMergedExternal <- data.frame(cbind(resultsTable, resultsExternal))
colnames(perfMergedExternal) <- c("singleDeep",  "LR", "SVM", "RF", "KNN", "NB", "LDA", "FNN", "DT", "CloudPred", "ProtoCell4P", "ScRAT")
rownames(perfMergedExternal) <- c("Accuracy", "Precision", "Recall", "F1", "MCC", "normMCC")
perfMergedExternal <- perfMergedExternal[c("MCC", "Accuracy", "Precision", "Recall", "F1"), 
                                         c("singleDeep", "CloudPred", "ProtoCell4P", "ScRAT", "LR", "SVM", "RF", "KNN", "NB", "LDA", "FNN", "DT")]

write.table(round(perfMergedExternal, 2), "SLE/Table1.tsv", sep = "\t", quote = F, col.names = NA)

# singleDeep performance by cell type - Figure 2b -------------------------

MCCClust <- read.delim("SLE/results_SLE/Status_clusterResults.tsv", row.names = 1)[,"MCC",drop=F]
cellTypes <- rownames(MCCClust)[order(MCCClust$MCC, decreasing = T)]
MCCClust <- MCCClust[cellTypes,,drop=F]

perfList <- list()
for (cluster in cellTypes) {
    perfList[[cluster]] <- read.delim(paste0("SLE/results_SLE/folds_performance/cluster_", cluster, ".tsv"), row.names = 1)[,"MCC"]
}

MCCClustSLE <- do.call("cbind", perfList)
rownames(MCCClustSLE) <- 1:(nrow(MCCClustSLE))
datFolds <- stack(as.data.frame(MCCClustSLE))

dat <- stack(as.data.frame(t(MCCClust)))

# Calculate mean and sd for each cell type
for (cluster in cellTypes) {
    values <- datFolds[datFolds$ind == cluster, "values"]
    sdValues <- sd(values)
    dat[dat$ind == cluster, "sd"] <- sdValues
}

p <- ggplot(dat, aes(x=ind, y=values, fill=values)) +
    geom_bar( stat = "identity") +
    geom_errorbar(aes(ymin=values-sd, ymax=values+sd), width=.2,
                  position=position_dodge(.9)) +
    scale_fill_continuous(low="#F05A3E", high="#0AA363") +
    xlab("Cell type") +
    ylab("MCC") +
    theme_classic() +
    theme(legend.position = "none", axis.text.x=element_text(angle=90, vjust=0.5, hjust = 0.95))

ggsave("figures/figure2b.pdf", p, scale=1.2, width = 3.5, height = 2.6)


# Heatmap of gene contributions across cell types - Figure 2c -------------

pseudobulkSLE <- read.delim("SLE/data/Science/pseudobulk_Science_raw_sum.tsv", check.names = F, row.names = 1)

phenoSLE <- read.delim("SLE/data/SLE_Science/Phenodata.tsv")
cellTypesTop <- cellTypes[1:7]
clusterContributions <- list()
for (cluster in cellTypesTop) {
    clustContrib <- read.delim(paste0("SLE/results_SLE/gene_contributions/geneContributions_cluster_", cluster, ".tsv"),
                               row.names = 1, check.names = F)
    clusterContributions[[cluster]] <- clustContrib
}

genesSLE <- lapply(clusterContributions, rowMeans)
genesSLE <- do.call(cbind, genesSLE)
rownames(genesSLE) <- rownames(clustContrib)
genesSLE <- genesSLE[,cellTypesTop]

genesSLERank <- apply(-genesSLE, 2, rank) # Negative to rank inversely
rownames(genesSLERank) = rownames(genesSLE)

selectedGenes <- c()

for (cluster in cellTypesTop) {
    topGenes <- rownames(genesSLE)[order(rank(-abs(genesSLE[,cluster])), decreasing = F)[1:10]]
    selectedGenes <- c(selectedGenes, topGenes)
}

selectedGenes <- names(sort(table(selectedGenes), decreasing = T))
genesSLERank[(genesSLERank > 100 & genesSLERank < (nrow(genesSLERank) - 100))] <- NA

for (cellType in cellTypesTop) {
    samplesCellType <- grep(paste0("_", cellType), rownames(pseudobulkSLE), value = T)
    exprCellType <- colMeans(pseudobulkSLE[samplesCellType,])
    cutValues = cut2(exprCellType, g=20, onlycuts = T)
    topExprGenes <- names(exprCellType)[exprCellType >= cutValues[20]]
    for (gene in topExprGenes) {
        if (gene %in% rownames(genesSLERank)) {
            if(!is.na(genesSLERank[gene, cellType])) {
                if (genesSLERank[gene, cellType] <= 100) {
                    genesSLERank[gene, cellType] <- 350
                }
                else if (genesSLERank[gene, cellType] > 1000){
                    genesSLERank[gene, cellType] <- (nrow(genesSLERank) - 350)
                }
            }
        }
    }
}

pheatmap(-genesSLERank[selectedGenes,], cluster_cols = F, cluster_rows = F,
         color = paletteer_d("beyonce::X39"), border_color = "black",
         angle_col = 90, width = 2.5, height = 4.93, legend = F, fontsize = 7,
         filename = "figures/figure2c.pdf")
