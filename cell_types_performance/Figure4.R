# Load and install missing packages ---------------------------------------

packages <- c("ggplot2", "matrixStats", "reshape2")
package.check <- lapply(
    packages,
    FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    }
)


# Correlation performance, # Cells and cell types ratio - Figure 4a -------

MCCClust <- read.delim("SLE/results_SLE/Status_clusterResults.tsv", row.names = 1)[,"MCC",drop=F]
cellTypes <- rownames(MCCClust)[order(MCCClust$MCC, decreasing = T)]
MCCClust <- MCCClust[cellTypes,,drop=F]
NCells <- c()
ratioCells <- c()
for (cellType in cellTypes) {
    metaCell <- read.delim(paste0("SLE/data/SLE_Science/Metadata_", cellType, ".tsv"), row.names = 1)
    NCells <- c(NCells, nrow(metaCell))
    counts <- table(metaCell$Status)
    ratioCells <- c(ratioCells, counts[1]/counts[2])
}

cor.test(ratioCells, MCCClust[,1], method = "spearman")
cor.test(NCells, MCCClust[,1], method = "spearman")

df <- data.frame(x = NCells, y = MCCClust[,1], N = ratioCells, name = cellTypes)

p <- ggplot(df, aes(x = x, y = y, fill = N)) +
    geom_point(size = 3, color="black", shape = 21) +
    geom_text(aes(x = x, y = y, label = name), hjust = 1.5, vjust = 0.5, size=3, color="black") +
    xlab("# Cells") +
    ylab("MCC") +
    ylim(min(df$y), 0.8)+
    theme_classic() +
    theme(legend.position = c(0, 1), 
          legend.justification = c("left", "top"),
          legend.key.size = unit(0.05, "npc"),
          text=element_text(size=16)) +
    scale_fill_gradient(low = "yellow", high = "red") +
    labs(fill = "Healthy/SLE ratio")

ggsave("figures/figure4a.pdf", p, scale=2.5, width = 2.05, height = 1.5)



# In silico analysis - Figure 4b ------------------------------------------

iterationsList <- as.character(seq(10))
simulationResults <- list()
for (iteration in iterationsList) {
    simulationResults[[iteration]] <- read.delim(paste0("cell_types_performance/results_simulation/results_simulation_", iteration, 
                                                      "/Condition_clusterResults.tsv"), row.names = 1)[,"MCC", drop = FALSE]
}

simulationResults <- do.call(cbind, simulationResults)
colnames(simulationResults) <- iterationsList
simulationResults <- simulationResults[c(1, 3:10, 2),]


meansimulation = rowMedians(as.matrix(simulationResults))
SDsimulation = rowSds(as.matrix(simulationResults))
rangeIQR <- colIQRs(as.matrix(simulationResults))
upperLim <- meansimulation + rangeIQR
lowerLim <- meansimulation - rangeIQR
datLineplot <- cbind(meansimulation, upperLim, lowerLim)
datLineplot <- data.frame(datLineplot)
datLineplot$ind <- factor(rev(seq(10, 100, 10)), levels=rev(seq(10, 100, 10)))

p <- ggplot(datLineplot) + aes(x=ind, y=meansimulation, ymin=lowerLim, ymax=upperLim) +
    geom_point(size=2) +
    geom_line(group=1, lty=2) +
    geom_ribbon(group=1,alpha=0.15, fill = "#E15759") +
    labs(x="Differential expression (%)", y = "MCC") +
    theme_classic() +
    theme(legend.title = element_blank(), text=element_text(size=16))

ggsave("figures/figure4b.pdf", p, scale=2.5, width = 2.05, height = 1.5)


# Ablation analysis - Figure 4c -------------------------------------------

NCellsList <- paste(c("30000", "25000", "20000", "15000", "10000", "5000"), rep(seq(10), each=6), sep = "_")
NCellsResults <- list()
for (NCells in NCellsList) {
    NCellsResults[[NCells]] <- read.delim(paste0("cell_types_performance/ablation/results_ablat_N", 
                                                 NCells, "/Status_clusterResults.tsv"), 
                                          row.names = 1)[,"MCC", drop = FALSE]
    
}

NCellsResults <- do.call(cbind, NCellsResults)
colnames(NCellsResults) <- NCellsList

NCellsResults <- NCellsResults[order(NCellsResults[,1], decreasing = T),]

mediansCellTypes <- list()
for (cellType in rownames(NCellsResults)) {
    for (NCells in c("30000", "25000", "20000", "15000", "10000", "5000")) {
        MCCs <- NCellsResults[cellType, grep(paste0("^",NCells, "_"), colnames(NCellsResults))]
        mediansCellTypes[[cellType]] <- c(mediansCellTypes[[cellType]], median(unlist(MCCs[1,,drop=T])))
    }
}

mediansCellTypes <- do.call(cbind, mediansCellTypes)
rownames(mediansCellTypes) <- c("30000", "25000", "20000", "15000", "10000", "5000")
melted_NCellsResults <- melt(mediansCellTypes)
melted_NCellsResults$Var1 <- factor(melted_NCellsResults$Var1, levels = c("30000", "25000", "20000", "15000", "10000", "5000"))

p <- ggplot(melted_NCellsResults, aes(x = Var1, y = value, group = Var2, color = Var2)) +
    geom_line() +
    geom_point() +
    labs(x = "# Cells",
         y = "MCC",
         color = "Cell type") +
    theme_classic() +
    guides(color="none")+
    theme(text=element_text(size=16))

ggsave("figures/figure4c.pdf", p, scale=2.5, width = 2.05, height = 1.5)
