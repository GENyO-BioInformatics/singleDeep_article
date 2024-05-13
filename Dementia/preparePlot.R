# singleDeep performance by cell type - Figure 5a -------------------------

MCCClust <- read.delim("Dementia/results_Dementia/Status_clusterResults.tsv", row.names = 1)[,"MCC",drop=F]
cellTypes <- rownames(MCCClust)[order(MCCClust$MCC, decreasing = T)]

phenoDementiaMTG <- read.delim("Dementia/dataMTG/Phenodata.tsv")
clusterContributions <- list()
clusterContributionsRank <- list()
dir.create("Dementia/results_Dementia/gene_contributions_rank_normal")
dir.create("Dementia/results_Dementia/gene_contributions_rank__dementia")

for (cluster in cellTypes) {
    clustContrib <- read.delim(paste0("Dementia/results_Dementia/gene_contributions/geneContributions_cluster_", cluster, ".tsv"),
                               row.names = 1, check.names = F)
    clusterContributions[[cluster]] <- clustContrib
    clusterContributionsRank[[cluster]] <- -apply(clustContrib, 2, rank)
    clusterContributionsRank[[cluster]][clusterContributionsRank[[cluster]] < -100] <- NA
    write.table(clusterContributionsRank[[cluster]], file = paste0("Dementia/results_Dementia/gene_contributions_rank_normal/", cluster, ".tsv"),
                sep = "\t")
    clusterContributionsRank[[cluster]] <- -apply(-clustContrib, 2, rank)
    clusterContributionsRank[[cluster]][clusterContributionsRank[[cluster]] < -100] <- NA
    write.table(clusterContributionsRank[[cluster]], file = paste0("Dementia/results_Dementia/gene_contributions_rank_dementia/", cluster, ".tsv"),
                sep = "\t")
}