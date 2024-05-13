# Load and install missing packages ---------------------------------------
packages <- c("devtools", "Seurat")
package.check <- lapply(
    packages,
    FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    }
)

if (!require("splatter", character.only = TRUE)) {
    devtools::install_github("Oshlack/splatter")
}
library(splatter)
source("cell_types_performance/custom_splatter_functions.R")


# Simulate data over 10 iterations ----------------------------------------
dir.create("./cell_types_performance/dataSimulation")
for (iteration in seq(10)) {
    fileName <- paste0("simulation_", iteration)
    sample_heterogeneity = 0.1
    nSamples <- 100
    nGroups <- 10

    # Groups ore the cell types
    # Random cell types proportions
    set.seed(iteration)
    prop <- runif(nGroups)
    prop = prop/(sum(prop))

    params <- newSplatPopParams(nGenes = 2000,
                                batchCells=100,
                                batch.size=nSamples,
                                condition.prob=c(0.5, 0.5),
                                cde.prob = 0.025,
                                cde.facLoc = 0.05,
                                group.prob = prop,
                                similarity.scale=1000, 
                                de.prob = 0.8, 
                                de.facLoc = 0.05,
                                seed=iteration)

    # Reduces the differential expression level for each group
    cde.groupSpecific = list(Group1 = 1, Group2 = 0.9, Group3 = 0.8, Group4 = 0.7, Group5 = 0.6, Group6 = 0.5,
                             Group7 = 0.4, Group8 = 0.3, Group9 = 0.2, Group10 = 0.1)

    sim.means <- splatPopSimulateMeans(params,
                                       vcf=mockVCF(n.samples=nSamples),
                                       cde.groupSpecific = cde.groupSpecific,
                                       sample_heterogeneity = sample_heterogeneity)
    sim <- splatPopSimulateSC(sim.means = sim.means$means,
                              params = params, key = sim.means$key,
                              conditions = sim.means$conditions, counts.only = F,
                              method = "groups", sparsify = FALSE, verbose = T)


    # Save Seurat object
    colnames(sim@assays@data@listData$counts) <- paste0(colData(sim)$Sample, "_", colData(sim)$Cell, "_", colData(sim)$Group)
    rownames(colData(sim)) = colnames(sim@assays@data@listData$counts)
    rownames(sim@assays@data@listData$counts) <- gsub("_", "", rownames(sim@assays@data@listData$counts))

    seuratobj = CreateSeuratObject(counts = sim@assays@data@listData$counts,
                                   project = "SCSimulation",
                                   meta.data = as.data.frame(colData(sim))[,c(3,5,6)],
                                   names.field = 1,
                                   names.delim = "NA")
    seuratobj <- NormalizeData(seuratobj, normalization.method = "LogNormalize")
    seuratobj <- ScaleData(seuratobj)
    seuratobj <- FindVariableFeatures(seuratobj, nfeatures = 2000)
    saveRDS(seuratobj, file = paste0("./cell_types_performance/dataSimulation/", fileName, ".rds"))
    
    commandPrepare <- paste0("Rscript singleDeep/PrepareData.R --inputPath cell_types_performance/dataSimulation/", fileName,
                             ".rds --fileType seurat --sampleColumn Sample --clusterColumn Group --clinicalColumns Condition --targetColumn Condition --outPath cell_types_performance/dataSimulation/", 
                             fileName)
    system(commandPrepare)
}


# Run singleDeep ----------------------------------------------------------
for (iteration in seq(10)) {
    fileName <- paste0("simulation_", iteration)
    commandSingleDeep <- paste0("python3.8 singleDeep/singleDeep.py --inPath cell_types_performance/dataSimulation/",
                                fileName, " --sampleColumn Sample --logPath cell_types_performance/log --resultsPath cell_types_performance/results_simulation/results_",
                                fileName, "/ --varColumn Condition --num_epochs 250 --resultsFilenames Condition --KOuter 5 --KInner 4 --lr 0.01")
    system(commandSingleDeep)
}