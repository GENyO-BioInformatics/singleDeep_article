# Load and install missing packages
packages = c("optparse")
if (!requireNamespace("optparse", quietly = TRUE)) {
    install.packages("optparse")
}

# Parameters definiton
option_list <- list(
    make_option(c("--inputPath"), type="character", default = NULL,
                help="Data folder"),
    make_option(c("--sampleColumn"), type="character", default= NULL,
                help="Metadata column with samples labels."))

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
inputPath <- opt$inputPath
sampleColumn <- opt$sampleColumn

# Pseudobulk
phenodata <- read.delim(paste0(inputPath, "/Phenodata.tsv"))
samples <- phenodata$Sample
cellFiles <- list.files(inputPath, pattern = "Metadata")
cellFiles <- strsplit(cellFiles, "Metadata_|.tsv")
cellTypes <- sapply(cellFiles, function(x)  {return(x[[2]])})

pseudobulk <- list()
for (cluster in cellTypes) {
    clusterTable <- read.delim(paste0(inputPath, "/", cluster, ".tsv"))
    clusterMeta <- read.delim(paste0(inputPath, "/Metadata_", cluster, ".tsv"))
    samplesCluster <- unique(clusterMeta[,sampleColumn])
    pseudobulkCluster <- matrix(NA, nrow = nrow(clusterTable), ncol = length(samplesCluster),
                                dimnames = list(rownames(clusterTable), samplesCluster))
    for (sample in samplesCluster) {
        cellsSample <- rownames(clusterMeta)[clusterMeta[,sampleColumn] == sample]
        cellsSample <- gsub("-", ".", cellsSample, fixed = T) # To match with data frame names
        exprSample <- clusterTable[,cellsSample, drop=F]
        pseudobulkCluster[,sample] <- rowMeans(exprSample)
    }
    pseudobulk[[cluster]] <- pseudobulkCluster
}

pseudobulk_whole <- matrix(0, nrow = nrow(pseudobulk[[1]]), ncol <- length(samples),
                           dimnames = list(rownames(pseudobulk[[1]]), samples))
nCells <- list()

for (cluster in cellTypes) {
    clusterTable <- read.delim(paste0(inputPath, "/", cluster, ".tsv"))
    clusterMeta <- read.delim(paste0(inputPath, "/Metadata_", cluster, ".tsv"))
    samplesCluster <- unique(clusterMeta[,sampleColumn])
    for (sample in samplesCluster) {
        cellsSample <- rownames(clusterMeta)[clusterMeta[,sampleColumn] == sample]
        cellsSample <- gsub("-", ".", cellsSample, fixed = T) # To match with data frame names
        exprSample <- clusterTable[,cellsSample, drop=F]
        pseudobulk_whole[,sample] <- pseudobulk_whole[,sample] + rowSums(exprSample)
        if (!sample %in% names(nCells)) {
            nCells[[sample]] <- length(cellsSample)
        }
        else {
            nCells[[sample]] <- nCells[[sample]] + length(cellsSample)
        }
    }
}

pseudobulk_whole2 <- sapply(colnames(pseudobulk_whole), function(x) {return(pseudobulk_whole[,x] / nCells[[x]])})

write.table(pseudobulk_whole2, 
            file = paste0(inputPath, "/pseudobulk_whole.tsv"),
            sep="\t")
