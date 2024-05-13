# Dependencies ------------------------------------------------------------

packagesCRAN <- c("remotes", "hdf5r", "glue", "R.utils", "Seurat", "parallel", 
                  "tidyverse", "BiocManager")

package.check <- lapply(
    packages,
    FUN = function(x) {
        if (!require(x, character.only = TRUE)) {
            install.packages(x, dependencies = TRUE)
            library(x, character.only = TRUE)
        }
    }
)

packagesGithub <- c("SeuratDisk", "SeuratData")
if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
    remotes::install_github("mojaveazure/seurat-disk", dependencies = TRUE)
}

if (!requireNamespace("SeuratData", quietly = TRUE)) {
    remotes::install_github("satijalab/seurat-data")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
    BiocManager::install("GEOquery")
}


# Download and process data -----------------------------------------------

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

dir.create("SLE/data/pediatrics/data_mtx")
dir.create("SLE/data/pediatrics/data")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135779/suppl/GSE135779_genes.tsv.gz",
              "SLE/data/pediatrics/data_mtx/features.tsv.gz")

untar("SLE/data/pediatrics/GSE135779_RAW.tar",exdir = "SLE/data/pediatrics/data_mtx")

patterns = list.files("SLE/data/pediatrics/data_mtx",pattern = "mtx.gz")
gsms = sapply(strsplit(patterns,"_"), `[`, 1)

st = getGEO("GSE135779")
phen = data.frame(phenoData(st$GSE135779_series_matrix.txt.gz)@data)

clin_info = read.csv("SLE/data/pediatrics/clinical_info.csv")
clin_info$Condition = ifelse(clin_info$Groups == "cHD" | clin_info$Groups == "aHD", "Healthy", "SLE")

cc = strsplit(phen$title," \\[JB")
phen$title <- unlist(cc)[2*(1:length(phen$title))-1]
phen = merge(phen,clin_info,by.x = "title",by.y = "Names")

phen <- phen[,c("title","geo_accession","Condition","Batch")]
gsms = phen$geo_accession
for (gsm in gsms){
  dir.create(glue("SLE/data/pediatrics/data/{gsm}"))
  files = list.files("SLE/data/pediatrics/data_mtx",gsm)
  file.copy(file.path("SLE/data/pediatrics/data_mtx", files), glue("SLE/data/pediatrics/data/{gsm}"), overwrite = TRUE)
  file.copy("SLE/data/pediatrics/data_mtx/features.tsv.gz",glue("SLE/data/pediatrics/data/{gsm}"))
  files = list.files(glue("SLE/data/pediatrics/data/{gsm}"))
  barfile = files[grep("barcodes",files)]
  file.rename(glue("SLE/data/pediatrics/data/{gsm}/{barfile}"),glue("SLE/data/pediatrics/data/{gsm}/barcodes.tsv.gz"))
  mtxfile = files[grep("mtx",files)]
  file.rename(glue("SLE/data/pediatrics/data/{gsm}/{mtxfile}"),glue("SLE/data/pediatrics/data/{gsm}/matrix.mtx.gz"))
}

write.table(phen,"SLE/data/pediatrics/clin_data.tsv",row.names =F, quote = F,sep = "\t")

rownames(phen) <- phen$geo_accession
datasets = phen$geo_accession
seurat_objects <- mclapply(datasets,function(dataset){
  path <- glue("SLE/data/pediatrics/data/{dataset}")
  seu.data <- Read10X(data.dir = path)
  seu <- CreateSeuratObject(counts = seu.data,
                            project = dataset)
  metadata <- phen[phen$geo_accession == dataset,]
  for(metadata_column in colnames(metadata)){
    metadata_value <- metadata[,c(metadata_column)]
    seu <- AddMetaData(seu,metadata_value,col.name = metadata_column)
  }
  return(seu)
}, mc.cores = 10)

names(seurat_objects) <- datasets
first_object <- seurat_objects[[1]]
others <- c()

for (i in 2:length(seurat_objects)){
  others <- c(others,seurat_objects[[i]])
}

alldata <- merge(first_object,others,add.cell.ids = datasets)

# Save object
SaveH5Seurat(alldata, filename = "SLE/data/pediatrics/pediatrics_raw.h5Seurat", overwrite = T)
Convert("SLE/data/pediatrics/pediatrics_raw.h5Seurat", dest = "h5ad", overwrite = T)

