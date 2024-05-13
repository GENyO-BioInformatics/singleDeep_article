
import scanpy as sc



# Data reading
input_file = "COVID19/data/COVID19_ALL.h5ad"
adata = sc.read_h5ad(input_file)

# Select PBMC samples
adata = adata[adata.obs['Sample type'].isin(['frozen PBMC', 'fresh PBMC']),:]

# Get highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, inplace=True, subset=True)

# Scale
sc.pp.scale(adata)

# Write the output file
adata.write_h5ad("COVID19/data/COVID19_processed.h5ad")
