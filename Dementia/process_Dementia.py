
import scanpy as sc



# Data reading
input_file = "Dementia/SEA_AD_MTG.h5ad"
adata = sc.read_h5ad(input_file)

# Get highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000, inplace=True, subset=True)

# Scale
sc.pp.scale(adata)

# Write the output file
adata.write_h5ad("Dementia/Dementia_processed.h5ad")
