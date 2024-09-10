
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# Data reading
input_file = "Dementia/Dementia_processed.h5ad"
adata = sc.read_h5ad(input_file)


# Add APOE contribution
adata.obs['APOEContr'] = 0.0
adata.obs['APOC1ContrDementia'] = 0.0
adata.obs['SPP1ContrDementia'] = 0.0
adata.obs['MCC'] = 0.0

cellTypes = adata.obs['cell_type'].unique()

for cellType in cellTypes:
    cellTypePath = cellType.replace("/", "_")
    PerformanceTable =  pd.read_table('Dementia/results_Dementia/Status_clusterResults.tsv', index_col=0)
    MCCCellType = PerformanceTable.loc[cellTypePath, 'MCC']
    clusterCells = (adata.obs['cell_type'] == cellType)
    adata.obs.loc[clusterCells, 'MCC'] = MCCCellType
    contributions = pd.read_table('Dementia/results_Dementia/gene_contributions_rank_normal/' + cellTypePath + '.tsv', index_col=0)
    for sample in contributions.columns:
        APOESample = contributions.loc["APOE"][sample]
        clusterCellsSample = (adata.obs['cell_type'] == cellType) & (adata.obs['donor_id'] == sample)
        adata.obs.loc[clusterCellsSample, 'APOEContr'] = APOESample

    contributions = pd.read_table('Dementia/results_Dementia/gene_contributions_rank_dementia/' + cellTypePath + '.tsv', index_col=0)
    for sample in contributions.columns:
        APOC1Sample = contributions.loc["APOC1"][sample]
        SPP1Sample = contributions.loc["SPP1"][sample]
        clusterCellsSample = (adata.obs['cell_type'] == cellType) & (adata.obs['donor_id'] == sample)
        adata.obs.loc[clusterCellsSample, 'APOC1ContrDementia'] = APOC1Sample
        adata.obs.loc[clusterCellsSample, 'SPP1ContrDementia'] = SPP1Sample


sc.set_figure_params(dpi_save=300, figsize=(6,6), format="tiff")
sc.pl.scatter(adata, basis="umap", color='cell_type', save="figures/figure4a",  palette='tab20', show=False, title = "Cell types", size=0.4)
sc.pl.umap(adata, color='MCC', color_map="viridis_r", save="figures/figure4b",  show=False, title='MCC', size=0.4)
sc.pl.umap(adata, color='APOEContr', color_map="plasma_r", save="figures/figure4c",  show=False, title='APOE Healthy', size=0.4, na_color="lightgray", colorbar_loc=None)
sc.pl.umap(adata, color='APOC1ContrDementia', color_map="plasma_r", save="figures/figure4d",  show=False, title='APOC1 Dementia', size=0.4, na_color="lightgray", colorbar_loc=None)
sc.pl.umap(adata, color='SPP1ContrDementia', color_map="plasma_r", save="figures/figure4e",  show=False, title='SPP1 Dementia', size=0.4, na_color="lightgray", colorbar_loc=None)
