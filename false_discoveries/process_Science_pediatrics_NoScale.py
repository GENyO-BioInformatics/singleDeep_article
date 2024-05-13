################
# Dependencies #
################
import scanpy as sc
import scanpy.external as sce
import anndata as ad
import pandas as pd
import decoupler as dc
import gc as py_gc

##############################
# Science dataset Processing #
##############################

# Data reading
input_file = "SLE/data/Science/GSE174188_CLUES1_adjusted.h5ad"

adata = sc.read_h5ad(input_file)

X = adata.raw.X
obs = pd.DataFrame()
obs['Status'] = adata.obs['SLE_status'].tolist()
obs['ind_cov'] = adata.obs['ind_cov'].tolist()
obs['batch_cov'] = adata.obs['batch_cov'].tolist()
obs['cg_cov'] = adata.obs['cg_cov'].tolist()
var_names = adata.raw.var_names.tolist()
var = pd.DataFrame(index=var_names)
cdata = ad.AnnData(X, obs=obs, var=var, dtype='int32')
cdata.obs_names = adata.obs_names

del(adata)
py_gc.collect()

# Quality control
## annotate the group of mitochondrial genes as 'mt'
cdata.var['mt'] = cdata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(
    cdata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

ncells = cdata.shape[0]

# Cells and genes filtering
sc.pp.filter_cells(cdata, min_genes=100)
sc.pp.filter_cells(cdata, max_genes=4000)
sc.pp.filter_genes(cdata, min_cells=ncells/1000)

# Normalization
sc.pp.normalize_total(cdata)
sc.pp.log1p(cdata)
py_gc.collect()

# Batch effect correction
sc.pp.combat(cdata, key="batch_cov")
py_gc.collect()

################################
# Pediatric dataset Processing #
################################

# Data reading
input_file = "SLE/data/pediatrics/pediatrics_raw.h5ad"

adata = sc.read_h5ad(input_file)
X = adata.raw.X
obs = pd.DataFrame()
obs['Status'] = adata.obs['Condition'].tolist()
obs['ind_cov'] = adata.obs['orig.ident'].tolist()
obs['batch_cov'] = adata.obs['Batch'].tolist()
var_names = adata.var_names.tolist()
var = pd.DataFrame(index=var_names)
pedata = ad.AnnData(X, obs=obs, var=var, dtype='int32')
pedata.obs_names = adata.obs_names

del(adata)
py_gc.collect()

# Quality control
## annotate the group of mitochondrial genes as 'mt'
pedata.var['mt'] = pedata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(
    pedata, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)

# Cells and genes filtering
ncells = pedata.shape[0]
sc.pp.filter_cells(pedata, min_genes=100)
sc.pp.filter_cells(pedata, max_genes=4000)
sc.pp.filter_genes(pedata, min_cells=ncells/1000)


##############
# Processing #
##############

# Subset to common genes
common_genes = list(set(cdata.var_names) & set(pedata.var_names))
cdata = cdata[:, common_genes]
py_gc.collect()

# PCA
initialization = 1
sc.tl.pca(cdata, random_state=initialization)
py_gc.collect()

# Select highly variable genes from Science dataset
sc.pp.highly_variable_genes(cdata, n_top_genes=2000, inplace=True, subset=True)

# Regress by total and mitochondrial expression
sc.pp.regress_out(cdata, ['total_counts', 'pct_counts_mt'])
sc.pp.regress_out(pedata, ['total_counts', 'pct_counts_mt'])

# Subset to highly variable genes
hv_genes = list(cdata.var_names)
cdata_filt = cdata[:, hv_genes]

# Write Science dataset
cdata_filt.write_h5ad("false_discoveries/science_NoScaled.h5ad")
