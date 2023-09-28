import numpy as np
import anndata as ad
import scanpy as sc
import sys, os


ad_filename = '/Users/mmd47/Google Drive/My Drive/DiStasio Lab Share/03 Data/Retina_SlideSeq_Curio/Curio_Pilot Study on Control Human Retina/A0052_029/OUTPUT/A0052_029/A0052_029_anndata.h5ad'

adata = ad.read_h5ad(ad_filename)


OUTPUT_IMG_RESOLUTION_X = 2048
OUTPUT_IMG_RESOLUTION_Y = 2048

# QC
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)


## Normalize
print('Normalizing...')
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## Filtering
print('Filtering...')
#sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

## Clustering
print('Clustering...')
sc.tl.pca(adata, svd_solver='arpack', n_comps=30)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)
sc.tl.leiden(adata)
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)
sc.tl.umap(adata, init_pos='paga')



minX = np.min(adata.obsm['X_spatial'][:,0])
maxX = np.max(adata.obsm['X_spatial'][:,0])
minY = np.min(adata.obsm['X_spatial'][:,1])
maxY = np.max(adata.obsm['X_spatial'][:,1])


outgrid = np.meshgrid(np.linspace(minX, maxX, num=OUTPUT_IMG_RESOLUTION_X), np.linspace(minY, maxY, num=OUTPUT_IMG_RESOLUTION_X))








