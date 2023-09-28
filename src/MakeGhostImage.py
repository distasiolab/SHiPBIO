import numpy as np
import anndata as ad
import scanpy as sc
import sys, os
from matplotlib import pyplot as plt
from PIL import Image
from matplotlib import cm



ad_filename = '/Users/mmd47/Google Drive/My Drive/DiStasio Lab Share/03 Data/Retina_SlideSeq_Curio/Curio_Pilot Study on Control Human Retina/A0052_029/OUTPUT/A0052_029/A0052_029_anndata.h5ad'
ad_filename = '/home/mdistasio/YaleGoogleDrive/DiStasio Lab Share/03 Data/Retina_SlideSeq_Curio/Curio_Pilot Study on Control Human Retina/A0052_029/OUTPUT/A0052_029/A0052_029_anndata.h5ad'

print('Loading...')
adata = ad.read_h5ad(ad_filename)

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


print('Building Output Image...')

#fig = sc.pl.embedding(adata, basis='spatial', color='leiden', return_fig=True)


minX = np.min(adata.obsm['X_spatial'][:,0])
maxX = np.max(adata.obsm['X_spatial'][:,0])
minY = np.min(adata.obsm['X_spatial'][:,1])
maxY = np.max(adata.obsm['X_spatial'][:,1])


downsample = 1

X,Y = np.meshgrid(np.arange(minX, maxX, downsample), np.arange(minY, maxY, downsample))


@np.vectorize
def BuildClusterImage(x,y):
    # Find all spots within distance 0.5 of the grid location
    spots = np.where( np.logical_and( np.abs(adata.obsm['X_spatial'][:,0] - x) < downsample, np.abs(adata.obsm['X_spatial'][:,1] - y) < downsample ))
    spots = spots[0]
    if spots.size > 0:
        # Get their cluster labels
        clusters = adata.obs['leiden'][spots]
        # Return the most common value
        values, counts = np.unique(clusters, return_counts=True); ind = np.argmax(counts); 
        return values[ind]
    else:
        return np.nan

OUTIMG = BuildClusterImage(X,Y)
OUTIMG = OUTIMG.astype(float)

OUTIMG_norm = (OUTIMG - np.nanmin(OUTIMG))/(np.nanmax(OUTIMG) - np.nanmin(OUTIMG))
im = Image.fromarray(np.uint8(cm.gist_earth(OUTIMG_norm)*255)).convert('RGB')
im.save('test.tif')

print('Done!')

