import argparse

import numpy as np
import anndata as ad
import scanpy as sc

from PIL import Image, ImageFilter, ImageEnhance
from matplotlib import cm

## --------------------------------------------------------------------------------
## Argument parsing
## --------------------------------------------------------------------------------

# Initialize parser
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--file", help="annData *.h5ad file")
parser.add_argument("-o", "--out", help="Output filename", nargs='?')
parser.add_argument("-g", "--genes", help="List of gene names", nargs='?')
parser.add_argument("-c", "--clusters", help="Make clusters image", nargs='?')
parser.add_argument("-n", "--counts", help="Make counts image", nargs='?')

# Read arguments from command line
args = parser.parse_args()

print("\n\n")
print("Input image file: {0}".format(args.file))
print("\n")

if args.out == None:
    outfile =  args.file + '.FullSizeImage_SINGLEGENES.tif'
else:
    outfile = args.out
print("\n\n")
print("Output image file: {0}".format(outfile))
print("\n")

## --------------------------------------------------------------------------------
## Run
## --------------------------------------------------------------------------------


print('Loading...')

ad_filename = args.file
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

def BuildCountsImage(x,y):
    # Find all spots within distance 0.5 of the grid location
    spots = np.where( np.logical_and( np.abs(adata.obsm['X_spatial'][:,0] - x) < downsample, np.abs(adata.obsm['X_spatial'][:,1] - y) < downsample ))
    spots = spots[0]
    if spots.size > 0:
        # Get their cluster labels
        return np.mean(adata.obs['total_counts'])
    else:
        return np.nan

    

OUTIMG = BuildClusterImage(X,Y)
OUTIMG = OUTIMG.astype(float)

OUTIMG_norm = (OUTIMG - np.nanmin(OUTIMG))/(np.nanmax(OUTIMG) - np.nanmin(OUTIMG))
im = Image.fromarray(np.uint8(cm.gist_rainbow(OUTIMG_norm)*255)).convert('RGB')
imc = ImageEnhance.Contrast(im.filter(ImageFilter.BoxBlur(radius=1)).filter(ImageFilter.GaussianBlur(radius=2))).enhance(50)

# Save output image file
imc.save(outfile)

print('Done!')

