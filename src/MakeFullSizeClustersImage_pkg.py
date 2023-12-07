import numpy as np
import anndata as ad
import scanpy as sc

from PIL import Image, ImageFilter, ImageEnhance
from matplotlib import cm


class Cluster:
    def __init__(self, file, out_file):
        self.ad_filename = file
        self.adata = ad.read_h5ad(self.ad_filename)
        self.outfile = out_file

    def process(self):
        print('Loading...')

        # QC
        self.adata.var["mt"] = self.adata.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(self.adata, qc_vars=["mt"], inplace=True)

        ## Normalize
        print('Normalizing...')
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)

        ## Filtering
        print('Filtering...')
        # sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)

        ## Clustering
        print('Clustering...')
        sc.tl.pca(self.adata, svd_solver='arpack', n_comps=30)
        sc.pp.neighbors(self.adata, n_neighbors=20, n_pcs=30)
        sc.tl.leiden(self.adata)
        sc.tl.paga(self.adata)
        sc.pl.paga(self.adata, plot=False)
        sc.tl.umap(self.adata, init_pos='paga')

    def build_image(self, X, Y, downsample):
        # Find all spots within distance 0.5 of the grid location
        spots = np.where(np.logical_and(np.abs(self.adata.obsm['X_spatial'][:, 0] - X) < downsample,
                                        np.abs(self.adata.obsm['X_spatial'][:, 1] - Y) < downsample))
        spots = spots[0]
        if spots.size > 0:
            # Get their cluster labels
            clusters = self.adata.obs['leiden'][spots]
            # Return the most common value
            values, counts = np.unique(clusters, return_counts=True);
            ind = np.argmax(counts);
            return values[ind]
        else:
            return np.nan

    def vectorized_method(self, downsample=1):
        print('Building Output Image...')
        minX = np.min(self.adata.obsm['X_spatial'][:, 0])
        maxX = np.max(self.adata.obsm['X_spatial'][:, 0])
        minY = np.min(self.adata.obsm['X_spatial'][:, 1])
        maxY = np.max(self.adata.obsm['X_spatial'][:, 1])
        X, Y = np.meshgrid(np.arange(minX, maxX, downsample), np.arange(minY, maxY, downsample))
        return np.vectorize(self.build_image)(X, Y, downsample)

    def write_image(self, outimg):
        outimg = outimg.astype(float)
        outimg_norm = (outimg - np.nanmin(outimg)) / (np.nanmax(outimg) - np.nanmin(outimg))
        im = Image.fromarray(np.uint8(cm.gist_rainbow(outimg_norm) * 255)).convert('RGB')
        imc = ImageEnhance.Contrast(
            im.filter(ImageFilter.BoxBlur(radius=1)).filter(ImageFilter.GaussianBlur(radius=2))).enhance(50)
        imc.save(self.outfile + ".tiff")
        print('Done!')