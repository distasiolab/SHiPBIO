import anndata as ad
import squidpy as sq
import scanpy as sc

BASEDIR = './data/'
samples = os.listdir(os.path.join(BASEDIR,'AnnData'))
samplepathf = os.path.join(BASEDIR,'AnnData','{}')

adata = ad.read_h5ad(samplepathf.format(samples[0]))
retina1 = adata[adata.obs['Retina_1']]
retina2 = adata[adata.obs['Retina_2']]

adata = ad.read_h5ad(samplepathf.format(samples[1]))
retina3 = adata[adata.obs['Retina']]
