#!/usr/bin/env python
## BatchCorrectAndCluster.py
#
# Integrate and preprocess data using scvi-tools
#
# Marcello DiStasio, Oct 2023


import scvi
import scanpy as sc
import anndata as ad

import sys, os
from tqdm import tqdm

from scvi.model.utils import mde
import pymde

###############################################################################
## BUILD OUTPUT SCAFFOLDING
###############################################################################

sc.settings.figdir = './img/integration'

os.makedirs('./calc', exist_ok=True)
os.makedirs('./img/integration', exist_ok=True)

###############################################################################
## LOAD INDIVIDUAL DATASETS
###############################################################################

BASEDIR = './data/'
samples = os.listdir(os.path.join(BASEDIR,'AnnData'))

print(samples)

sys.exit(0)

adata_objects = []
sample_names = []
for s in tqdm(samples):
    adata = adata = ad.read_h5ad(s)
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'],
                               percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 1, :]
    adata_objects.append(adata)
    sample_names.append(s)
    
adata = ad.concat(adata_objects, label="batch", keys=sample_names)
adata.raw = adata
adata.layers["counts"] = adata.X.copy()

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat_v3",
    n_top_genes=2000,
    layer="counts",
    batch_key="batch",
    subset=True
)

# setup and train autoencoder
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()


adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")

sc.tl.leiden(adata, resolution=0.7)

pymde.seed(0)
adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])

fig = sc.pl.embedding(
    adata,
    basis="X_mde",
    color=["batch"],
    frameon=False,
    ncols=1,
    show=False,
    save='_all_unlabeled_batch.png'
)

fig = sc.pl.embedding(
    adata,
    basis="X_mde",
    color=["leiden"],
    frameon=False,
    ncols=1,
    show=False,
    save='_all_unlabeled_leiden.png'
)
