#!/usr/bin/env python
## BatchCorrectAndCluster.py
#
# Integrate and preprocess data using scvi-tools
#
# Marcello DiStasio, Oct 2023


import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import scvi
import scanpy as sc
import anndata as ad

import numpy as np
from matplotlib import pyplot as plt
import sys, os
from tqdm import tqdm

from scvi.model.utils import mde
import pymde

###############################################################################
## BUILD OUTPUT SCAFFOLDING
###############################################################################

sc.settings.figdir = './img/integration'

os.makedirs('./calc', exist_ok=True)
outfilename = os.path.join('.','calc','anndata_integrated.h5ad')

IMGDIR = './img/integration'
os.makedirs(IMGDIR, exist_ok=True)


###############################################################################
## LOAD INDIVIDUAL DATASETS
###############################################################################

BASEDIR = './data/'
samples = os.listdir(os.path.join(BASEDIR,'AnnData'))
samplepathf = os.path.join(BASEDIR,'AnnData','{}')

print(samples)

samplenames = []
retinas = []
adata = ad.read_h5ad(samplepathf.format(samples[0]))
retinas.append(adata[adata.obs['Retina_1']].copy())
samplenames.append(os.path.splitext(samples[0])[0] + '_retina1')
retinas.append(adata[adata.obs['Retina_2']].copy())
samplenames.append(os.path.splitext(samples[0])[0] + '_retina2')

adata = ad.read_h5ad(samplepathf.format(samples[1]))
retinas.append(adata[adata.obs['Retina']].copy())
samplenames.append(os.path.splitext(samples[1])[0] + '_retina1')

adata_objects = []
sn = 0
for s in tqdm(retinas):
    adata = s
    sc.pp.filter_genes(adata, min_cells=5)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'],
                               percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 200, :]
    adata = adata[adata.obs.total_counts < 300, :]
    adata = adata[adata.obs.pct_counts_mt < 10, :]
    adata_objects.append(adata)

#    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
#                 jitter=0.4, multi_panel=True, save='QC_'+samplenames[sn], show=False)

    sn=sn+1
    
adata = ad.concat(adata_objects, label="batch", join="outer")

# Clean up the NAs in manual annotation columns in adata.obs, which should be boolean
cs = adata.obs.select_dtypes(include='object').columns
adata.obs[cs] = adata.obs[cs].astype('boolean').fillna(False)


adata.raw = adata
adata.layers["counts"] = adata.X.copy()


# setup and train autoencoder
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()


adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")

sc.tl.leiden(adata, resolution=0.7)

pymde.seed(0)
adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])


fig = sc.pl.embedding(adata, basis="X_mde", color=["batch"], frameon=False, ncols=1, show=False, save=False, return_fig=True)
fig.savefig(os.path.join(IMGDIR,'MDE_batch.png'), dpi=300)

fig = sc.pl.embedding(adata, basis="X_mde", color=["leiden"], frameon=False, ncols=1, show=False, save=False, return_fig=True)
fig.savefig(os.path.join(IMGDIR,'MDE_leiden_cluster.png'), dpi=300)


adata.write(outfilename)
