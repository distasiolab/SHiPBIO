import warnings
warnings.filterwarnings('ignore')

import anndata as ad
import scanpy as sc

import scvi
from scvi.external import GIMVI
from scipy import sparse

import numpy as np
import pandas as pd
import os
from pathlib import Path
import argparse


import re
import json

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from cycler import cycler


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
parser.add_argument('-c', '--singlecell', type=str, help='Path to single cell seq *.h5ad file to use for integration')
parser.add_argument('-o', '--output', type=str, help='Path to output *.h5ad file to create')
args = parser.parse_args()


# --------------------------------------------------------------------------------
# File I/O Setup
# --------------------------------------------------------------------------------
FILEPATHBASE = args.basepath

SAVEDATA = True
SAVEFIGS = True
if SAVEFIGS:
    IMGDIR = os.path.join(FILEPATHBASE, 'img', 'out')
    Path(IMGDIR).mkdir(parents=True, exist_ok=True)
    
# --------------------------------------------------------------------------------
# Load datasets 
# --------------------------------------------------------------------------------
# Load the spatial data
filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_harmony_unfiltered.h5ad')
print("Loading Data from: " + filename + ' ...')
samples_all = ad.read_h5ad(filename)

SampleKey = samples_all.uns["SampleKey"]
Samples = list(samples_all.obs['dataset'].cat.categories)


# Load the single nucleus RNA Seq data
if args.singlecell is None:
    filename_sn = os.path.join(FILEPATHBASE, 'data', 'sn_combined.h5ad')
else:
    filename_sn = args.singlecell

print('Loading snRNAseq data from' + filename_sn + ' ...')
samples_sn = ad.read_h5ad(filename_sn)
print('Done')

SampleKey = samples_all.uns["SampleKey"]
Samples = list(samples_all.obs['dataset'].cat.categories)

# --------------------------------------------------------------------------------
# Imputation from snRNAseq data
# --------------------------------------------------------------------------------


# Fit an scVI model to the spatial data
print('Fitting an scVI model to the spatial data...')
scvi.model.SCVI.setup_anndata(samples_all, layer="counts", batch_key="dataset")
model = scvi.model.SCVI(samples_all, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train(early_stopping=True, enable_progress_bar=True)
print('Done!')

print('Getting latent representation')
SCVI_LATENT_KEY = "X_scVI"

latent = model.get_latent_representation()
if hasattr(latent, "todense"):
    latent = latent.todense()
latent = np.array(latent).astype(np.float32)
samples_all.obsm[SCVI_LATENT_KEY] = latent

SCVI_NEIGHBORS_KEY = "neighbors_scVI"
sc.pp.neighbors(samples_all, use_rep=SCVI_LATENT_KEY, key_added=SCVI_NEIGHBORS_KEY)
sc.tl.leiden(samples_all, key_added="leiden_scVI", neighbors_key=SCVI_NEIGHBORS_KEY)

counts_COO = model.posterior_predictive_sample().tocoo()

samples_all.layers['counts_scvi'] = sparse.csr_matrix((counts_COO.data, counts_COO.coords), shape=counts_COO.shape)
print('Done fitting an scVI model to the spatial data. Added samples_all.layers[\'counts_scvi\']')


print('Fitting GIMVI model to the single cell and spatial data...')
# filter genes to be the same on the spatial data and copy the data for the imputation
intersect = np.intersect1d(samples_sn.var_names, samples_all.var_names)
st_adata = samples_all[:, intersect].copy()
sc_adata = samples_sn[:, intersect].copy()


# setup_anndata for spatial and sequencing data
GIMVI.setup_anndata(st_adata, layer="counts", batch_key="dataset")
GIMVI.setup_anndata(sc_adata, layer="counts", labels_key="CellClass")

model = GIMVI(sc_adata, st_adata)
model.train(max_epochs=200, early_stopping=True, enable_progress_bar=True)

# get the latent representations for the sequencing and spatial data
_, imputed = model.get_imputed_values(normalized=True)

def transform(data):
    return np.log(1 + 100 * data)

samples_all[:,st_adata.var_names] = transform(imputed)
samples_all.layers['GIMVI_imputed'] = samples_all.X
print('Done!')

# --------------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------------
if SAVEDATA:
    # Save
    if args.output is None:
        out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_snRNAseq_imputed.h5ad')
    else:
        out_filename = args.output
    
    samples_all.write_h5ad(out_filename)
    print('Saved ' + out_filename)

