# Cluster_CellCharter_Subcluster.py
# Part of SHiPBIO
# Marcello DiStasio
# Sept 2024
##################################################

import warnings
warnings.filterwarnings('ignore')

from pytorch_lightning import Trainer
from pytorch_lightning.tuner import batch_size_scaling

import cellcharter as cc

import anndata as ad
import squidpy as sq
import scanpy as sc

import numpy as np
import pandas as pd
import os
from pathlib import Path
import argparse


import re
import json

#import seaborn as sns
#from matplotlib import pyplot as plt
#from matplotlib.colors import ListedColormap
#from cycler import cycler


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
parser.add_argument('-i', '--inputfile', type=str, help='Path to input AnnData *.h5ad file')
parser.add_argument('-n', '--n_clusters', type=int, default=11, help='Number of clusters for CellCharter to find')
parser.add_argument('-d', '--distance', type=int, default=3, help='Distance; Number of hops to use to build neighborhood graph')
#parser.add_argument('-m', '--markers', type=str, default=3, help='Path to marker gene file')
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
if args.inputfile is None:
    filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_cellcharter_3hop_clustered.h5ad')
else:
    filename = args.inputfile
print("Loading Data from: " + filename + '...')
samples_all = ad.read_h5ad(filename)

SampleKey = samples_all.uns["SampleKey"]
Samples = list(samples_all.obs['dataset'].cat.categories)


# Arrange samples in space
print('Arrange samples in space...')
samples_all.obsm['X_spatial_fov'] = samples_all.obsm['X_spatial'].copy()
y_max_p = [np.max(samples_all[samples_all.obs['dataset']==Samples[r]].obsm['X_spatial'][:,1]) for r in np.arange(1,len(Samples))]
y_offsets = np.append(0,np.cumsum(y_max_p)+10000)
for r in np.arange(1,len(Samples)):
    samples_all.obsm['X_spatial_fov'][samples_all.obs['dataset']==Samples[r],1] = samples_all.obsm['X_spatial_fov'][samples_all.obs['dataset']==Samples[r],1] + y_offsets[r]

n_hops = args.distance
sq.gr.spatial_neighbors(samples_all, coord_type='generic', delaunay=True, spatial_key='X_spatial_fov')
cc.gr.remove_long_links(samples_all)
print(f"Aggregating neighbors with {n_hops}-hops...")
cc.gr.aggregate_neighbors(samples_all, n_layers=n_hops, use_rep='X_scVI', out_key='X_cellcharter_subcluster', sample_key='batch') #n_layers = 3 means 1,2,3-hop neighbors

# --------------------------------------------------------------------------------
# Subclustering
# --------------------------------------------------------------------------------
groups = np.array(sorted(np.unique(samples_all.obs['spatial_cluster_label'])))

n_clusters = args.n_clusters
    
s_all = {}
for group in groups:
    samples_all_group = samples_all[samples_all.obs['spatial_cluster_label'] == group]
    print(f"Fitting model with AutoK for number of clusters to {group} data...")
    
    autok = cc.tl.ClusterAutoK(
        n_clusters=(2,10), 
        max_runs=10,
        convergence_tol=0.001
    )

    autok.fit(samples_all_group, use_rep='X_cellcharter_subcluster')
    samples_all_group.obs['spatial_subcluster'] = autok.predict(samples_all_group, use_rep='X_cellcharter_subcluster')
    
    # gmm = cc.tl.Cluster(n_clusters=n_clusters,
    #                     random_state=12345,
    #                     covariance_type='full',
    #                     batch_size=256,
    #                     trainer_params=dict(accelerator='gpu', devices=1, default_root_dir=os.path.join(FILEPATHBASE, 'tmp')))
    
    # gmm.fit(samples_all_group, use_rep='X_cellcharter_subcluster')

    # samples_all_group.obs['spatial_subcluster'] = gmm.predict(samples_all_group, use_rep='X_cellcharter_subcluster')

    s_all[group] = samples_all_group

samples_all = ad.concat(s_all, label="spatial_cluster_label", uns_merge="first", join='outer')

    
# --------------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------------
if SAVEDATA:
    # Save
    try:
        del samples_all.layers['counts_magic']
    except:
        pass
    
    if args.output is None:
        out_filename = os.path.join(FILEPATHBASE, 'calc', f'samples_all_integrated_imputed_cellcharter_clustered_subclustered_{n_hops}hops.h5ad')
    else:
        out_filename = args.output

    samples_all.write_h5ad(out_filename)
    print('Saved ' + out_filename)

print('Done. Exiting!')
