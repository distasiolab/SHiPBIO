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

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from cycler import cycler

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
parser.add_argument('-n', '--n_clusters', type=int, default=11, help='Number of clusters for CellCharter to find')
parser.add_argument('-d', '--distance', type=int, default=3, help='Distance; Number of hops to use to build neighborhood graph')
parser.add_argument('-m', '--markers', type=str, default=3, help='Path to marker gene file')
parser.add_argument('-o', '--output', type=str, help='Path to output *.h5ad file to create')
args = parser.parse_args()



# --------------------------------------------------------------------------------
# File I/O Setup
# --------------------------------------------------------------------------------
FILEPATHBASE = args.basepath

SAVEDATA = False
SAVEFIGS = True
if SAVEFIGS:
    IMGDIR = os.path.join(FILEPATHBASE, 'img', 'out')
    Path(IMGDIR).mkdir(parents=True, exist_ok=True)
    
# --------------------------------------------------------------------------------
# Load datasets 
# --------------------------------------------------------------------------------
filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_snRNAseq_imputed.h5ad')
print("Loading Data from: " + filename + '...')
samples_all = ad.read_h5ad(filename)

SampleKey = samples_all.uns["SampleKey"]
Samples = list(samples_all.obs['dataset'].cat.categories)



# Arrange samples in space
samples_all.obsm['X_spatial_fov'] = samples_all.obsm['X_spatial'].copy()
y_max_p = [np.max(samples_all[samples_all.obs['dataset']==Samples[r]].obsm['X_spatial'][:,1]) for r in np.arange(1,len(Samples))]
y_offsets = np.append(0,np.cumsum(y_max_p))
for r in np.arange(1,len(Samples)):
    samples_all.obsm['X_spatial_fov'][samples_all.obs['dataset']==Samples[r],1] = samples_all.obsm['X_spatial_fov'][samples_all.obs['dataset']==Samples[r],1] + y_offsets[r]

# --------------------------------------------------------------------------------
# Compute Neighborhood Graph
# --------------------------------------------------------------------------------
n_hops = args.distance
sq.gr.spatial_neighbors(samples_all, coord_type='generic', delaunay=True, spatial_key='X_spatial_fov')
cc.gr.remove_long_links(samples_all)
cc.gr.aggregate_neighbors(samples_all, n_layers=n_hops, use_rep='X_scVI', out_key='X_cellcharter', sample_key='batch') #n_layers = 3 means 1,2,3-hop neighbors


# --------------------------------------------------------------------------------
# Cluster
# --------------------------------------------------------------------------------
n_clusters = args.n_clusters
print(f"Fitting Gaussian Mixture model with {n_clusters} clusters...")


gmm = cc.tl.Cluster(
    n_clusters=n_clusters, 
    random_state=12345,
    covariance_type='full',
    batch_size=2,
    # If running on GPU
    #trainer_params=dict(accelerator='gpu', devices=1, auto_scale_batch_size='binsearch')
    trainer_params=dict(accelerator='gpu', devices=1, default_root_dir=os.path.join(FILEPATHBASE, 'tmp'))
)


import logging


# configure logging on module level, redirect to file
gmm_logger = logging.getLogger(type(gmm.trainer()).__name__)
gmm_logger.setLevel(logging.DEBUG)
gmm_logger.addHandler(logging.FileHandler(os.path.join(FILEPATHBASE, 'tmp', 'gmm.log')))

                      
gmm.fit(samples_all, use_rep='X_cellcharter')
samples_all.obs['spatial_cluster'] = gmm.predict(samples_all, use_rep='X_cellcharter')

# --------------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------------
if SAVEDATA:
    # Save
    if args.output is None:
        out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_cellcharter_clustered_' + str(n_hops) + '_hops_ ' + str(n_clusters) + '_clusters.h5ad')
    else:
        out_filename = args.output

    samples_all.write_h5ad(out_filename)
    print('Saved ' + out_filename)



# --------------------------------------------------------------------------------
# Plots of Gene Marker vs. Cluster ID Matrix
# --------------------------------------------------------------------------------

with open(args.markers) as f:
    gates = json.load(f)

marker_genes = {}
for item in gates:
    cell_type = item['cell']
    gene_list = [gene['gene'] for gene in item['gates']]
    marker_genes[cell_type] = gene_list
  
ov = None

try:

    sc.tl.rank_genes_groups(samples_all, groupby='spatial_cluster', use_raw=False, layer='counts_scvi')
    ov = sc.tl.marker_gene_overlap(samples_all, marker_genes)
    
    fig, ax = plt.subplots(1, 1, figsize=(10,8))
    sns.heatmap(ov, ax=ax, annot=True)
    
    if SAVEFIGS:
        filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_ ' + str(n_clusters) + '_clusters_AllSamples_clusters_markermatrix.png')
        fig.savefig(filename_out, dpi=300)
        print('Saved: ' + filename_out)
        filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_ ' + str(n_clusters) + '_clusters_AllSamples_clusters_markermatrix.svg')
        fig.savefig(filename_out, dpi=300)
        print('Saved: ' + filename_out)

except:
    print("Failure of sc.tl.rank_genes_groups.")




# --------------------------------------------------
# Plot of spatial scatter montage colored by cluster label
groups = np.array(sorted(np.unique(samples_all.obs['spatial_cluster'])))
nGroupsToColor = len(groups) 
spect = plt.cm.gist_rainbow.resampled(nGroupsToColor)
newcolors = spect(np.linspace(0,1,nGroupsToColor))
newpalette = ListedColormap(newcolors)
color_cycler = cycler(color=newpalette.colors)

nRow = 3
nCol = int(np.ceil(len(Samples)/3))
fig, axs = plt.subplots(nRow, nCol, figsize=(40,30))
for r in np.arange(len(Samples)):
    ax = axs.reshape(-1)[r]
            
    locX = samples_all.obsm['X_spatial'][samples_all.obs['dataset']==Samples[r]][:,0]
    locY = samples_all.obsm['X_spatial'][samples_all.obs['dataset']==Samples[r]][:,1]
    clust = samples_all.obs['spatial_cluster'][samples_all.obs['dataset']==Samples[r]]

    for g in groups:
        mask = (clust == g)
        ax.scatter(locX[mask],
                   locY[mask],
                   c=newcolors[g],
                   label=g,
                   cmap=newpalette)

    
    ax.set_title(SampleKey[Samples[r]])

for ax in axs.reshape(-1):
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.axis('off')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
#    try:
#        ax.get_legend().remove()
#    except:
#        pass

fig.tight_layout()

fig.set_facecolor('k')
for text in fig.findobj(match=lambda x: isinstance(x, plt.Text)):
    if hasattr(text, 'set_color'):
        text.set_color('white')

#labels_handles = {  label: handle for ax in fig.axes for handle, label in zip(*ax.get_legend_handles_labels())    }
#fig.legend( labels_handles.values(), labels_handles.keys(), loc = "lower right", ncol=3)


if SAVEFIGS:
    n_clusters = len(samples_all.obs['spatial_cluster'].cat.categories)
    filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + '_hops_ ' + str(n_clusters) + '_clusters_AllSamples_clustered_spatial.png')
    fig.savefig(filename_out, dpi=300)
    print('Saved: ' + filename_out)
    filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + '_hops_ ' + str(n_clusters) + '_clusters_AllSamples_clustered_spatial.svg')
    fig.savefig(filename_out, dpi=300)
    print('Saved: ' + filename_out)
# --------------------------------------------------
