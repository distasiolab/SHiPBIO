# Cluster_CellCharter_IndividualSamples.py
# Part of SHiPBIO
# Marcello DiStasio
# July 2024
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

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from cycler import cycler


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
parser.add_argument('-i', '--inputfile', type=str, help='Path to input AnnData *.h5ad file')
parser.add_argument('-n', '--n_clusters', type=int, default=11, help='Number of clusters for CellCharter to find')
parser.add_argument('-d', '--distance', type=int, default=3, help='Distance; Number of hops to use to build neighborhood graph')
parser.add_argument('-m', '--markers', type=str, default=3, help='Path to marker gene file')
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
    filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_snRNAseq_imputed.h5ad')
else:
    filename = filename = args.inputfile
print("Loading Data from: " + filename + '...')
samples_all = ad.read_h5ad(filename)

SampleKey = samples_all.uns["SampleKey"]
Samples = list(samples_all.obs['dataset'].cat.categories)


# --------------------------------------------------------------------------------
# Preprocess
# --------------------------------------------------------------------------------
print('Filtering outliers in X_scVI projected into MDE space...')
# Filter outliers in X_scVI projected into MDE space
MDE_min_x = -4
MDE_max_x = 4
MDE_min_y = -4
MDE_max_y = 4
m0 = (MDE_min_x < samples_all.obsm['X_scVI_MDE'][:,0]) & (samples_all.obsm['X_scVI_MDE'][:,0] < MDE_max_x)
m1 = (MDE_min_y < samples_all.obsm['X_scVI_MDE'][:,1]) & (samples_all.obsm['X_scVI_MDE'][:,1] < MDE_max_y)
mask = m0 & m1
samples_all = samples_all[mask,:]


print('Log transform scvi_counts')
sc.pp.log1p(samples_all, layer='counts_scvi')





# --------------------------------------------------------------------------------
# Cluster
# --------------------------------------------------------------------------------
n_clusters = args.n_clusters
print(f"Fitting Gaussian Mixture model with {n_clusters} clusters...")

gmm = cc.tl.Cluster(
    n_clusters=n_clusters, 
    random_state=12346,
    covariance_type='full',
    batch_size=100000,
    # If running on GPU
    #trainer_params=dict(accelerator='gpu', devices=1, auto_scale_batch_size='binsearch')
    trainer_params=dict(accelerator='gpu', devices=1, default_root_dir=os.path.join(FILEPATHBASE, 'tmp'))
)

Samples = list(samples_all.obs['dataset'].cat.categories)

s_all = {}
for r in np.arange(len(Samples)):
    # Select each sample individually
    sample = samples_all[samples_all.obs['dataset']==Samples[r]]

    # --------------------------------------------------------------------------------
    # Compute Neighborhood Graph
    # --------------------------------------------------------------------------------
    n_hops = args.distance
    sq.gr.spatial_neighbors(sample, coord_type='generic', delaunay=True, spatial_key='X_spatial')
    cc.gr.remove_long_links(sample)
    cc.gr.aggregate_neighbors(sample, n_layers=n_hops, use_rep='X_scVI', out_key='X_cellcharter', sample_key='batch') #n_layers = 3 means 1,2,3-hop neighbors
    
    # Fit the model for this sample
    gmm.fit(sample, use_rep='X_cellcharter')
    # Cluster on this sample
    sample.obs['spatial_cluster'] = gmm.predict(sample, use_rep='X_cellcharter')
    
    #samples_all[samples_all.obs['dataset']==Samples[r]] = sample
    s_all[Samples[r]] = sample

samples_all = ad.concat(s_all, label="dataset", uns_merge="first", join='outer')
    
# --------------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------------
if SAVEDATA:
    # Save
    if args.output is None:
        out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_cellcharter_clustered_' + str(n_hops) + '_hops_' + str(n_clusters) + '_clusters_' + '_individual.h5ad')
    else:
        out_filename = args.output

    samples_all.write_h5ad(out_filename)
    print('Saved ' + out_filename)


# --------------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------------


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

   
for r in np.arange(len(Samples)):
    # Select each sample individually
    sample = samples_all[samples_all.obs['dataset']==Samples[r]]

    if (len(sample.obs['spatial_cluster'].cat.categories) > 1):

        ov = None
        
        try:
        
            sc.tl.rank_genes_groups(sample, groupby='spatial_cluster', use_raw=False, layer='counts_scvi')
            ov = sc.tl.marker_gene_overlap(sample, marker_genes)
            
            fig, ax = plt.subplots(1, 1, figsize=(10,8))
            sns.heatmap(ov, ax=ax, annot=True)
            
            if SAVEFIGS:
                filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_' + str(n_clusters) + '_clusters_' + SampleKey[Samples[r]] +  '_clustered_markermatrix.png')
                fig.savefig(filename_out, dpi=300)
                print('Saved: ' + filename_out)
                filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_' + str(n_clusters) + '_clusters_' + SampleKey[Samples[r]] +  '_clustered_markermatrix.svg')
                fig.savefig(filename_out, dpi=300)
                print('Saved: ' + filename_out)

        except:
            print("For sample " + samples_all.uns['SampleKey'][Samples[r]] + ": Failure of sc.tl.rank_genes_groups.")

    else:
        print("For sample " + samples_all.uns['SampleKey'][Samples[r]] + ": Only one cluster found.")

# ----------------------------------------------------------
# Plots of spatial scatter montage colored by cluster label
# ----------------------------------------------------------

groups = np.array(sorted(np.unique(samples_all.obs['spatial_cluster'])))
nGroupsToColor = len(groups) 
spect = plt.cm.gist_rainbow.resampled(nGroupsToColor)
newcolors = spect(np.linspace(0,1,nGroupsToColor))
newpalette = ListedColormap(newcolors)
color_cycler = cycler(color=newpalette.colors)
    
for r in np.arange(len(Samples)):
        
    locX = samples_all.obsm['X_spatial'][samples_all.obs['dataset']==Samples[r]][:,0]
    locY = samples_all.obsm['X_spatial'][samples_all.obs['dataset']==Samples[r]][:,1]
    clust = samples_all.obs['spatial_cluster'][samples_all.obs['dataset']==Samples[r]]

    fig, ax = plt.subplots(1, 1, figsize=(40,30))
    for g in groups:
        mask = (clust == g)
        ax.scatter(locX[mask],
                   locY[mask],
                   c=newcolors[g],
                   label=g,
                   cmap=newpalette)
        
    ax.set_title(SampleKey[Samples[r]])
    ax.set_aspect('equal')

    ax.legend(title='Cluster', facecolor='black', bbox_to_anchor=(1.01,0.5), fontsize='xx-large')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.axis('off')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    fig.tight_layout()
    fig.set_facecolor('k')
    for text in fig.findobj(match=lambda x: isinstance(x, plt.Text)):
        if hasattr(text, 'set_color'):
            text.set_color('white')

    if SAVEFIGS:
        n_clusters = len(samples_all.obs['spatial_cluster'].cat.categories)
        filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_' + str(n_clusters) + '_clusters_' + SampleKey[Samples[r]] + '_clustered_spatial.png')
        fig.savefig(filename_out, dpi=300)
        print('Saved: ' + filename_out)
        filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_' + str(n_clusters) + '_clusters_' + SampleKey[Samples[r]] + '_clustered_spatial.svg')
        fig.savefig(filename_out, dpi=300)
        print('Saved: ' + filename_out)
        
print('Done. Exiting!')
