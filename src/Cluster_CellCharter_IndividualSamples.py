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
parser.add_argument('-n', '--n_clusters', type=int, default=11, help='Number of clusters for CellCharter to find')
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
#filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed.h5ad')
filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_snRNAseq_imputed.h5ad')
print("Loading Data from: " + filename + '...')
samples_all = ad.read_h5ad(filename)

SampleKey = samples_all.uns["SampleKey"]
Samples = list(samples_all.obs['dataset'].cat.categories)


# --------------------------------------------------------------------------------
# Compute Neighborhood Graph
# --------------------------------------------------------------------------------
n_hops = 5
sq.gr.spatial_neighbors(samples_all, coord_type='generic', delaunay=True, spatial_key='X_spatial')
cc.gr.remove_long_links(samples_all)
cc.gr.aggregate_neighbors(samples_all, n_layers=5, use_rep='X_scVI', out_key='X_cellcharter', sample_key='batch') #n_layers = 3 means 1,2,3-hop neighbors



# --------------------------------------------------------------------------------
# Cluster
# --------------------------------------------------------------------------------
n_clusters = args.n_clusters
print(f"Fitting Gaussian Mixture model with {n_clusters} clusters...")

gmm = cc.tl.Cluster(
    n_clusters=n_clusters, 
    random_state=12345,
    covariance_type='full',
    batch_size=None,
    # If running on GPU
    #trainer_params=dict(accelerator='gpu', devices=1, auto_scale_batch_size='binsearch')
    trainer_params=dict(accelerator='gpu', devices=1, default_root_dir=os.path.join(FILEPATHBASE, 'tmp'))
)

import logging

# configure logging on module level, redirect to file
gmm_logger = logging.getLogger(type(gmm.trainer()).__name__)
gmm_logger.setLevel(logging.DEBUG)
Path(os.path.join(FILEPATHBASE, 'tmp')).mkdir(parents=True, exist_ok=True)
gmm_logger.addHandler(logging.FileHandler(os.path.join(FILEPATHBASE, 'tmp', 'gmm.log')))

Samples = list(samples_all.obs['dataset'].cat.categories)


s_all = {}
for r in np.arange(len(Samples)):
    # Select each sample individually
    sample = samples_all[samples_all.obs['dataset']==Samples[r]]
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
        out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_cellcharter_clustered_' + str(n_hops) + 'hops_individual.h5ad')
    else:
        out_filename = args.output

    samples_all.write_h5ad(out_filename)
    print('Saved ' + out_filename)


# --------------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------------

    
with open(os.path.join(FILEPATHBASE,'data', 'retinal_celltype_gates.json')) as f:
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

        try:
        
            sc.tl.rank_genes_groups(sample, groupby='spatial_cluster', use_raw=False, layer='counts_scvi')
            ov = sc.tl.marker_gene_overlap(sample, marker_genes)
            
            fig, ax = plt.subplots(1, 1, figsize=(10,8))
            sns.heatmap(ov, ax=ax, annot=True)
            
            if SAVEFIGS:
                filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_clustered_' + SampleKey[Samples[r]] + '_' + str(n_clusters) + '_clusters_markermatrix.png')
                fig.savefig(filename_out, dpi=300)
                print('Saved: ' + filename_out)
                filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_clustered_' + SampleKey[Samples[r]] + '_' + str(n_clusters) + '_clusters_markermatrix.svg')
                fig.savefig(filename_out, dpi=300)
                print('Saved: ' + filename_out)

        except:
            print("For sample " + samples_all.uns['SampleKey'][Samples[r]] + ": Failure of sc.tl.rank_genes_groups.")

    else:
        print("For sample " + samples_all.uns['SampleKey'][Samples[r]] + ": Only one cluster found.")

    # --------------------------------------------------------------------------------
    # After inspection of the figures, establishment of cluster label for each sample can be performed like this:
    # Use marker genes to establish cluster labels for each sample
    # --------------------------------------------------------------------------------
    #spatial_cluster_label_key = {0 : 'Photoreceptor',
    #                             1 : 'RGC',
    #                             2 : 'Vascular',
    #                             3 : 'Other',
    #                             4 : 'Astrocyte',
    #                             5 : 'Other',
    #                             6 : 'Photoreceptor',
    #                             7 : 'Muller glia',
    #                             8 : 'Muller glia',
    #                             9 : 'Bipolar',
    #                             10: 'Vascular',
    #                             11: 'RGC',
    #                             12: 'Muller glia',
    #                             13: 'Other',
    #                             14: 'Photoreceptor',
    #                             15: 'Other',
    #                             16: 'RPE',
    #                             17: 'Photoreceptor',
    #                             18: 'RPE',
    #                             19: 'Photoreceptor',
    #                             20: 'Other'}
    #
    #samples_all.obs['spatial_cluster_label'] = samples_all.obs['spatial_cluster'].map(spatial_cluster_label_key)
    

    # --------------------------------------------------
    # Plot of spatial scatter montage colored by cluster label
    groups = np.array(sorted(np.unique(sample.obs['spatial_cluster'])))
    nGroupsToColor = len(groups) 
    spect = plt.cm.tab10.resampled(nGroupsToColor)
    newcolors = spect(np.linspace(0,0.5,nGroupsToColor))
    newpalette = ListedColormap(newcolors)
    color_cycler = cycler(color=newpalette.colors)
    

    fig, ax = plt.subplots(1, 1, figsize=(40,30))

    ss = 40
    sq.pl.spatial_scatter(samples_all[samples_all.obs['dataset']==Samples[r]],
                          color='spatial_cluster',
                          size=ss,
                          shape=None,
                          groups=groups,
                          ax=ax)
    #palette=newpalette)
    ax.set_title(SampleKey[Samples[r]])

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
        filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_clustered_' + SampleKey[Samples[r]] + '_' + str(n_clusters) + '_clusters_spatial.png')
        fig.savefig(filename_out, dpi=300)
        print('Saved: ' + filename_out)
        filename_out = os.path.join(IMGDIR, 'Clusters_integrated_imputed_cellcharter_' + str(n_hops) + 'hops_clustered_' + SampleKey[Samples[r]] + '_' + str(n_clusters) + '_clusters_spatial.svg')
        fig.savefig(filename_out, dpi=300)
        print('Saved: ' + filename_out)
        
print('Done. Exiting!')
