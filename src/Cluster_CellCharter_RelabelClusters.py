# Cluster_CellCharter_RelabelClusters.py
# Part of SHiPBIO
# Marcello DiStasio
# July 2024
##################################################

import warnings
warnings.filterwarnings('ignore')

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
parser.add_argument('-c', '--clusterlabels', type=str, help='Path to *.json file with labels for clusters of each individual sample')
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
filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_cellcharter_clustered_individual.h5ad')
print("Loading Data from: " + filename + '...')
samples_all = ad.read_h5ad(filename)

SampleKey = samples_all.uns["SampleKey"]
Samples = list(samples_all.obs['dataset'].cat.categories)

filename = args.clusterlabels
print("Loading cluster labels from: " + filename + '...')
with open(args.clusterlabels) as f:
    cldata = json.load(f)

samples_all.obs['spatial_cluster_label'] = ''
for r in np.arange(len(Samples)):
    label_dict = cldata[np.where(np.array([d.get('SampleName', None) for d in cldata]) == SampleKey[Samples[r]])[0][0]]['cluster_labels']
    #Repack label dictionary
    label_dict = {int(key): value for d in label_dict for key, value in d.items()}

    samples_all.obs.loc[(samples_all.obs['dataset']==Samples[r]),'spatial_cluster_label'] = samples_all.obs.loc[(samples_all.obs['dataset']==Samples[r]),'spatial_cluster'].map(label_dict)
    
# --------------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------------
if SAVEDATA:
    # Save
    filename_out = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_cellcharter_clustered_individual_labeled.h5ad')
    samples_all.write_h5ad(filename_out)
    print('Saved ' + filename_out)


# --------------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------------

# --------------------------------------------------
# Plot of spatial scatter montage colored by cluster label
groups = np.array(sorted(np.unique(samples_all.obs['spatial_cluster_label'])))
nGroupsToColor = len(groups) 
spect = plt.cm.tab10.resampled(nGroupsToColor)
newcolors = spect(np.linspace(0,0.5,nGroupsToColor))
newpalette = ListedColormap(newcolors)
color_cycler = cycler(color=newpalette.colors)
    
nRow = 3
nCol = int(np.ceil(len(Samples)/3))
fig, axs = plt.subplots(nRow, nCol, figsize=(40,30))
for r in np.arange(len(Samples)):
    ax = axs.reshape(-1)[r]
    ss = 1
    sq.pl.spatial_scatter(samples_all[samples_all.obs['dataset']==Samples[r]],
                          color='spatial_cluster_label',
                          size=ss,
                          shape=None,
                          groups=groups,
                          ax=ax)#,
#                          palette=newpalette)
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
    filename_out = os.path.join(IMGDIR, 'AllSamples_integrated_imputed_cellcharter_3hop_clustered_individually_clusters_labeled.png')
    fig.savefig(filename_out, dpi=300)
    print('Saved: ' + filename_out)
    filename_out = os.path.join(IMGDIR, 'AllSamples_integrated_imputed_cellcharter_3hop_clustered_individually_clusters_labeled.svg')
    fig.savefig(filename_out, dpi=300)
    print('Saved: ' + filename_out)
# --------------------------------------------------




        
print('Done. Exiting!')
