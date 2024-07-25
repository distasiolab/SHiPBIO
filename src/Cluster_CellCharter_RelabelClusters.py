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
parser.add_argument('-i', '--inputfile', type=str, help='Path to input AnnData *.h5ad file')
parser.add_argument('-c', '--clusterlabels', type=str, help='Path to *.json file with labels for clusters of each individual sample')
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
filename = args.inputfile
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
    if args.output is None:
        out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_cellcharter_clustered_individual_labeled.h5ad')
    else:
        out_filename = args.output

    samples_all.write_h5ad(out_filename)
    print('Saved ' + out_filename)


# --------------------------------------------------------------------------------
# Plotting
# --------------------------------------------------------------------------------

# --------------------------------------------------
# Plot of spatial scatter montage colored by cluster label
groups = np.array(sorted(np.unique(samples_all.obs['spatial_cluster_label'])))
nGroupsToColor = len(groups) 
spect = plt.cm.gist_rainbow.resampled(nGroupsToColor)
newcolors = spect(np.linspace(0,1,nGroupsToColor))
try:
    newcolors[np.where(groups=='Other')[0][0],:] = [0.4,0.4,0.4,1] # Make 'Other' gray
    newcolors[np.where(groups=='Photoreceptor')[0][0],:] = np.array([50, 136, 189, 255])/255
    newcolors[np.where(groups=='RPE')[0][0],:]           = np.array([254, 224, 139, 255])/255
    newcolors[np.where(groups=='MullerGlia')[0][0],:]    = np.array([158, 1, 66, 255])/255
    newcolors[np.where(groups=='RGC')[0][0],:]           = np.array([244, 109, 67, 255])/255
    newcolors[np.where(groups=='Astrocyte')[0][0],:]     = np.array([102, 194, 165, 255])/255
    newcolors[np.where(groups=='Vascular')[0][0],:]      = np.array([94, 79, 162, 255])/255
    newcolors[np.where(groups=='Amacrine')[0][0],:]      = np.array([213, 62, 79, 255])/255
except:
    pass

newpalette = ListedColormap(newcolors)
color_cycler = cycler(color=newpalette.colors)



nRow = 3
nCol = int(np.ceil(len(Samples)/3))
fig, axs = plt.subplots(nRow, nCol, figsize=(60,30))

for r in np.arange(len(Samples)):
    ax = axs.reshape(-1)[r]
            
    locX = samples_all.obsm['X_spatial'][samples_all.obs['dataset']==Samples[r]][:,0]
    locY = samples_all.obsm['X_spatial'][samples_all.obs['dataset']==Samples[r]][:,1]
    clust = samples_all.obs['spatial_cluster_label'][samples_all.obs['dataset']==Samples[r]]

    for g in np.arange(len(groups)):
        mask = (clust == groups[g])
        ax.scatter(locX[mask],
                   locY[mask],
                   c=newcolors[g],
                   label=groups[g],
                   s=1,
                   cmap=newpalette)

    
    ax.set_title(SampleKey[Samples[r]])

    if (r == len(Samples)-1):
        lgnd = ax.legend(title='Cluster', facecolor='black', bbox_to_anchor=(1.01,0.5), fontsize='xx-large')
        for h in lgnd.legendHandles:
            h._sizes = [30]
    ax.set_title(SampleKey[Samples[r]])
    ax.set_aspect('equal')
    
for ax in axs.reshape(-1):
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
    filename_out = os.path.join(IMGDIR, 'AllSamples_integrated_imputed_cellcharter_clustered_individually_clusters_labeled.png')
    fig.savefig(filename_out, dpi=300)
    print('Saved: ' + filename_out)
    filename_out = os.path.join(IMGDIR, 'AllSamples_integrated_imputed_cellcharter_clustered_individually_clusters_labeled.svg')
    fig.savefig(filename_out, dpi=300)
    print('Saved: ' + filename_out)
# --------------------------------------------------




        
print('Done. Exiting!')
