# GWAS_Analysis.py
# Part of SHiPBIO
# Marcello DiStasio
# July 2024
##################################################

import warnings
warnings.filterwarnings('ignore')

import scdrs

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
import csv

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from cycler import cycler

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
parser.add_argument('-i', '--inputfile', type=str, help='Path to input AnnData *.h5ad file')
parser.add_argument('-g', '--gwas_worksheet', type=str, default=3, help='Path to *.csv file with list of GWAS datafiles')
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
    filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_magic.h5ad')
else:
    filename = filename = args.inputfile
print("Loading Data from: " + filename + '...')
samples_all = ad.read_h5ad(filename)


# Delete duplicate indices
samples_all = samples_all[~samples_all.obs.index.duplicated(keep='first')]
SampleKey = samples_all.uns["SampleKey"]
Samples = list(samples_all.obs['dataset'].cat.categories)
print('Done')

        
# --------------------------------------------------------------------------------------
# Load GWAS data 
# --------------------------------------------------------------------------------------
def read_csv_into_dict(filename, known_columns):
    result = []
    with open(filename, mode='r', newline='\n') as file:
        reader = csv.DictReader(file, delimiter=',')
        for row in reader:
            entry = {}
            unknown_cols = []
            for key, value in row.items():
                if key.strip() in known_columns:
                    entry[key.strip()] = value.strip()
                else:
                    if value is not None:
                        unknown_cols.append(value.strip())
            entry['Annotations'] = unknown_cols
            result.append(entry)
    return result

known_columns =['Condition', 'filename']
data = read_csv_into_dict(args.gwas_worksheet, known_columns)

GWAS_data = {}
cnt = 1
for c in data:
    condition = c['Condition']
    filename = c['filename']
    print('Loading GWAS data from' + filename + '(' + str(cnt) + '/' + str(len(data)) + ')...')
    GWAS_data[condition] = scdrs.util.load_gs(filename, 'hsapiens', 'hsapiens', samples_all.var_names.tolist())
    cnt=cnt+1
print("Done loading data")


# --------------------------------------------------------------------------------
# scDRS processing
# https://www.nature.com/articles/s41588-022-01167-z
# --------------------------------------------------------------------------------

samples_all_gwas = samples_all.copy()
#samples_all_gwas.raw = samples_all_gwas.layers['counts_scvi']
sc.external.pp.magic(samples_all_gwas)
#samples_all_gwas.obsp = samples_all.obsp
#samples_all_gwas.uns = samples_all.uns

print('QC and filtering')
sc.pp.calculate_qc_metrics(samples_all_gwas)
sc.pp.filter_cells(samples_all_gwas, min_counts=100)
sc.pp.filter_genes(samples_all_gwas, min_cells=10)
print('Normalization...')
sc.pp.normalize_total(samples_all_gwas, target_sum=1e4)
sc.pp.log1p(samples_all_gwas)
print('Computing nearest neighbors distance matrix')
sc.pp.neighbors(samples_all_gwas, use_rep='X', n_neighbors=32)
print('Done!')

df_cov = pd.DataFrame(index=samples_all_gwas.obs.index)
df_cov["const"] = 1
df_cov["n_counts"] = samples_all_gwas.obs["n_counts"]

print('scDRS preproccessing')
scdrs.preprocess(samples_all_gwas, adj_prop='spatial_cluster_label', cov=df_cov)
print('Done!')

for disease in GWAS_data.keys():
    d = disease
    df = scdrs.score_cell(samples_all_gwas, GWAS_data[d][d][0], GWAS_data[d][d][1])
        
    samples_all_gwas.uns['scdrs_' + disease + '_celltypes'] = scdrs.method.downstream_group_analysis(samples_all_gwas, df, group_cols=['spatial_cluster_label'])
    print(samples_all_gwas.uns['scdrs_' + disease + '_celltypes'])

    df.columns = ['scdrs_' + disease + '_' + c for c in df.columns]
    samples_all_gwas.obs = pd.concat([samples_all_gwas.obs, df], axis=1)

if SAVEDATA:

    # Delete the MAGIC imputation to save space
    samples_all_gwas.X = samples_all_gwas.layers['counts']
    try:
        del samples_all.layers['counts_magic']
    except:
        pass
    try:
        del samples_all_gwas.obsm['X_magic']
    except:
        pass
        
    # Save
    if args.output is None:
        out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_magic_GWAS_SCDRS.h5ad') 
    else:
        out_filename = args.output

    samples_all_gwas.write_h5ad(out_filename)
    print('Saved ' + out_filename)


# ----------------------------------------------------------------------------------------------
# 3d Plots
# ----------------------------------------------------------------------------------------------
print('Generating Plots...')

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
newcmap = dict(map(lambda i,j : (i,j) , groups, newcolors))

for disease in GWAS_data.keys():
    for r in np.arange(len(Samples)):

        sample = samples_all_gwas[samples_all_gwas.obs['dataset']==Samples[r]]

        fig, ax = plt.subplots(1, 1, figsize=(60,30), subplot_kw=dict(projection='3d'))

        x = np.array(sample.obsm['X_spatial'][:,0])
        y = np.array(sample.obsm['X_spatial'][:,1])
        z0 = np.zeros_like(x)
        z1 = np.array(sample.obs['scdrs_' + disease + '_raw_score'])
        
        if (len(z1) > 0):
            threshold = np.percentile(z1, 95)
            z1[z1<threshold]=0

            z1 = 0.1*(z1-np.min(z1))/(np.max(z1)-np.min(z1))

            c = list(map(newcmap.get, sample.obs['spatial_cluster_label'].tolist()))
            graycolor = np.array([0.4, 0.4, 0.4, 0.2])
            c = [graycolor if v is None else v for v in c]

            ax.bar3d(x, y, z0, 30, 30, z1, color=c, shade=True)
            ax.set_aspect('equalxy')
            ax.set_zlim(0,1)
            ax.view_init(elev=60, azim=120)
            ax.axis('off')
            ax.set_facecolor('k')

            # Create custom legend entries with colors and labels
            legend_colors = newcolors
            legend_labels = [g + ' enriched' for g in groups]
            
            # Add the custom legend
            from matplotlib.lines import Line2D
            
            legend_patches = [Line2D([0], [0], color='w', marker='s', markersize=10, markerfacecolor=color) for color in legend_colors]
            
            # Add the custom legend
            ax.legend(legend_patches, legend_labels, facecolor='k', fontsize="18", loc = "lower right")
            
            for text in fig.findobj(match=lambda x: isinstance(x, plt.Text)):
                if hasattr(text, 'set_color'):
                    text.set_color('white')
                    
            if SAVEFIGS:

                Path(os.path.join(IMGDIR,'GWAS_3D')).mkdir(parents=True, exist_ok=True)
                filename_base = os.path.join(IMGDIR, 'GWAS_3D', 'GIMVI_-_CellCharter_Clusters_-_scDRS_Scores_3D_-_' + disease + '_-_Sample_' + SampleKey[Samples[r]])
                
                filename_out = filename_base + '.png'
                fig.savefig(filename_out, dpi=300)
                print('Saved: ' + filename_out)
                filename_out = filename_base + '.svg'
                fig.savefig(filename_out, dpi=300)
                print('Saved: ' + filename_out)

print('Done')
