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

import re
import json

import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from cycler import cycler


# --------------------------------------------------------------------------------
# File I/O Setup
# --------------------------------------------------------------------------------

# OS X Laptop
# FILEPATHBASE = '/Users/mmd47/Library/CloudStorage/GoogleDrive-mmd47@yale.edu/My Drive/DiStasio Lab/DiStasio Lab Share/'

# Fedora Desktop
# conda activate cellcharter-env
#FILEPATHBASE = '/home/mdistasio/YaleGoogleDrive/DiStasio Lab/DiStasio Lab Share/'
FILEPATHBASE = '/home/mdistasio/Workspace/'

# --------------------------------------------------------------------------------

SAVEDATA = True
SAVEFIGS = True
if SAVEFIGS:
    IMGDIR = os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate', 'img')
    Path(IMGDIR).mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------------
# Load datasets 
# --------------------------------------------------------------------------------
filename_retinas_all_magic = os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate', 'calc', 'retinas_all_unfiltered_snRNAseq_imputed_GIMVI_MAGIC.h5ad')
if os.path.isfile(filename_retinas_all_magic):
    print("Loading Data from: " + filename_retinas_all_magic + ' ...')
    retinas_all_magic = ad.read_h5ad(filename_retinas_all_magic)
else:
    filename = os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate', 'calc', 'retinas_all_unfiltered_snRNAseq_imputed_GIMVI_ClustersLabeled.h5ad')
    print("Loading Data from: " + filename + ' ...')
    retinas_all = ad.read_h5ad(filename)

    # Delete duplicate indices
    retinas_all = retinas_all[~retinas_all.obs.index.duplicated(keep='first')]
    SampleKey = retinas_all.uns["SampleKey"]
    Samples = list(retinas_all.obs['dataset'].cat.categories)
    print('Done')

    # --------------------------------------------------------------------------------
    ## MAGIC imputation
    # --------------------------------------------------------------------------------
    print('Imputation with MAGIC')
    retinas_all_magic = sc.external.pp.magic(retinas_all, copy=True)
    retinas_all_magic.obsp = retinas_all.obsp
    retinas_all_magic.uns = retinas_all.uns
    print('Done!')
    
    if SAVEDATA:
        # Save
        filename_out = os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate', 'calc', 'retinas_all_unfiltered_snRNAseq_imputed_GIMVI_MAGIC.h5ad')
        retinas_all_magic.write_h5ad(filename_out)
        print('Saved ' + filename_out)
        


# --------------------------------------------------------------------------------------
# Load GWAS data 
# --------------------------------------------------------------------------------------

GWAS_files = {'AMD'   : os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate','data', 'gwas', 'AMD_Fritsche-26691988.gs'),
              'MS'    : os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate','data', 'gwas', 'MS_IMSGC-24076602.gs'),
              'Glaucoma'   : os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate','data', 'gwas', 'Glaucoma_GlobalBiobank_36777996_GCST90399726.gs'),
              'Alzheimers' : os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate','data', 'gwas', 'Alzheimers_35379992-GCST90027158-MONDO_0004975_PMID_35379992.gs')
             }

GWAS_data = {}
for disease in GWAS_files.keys():
    filename = GWAS_files[disease]
    print('Loading GWAS data from' + filename + '...')
    GWAS_data[disease] = scdrs.util.load_gs(filename, 'hsapiens', 'hsapiens', retinas_all_magic.var_names.tolist())
    print('Done.')


# --------------------------------------------------------------------------------
# scDRS processing
# https://www.nature.com/articles/s41588-022-01167-z
# --------------------------------------------------------------------------------

#retinas_all_gwas = retinas_all.copy()
retinas_all_gwas = retinas_all_magic.copy()
print('QC and filtering')
sc.pp.calculate_qc_metrics(retinas_all_gwas)
sc.pp.filter_cells(retinas_all_gwas, min_counts=100)
sc.pp.filter_genes(retinas_all_gwas, min_cells=10)
print('Normalization...')
sc.pp.normalize_total(retinas_all_gwas, target_sum=1e4)
sc.pp.log1p(retinas_all_gwas)
print('Done!')

df_cov = pd.DataFrame(index=retinas_all_gwas.obs.index)
df_cov["const"] = 1
df_cov["n_counts"] = retinas_all_gwas.obs["n_counts"]

print('scDRS preproccessing')
scdrs.preprocess(retinas_all_gwas, adj_prop='spatial_cluster_label', cov=df_cov)
print('Done!')

for disease in GWAS_files.keys():
    d = disease
    df = scdrs.score_cell(retinas_all_gwas, GWAS_data[d][d][0], GWAS_data[d][d][1])
        
    retinas_all_gwas.uns['scdrs_' + disease + '_celltypes'] = scdrs.method.downstream_group_analysis(retinas_all_gwas, df, group_cols=['spatial_cluster_label'])
    print(retinas_all_gwas.uns['scdrs_' + disease + '_celltypes'])

    df.columns = ['scdrs_' + disease + '_' + c for c in df.columns]
    retinas_all_gwas.obs = pd.concat([retinas_all_gwas.obs, df], axis=1)
    
    if SAVEDATA:
        # Save
        filename_out = os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate', 'calc', 'retinas_all_unfiltered_snRNAseq_imputed_GIMVI_ClustersLabeled_scDRS_GWAS_scores.h5ad')
        retinas_all_gwas.write_h5ad(filename_out)
        print('Saved ' + filename_out)    
print('Done!')




# ----------------------------------------------------------------------------------------------
# 3d Plots
# ----------------------------------------------------------------------------------------------
print('Generating Plots')
for disease in GWAS_files.keys():
    for r in np.arange(len(Samples)):

        retina = retinas_all[retinas_all.obs['dataset']==Samples[r]]


        # Selection of groups to plot
        groups = np.array(['RGC',
                           'Astrocyte',
                           'Muller glia',
                           'Photoreceptor',
                           'RPE',
                           'Vascular'])
        #groups = np.array(retinas_all.obs['spatial_cluster_label'].cat.categories.tolist())
        nGroupsToColor = len(groups)
        spect = plt.cm.gist_rainbow.resampled(nGroupsToColor)
        newcolors = spect(np.linspace(0,1,nGroupsToColor))
        newcolors = np.flipud(newcolors * np.array([0.8, 0.8, 0.8, 1]))
        newcmap = dict(map(lambda i,j : (i,j) , groups, newcolors))




        fig, ax = plt.subplots(1, 1, figsize=(60,30), subplot_kw=dict(projection='3d'))

        x = np.array(retina.obsm['X_spatial'][:,0])
        y = np.array(retina.obsm['X_spatial'][:,1])
        z0 = np.zeros_like(x)
        z1 = np.array(retina.obs['scdrs_' + disease + '_raw_score'])
        threshold = np.percentile(z1, 95)
        z1[z1<threshold]=0

        z1 = 0.1*(z1-np.min(z1))/(np.max(z1)-np.min(z1))

        c = list(map(newcmap.get, retina.obs['spatial_cluster_label'].tolist()))
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

            Path(os.path.join(IMGDIR,'GWAS_3D_raw')).mkdir(parents=True, exist_ok=True)
            filename_base = os.path.join(IMGDIR, 'GWAS_3D', 'GIMVI_-_CellCharter_Clusters_-_scDRS_Scores_3D_-_' + disease + '_-_Sample_' + SampleKey[Samples[r]])

            filename_out = filename_base + '.png'
            fig.savefig(filename_out, dpi=300)
            print('Saved: ' + filename_out)
            filename_out = filename_base + '.svg'
            fig.savefig(filename_out, dpi=300)
            print('Saved: ' + filename_out)

print('Done')
