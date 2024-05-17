# Preprocess.py
# Part of anndata_manualAnnotate for processing spatial omics data
# Marcello DiStasio
# Oct 2023


import argparse

import anndata as ad
import squidpy as sq
import scanpy as sc

import os
from pathlib import Path
    
parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
args = parser.parse_args()

FILEPATHBASE = args.basepath

SAVEFIGS = True
if SAVEFIGS:
    IMGDIR = os.path.join(FILEPATHBASE,'img')


# --------------------------------------------------------------------------------
# If manual spatial annotations are desired, then Prior to this loading, each
# *_anndata.h5ad file should be annotated with
#
#   1. MakeFullSizeClustersImage.py
#   2. QuPath (https://qupath.github.io/), to draw ROIs and export with Export_Annotations_GeoJSON.groovy
#   3. ReadAnnotationsToAnnData.py
#
# --------------------------------------------------------------------------------

print('Loading data files and selecting annotated regions...')

## Muscle data

datadir = os.path.join(FILEPATHBASE, 'data')
files = os.listdir(datadir)
h5ad_files = [file for file in files if file.endswith('.h5ad')]

sample = 1
samples_all = dict()
for filename in h5ad_files:
    filepath = os.path.join(datadir, filename)
    print("Loading data from " + filepath)
    adata = ad.read_h5ad(filepath)
    library_id = 'IBM_Muscle_bx'
    adata.uns["spatial"] = dict()
    adata.uns["spatial"][library_id] = dict()
    samples_all["Muscle_"+str(sample)] = adata
    print('Loaded')
    sample = sample + 1
print('Done loading all data sets')

### Done Loading 

# --------------------------------------------------------------------------------
# Defining Sample-Level Characteristics
# --------------------------------------------------------------------------------

SampleKey = {"Muscle_1": "DM IBM 1",
             "Muscle_2": "PTV IBM 2"
             }

# --------------------------------------------------------------------------------
# Concatenation of all datasets into one
# --------------------------------------------------------------------------------

print('Concatenating...')

samples_all = ad.concat(samples_all, label="dataset", uns_merge="first", join='outer')
samples_all.obs_names_make_unique()

# Clean up the NAs in manual annotation columns in adata.obs, which should be boolean
cs = samples_all.obs.select_dtypes(include='object').columns
samples_all.obs[cs] = samples_all.obs[cs].astype('boolean').fillna(False)

samples_all.raw = samples_all
samples_all.layers["counts"] = samples_all.X.copy()   

library_id = 'Curio_Seeker_v1.1_AllRetinas'
samples_all.uns["spatial"] = dict()
samples_all.uns["spatial"][library_id] = dict()


## Sample-level info
samples_all.uns["SampleKey"] = SampleKey



# --------------------------------------------------------------------------------
# Save concatenated data
# --------------------------------------------------------------------------------
Path(os.path.join(FILEPATHBASE, 'calc')).mkdir(parents=True, exist_ok=True)
out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all.h5ad')
samples_all.write(out_filename)

print('Saved concatenated data to: ' + out_filename)
print('Done!')
