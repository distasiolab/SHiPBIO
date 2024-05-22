import warnings
warnings.filterwarnings('ignore')

import anndata as ad
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

import magic


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
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
# Load dataset
# --------------------------------------------------------------------------------
filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_harmony_unfiltered.h5ad')
print("Loading Data from: " + filename + ' ...')
samples_all = ad.read_h5ad(filename)


# --------------------------------------------------------------------------------
# Imputation
# --------------------------------------------------------------------------------
magic_operator = magic.MAGIC()
samples_all.layers["X_magic"] = magic_operator.fit_transform(samples_all.X)

# --------------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------------
if SAVEDATA:
    # Save
    filename_out = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed.h5ad')
    samples_all.write_h5ad(filename_out)
    print('Saved ' + filename_out)

