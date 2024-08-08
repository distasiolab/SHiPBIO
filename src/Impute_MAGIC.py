# Impute_MAGIC.py
# Part of SHiPBIO
# Marcello DiStasio
# July 2024
#
# Uses MAGIC [GitHub](https://github.com/KrishnaswamyLab/MAGIC) [Cell Paper](https://www.cell.com/cell/abstract/S0092-8674(18)30724-4)
# 
##################################################

import warnings
warnings.filterwarnings('ignore')

import anndata as ad

import os
from pathlib import Path
import argparse

import re

import magic


parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
parser.add_argument('-i', '--inputfile', type=str, help='Path to input AnnData *.h5ad file')
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

# --------------------------------------------------------------------------------
# Imputation
# --------------------------------------------------------------------------------
magic_operator = magic.MAGIC()
samples_all.layers['counts_magic'] = magic_operator.fit_transform(samples_all.layers['counts_scvi'])

# --------------------------------------------------------------------------------
# Save
# --------------------------------------------------------------------------------
if SAVEDATA:
    # Save
    if args.output is None:
        out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all_integrated_imputed_magic.h5ad') 
    else:
        out_filename = args.output

    samples_all.write_h5ad(out_filename)
    print('Saved ' + out_filename)

