# Preprocess.py
# Part of SHiPBIO for processing spatial omics data
# Marcello DiStasio
# Jul 2024


import anndata as ad
import squidpy as sq
import scanpy as sc

import os
import csv
import pprint

from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
parser.add_argument('-s', '--worksheet', type=str, help='Path to *.csv file with sample file paths and information')
parser.add_argument('-o', '--output', type=str, help='Path to output *.h5ad file to create')
args = parser.parse_args()

FILEPATHBASE = args.basepath


# --------------------------------------------------------------------------------
# Prior to this loading, each *_anndata.h5ad should be annotated with
#
#   1. MakeFullSizeClustersImage.py
#   2. QuPath (https://qupath.github.io/), to draw ROIs and export with Export_Annotations_GeoJSON.groovy
#   3. ReadAnnotationsToAnnData.py
#
# --------------------------------------------------------------------------------
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


known_columns =['sampleID', 'sampleName', 'Condition', 'AnatomicLocation', 'GoodSample', 'filename']
data = read_csv_into_dict(args.worksheet, known_columns)

pp = pprint.PrettyPrinter(indent=4, width=200)
print("Sample worksheet info from " + args.worksheet)
pp.pprint(data)

print('Loading data files and selecting annotated regions...')

cnt = 1
for sample in data:
    print("Loading:" + sample['filename'] + "...")
    adata = ad.read_h5ad(sample['filename'])
    adata.uns["spatial"] = dict()
    adata.uns["spatial"][sample['sampleID']] = dict()
    sample['data'] = adata[adata.obs[sample['Annotations']].any(axis=1)]
    print("Done loading sample " + str(cnt) + "/" + str(len(data)) + ".")
    cnt=cnt+1
print("Done loading data")

# --------------------------------------------------------------------------------
# Defining Sample-Level Characteristics
# --------------------------------------------------------------------------------

r_all =     dict([[sample['sampleID'],sample['data']] for sample in data])
SampleKey = dict([[sample['sampleID'],sample['sampleName']] for sample in data])
AnatLoc =    dict([[sample['sampleID'],sample['AnatomicLocation']] for sample in data])
ConditionKey =  dict([[sample['sampleID'],sample['Condition']] for sample in data])
GoodSampleKey = dict([[sample['sampleID'],bool(sample['GoodSample'])] for sample in data])

# --------------------------------------------------------------------------------
# Concatenation of all datasets into one
# --------------------------------------------------------------------------------

print('Concatenating...')

samples_all = ad.concat(r_all, label="dataset", uns_merge="first", join='outer')

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
samples_all.uns["Anatomic_Location"] = AnatLoc
samples_all.uns["ConditionKey"] = ConditionKey
samples_all.uns["GoodSampleKey"] = GoodSampleKey


# --------------------------------------------------------------------------------
# Subset only good samples
# --------------------------------------------------------------------------------
SUBSET_GOOD = True
if SUBSET_GOOD:
    samples_all = samples_all[samples_all.obs['dataset'].map(samples_all.uns['GoodSampleKey']),:]

# --------------------------------------------------------------------------------
# Quality Control Tests
# --------------------------------------------------------------------------------
sc.pp.calculate_qc_metrics(samples_all, 
                           expr_type='counts', 
                           var_type='genes',
                           percent_top=(5, 10, 20, 50), 
                           layer='counts',
                           inplace=True)

# Gene thresholds
F_genes_counts = sc.pp.filter_genes(samples_all,
                       min_counts = 3,
                       inplace=False)
F_genes_cells = sc.pp.filter_genes(samples_all,
                       min_cells = 3,
                       inplace=False)

# SPOT thresholds
F_cells_counts = sc.pp.filter_cells(samples_all,
                       min_counts = 40,
                       inplace=False)
F_cells_genes = sc.pp.filter_cells(samples_all,
                       min_genes = 3,
                       inplace=False)

samples_all_filtered = samples_all[F_cells_counts[0] & F_cells_genes[0], F_genes_counts[0] & F_genes_cells[0]].copy()


# --------------------------------------------------------------------------------
# Save concatenated data
# --------------------------------------------------------------------------------
if not os.path.exists(os.path.join(FILEPATHBASE, 'calc')):
    os.makedirs(os.path.join(FILEPATHBASE, 'calc'))

if args.output is None:
    out_filename = os.path.join(FILEPATHBASE, 'calc', 'samples_all.h5ad')
else:
    out_filename = args.output
    
samples_all.write(out_filename)

print('Saved concatenated data to: ' + out_filename)
print('Done!')
