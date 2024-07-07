# Preprocess.py
# Part of SHiPBIO for processing spatial omics data
# Marcello DiStasio
# Jul 2024


import anndata as ad
import squidpy as sq
import scanpy as sc

import os
import csv

from pathlib import Path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-b', '--basepath', type=str, help='Path to base directory for the project; should contain directories \'data\' and \'calc\'')
parser.add_argument('-s', '--worksheet', type=str, help='Path to *.csv file with sample file paths and information')
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


known_columns =['sampleID', 'sampleName', 'AnatomicLocation', 'filename']
data = read_csv_into_dict(args.worksheet, known_columns)
    
import pprint
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
print("Done loading data")

# --------------------------------------------------------------------------------
# Defining Sample-Level Characteristics
# --------------------------------------------------------------------------------

r_all =     dict([[sample['sampleID'],sample['data']] for sample in data])
SampleKey = dict([[sample['sampleID'],sample['sampleName']] for sample in data])
AnaLoc =    dict([[sample['sampleID'],sample['AnatomicLocation']] for sample in data])

# --------------------------------------------------------------------------------
# Concatenation of all datasets into one
# --------------------------------------------------------------------------------

print('Concatenating...')

retinas_all = ad.concat(r_all, label="dataset", uns_merge="first", join='outer')

# Clean up the NAs in manual annotation columns in adata.obs, which should be boolean
cs = retinas_all.obs.select_dtypes(include='object').columns
retinas_all.obs[cs] = retinas_all.obs[cs].astype('boolean').fillna(False)

retinas_all.raw = retinas_all
retinas_all.layers["counts"] = retinas_all.X.copy()   

library_id = 'Curio_Seeker_v1.1_AllRetinas'
retinas_all.uns["spatial"] = dict()
retinas_all.uns["spatial"][library_id] = dict()


## Sample-level info
retinas_all.uns["SampleKey"] = SampleKey
retinas_all.uns["Anatomic_Location"] = AnatLoc



# --------------------------------------------------------------------------------
# Save concatenated data
# --------------------------------------------------------------------------------
out_filename = os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate', 'calc', 'retinas_all.h5ad')
retinas_all.write(out_filename)

print('Saved concatenated data to: ' + out_filename)
print('Done!')
