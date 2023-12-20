# Preprocess.py
# Part of anndata_manualAnnotate for processing spatial omics data
# Marcello DiStasio
# Oct 2023


import anndata as ad
import squidpy as sq
import scanpy as sc

import os

#FILEPATHBASE = '/Users/mmd47/Library/CloudStorage/GoogleDrive-mmd47@yale.edu/My Drive/DiStasio Lab/DiStasio Lab Share/'
FILEPATHBASE = '/home/mdistasio/YaleGoogleDrive/DiStasio Lab/DiStasio Lab Share/'

SAVEFIGS = True
if SAVEFIGS:
    IMGDIR = os.path.join(FILEPATHBASE,'02 Analysis/annData_ManualAnnotate/img/')


# --------------------------------------------------------------------------------
# Prior to this loading, each *_anndata.h5ad should be annotated with
#
#   1. MakeFullSizeClustersImage.py
#   2. QuPath (https://qupath.github.io/), to draw ROIs and export with Export_Annotations_GeoJSON.groovy
#   3. ReadAnnotationsToAnnData.py
#
# --------------------------------------------------------------------------------

print('Loading data files and selecting annotated regions...')

# File 1    
filename = os.path.join(FILEPATHBASE,'03 Data','Retina_SlideSeq_Curio','A22_3781_AMD_SlideSeq_001','A0052_029_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A22_3781'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina1 = adata[adata.obs['Retina_1']]
retina2 = adata[adata.obs['Retina_2']]

# File 2
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-914_OS_Ctrl_SlideSeq_001/A0052_030_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-914'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina3 = adata[adata.obs['Retina']]

# File 3
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1277_AMD_SlideSeq_001/OUTPUT/A23-1277-OS_macula/A23-1277-OS_macula_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1277'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina4 = adata[adata.obs['Retina']]

# File 4
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1279-OS_SlideSeq_001/23-1279-OS_OUTPUT/OUTPUT/A23-1279-OS_macula/A23-1279-OS_macula_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1279'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina5 = adata[adata.obs['Retina']]

# File 5
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1422_Ctrl_Macula_SlideSeq_001/OUTPUT/A23-1422_-_macula/A23-1422_-_macula_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1422'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina6 = adata[adata.obs['Retina']]

# File 6
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1425_Ctrl_Macula_SlideSeq_001/OUTPUT/A23-1425_-_macula1/A23-1425_-_macula1_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1425-M1'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina7 = adata[adata.obs['Retina']]


# File 7
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1152_OD_AMD_SlideSeq_002/OUTPUT/A23-1152_-_Mac_AMD/A23-1152_-_Mac_AMD_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1152_OD'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina8 = adata[adata.obs['Retina']]

r_all = {"R1": retina1, "R2": retina2, "R3": retina3, "R4": retina4, "R5": retina5, "R6": retina6, "R7": retina7, "R8": retina8}
SampleKey = {"R1": "A22_3781", "R2": "A22_3781", "R3": "A23-914", "R4": "A23-1277", "R5": "A23-1279", "R6": "A23-1422", "R7": "A23-1425-M1", "R8": "A23-1152_OD"}

# --------------------------------------------------------------------------------
# Done Loading 
# --------------------------------------------------------------------------------



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

retinas_all.uns["SampleKey"] = SampleKey

# --------------------------------------------------------------------------------
# Save concatenated data
# --------------------------------------------------------------------------------
out_filename = os.path.join(FILEPATHBASE,'02 Analysis', 'annData_ManualAnnotate', 'calc', 'retinas_all.h5ad')
retinas_all.write(out_filename)

print('Saved concatenated data to: ' + out_filename)
print('Done!')
