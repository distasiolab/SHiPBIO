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

## MACULAR RETINA

# File 1    
filename = os.path.join(FILEPATHBASE,'03 Data','Retina_SlideSeq_Curio','A22_3781_AMD_SlideSeq_001','A0052_029_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A22_3781'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina1 = adata[(adata.obs['Retina_1'] == True) | (adata.obs['Retina_2'] == True)]

# File 2
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-914_OS_Ctrl_SlideSeq_001/A0052_030_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-914'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina2 = adata[adata.obs['Retina']]

# File 3
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1277_AMD_SlideSeq_001/OUTPUT/A23-1277-OS_macula/A23-1277-OS_macula_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1277'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina3 = adata[adata.obs['Retina']]

# File 4
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1279-OS_SlideSeq_001/23-1279-OS_OUTPUT/OUTPUT/A23-1279-OS_macula/A23-1279-OS_macula_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1279'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina4 = adata[adata.obs['Retina']]

# File 5
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1422_Ctrl_Macula_SlideSeq_001/OUTPUT/A23-1422_-_macula/A23-1422_-_macula_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1422'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina5 = adata[adata.obs['Retina']]

# File 6
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1425_Ctrl_Macula_SlideSeq_001/OUTPUT/A23-1425_-_macula1/A23-1425_-_macula1_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1425-M1'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina6 = adata[adata.obs['Retina']]


# File 7
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/A23-1152_OD_AMD_SlideSeq_002/OUTPUT/A23-1152_-_Mac_AMD/A23-1152_-_Mac_AMD_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1152_OD'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina7 = adata[adata.obs['Retina']]


## PERIPHERAL RETINA

# File 8
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/Peripheral_Retina/1279_Peripheral_SlideSeq_001/OUTPUT/1279_P/1279_P_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1279_P'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina8 = adata[adata.obs['Retina']]


# File 9
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/Peripheral_Retina/1341-1_Peripheral_SlideSeq_001/OUTPUT/1341-1_P/1341-1_P_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1341_P'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina9 = adata[adata.obs['Retina']]


# File 10
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/Peripheral_Retina/1422_P_Peripheral_SlideSeq_001/OUTPUT/1422_P/1422_P_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1422_P'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina10 = adata[adata.obs['Retina']]

# File 11
filename = os.path.join(FILEPATHBASE,'03 Data/Retina_SlideSeq_Curio/Peripheral_Retina/1425-3_Peripheral_SlideSeq_001/OUTPUT/1425-3_P/1425-3_P_anndata_annotated.h5ad')
adata = ad.read_h5ad(filename)

library_id = 'A23-1425-3_P'
adata.uns["spatial"] = dict()
adata.uns["spatial"][library_id] = dict()
retina11 = adata[adata.obs['Retina']]



### Done Loading 


# --------------------------------------------------------------------------------
# Defining Sample-Level Characteristics
# --------------------------------------------------------------------------------

r_all = {"R1": retina1,
         "R2": retina2,
         "R3": retina3,
         "R4": retina4,
         "R5": retina5,
         "R6": retina6,
         "R7": retina7,
         "R8": retina8,
         "R9": retina9,
         "R10": retina10,
         "R11": retina11}

SampleKey = {"R1": "A22_3781",
             "R2": "A23-914",
             "R3": "A23-1277",
             "R4": "A23-1279",
             "R5": "A23-1422",
             "R6": "A23-1425-M1",
             "R7": "A23-1152_OD",
             "R8": "A23-1279_P",
             "R9": "A23-1341_P",
             "R10": "A23-1422_P",
             "R11": "A23-1425-3_P"}

AnatLoc = {"R1": "macula_retina",
           "R2": "macula_retina",
           "R3": "macula_retina",
           "R4": "macula_retina",
           "R5": "macula_retina",
           "R6": "macula_retina",
           "R7": "macula_retina",
           "R8": "peripheral_retina",
           "R9": "peripheral_retina",
           "R10": "peripheral_retina",
           "R11": "peripheral_retina"}






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
