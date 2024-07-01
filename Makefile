########################################################################################################################
## Environment variables for the project; BASEDIR should contain directories 'src', 'data', and 'calc'
##
##
### Mac OS X 
#BASEDIR := "/Users/mmd47/Google Drive/My Drive/DiStasio Lab/DiStasio Lab Share/02 Analysis/Muscle_IBM/SHiPBIO"
########################################################################################################################
### Fedora
#BASEDIR := "/home/mdistasio/YaleGoogleDrive/DiStasio Lab/DiStasio Lab Share/02 Analysis/Muscle_IBM/SHiPBIO"
##
########################################################################################################################
### mccleary.ycrc.yale.edu
SHELL=/bin/bash -i
BASEDIR := "/home/mmd47/project/retina_curio/SHiPBIO"
CONDA_ENV_CELLCHARTER := "/gpfs/gibbs/project/distasio/mmd47/envs/cellcharter-env"
CONDA_ENV_PHATE := "/gpfs/gibbs/project/distasio/mmd47/envs/phate"

##
########################################################################################################################

SOURCE = ./src/
CALC = ./calc/

SINGLE_CELL_DATA := $(CALC)retina_sn_combined.h5ad

PREPROCESS_RESULT := $(CALC)samples_all.h5ad
BATCH_INTEGRATE_RESULT := $(CALC)samples_all_integrated_harmony_unfiltered.h5ad
SINGLECELL_INTEGRATE_RESULT := $(CALC)samples_all_integrated_snRNAseq_imputed.h5ad
IMPUTATION_RESULT := $(CALC)samples_all_integrated_imputed.h5ad
CLUSTER_RESULT := $(CALC)samples_all_integrated_imputed_cellcharter_clustered.h5ad


.dummy: preprocess integrate impute cluster

all: cluster

preprocess: $(PREPROCESS_RESULT)
	@echo "Preprocessing completed."
	@echo $(PREPROCESS_RESULT) " exists."

batch_integrate: $(BATCH_INTEGRATE_RESULT)
	@echo "Batch integration completed."
	@echo $(BATCH_INTEGRATE_RESULT) " exists."

singlecell_integrate: $(SINGLECELL_INTEGRATE_RESULT)
	@echo "Single cell integration completed."
	@echo $(SINGLECELL_INTEGRATE_RESULT) " exists."

impute: $(IMPUTATION_RESULT)
	@echo "Imputation completed."
	@echo $(IMPUTATION_RESULT) " exists."

cluster: $(CLUSTER_RESULT)
	@echo "Clustering completed."
	@echo $(CLUSTER_RESULT) " exists."




$(PREPROCESS_RESULT): 
	@echo "Preprocessing..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}/lib/; python ${SOURCE}/Preprocess.py -b ${BASEDIR}' | bash -i

$(BATCH_INTEGRATE_RESULT): $(PREPROCESS_RESULT)
	@echo "Batch integration..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}/lib/; python ${SOURCE}/IntegrateHarmony.py -b ${BASEDIR}' | bash -i

$(SINGLECELL_INTEGRATE_RESULT): $(BATCH_INTEGRATE_RESULT) $(SINGLE_CELL_DATA)
	@echo "Single cell integration..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}/lib/; python ${SOURCE}/Integrate_scRNAseq_GIMVI.py -b ${BASEDIR}' | bash -i 

$(IMPUTATION_RESULT): $(SINGLECELL_INTEGRATE_RESULT)
	@echo "Imputation (MAGIC)..."
	echo 'conda activate ${CONDA_ENV_PHATE}; export LD_LIBRARY_PATH=${CONDA_ENV_PHATE}/lib/; python ${SOURCE}/Impute_MAGIC.py -b ${BASEDIR}' | bash -i

$(CLUSTER_RESULT): $(IMPUTATION_RESULT)
	@echo "Clustering..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}/lib/; python ${SOURCE}/Cluster_CellCharter.py -b ${BASEDIR}' | bash -i


clean:
	rm -rf calc
	rm -rf img/out
