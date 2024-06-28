.ONESHELL:
SHELL = /bin/bash
CONDA_ACTIVATE = source $(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate


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
BASEDIR := "/home/mmd47/project/retina_curio/SHiPBIO"
CONDA_ENV_CELLCHARTER := "/gpfs/gibbs/project/distasio/mmd47/envs/cellcharter-env"
CONDA_ENV_PHATE := "/gpfs/gibbs/project/distasio/mmd47/envs/phate"

.ONESHELL:
SHELL = /bin/bash
CONDA_ACTIVATE = source $$(conda info --base)/etc/profile.d/conda.sh ; conda activate ; conda activate
##
########################################################################################################################

SOURCE = ./src/
CALC = ./calc/

SINGLE_CELL_DATA := $(CALC)retina_sn_combined.h5ad

PREPROCESS_RESULT := $(CALC)samples_all.h5ad
BATCH_INTEGRATE_RESULT := $(CALC)samples_all_integrated_harmony_unfiltered.h5ad
SINGLECELL_RESULT := $(CALC)samples_all_integrated_snRNAseq_imputed.h5ad
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
	conda init
	conda deactivate
	$(CONDA_ACTIVATE) $(CONDA_ENV_CELLCHARTER)
	python $(SOURCE)Preprocess.py -b $(BASEDIR)
	conda deactivate

$(BATCH_INTEGRATE_RESULT): $(PREPROCESS_RESULT)
	@echo "Batch integration..."
	bash --init-file <(echo ". "$(HOME)/.bashrc"; conda activate $(CONDA_ENV_CELLCHARTER); export LD_LIBRARY_PATH=$(CONDA_ENV_CELLCHARTER)/lib/; python $(SOURCE)/IntegrateHarmony.py -b $(BASEDIR)")

$(SINGLECELL_INTEGRATE_RESULT): $(BATCH_INTEGRATE_RESULT) $(SINGLE_CELL_DATA)
	@echo "Single cell integration..."
	bash --init-file <(echo ". "$(HOME)/.bashrc"; conda activate $(CONDA_ENV_CELLCHARTER); export LD_LIBRARY_PATH=$(CONDA_ENV_CELLCHARTER)/lib/; python $(SOURCE)/Integrate_scRNAseq_GIMVI.py -b $(BASEDIR)")

$(IMPUTATION_RESULT): $(SINGLECELL_INTEGRATE_RESULT)
	@echo "Imputation (MAGIC)..."
	bash --init-file <(echo ". "$(HOME)/.bashrc"; conda activate $(CONDA_ENV_PHATE); export LD_LIBRARY_PATH=$(CONDA_ENV_PHATE)/lib/; python $(SOURCE)/Impute_MAGIC.py -b $(BASEDIR)")

$(CLUSTER_RESULT): $(IMPUTATION_RESULT)
	@echo "Clustering..."
	bash --init-file <(echo ". "$(HOME)/.bashrc"; conda activate $(CONDA_ENV_CELLCHARTER); export LD_LIBRARY_PATH=$(CONDA_ENV_CELLCHARTER)/lib/; python $(SOURCE)/Cluster_CellCharter.py -b $(BASEDIR)")


clean:
	rm -rf calc
	rm -rf img/out
