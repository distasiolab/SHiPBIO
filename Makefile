
## Path to base directory for the project; should contain directories 'src', 'data', and 'calc'
##
##
# Mac OS X 
#BASEDIR := "/Users/mmd47/Google Drive/My Drive/DiStasio Lab/DiStasio Lab Share/02 Analysis/Muscle_IBM/SHiPBIO"
BASEDIR := "/home/mdistasio/YaleGoogleDrive/DiStasio Lab/DiStasio Lab Share/02 Analysis/Muscle_IBM/SHiPBIO"
##

SOURCE = ./src/
CALC = ./calc/


PREPROCESS_RESULT := $(CALC)samples_all.h5ad
INTEGRATE_RESULT := $(CALC)samples_all_integrated_harmony_unfiltered.h5ad
IMPUTATION_RESULT := $(CALC)samples_all_integrated_imputed.h5ad
CLUSTER_RESULT := $(CALC)samples_all_integrated_imputed_cellcharter_clustered.h5ad


.dummy: preprocess integrate impute cluster

all: cluster

preprocess: $(PREPROCESS_RESULT)
	@echo "Preprocessing completed."
	@echo $(PREPROCESS_RESULT) " exists."

integrate: $(INTEGRATE_RESULT)
	@echo "Integration completed."
	@echo $(INTEGRATE_RESULT) " exists."

impute: $(IMPUTATION_RESULT)
	@echo "Imputation completed."
	@echo $(IMPUTATION_RESULT) " exists."

cluster: $(CLUSTER_RESULT)
	@echo "Clustering completed."
	@echo $(CLUSTER_RESULT) " exists."



$(PREPROCESS_RESULT): 
	@echo "Preprocessing..."
	python $(SOURCE)Preprocess.py -b $(BASEDIR)

$(INTEGRATE_RESULT): $(PREPROCESS_RESULT)
	@echo "Batch Integration..."
	python $(SOURCE)IntegrateHarmony.py -b $(BASEDIR)


$(IMPUTATION_RESULT): $(INTEGRATE_RESULT)
	@echo "Imputation (MAGIC)..."
	python $(SOURCE)Impute_MAGIC.py -b $(BASEDIR)

$(CLUSTER_RESULT): $(IMPUTATION_RESULT)
	@echo "Clustering..."
	python $(SOURCE)Cluster_CellCharter.py -b $(BASEDIR)


clean:
	rm -rf calc
	rm -rf img/out
