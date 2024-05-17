
## Path to base directory for the project; should contain directories 'src', 'data', and 'calc'
##
##
BASEDIR := "/Users/mmd47/Google Drive/My Drive/DiStasio Lab/DiStasio Lab Share/02 Analysis/Muscle_IBM/SHiPBIO"
##

SOURCE = ./src/
CALC = ./calc/


PREPROCESS_RESULT := $(CALC)IBM_muscle_all.h5ad
INTEGRATE_RESULT := integrate_result.txt
CLUSTER_RESULT := cluster_result.txt


.dummy: preprocess integrate cluster

all: preprocess

preprocess: $(PREPROCESS_RESULT)
	@echo "Preprocessing completed."
	@echo $(PREPROCESS_RESULT) "exists."

integrate: $(INTEGRATE_RESULT)
	@echo "Integration completed."

cluster: $(CLUSTER_RESULT)
	@echo "Clustering completed."


$(PREPROCESS_RESULT):
	mkdir -p $(CALC)
	@echo "Preprocessing..."
	python $(SOURCE)Preprocess.py -b $(BASEDIR)

$(INTEGRATE_RESULT):
	@echo "Batch Integration..."
	#python $(SOURCE)IntegrateHarmony.py -b $(BASEDIR)

$(CLUSTER_RESULT):
	@echo "Clustering..."
	#python $(SOURCE)Cluster_CellCharter.py -b $(BASEDIR)

