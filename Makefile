########################################################################################################################
## Environment variables for the project; BASEDIR should contain directories 'src', 'data', and 'calc'
##
##
### Mac OS X
##
#BASEDIR := "/Users/mmd47/Google Drive/My Drive/DiStasio Lab/DiStasio Lab Share/02 Analysis/Muscle_IBM/SHiPBIO"
##
##
##
##
########################################################################################################################
### Fedora
##
BASEDIR := "/home/mdistasio/YaleGoogleDrive/DiStasio Lab/DiStasio Lab Share/02 Analysis/Muscle_IBM/SHiPBIO"
SHELL=/bin/bash -i
BASEDIR := "/home/mdistasio/Workspace/SHiPBIO/"
CONDA_ENV_CELLCHARTER := "/home/mdistasio/miniconda3/envs/cellcharter-env/"
CONDA_ENV_PHATE := "/home/mdistasio/miniconda3/envs/phate-env/"
SAMPLE_WORKSHEET := ${BASEDIR}/sample_worksheet.csv
GWAS_WORKSHEET := ${BASEDIR}/gwas_worksheet.csv
##
##
##
##
########################################################################################################################
### mccleary.ycrc.yale.edu
##
#SHELL=/bin/bash -i
#BASEDIR := "/home/mmd47/project/retina_curio/SHiPBIO"
#CONDA_ENV_CELLCHARTER := "/gpfs/gibbs/project/distasio/mmd47/envs/cellcharter-env/"
#CONDA_ENV_PHATE := "/gpfs/gibbs/project/distasio/mmd47/envs/phate/"
#SAMPLE_WORKSHEET := ${BASEDIR}/sample_worksheet.csv
#GWAS_WORKSHEET := ${BASEDIR}/gwas_worksheet.csv
##
##
##
##
########################################################################################################################

SOURCE = ./src/
CALC = ./calc/
DATA = ./data/

N_CLUSTERS ?=23
N_HOPS ?=3

SINGLE_CELL_DATA := $(DATA)retina_sn_combined.h5ad
MARKER_GENE_FILE := $(DATA)retinal_celltype_gates.json
CLUSTER_LABELS_FILE := $(DATA)samples_all_integrated_imputed_cellcharter_3hops_23_clusters_individual_clusterlabels.json

PREPROCESS_RESULT := $(CALC)samples_all.h5ad
BATCH_INTEGRATE_RESULT := $(CALC)samples_all_integrated_harmony_unfiltered.h5ad
SINGLECELL_INTEGRATE_RESULT := $(CALC)samples_all_integrated_snRNAseq_imputed.h5ad
CLUSTER_RESULT := $(CALC)samples_all_integrated_imputed_cellcharter_clustered.h5ad
CLUSTER_INDIVIDUAL_RESULT := $(CALC)samples_all_integrated_imputed_cellcharter_clustered_individual_$(N_HOPS)_hops_$(N_CLUSTERS)_clusters.h5ad
CLUSTER_INDIVIDUAL_RELABELED_RESULT := $(CALC)samples_all_integrated_imputed_cellcharter_clustered_individual_$(N_HOPS)_hops_$(N_CLUSTERS)_clusters_relabeled.h5ad
IMPUTATION_RESULT := $(CALC)samples_all_integrated_imputed_cellcharter_clustered_individual_$(N_HOPS)_hops_$(N_CLUSTERS)_clusters_relabeled_magic.h5ad
GWAS_RESULT := $(CALC)samples_all_integrated_imputed_cellcharter_clustered_individual_$(N_HOPS)_hops_$(N_CLUSTERS)_clusters_relabeled_magic_GWAS_SCDRS.h5ad

#.dummy: preprocess integrate cluster_individual

all: cluster_individual

preprocess: $(PREPROCESS_RESULT)
	@echo "Preprocessing completed."
	@echo $(PREPROCESS_RESULT) " exists."

batch_integrate: $(BATCH_INTEGRATE_RESULT)
	@echo "Batch integration completed."
	@echo $(BATCH_INTEGRATE_RESULT) " exists."

singlecell_integrate: $(SINGLECELL_INTEGRATE_RESULT)
	@echo "Single cell integration completed."
	@echo $(SINGLECELL_INTEGRATE_RESULT) " exists."

cluster: $(CLUSTER_RESULT)
	@echo "Clustering completed."
	@echo $(CLUSTER_RESULT) " exists."

cluster_individual: $(CLUSTER_INDIVIDUAL_RESULT)
	@echo "Clustering completed."
	@echo $(CLUSTER_INDIVIDUAL_RESULT) " exists."

relabel: $(CLUSTER_INDIVIDUAL_RELABELED_RESULT)
	@echo "Cluster relabeling completed."
	@echo $(CLUSTER_INDIVIDUAL_RELABELED_RESULT) " exists."

impute: $(IMPUTATION_RESULT)
	@echo "Imputation completed."
	@echo $(IMPUTATION_RESULT) " exists."

gwas: $(GWAS_RESULT)
	@echo "Here!!"
	@echo "GWAS analysis completed."
	@echo $(GWAS_RESULT) " exists."


$(PREPROCESS_RESULT):
	@echo "$(PREPROCESS_RESULT) does not exist. Building."
	@echo "Preprocessing..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}lib/; python ${SOURCE}Preprocess.py -b ${BASEDIR} -o ${PREPROCESS_RESULT}' -s ${SAMPLE_WORKSHEET} | bash -i

$(BATCH_INTEGRATE_RESULT): $(PREPROCESS_RESULT)
	@echo "Batch integration..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}lib/; python ${SOURCE}IntegrateHarmony.py -b ${BASEDIR} -o ${BATCH_INTEGRATE_RESULT}' | bash -i

$(SINGLECELL_INTEGRATE_RESULT): $(BATCH_INTEGRATE_RESULT) $(SINGLE_CELL_DATA)
	@echo "Single cell integration..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}lib/; python ${SOURCE}Integrate_scRNAseq_GIMVI.py -b ${BASEDIR} -c ${SINGLE_CELL_DATA} -o ${SINGLECELL_INTEGRATE_RESULT}' | bash -i 


$(CLUSTER_RESULT): $(SINGLECELL_INTEGRATE_RESULT)
	@echo "Clustering..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}lib/; export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:32; python ${SOURCE}Cluster_CellCharter.py -b ${BASEDIR} -n ${N_CLUSTERS} -d ${N_HOPS} -m ${MARKER_GENE_FILE} -o ${CLUSTER_RESULT}' | bash -i

$(CLUSTER_INDIVIDUAL_RESULT): $(SINGLECELL_INTEGRATE_RESULT)
	@echo "Clustering..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}lib/; export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:32; python ${SOURCE}Cluster_CellCharter_IndividualSamples.py -b ${BASEDIR} -i ${SINGLECELL_INTEGRATE_RESULT} -n ${N_CLUSTERS} -d ${N_HOPS} -m ${MARKER_GENE_FILE} -o ${CLUSTER_INDIVIDUAL_RESULT}' | bash -i

$(CLUSTER_INDIVIDUAL_RELABELED_RESULT): $(CLUSTER_INDIVIDUAL_RESULT) $(CLUSTER_LABELS_FILE)
	@echo "Relabeling clusters"
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}lib/; export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:32; python ${SOURCE}Cluster_CellCharter_RelabelClusters.py -b ${BASEDIR} -i ${CLUSTER_INDIVIDUAL_RESULT} -c ${CLUSTER_LABELS_FILE} -o ${CLUSTER_INDIVIDUAL_RELABELED_RESULT}' | bash -i

$(IMPUTATION_RESULT): $(CLUSTER_INDIVIDUAL_RELABELED_RESULT)
	@echo "Imputation (MAGIC)..."
	echo 'conda activate ${CONDA_ENV_PHATE}; export LD_LIBRARY_PATH=${CONDA_ENV_PHATE}lib/; python ${SOURCE}Impute_MAGIC.py -b ${BASEDIR} -i ${CLUSTER_INDIVIDUAL_RELABELED_RESULT} -o ${IMPUTATION_RESULT}' | bash -i

$(GWAS_RESULT): $(CLUSTER_INDIVIDUAL_RELABELED_RESULT)
	@echo "Here!"
	@echo "GWAS analysis (SCDRS)..."
	echo 'conda activate ${CONDA_ENV_CELLCHARTER}; export LD_LIBRARY_PATH=${CONDA_ENV_CELLCHARTER}lib/; python ${SOURCE}GWAS_Analysis.py -b ${BASEDIR} -i $(CLUSTER_INDIVIDUAL_RELABELED_RESULT) -g ${GWAS_WORKSHEET} -o ${GWAS_RESULT}' | bash -i

clean:
	rm -rf calc
	rm -rf img/out
