SOURCE = ./src/
CALC = ./calc/

all : $(CALC)retinas_all_unfiltered_snRNAseq_imputed_GIMVI.h5ad

$(CALC)retinas_all.h5ad :
	python $(SOURCE)Preprocess.py

$(CALC)retinas_all_integrated_harmony_unfiltered.h5ad : $(CALC)retinas_all.h5ad
	python $(SOURCE)IntegrateHarmony.py

$(CALC)retinas_all_unfiltered_snRNAseq_imputed_GIMVI.h5ad : $(CALC)retinas_all_integrated_harmony_unfiltered.h5ad
	python $(SOURCE)Cluster_CellCharter.py

$(CALC)retinas_all_unfiltered_snRNAseq_imputed_GIMVI_ClustersLabeled_scDRS_GWAS_scores.h5ad' : $(CALC)retinas_all_unfiltered_snRNAseq_imputed_GIMVI_ClustersLabeled.h5ad
	python $(SOURCE)GWAS_Analysis.py
