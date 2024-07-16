# SHiPBIO
## Spatial High-Plex Biological Analysis

Tools for 

* Creating ROIs / Annotations manually for annData objects in QuPath
* Clustering in Molecular Space
* Analyzing Patterns in anatomically defined compartments

```
SHiPBIO
|	README.md
|_______calc
|	|_______Used for intermediate results (will not be tracked by git, per .gitignore)
|_______data
|	|_______Put your data here (will not be tracked by git, per .gitignore)
|_______doc
|_______src
|	|_______Analysis source code
|_______img
|	|_______Output images go here
|_______local
|	|_______To keep your own local config files (e.g. sample_worksheet.csv, cluster_labels.json)
|_______Scripts
|	|_______Some helpful shell scripts
|_______Docs
	|_______InfoAboutProject.txt
```

## Installation of prerequisites
1. Create a conda environment
```
conda create -n cellcharter-env -c conda-forge python=3.10 mamba
conda activate cellcharter-env
```
2. Install cellcharter (which should bring along pytorch)
```
pip install cellcharter
```
3. Install [scvi-tools](https://scvi-tools.org/)
```
pip install --upgrade --upgrade-strategy only-if-needed scvi-tools
```
4. Install [harmonypy](https://github.com/slowkow/harmonypy) and [pymde](https://pymde.org/)
```
pip install --upgrade --upgrade-strategy only-if-needed harmonypy
pip install --upgrade --upgrade-strategy only-if-needed pymde
pip install --upgrade --upgrade-strategy only-if-needed geojson
```
5. Install [MAGIC](https://github.com/KrishnaswamyLab/MAGIC)
```
conda deactivate
conda create -n magic-env -c conda-forge python=3.10 mamba
conda activate magic-env
pip install magic-impute
```

## Set up
* If needed, you may want to add subdirectories for raw data (like *.fastq) and processed data (like matrix files, etc)
```
mkdir data/raw
mkdir data/processed
mkdir ReferenceFiles
```

# Usage
1. For each sample with data in anndata format (`*h5ad`):
   1. Run `MakeFullSizeClustersImage.py -f /path/to/*.h5ad` which will output a rough (KNN in UMAP space) clustered image (into the same directory, unless you specify another location with the `-o` flag).
   2. Open each image in [QuPath](https://qupath.github.io/). Draw annotations around your ROIs and set the classes of each anootation to your desired annotation labels. Then, in QuPath, run the script `src/QuPath/Export_Annotations_GeoJSON.groovy`. This will create a *.json file in that same directory
   3. Run `ReadAnnotationsToAnnData.py -f /path/to/*.h5ad -c /path/to/*.h5ad_annotations.json`.  This will create a file `/path/to/*_annotated.h5ad`
2. Make a copy of `sample_worksheet.csv` and edit so that each line contains information about the sample (including path to `/path/to/*_annotated.h5ad`)
3. Edit `Makefile`:
   1. Set paths for your system:
      1. `BASEDIR` is the full path to the SHiPBIO directory itself
      2. `CONDA_ENV_CELLCHARTER` is the full path to the conda environment you installe cellcharter in (see `conda activate cellcharter-env; conda info`)
      3. `SAMPLE_WORKSHEET` is the path to the `sample_worksheet.csv` file
      4. `MARKER_GENE_FILE` is the path to the marker gene `*.json` file
   2. `SINGLE_CELL_DATA` is the path to an AnnData `*.h5ad` file containing single cell sequencing data to integrate with the spatial data (using GIMVI)
   3. Edit `N_CLUSTERS` and `N_HOPS` to set the parameters for cellcharter clustering (number of GMM clusters, and hops along neighborhood graph, respectively)
4. Run analysis with `make` or `make cluster_individual`. This will create a new file with all data concatenated in `calc/`, as well as figures in `img/out/` for each sample
5. Inspect the spatial cluster images and marker gene matrix images generated in `img/out/` to determine the appropriate cluster label for each cluster numeric ID for each sample.  Edit `example_cluster_labels.json` to reflect the chosen labels for each sample. Then run `Cluster_CellCharter_RelabelClusters.py -b /path/to/SHiPBIO/ -c /path/to/*_cluster_labels.json`.  This will create the relabeled output figures in `img/out`.