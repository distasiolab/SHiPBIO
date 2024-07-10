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

* Install [CellCharter](https://github.com/CSOgroup/cellcharter)
1. Create a conda environment
`conda create -n cellcharter-env -c conda-forge python=3.10 mamba`
`conda activate cellcharter-env`
2. Install cellcharter (which will bring along pytorch)
`pip install cellcharter`
3. Install [scvi-tools](https://scvi-tools.org/)
`pip install --upgrade --upgrade-strategy only-if-needed scvi-tools`
4. Install [harmonypy](https://github.com/slowkow/harmonypy)
`pip install harmonypy`

# Set up
* If needed, you may want to add subdirectories for raw data (like *.fastq) and processed data (like matrix files, etc)
```
mkdir data/raw
mkdir data/processed
mkdir ReferenceFiles
```
# Usage
1. For each sample with data in anndata format (`*h5ad`):
  1. Run `MakeFullSizeClustersImage.py -f /path/to/*.h5ad` which will output a rough (KNN in UMAP space) clustered image (into the same directory, unless you specify another location with the `-o` flag).
  2. Open each image in [QuPath](https://qupath.github.io/). Draw annotations around your ROIs and set the classes of each anootation to your desired annotation labels. Then, in QuPath, run the script `Export_Annotations_GeoJSON.groovy`. This will create a *.json file in that same directory
  3. Run `ReadAnnotationsToAnnData.py -f /path/to/*.h5ad -c /path/to/*.h5ad_annotations.json`.  This will create a file `/path/to/*_annotated.h5ad`
2. Make a copy of `sample_worksheet.csv` and edit so that each line contains information about the sample (including path to `/path/to/*_annotated.h5ad`)
3. Run analysis with `make` or `make cluster_individual`
