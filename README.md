# annData_ManualAnnotate

Tools for creating ROIs / Annotations manually for annData objects containing spatial transcriptomic data in QuPath. Also has functions for batch integration and plotting.

```
annData_ManualAnnotate
|	README.md
|	.gitignore
|
|_______data
|	|_______Put your data here (will not be tracked by git, per .gitignore)
|	
|_______src
|	|_______Analysis source code
|	
|_______img
|	|_______Output images go here
|
|_______Scripts
|	|_______Some helpful shell scripts
|
|
|_______Docs
	|_______InfoAboutProject.txt

```

Code and Workbooks are in `src` directory.

Current workflow is:

1) Move individual sample *.h5ad files containing cropped & annotated (e.g. 'Retina' ROIs) data into a known location
2) Preprocess.py (basically just loads and concatenates all datasets, saves resulting object in `calc` directory)
3) IntegrateHarmony.py (batch correction)
4) GeneratePlots.py
