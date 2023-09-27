# ProjectTemplate
An empty project template to start with. Contains the directory structure:


```
ProjectName_version
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

If needed, you may want to add subdirectories for raw data (like *.fastq) and processed data (like matrix files, etc)

```
mkdir data/raw
mkdir data/processed
mkdir ReferenceFiles
```
