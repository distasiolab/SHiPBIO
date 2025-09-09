#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=GeneLookup
#SBATCH --cpus-per-task=4
#SBATCH --mem=264G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu
#SBATCH -o tmp/slurm-%j.out

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment
module purge
module load miniconda

conda activate cellcharter-env

python rsid_to_gene.py --map  rsid_to_gene.full.tsv --annotate *.sumstats

conda deactivate



