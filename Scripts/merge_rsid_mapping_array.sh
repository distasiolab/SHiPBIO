#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=merge
#SBATCH --cpus-per-task=1
#SBATCH --array=1-698
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu
#SBATCH -o tmp/slurm-%j.out

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment
module purge

cat output/rsid_to_gene.*.tsv | sort -u > rsid_to_gene.full.tsv
