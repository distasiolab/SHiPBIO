#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --job-name=bigmemRetina
#SBATCH --cpus-per-task=12
#SBATCH --mem=2400G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu
#SBATCH -o tmp/slurm-%j.out

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment
module load miniconda

make cluster



