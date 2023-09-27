#!/bin/bash
#SBATCH --partition=bigmem
#SBATCH --job-name=MyJob
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=400G
#SBATCH --time=8:00:00
#SBATCH --mail-type=ALL
#SBATCH -o tmp/slurm-%j.out


# NOTE: --partition bigmem allows up to 1050G of mem

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment

# (uncomment if needed)
# module load miniconda
# conda activate MyCondaEnvironment


# NOTE: now issue commands to run your code

#python -u ../Python/MyPythonScript.py

