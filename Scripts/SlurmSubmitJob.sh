#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --job-name=Retina
#SBATCH --cpus-per-task=12
#SBATCH --gpus 1
#SBATCH --mem=512G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu
#SBATCH -o tmp/slurm-%j.out

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment
module load miniconda

make -f /home/mmd47/project/retina_curio/SHiPBIO/local/Makefile gwas



