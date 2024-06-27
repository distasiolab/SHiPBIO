#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --job-name=Muscle
#SBATCH --cpus-per-task=12
#SBATCH --gpus 1 
#SBATCH --mem=900G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment
module load miniconda

conda init

conda activate /gpfs/gibbs/project/distasio/mmd47/envs/cellcharter-env

#python src/MakePlots.py -b /gpfs/gibbs/project/distasio/mmd47/muscle_curio/SHiPBIO
python src/Cluster_CellCharter.py -b /gpfs/gibbs/project/distasio/mmd47/muscle_curio/SHiPBIO


