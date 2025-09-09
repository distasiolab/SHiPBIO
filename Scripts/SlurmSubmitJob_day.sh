#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=rsid_map
#SBATCH --cpus-per-task=1
#SBATCH --array=1-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu
#SBATCH -o tmp/slurm-%j.out

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment
module purge
module load miniconda
module load tabix

# --------------------------
# Configuration
# --------------------------
VCF="/home/mmd47/ycga_work/References/RSID_VCF_GeneInfoMapping/GCF_000001405.40.gz"
MAP_DIR="chrom_maps"
MERGED_MAP="rsid_to_gene_full.tsv"
mkdir -p $MAP_DIR
mkdir -p logs

# Map SLURM_ARRAY_TASK_ID to chromosome
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)
CHR=${CHROMS[$SLURM_ARRAY_TASK_ID-1]}

OUT_MAP="$MAP_DIR/rsid_to_gene_chr${CHR}.tsv"

echo "Processing chromosome $CHR -> $OUT_MAP"

# Extract chromosome-specific VCF with tabix
TMP_VCF="$MAP_DIR/dbsnp_chr${CHR}.vcf.gz"
tabix -h $VCF $CHR > $TMP_VCF

# Run the Python mapping script
python rsid_to_gene.py --vcf $TMP_VCF --map $OUT_MAP

# Remove temporary VCF to save space
rm $TMP_VCF
