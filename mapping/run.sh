#!/bin/bash
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH -p short
#SBATCH --mem=64G
#SBATCH --job-name=run


module load gcc/6.2.0
module load python/3.7.4
module load perl/5.30.0
module load R/3.6.1

snakemake -j 8 -p



