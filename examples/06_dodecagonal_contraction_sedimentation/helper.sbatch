#!/bin/bash

##SBATCH --account=lilly-nan6
#SBATCH --account=lilly-roa3
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=14-00:00:0
#SBATCH --mem-per-cpu=4G

module purge
module load gcc/14.1.0
module load python/3.11.9

python3 force_top_only.py log.pimpleLYJHFDIBFoam --top-n 100 > force_summary.txt