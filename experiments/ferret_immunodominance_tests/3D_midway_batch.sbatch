#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --job-name=2D-fits
#SBATCH --account=pi-cobey
#SBATCH --partition=cobey
#SBATCH --ntasks-per-node=1  # number of tasks
#SBATCH --cpus-per-task=4    # number of threads per task
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=kgostic

# LOAD MODULES
module load R

# DO COMPUTE WORK
Rscript infer_maps_3D.R