#!/bin/bash

#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad
#SBATCH --job-name=IBMR
#SBATCH --mem-per-cpu=24gb
#SBATCH -t 4320
#SBATCH --output=logs/simulation_real_X_and_Beta_prepare_data.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

export OMP_NUM_THREADS=1

ml R/4

Rscript scripts/simulation_real_X_and_Beta_prepare_data.R
