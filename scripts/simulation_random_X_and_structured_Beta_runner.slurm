#!/bin/bash

#SBATCH --array=1-1400
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=IBMR
#SBATCH --mem-per-cpu=6gb
#SBATCH -t 4320
#SBATCH --output=logs/simulation_random_X_and_structured_Beta_%a.log

export OMP_NUM_THREADS=1

ml R/4.0

Rscript scripts/simulation_random_X_and_structured_Beta_runner.R
