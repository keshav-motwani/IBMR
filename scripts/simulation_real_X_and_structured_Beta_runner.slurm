#!/bin/bash

#SBATCH --array=1-900%800
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=IBMR
#SBATCH --mem-per-cpu=6gb
#SBATCH -t 1320
#SBATCH --output=logs/simulation_real_X_and_structured_Beta_%a.log

ml R/4

Rscript scripts/simulation_real_X_and_structured_Beta_runner.R
