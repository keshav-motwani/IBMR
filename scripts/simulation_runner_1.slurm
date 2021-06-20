#!/bin/bash

#SBATCH --array=1-2000
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=IBMR
#SBATCH --mem-per-cpu=2gb
#SBATCH -t 4320
#SBATCH --output=logs_1/IBMR_%a.log

ml R/4

Rscript scripts/simulation_runner_1.R
