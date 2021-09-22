#!/bin/bash

#SBATCH --array=1-2520
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=IBMR
#SBATCH --mem-per-cpu=14gb
#SBATCH -t 4320
#SBATCH --output=logs/application_%a.log

export OMP_NUM_THREADS=1

ml R/4

Rscript scripts/application_runner.R
