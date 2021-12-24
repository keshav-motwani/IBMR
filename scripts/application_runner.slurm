#!/bin/bash

#SBATCH --array=304
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad
#SBATCH --job-name=IBMR
#SBATCH --mem-per-cpu=16gb
#SBATCH -t 4320
#SBATCH --output=logs/application_%a.log

export OMP_NUM_THREADS=1

ml R/4.0

Rscript scripts/application_runner.R
