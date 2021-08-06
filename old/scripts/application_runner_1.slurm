#!/bin/bash

#SBATCH --array=1-8
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=IBMR_application_1
#SBATCH --mem-per-cpu=6gb
#SBATCH -t 4320
#SBATCH --output=logs/IBMR_application_1_%a.log
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20

export OMP_NUM_THREADS=1

ml R/4

Rscript scripts/application_runner_1.R
