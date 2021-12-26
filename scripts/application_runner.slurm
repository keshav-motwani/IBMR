#!/bin/bash

#SBATCH --constraint=milan
#SBATCH --array=1-2520%400
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=IBMR
#SBATCH --mem-per-cpu=16gb
#SBATCH -t 4320
#SBATCH --output=logs/application_%a.log

export OMP_NUM_THREADS=1

ml R/4.0

PER_TASK=4

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

for (( PARAMETER_ID=$START_NUM; PARAMETER_ID<=END_NUM; PARAMETER_ID++ )); do
  Rscript --vanilla scripts/application_runner.R $PARAMETER_ID
done
