#!/bin/bash

#SBATCH --constraint=milan
#SBATCH --array=2167,2175,2199,2207,2215,2217,2218,2221,2222,2239,2247,2271,2287,2290,2293,2319,2343,2362,2365,2383,2391,2415,2433,2434,2437,2463,2487,2506,2509
#SBATCH --mail-user=keshav.motwani@ufl.edu
#SBATCH --mail-type=END
#SBATCH --account=amolstad
#SBATCH --qos=amolstad-b
#SBATCH --job-name=IBMR
#SBATCH --mem-per-cpu=128gb
#SBATCH -t 4320
#SBATCH --output=logs/application_Seurat_%a.log

export OMP_NUM_THREADS=1

ml R/4.0

PER_TASK=1

START_NUM=$(( ($SLURM_ARRAY_TASK_ID - 1) * $PER_TASK + 1 ))
END_NUM=$(( $SLURM_ARRAY_TASK_ID * $PER_TASK ))

for (( PARAMETER_ID=$START_NUM; PARAMETER_ID<=END_NUM; PARAMETER_ID++ )); do
  Rscript --vanilla scripts/application_runner.R $PARAMETER_ID Seurat
done
