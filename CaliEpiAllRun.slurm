#!/bin/bash

#SBATCH -J CaliEpiScript

#SBATCH -t 100:00:00
#SBATCH -N 1 # node
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=4096M


#SBATCH --mail-user=rbirger@princeton.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL



Rscript CaliSimulations.R  $SLURM_ARRAY_TASK_ID



# submit using sbatch
# scontrol show job 1373263
# squeue -u rbirger
# scancel 1423435