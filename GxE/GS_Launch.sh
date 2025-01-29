#!/bin/bash

#SBATCH --export=NONE               # do not export current env to the job
#SBATCH --job-name=GblupparGXE.fixed.brr.fixed.LOO      # job name
#SBATCH --time=01-00:00:00             # max job run time dd-hh:mm:ss
#SBATCH --ntasks-per-node=1         # tasks (commands) per compute node
#SBATCH --cpus-per-task=48          # CPUs (threads) per command
#SBATCH --mem=360GB                  # total memory per node
#SBATCH --output=stdout.%x.%j       # save stdout to file
#SBATCH --error=stderr.%x.%j        # save stderr to file
#SBATCH --mail-user=
#SBATCH --mail-type=all


module purge
ml GCC/11.2.0  OpenMPI/4.1.1 R_tamu/4.2.0 
Rscript Gxe.R
