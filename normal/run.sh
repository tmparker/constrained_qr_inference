#!/bin/bash
#SBATCH --account=def-tmparker
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00-00:20 # time (DD-HH:MM)

module load gcc r
Rscript sim.R 

