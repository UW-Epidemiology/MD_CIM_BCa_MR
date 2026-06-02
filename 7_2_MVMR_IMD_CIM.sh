#!/bin/bash
#SBATCH --job-name=MVMR
#SBATCH --output=MVMR_%j.log
#SBATCH --error=MVMR_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

Rscript 7_2_MVMR_IMD_CIMs.R