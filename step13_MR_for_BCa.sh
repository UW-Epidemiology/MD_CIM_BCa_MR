#!/bin/bash
#SBATCH --job-name=MR
#SBATCH --output=mr_%j.log
#SBATCH --error=mr_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

Rscript step13_MR_for_BCa.R