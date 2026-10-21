#!/bin/bash
#SBATCH --job-name=MR
#SBATCH --output=mr_%j.log
#SBATCH --error=mr_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16

Rscript step8_MVMR_OpenGWAS_project_CIMs.R