#!/bin/bash
#SBATCH --job-name=MR
#SBATCH --output=mr_%j.log
#SBATCH --error=mr_%j.err
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

Rscript step7_MR_OpenGWAS_project_CIMs.R