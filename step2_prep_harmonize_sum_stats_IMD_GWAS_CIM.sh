#!/bin/bash
#SBATCH --job-name=prep_sum_stats
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

Rscript step2_prep_harmonize_sum_stats_IMD_GWAS_CIMs.R