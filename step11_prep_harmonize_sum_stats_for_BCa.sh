#!/bin/bash
#SBATCH --job-name=prep_sum_stats
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G

Rscript step11_prep_harmonize_sum_stats_for_BCa.R