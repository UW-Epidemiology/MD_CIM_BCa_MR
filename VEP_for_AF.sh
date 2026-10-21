#!/bin/bash
#SBATCH --job-name=VEP_MD_MR
#SBATCH --time=4-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G

# Change to the directory where files are located
cd /projects/lindstroem/MD_MR_Inflammatory_Biomarkers/instrument_data/

# List of root names for input files
root_names=("IL6" "IL8" "TNFR2" "VEGF" "CRP" "CRP_and_BMI")  

# Loop over each root name
for input_name in "${root_names[@]}"; do
    # Set input and output file paths
    input="./${input_name}_snps.txt"
    output="./${input_name}_VEP.tsv"
    
    # Print paths for debugging
    echo "Processing: $input -> $output"
    
    # Run VEP
    /projects/lindstroem/programs/ensembl-vep/vep \
        -i ${input} \
        --format "id" \
        --assembly GRCh38 \
        --cache --dir /projects/lindstroem/programs/ensembl-vep/caches/ \
        --pick \
        --gene_phenotype \
        --symbol \
        --af \
        --af_1kg \
        --af_gnomade \
        --af_gnomadg \
        --no_stats \
        --fork 4 \
        --buffer_size 25000 \
        --tab \
        -o ${output} 

done