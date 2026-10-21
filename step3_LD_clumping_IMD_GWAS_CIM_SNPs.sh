#!/bin/bash
#SBATCH --job-name=ld_clump
#SBATCH --output=ld_clump_%j.log
#SBATCH --error=ld_clump_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Directory with your summary stats files
SUMSTATS_DIR="/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/harmonized_p_filtered/pre_clumping_SNPs/"

# Directory with chromosome-specific PLINK reference files
REF_DIR="/projects/lindstroem/1000Genomes/hg38/PLINK/eur_subset/plink1_files/"

# Output directory
OUT_DIR="/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/harmonized_fully_filtered/post_clumping_SNPs/"
mkdir -p "${OUT_DIR}"

# Loop over all .txt files
for file in "${SUMSTATS_DIR}"/*.txt; do
    # Get filename without extension
    filename=$(basename "${file}" .txt)
    
    # Extract substring before first underscore
    base="${filename%%_*}"
    
    echo "Processing ${file} (base=${base})..."
    
    # Read the file to get unique chromosomes
    chromosomes=$(awk 'NR>1 {print $2}' "${file}" | sort -u)
    
    # Loop over each chromosome
    for chr in ${chromosomes}; do
        echo "Processing chromosome ${chr}..."
        
        # Count number of SNPs for this chromosome (excluding header)
        num_snps=$(awk -v c="${chr}" 'NR>1 && $2==c' "${file}" | wc -l)
        
        final_out="${OUT_DIR}/${base}_chr${chr}.clumped"
        
        if [ "$num_snps" -le 1 ]; then
            echo "Only $num_snps SNP(s) on chromosome ${chr}, skipping PLINK clumping."
            # Copy the original rows (header + SNP row) as-is
            awk -v c="${chr}" 'NR==1 || $2==c' "${file}" > "${final_out}"
        else
            # Use the chromosome-specific rows directly, no deduplication
            tmp_chr_file=$(mktemp)
            awk -v c="${chr}" 'NR==1 || $2==c' "${file}" > "${tmp_chr_file}"
            
            echo "Clumping and MAF filtering ${num_snps} SNPs on chromosome ${chr}..."
            plink2 \
                --bfile "${REF_DIR}eur_hg38" \
                --maf 0.05 \
                --clump "${tmp_chr_file}" \
                --clump-snp-field SNP \
                --clump-field P \
                --clump-p1 1e-8 \
                --clump-p2 1e-5 \
                --clump-r2 0.001 \
                --clump-kb 10000 \
                --allow-no-sex \
                --out "${OUT_DIR}/${base}_chr${chr}"
            
            rm -f "${tmp_chr_file}"
        fi
    done
done

echo "All clumping jobs finished."