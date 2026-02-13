# Setup
library(data.table)
library(R.utils)
library(tidyverse)
setDTthreads(percent = 100)

# Generate list of filepaths and corresponding proteins

dir_path <- "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/sum_stats/"

filepaths <- list.files(
  path = dir_path,
  full.names = TRUE
)

filepaths <- filepaths[!grepl("\\.txt$", filepaths)]

proteins <- sub("\\.tbl\\.gz$", "", basename(filepaths))


# Load MD GWAS outcome data (no adjustment for BMI) 
## Allele1 is the effect allele

## PMD
PMD_outcome_data = fread("/projects/lindstroem/MD_GWAS/FilesforCameron/Total_NoBMI/PMD_NoBMI_Stderr_MetaAnalysis_Cleaned_Final.txt")

## DA 
DA_outcome_data = fread("/projects/lindstroem/MD_GWAS/FilesforCameron/Total_NoBMI/DA_NoBMI_Stderr_MetaAnalysis_Cleaned_Final.txt.gz")

## NDA
NDA_outcome_data = fread("/projects/lindstroem/MD_GWAS/FilesforCameron/Total_NoBMI/Rerun_20210325/NDA_NoBMI_Stderr_MetaAnalysis_Cleaned_Final.txt")


# Get CIM GWAS exposure data 
## Allele1 is the effect allele
for (i in 1:length(proteins)) {
  # Load GWAS summary statistics and p-value (p < 5e-8) filter 
  sum_stats = fread(filepaths[i]) %>%
    filter(`log(P)` < log10(5E-8)) %>%
    select("Chromosome","Position","Allele2","Allele1","Effect","StdErr","log(P)")
  
  sum_stats = sum_stats %>%
    mutate(P_exposure = 10^`log(P)`,
           exposure = proteins[i]) %>%
    select(-`log(P)`)
  
  colnames(sum_stats) = c("Chr", "Pos", "A1", "A2",
                             "Effect_exposure", "StdErr_exposure",
                             "P_exposure", "exposure")
  # Create sorted allele representations
  sum_stats <- sum_stats %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  # Harmonize with outcome data
  
  ## PMD
  PMD_outcome = PMD_outcome_data %>%
    filter(position %in% sum_stats$Pos) %>%
    mutate(outcome = "PMD")
  
  PMD_outcome = PMD_outcome %>%
    select(c("chr","position","rsid","P-value","Effect",
             "StdErr","Allele2","Allele1","outcome"))
  
  colnames(PMD_outcome) = c("Chr", "Pos", "RSID", "P_outcome", 
                          "Effect_outcome", "StdErr_outcome",
                          "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  PMD_outcome <- PMD_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  PMD_harmonized <- merge(
    sum_stats,
    PMD_outcome,
    by = c("Chr", "Pos", "A_min", "A_max")
  ) %>%
    # Determine whether alleles are flipped relative to exposure
    mutate(
      flipped = (A1.x == A2.y & A2.x == A1.y),
      
      # Harmonize outcome effect
      Effect_outcome = ifelse(
        flipped,
        -Effect_outcome,
        Effect_outcome
      )
    ) %>%
    select(-c("A_min", "A_max"))
  
  ### Cleanup 
  PMD_harmonized = PMD_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  ## DA
  DA_outcome = DA_outcome_data %>%
    filter(position %in% sum_stats$Pos) %>%
    mutate(outcome = "DA")
  
  DA_outcome = DA_outcome %>%
    select(c("chr","position","rsid","P-value","Effect",
             "StdErr","Allele2","Allele1","outcome"))
  
  colnames(DA_outcome) = c("Chr", "Pos", "RSID", "P_outcome", 
                           "Effect_outcome", "StdErr_outcome",
                           "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  DA_outcome <- DA_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  DA_harmonized <- merge(
    sum_stats,
    DA_outcome,
    by = c("Chr", "Pos", "A_min", "A_max")
  ) %>%
    # Determine whether alleles are flipped relative to exposure
    mutate(
      flipped = (A1.x == A2.y & A2.x == A1.y),
      
      # Harmonize outcome effect
      Effect_outcome = ifelse(
        flipped,
        -Effect_outcome,
        Effect_outcome
      )
    ) %>%
    select(-c("A_min", "A_max"))
  
  ### Cleanup 
  DA_harmonized = DA_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  ## NDA
  NDA_outcome = NDA_outcome_data %>%
    filter(position %in% sum_stats$Pos) %>%
    mutate(outcome = "NDA")
  
  NDA_outcome = NDA_outcome %>%
    select(c("chr","position","rsid","P-value","Effect",
             "StdErr","Allele2","Allele1","outcome"))
  
  colnames(NDA_outcome) = c("Chr", "Pos", "RSID", "P_outcome", 
                           "Effect_outcome", "StdErr_outcome",
                           "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  NDA_outcome <- NDA_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  NDA_harmonized <- merge(
    sum_stats,
    NDA_outcome,
    by = c("Chr", "Pos", "A_min", "A_max")
  ) %>%
    # Determine whether alleles are flipped relative to exposure
    mutate(
      flipped = (A1.x == A2.y & A2.x == A1.y),
      
      # Harmonize outcome effect
      Effect_outcome = ifelse(
        flipped,
        -Effect_outcome,
        Effect_outcome
      )
    ) %>%
    select(-c("A_min", "A_max"))
  
  ### Cleanup 
  NDA_harmonized = NDA_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  # Join harmonized data
  harmonized = rbind.data.frame(PMD_harmonized,DA_harmonized,NDA_harmonized)
  
  # Save data
  write.table(
    harmonized,
    file = paste0("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/harmonized_p_filtered/",proteins[i],"_harmonized_data.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  ### Directory for harmonized SNP files for PLINK clumping
  out_dir <- "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/harmonized_p_filtered/pre_clumping_SNPs/"

  ### Create SNP sets by chromosome
  snps_by_chr <- harmonized %>%
    distinct(Chr, RSID) %>%
    group_by(Chr) %>%
    summarise(
      SNPs = list(unique(RSID)),
      .groups = "drop"
    )
  
  snps_by_chr_list <- setNames(snps_by_chr$SNPs, snps_by_chr$Chr)
  
  ### Write one file per protein per chromosome for PLINK clumping
  invisible(
    lapply(names(snps_by_chr_list), function(chr) {
      
      chr_numeric <- as.numeric(chr)
      
      # Subset harmonized data to this chromosome and these SNPs
      df_chr <- harmonized %>%
        filter(Chr == chr_numeric & RSID %in% snps_by_chr_list[[chr]]) %>%
        select(RSID, Chr, P_exposure) %>%
        distinct()
      
      # Build filename
      out_file <- paste0(out_dir, proteins[i], "_chr", chr_numeric, "_for_clumping.txt")
      
      # Write tab-delimited file with header for PLINK
      write.table(
        df_chr,
        file = out_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = c("SNP","CHR","P")
      )
      
    })
  )
  
}