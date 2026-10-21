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


# Load BCa GWAS outcome data
## Effect allele is labled

## Overall
Overall_BCa_outcome_data = fread("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/BCa_summary_statistics/Overall_BCa.tsv")

## HER2_pos
HER2_pos_BCa_outcome_data = fread("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/BCa_summary_statistics/HER2_pos_BCa.tsv")

## Luminal A 
Luminal_A_BCa_outcome_data = fread("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/BCa_summary_statistics/Luminal_A_BCa.tsv")

## Luminal B 
Luminal_B_BCa_outcome_data = fread("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/BCa_summary_statistics/Luminal_B_BCa.tsv")

## ER_neg
ER_neg_BCa_outcome_data = fread("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/BCa_summary_statistics/ER_neg_BCa.tsv")

## TNBC
TNBCa_outcome_data = fread("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/BCa_summary_statistics/TNBCa.tsv")

# 1000G .bim
bim_dir <- "/projects/lindstroem/1000Genomes/hg19_auto_all/EUR/"

# Find all .bim files
bim_files <- list.files(
  path = bim_dir,
  pattern = "\\.bim$",
  full.names = TRUE
)

# Read and merge all BIM files
bim_data <- do.call(
  rbind,
  lapply(bim_files, function(file) {
    fread(file)
  })
)

colnames(bim_data) = c("Chr","RSID","cM","Pos","A1","A2")
bim <- bim_data %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )


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
      A1 = A1 %>% toupper(),
      A2 = A2 %>% toupper(),
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  # Harmonize with outcome data
  
  ## Overall_BCa
  Overall_BCa_outcome = Overall_BCa_outcome_data %>%
    filter(base_pair_location %in% sum_stats$Pos) %>%
    mutate(outcome = "Overall_BCa")
  
  Overall_BCa_outcome = Overall_BCa_outcome %>%
    select(c("chromosome","base_pair_location","p_value","beta",
             "standard_error","effect_allele","other_allele","outcome"))
  
  colnames(Overall_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                          "Effect_outcome", "StdErr_outcome",
                          "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  Overall_BCa_outcome <- Overall_BCa_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  Overall_BCa_outcome <- merge(
    Overall_BCa_outcome,
    bim,
    by = c("Chr", "Pos", "A_min", "A_max"),
    all.x = TRUE,
    all.y = FALSE
  )
  
  Overall_BCa_harmonized <- merge(
    sum_stats,
    Overall_BCa_outcome,
    by = c("Chr", "Pos", "A_min", "A_max"),
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
  Overall_BCa_harmonized = Overall_BCa_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  ## HER2_pos_BCa
  HER2_pos_BCa_outcome = HER2_pos_BCa_outcome_data %>%
    filter(base_pair_location %in% sum_stats$Pos) %>%
    mutate(outcome = "HER2_pos_BCa",
           p_value =  10^(-neg_log_10_p_value))
  
  HER2_pos_BCa_outcome = HER2_pos_BCa_outcome %>%
    select(c("chromosome","base_pair_location","p_value","beta",
             "standard_error","effect_allele","other_allele","outcome"))
  
  colnames(HER2_pos_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  HER2_pos_BCa_outcome <- HER2_pos_BCa_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  HER2_pos_BCa_outcome <- merge(
    HER2_pos_BCa_outcome,
    bim,
    by = c("Chr", "Pos", "A_min", "A_max"),
    all.x = TRUE,
    all.y = FALSE
  )
  
  HER2_pos_BCa_harmonized <- merge(
    sum_stats,
    HER2_pos_BCa_outcome,
    by = c("Chr", "Pos", "A_min", "A_max"),
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
  HER2_pos_BCa_harmonized = HER2_pos_BCa_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  ## Luminal_A_BCa
  Luminal_A_BCa_outcome = Luminal_A_BCa_outcome_data %>%
    filter(base_pair_location %in% sum_stats$Pos) %>%
    mutate(outcome = "Luminal_A_BCa",
           p_value = 10^(-neg_log_10_p_value))
  
  Luminal_A_BCa_outcome = Luminal_A_BCa_outcome %>%
    select(c("chromosome","base_pair_location","p_value","beta",
             "standard_error","effect_allele","other_allele","outcome"))
  
  colnames(Luminal_A_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  Luminal_A_BCa_outcome <- Luminal_A_BCa_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  Luminal_A_BCa_outcome <- merge(
    Luminal_A_BCa_outcome,
    bim,
    by = c("Chr", "Pos", "A_min", "A_max"),
    all.x = TRUE,
    all.y = FALSE
  )
  
  Luminal_A_BCa_harmonized <- merge(
    sum_stats,
    Luminal_A_BCa_outcome,
    by = c("Chr", "Pos", "A_min", "A_max"),
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
  Luminal_A_BCa_harmonized = Luminal_A_BCa_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  
  ## Luminal_B_BCa
  Luminal_B_BCa_outcome = Luminal_B_BCa_outcome_data %>%
    filter(base_pair_location %in% sum_stats$Pos) %>%
    mutate(outcome = "Luminal_B_BCa",
           p_value = 10^(-neg_log_10_p_value))
  
  Luminal_B_BCa_outcome = Luminal_B_BCa_outcome %>%
    select(c("chromosome","base_pair_location","p_value","beta",
             "standard_error","effect_allele","other_allele","outcome"))
  
  colnames(Luminal_B_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  Luminal_B_BCa_outcome <- Luminal_B_BCa_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  Luminal_B_BCa_outcome <- merge(
    Luminal_B_BCa_outcome,
    bim,
    by = c("Chr", "Pos", "A_min", "A_max"),
    all.x = TRUE,
    all.y = FALSE
  )
  
  Luminal_B_BCa_harmonized <- merge(
    sum_stats,
    Luminal_B_BCa_outcome,
    by = c("Chr", "Pos", "A_min", "A_max"),
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
  Luminal_B_BCa_harmonized = Luminal_B_BCa_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  ## ER_neg_BCa
  ER_neg_BCa_outcome = ER_neg_BCa_outcome_data %>%
    filter(base_pair_location %in% sum_stats$Pos) %>%
    mutate(outcome = "ER_neg_BCa")

  ER_neg_BCa_outcome = ER_neg_BCa_outcome %>%
    select(c("chromosome","base_pair_location","p_value","beta",
             "standard_error","effect_allele","other_allele","outcome"))
  
  colnames(ER_neg_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                      "Effect_outcome", "StdErr_outcome",
                                      "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  ER_neg_BCa_outcome <- ER_neg_BCa_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  ER_neg_BCa_outcome <- merge(
    ER_neg_BCa_outcome,
    bim,
    by = c("Chr", "Pos", "A_min", "A_max"),
    all.x = TRUE,
    all.y = FALSE
  )
  
  ER_neg_BCa_harmonized <- merge(
    sum_stats,
    ER_neg_BCa_outcome,
    by = c("Chr", "Pos", "A_min", "A_max"),
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
  ER_neg_BCa_harmonized = ER_neg_BCa_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  ## Triple negative
  TNBCa_outcome = TNBCa_outcome_data %>%
    filter(base_pair_location %in% sum_stats$Pos) %>%
    mutate(outcome = "TNBCa",
           p_value = 10^(-neg_log_10_p_value))
  
  TNBCa_outcome = TNBCa_outcome %>%
    select(c("chromosome","base_pair_location","p_value","beta",
             "standard_error","effect_allele","other_allele","outcome"))
  
  colnames(TNBCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")
  
  ### Create sorted allele representations
  TNBCa_outcome <- TNBCa_outcome %>%
    mutate(
      A_min = pmin(A1, A2),
      A_max = pmax(A1, A2)
    )
  
  ### Merge and harmonize
  TNBCa_outcome <- merge(
    TNBCa_outcome,
    bim,
    by = c("Chr", "Pos", "A_min", "A_max"),
    all.x = TRUE,
    all.y = FALSE
  )
  
  TNBCa_harmonized <- merge(
    sum_stats,
    TNBCa_outcome,
    by = c("Chr", "Pos", "A_min", "A_max"),
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
  TNBCa_harmonized = TNBCa_harmonized %>%
    mutate(A1 = toupper(A1.x),
           A2 = toupper(A2.x)) %>%
    select(-c(A1.x,A1.y,A2.x,A2.y))
  
  
  # Join harmonized data
  harmonized = rbind.data.frame(Overall_BCa_harmonized,
                                HER2_pos_BCa_harmonized,
                                Luminal_A_BCa_harmonized,
                                Luminal_B_BCa_harmonized,
                                ER_neg_BCa_harmonized,
                                TNBCa_harmonized)
  # Save data
  write.table(
    harmonized,
    file = paste0("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/BCa/harmonized_p_filtered/",proteins[i],"_harmonized_data.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  ### Directory for harmonized SNP files for PLINK clumping
  out_dir <- "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/BCa/harmonized_p_filtered/pre_clumping_SNPs/"

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

# Load GWAS summary statistics for CRP 

## Define CRP instrument including BMI-associated SNPs 

CRP_sum_stats = fread("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/raw_summary_statistics/CRP.tsv.gz")

read_VEP <- function(x) {
  VEP_text <- readLines(x)
  VEP_text <- VEP_text[!startsWith(VEP_text, "##")]
  df <- fread(text = VEP_text, sep = "\t", na.strings = "-", data.table = TRUE)
  if (startsWith(names(df)[1], "#")) {
    setnames(df, old = names(df)[1], new = sub("^#", "", names(df)[1]))
  }
  
  return(df)
}

CRP_VEP_data = read_VEP("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/CPR_VEP.tsv") %>%
  filter(EUR_AF > 0.05) %>%
  select(Uploaded_variation) %>%
  pull()

SNPs = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_CRP_noBMIadj_harmonized_data.rds") %>%
  mutate(exposure = "CRP", outcome = "PMD") %>%
  filter(SNP %in% CRP_VEP_data, mr_keep == TRUE) %>%
  select(SNP) %>% pull()

CRP_sum_stats = CRP_sum_stats %>%
  filter(variant_id %in% SNPs)


### Define CRP instrument with no BMI-associated SNPs 
PMD_CRP_BMI = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_CRP_and_BMI_harmonized_data.rds")
BMI_associated_SNPs = PMD_CRP_BMI$exposure_pval[,2] < 0.05
BMI_CRP_SNPs = dimnames(PMD_CRP_BMI$exposure_pval)[[1]]
BMI_associated_SNPs = BMI_CRP_SNPs[BMI_associated_SNPs]

# Other BMI-associated SNPs identified by checking LD trait #
other_BMI_associated_SNPs <- c("rs10498240","rs10498240","rs10797119","rs10871777",
                               "rs112875651","rs1222209", "rs1260326", "rs12716972",
                               "rs12929503","rs13066686","rs13294945","rs13416992",
                               "rs1905505","rs2726033","rs28429148","rs288183",
                               "rs34092660","rs2280406","rs197422","rs35695283",
                               "rs3768321","rs3935032","rs429358","rs4655802",
                               "rs630372","rs7794705","rs78428995","rs994596",
                               "rs73137144","rs7103411","rs6992315","rs12251016",
                               "rs10760691","rs10969334","rs11078597","rs11118625",
                               "rs11242113","rs12118890","rs13216424","rs141783576",
                               "rs141936164","rs1470560","rs1935237","rs1969066",
                               "rs2657896","rs2786178","rs301798","rs3808348",
                               "rs4018180","rs424539","rs4690098","rs58936320",
                               "rs60669734","rs6435156","rs6443429","rs663015",
                               "rs6695572","rs6734238","rs6801781","rs685218",
                               "rs6905544","rs7127808","rs7528419","rs7970695")

BMI_associated_SNPs = c(BMI_associated_SNPs,other_BMI_associated_SNPs)

noBMI_CRP_sum_stats = CRP_sum_stats %>%
  filter(!variant_id %in% BMI_associated_SNPs)

# Prepare CRP summary statistics

CRP_sum_stats = CRP_sum_stats %>%
  mutate(exposure = "CRP")

noBMI_CRP_sum_stats = noBMI_CRP_sum_stats %>%
  mutate(exposure = "CRP_noBMI")


colnames(CRP_sum_stats) = c("RSID","Chr", "Pos", "A1", "A2",
                        "Effect_exposure", "StdErr_exposure",
                        "P_exposure", "exposure")

colnames(noBMI_CRP_sum_stats) = c("RSID","Chr", "Pos", "A1", "A2",
                            "Effect_exposure", "StdErr_exposure",
                            "P_exposure", "exposure")

## Create sorted allele representations
CRP_sum_stats <- CRP_sum_stats %>%
  mutate(
    A1 = A1 %>% toupper(),
    A2 = A2 %>% toupper(),
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

noBMI_CRP_sum_stats <- noBMI_CRP_sum_stats %>%
  mutate(
    A1 = A1 %>% toupper(),
    A2 = A2 %>% toupper(),
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

# CRP merge and harmonize

## Overall_BCa
Overall_BCa_outcome = Overall_BCa_outcome_data %>%
  filter(base_pair_location %in% CRP_sum_stats$Pos) %>%
  mutate(outcome = "Overall_BCa")

Overall_BCa_outcome = Overall_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Overall_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
Overall_BCa_outcome <- Overall_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

# Merge and harmonize
Overall_BCa_harmonized <- merge(
  CRP_sum_stats,
  Overall_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Overall_BCa_harmonized = Overall_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## HER2_pos_BCa
HER2_pos_BCa_outcome = HER2_pos_BCa_outcome_data %>%
  filter(base_pair_location %in% CRP_sum_stats$Pos) %>%
  mutate(outcome = "HER2_pos_BCa",
         p_value = 10^(-neg_log_10_p_value))

HER2_pos_BCa_outcome = HER2_pos_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(HER2_pos_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
HER2_pos_BCa_outcome <- HER2_pos_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

# Merge and harmonize
HER2_pos_BCa_harmonized <- merge(
  CRP_sum_stats,
  HER2_pos_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
HER2_pos_BCa_harmonized = HER2_pos_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))



## Luminal_A_BCa
Luminal_A_BCa_outcome = Luminal_A_BCa_outcome_data %>%
  filter(base_pair_location %in% CRP_sum_stats$Pos) %>%
  mutate(outcome = "Luminal_A_BCa",
         p_value = 10^(-neg_log_10_p_value))

Luminal_A_BCa_outcome = Luminal_A_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Luminal_A_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")

### Create sorted allele representations
Luminal_A_BCa_outcome <- Luminal_A_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

Luminal_A_BCa_harmonized <- merge(
  CRP_sum_stats,
  Luminal_A_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Luminal_A_BCa_harmonized = Luminal_A_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## Luminal_B_BCa
Luminal_B_BCa_outcome = Luminal_B_BCa_outcome_data %>%
  filter(base_pair_location %in% CRP_sum_stats$Pos) %>%
  mutate(outcome = "Luminal_B_BCa",
         p_value = 10^(-neg_log_10_p_value))

Luminal_B_BCa_outcome = Luminal_B_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Luminal_B_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")

### Create sorted allele representations
Luminal_B_BCa_outcome <- Luminal_B_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

Luminal_B_BCa_harmonized <- merge(
  CRP_sum_stats,
  Luminal_B_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Luminal_B_BCa_harmonized = Luminal_B_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## ER_neg_BCa
ER_neg_BCa_outcome = ER_neg_BCa_outcome_data %>%
  filter(base_pair_location %in% CRP_sum_stats$Pos) %>%
  mutate(outcome = "ER_neg_BCa")

ER_neg_BCa_outcome = ER_neg_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(ER_neg_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                 "Effect_outcome", "StdErr_outcome",
                                 "A1", "A2", "outcome")

### Create sorted allele representations
ER_neg_BCa_outcome <- ER_neg_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

ER_neg_BCa_harmonized <- merge(
  CRP_sum_stats,
  ER_neg_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
ER_neg_BCa_harmonized = ER_neg_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## TNBCa
TNBCa_outcome = TNBCa_outcome_data %>%
  filter(base_pair_location %in% CRP_sum_stats$Pos) %>%
  mutate(outcome = "TNBCa",
         p_value = 10^(-neg_log_10_p_value))

TNBCa_outcome = TNBCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(TNBCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
TNBCa_outcome <- TNBCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

# Merge and harmonize
TNBCa_harmonized <- merge(
  CRP_sum_stats,
  TNBCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
TNBCa_harmonized = TNBCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))


# Join harmonized data
CRP_harmonized = rbind.data.frame(Overall_BCa_harmonized,
                                  HER2_pos_BCa_harmonized,
                                  Luminal_A_BCa_harmonized,
                              Luminal_B_BCa_harmonized,
                              ER_neg_BCa_harmonized,
                              TNBCa_harmonized)

# No BMI SNP CRP merge and harmonize

## Overall_BCa
Overall_BCa_outcome = Overall_BCa_outcome_data %>%
  filter(base_pair_location %in% noBMI_CRP_sum_stats$Pos) %>%
  mutate(outcome = "Overall_BCa")

Overall_BCa_outcome = Overall_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Overall_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
Overall_BCa_outcome <- Overall_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

# Merge and harmonize
Overall_BCa_harmonized <- merge(
  noBMI_CRP_sum_stats,
  Overall_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Overall_BCa_harmonized = Overall_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## HER2_pos_BCa
HER2_pos_BCa_outcome = HER2_pos_BCa_outcome_data %>%
  filter(base_pair_location %in% noBMI_CRP_sum_stats$Pos) %>%
  mutate(outcome = "HER2_pos_BCa",
         p_value = 10^(-neg_log_10_p_value))

HER2_pos_BCa_outcome = HER2_pos_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(HER2_pos_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
HER2_pos_BCa_outcome <- HER2_pos_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

# Merge and harmonize
HER2_pos_BCa_harmonized <- merge(
  noBMI_CRP_sum_stats,
  HER2_pos_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
HER2_pos_BCa_harmonized = HER2_pos_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## Luminal_A_BCa
Luminal_A_BCa_outcome = Luminal_A_BCa_outcome_data %>%
  filter(base_pair_location %in% noBMI_CRP_sum_stats$Pos) %>%
  mutate(outcome = "Luminal_A_BCa",
         p_value = 10^(-neg_log_10_p_value))

Luminal_A_BCa_outcome = Luminal_A_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Luminal_A_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")

### Create sorted allele representations
Luminal_A_BCa_outcome <- Luminal_A_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

Luminal_A_BCa_harmonized <- merge(
  noBMI_CRP_sum_stats,
  Luminal_A_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Luminal_A_BCa_harmonized = Luminal_A_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## Luminal_B_BCa
Luminal_B_BCa_outcome = Luminal_B_BCa_outcome_data %>%
  filter(base_pair_location %in% noBMI_CRP_sum_stats$Pos) %>%
  mutate(outcome = "Luminal_B_BCa",
         p_value = 10^(-neg_log_10_p_value))

Luminal_B_BCa_outcome = Luminal_B_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Luminal_B_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")

### Create sorted allele representations
Luminal_B_BCa_outcome <- Luminal_B_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

Luminal_B_BCa_harmonized <- merge(
  noBMI_CRP_sum_stats,
  Luminal_B_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Luminal_B_BCa_harmonized = Luminal_B_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## ER_neg_BCa
ER_neg_BCa_outcome = ER_neg_BCa_outcome_data %>%
  filter(base_pair_location %in% noBMI_CRP_sum_stats$Pos) %>%
  mutate(outcome = "ER_neg_BCa")

ER_neg_BCa_outcome = ER_neg_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(ER_neg_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                 "Effect_outcome", "StdErr_outcome",
                                 "A1", "A2", "outcome")

### Create sorted allele representations
ER_neg_BCa_outcome <- ER_neg_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

ER_neg_BCa_harmonized <- merge(
  noBMI_CRP_sum_stats,
  ER_neg_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
ER_neg_BCa_harmonized = ER_neg_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## TNBCa
TNBCa_outcome = TNBCa_outcome_data %>%
  filter(base_pair_location %in% noBMI_CRP_sum_stats$Pos) %>%
  mutate(outcome = "TNBCa",
         p_value = 10^(-neg_log_10_p_value))

TNBCa_outcome = TNBCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(TNBCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
TNBCa_outcome <- TNBCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

# Merge and harmonize
TNBCa_harmonized <- merge(
  noBMI_CRP_sum_stats,
  TNBCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
TNBCa_harmonized = TNBCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

# Join harmonized data
noBMI_CRP_harmonized = rbind.data.frame(Overall_BCa_harmonized,
                                        HER2_pos_BCa_harmonized,
                                        Luminal_A_BCa_harmonized,
                                  Luminal_B_BCa_harmonized,
                                  ER_neg_BCa_harmonized,
                                  TNBCa_harmonized)

# Load GWAS summary statistics for TNFR2 

## Define CRP instrument including BMI-associated SNPs 

TNFR2_sum_stats = fread("MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/raw_summary_statistics/TNFR2.tsv")

TNFR2_SNPs = read_VEP("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/TNFR2_VEP.tsv") %>%
  filter(EUR_AF > 0.05) %>%
  select(Uploaded_variation) %>%
  filter(Uploaded_variation %in% c("rs2256974", "rs28399563"))%>%
  pull()

TNFR2_sum_stats = TNFR2_sum_stats %>%
  filter(p_value <5E-8)


# Prepare TNFR2 summary statistics

TNFR2_sum_stats = TNFR2_sum_stats %>%
  mutate(exposure = "TNFR2") %>%
  select(chromosome, base_pair_location, effect_allele, other_allele,
         beta, standard_error, p_value, exposure)

colnames(TNFR2_sum_stats) = c("Chr", "Pos", "A1", "A2",
                            "Effect_exposure", "StdErr_exposure",
                            "P_exposure", "exposure")

## Create sorted allele representations
TNFR2_sum_stats <- TNFR2_sum_stats %>%
  mutate(
    A1 = A1 %>% toupper(),
    A2 = A2 %>% toupper(),
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

# Harmonize TNFR2 summary statistics with outcome data

## Overall_BCa
Overall_BCa_outcome = Overall_BCa_outcome_data %>%
  filter(base_pair_location %in% TNFR2_sum_stats$Pos) %>%
  mutate(outcome = "Overall_BCa")

Overall_BCa_outcome = Overall_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Overall_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
Overall_BCa_outcome <- Overall_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

### Merge and harmonize
Overall_BCa_outcome <- merge(
  Overall_BCa_outcome,
  bim,
  by = c("Chr", "Pos", "A_min", "A_max"),
  all.x = TRUE,
  all.y = FALSE
)

Overall_BCa_harmonized <- merge(
  TNFR2_sum_stats,
  Overall_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Overall_BCa_harmonized = Overall_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## HER2_pos_BCa
HER2_pos_BCa_outcome = HER2_pos_BCa_outcome_data %>%
  filter(base_pair_location %in% TNFR2_sum_stats$Pos) %>%
  mutate(outcome = "HER2_pos_BCa",
         p_value = 10^(-neg_log_10_p_value))

HER2_pos_BCa_outcome = HER2_pos_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(HER2_pos_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
HER2_pos_BCa_outcome <- HER2_pos_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

### Merge and harmonize
HER2_pos_BCa_outcome <- merge(
  HER2_pos_BCa_outcome,
  bim,
  by = c("Chr", "Pos", "A_min", "A_max"),
  all.x = TRUE,
  all.y = FALSE
)

HER2_pos_BCa_harmonized <- merge(
  TNFR2_sum_stats,
  HER2_pos_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
HER2_pos_BCa_harmonized = HER2_pos_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## Luminal_A_BCa
Luminal_A_BCa_outcome = Luminal_A_BCa_outcome_data %>%
  filter(base_pair_location %in% TNFR2_sum_stats$Pos) %>%
  mutate(outcome = "Luminal_A_BCa",
         p_value = 10^(-neg_log_10_p_value))

Luminal_A_BCa_outcome = Luminal_A_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Luminal_A_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
Luminal_A_BCa_outcome <- Luminal_A_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

### Merge and harmonize
Luminal_A_BCa_outcome <- merge(
  Luminal_A_BCa_outcome,
  bim,
  by = c("Chr", "Pos", "A_min", "A_max"),
  all.x = TRUE,
  all.y = FALSE
)

Luminal_A_BCa_harmonized <- merge(
  TNFR2_sum_stats,
  Luminal_A_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Luminal_A_BCa_harmonized = Luminal_A_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## Luminal_B_BCa
Luminal_B_BCa_outcome = Luminal_B_BCa_outcome_data %>%
  filter(base_pair_location %in% TNFR2_sum_stats$Pos) %>%
  mutate(outcome = "Luminal_B_BCa",
         p_value = 10^(-neg_log_10_p_value))

Luminal_B_BCa_outcome = Luminal_B_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(Luminal_B_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                    "Effect_outcome", "StdErr_outcome",
                                    "A1", "A2", "outcome")

### Create sorted allele representations
Luminal_B_BCa_outcome <- Luminal_B_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

### Merge and harmonize
Luminal_B_BCa_outcome <- merge(
  Luminal_B_BCa_outcome,
  bim,
  by = c("Chr", "Pos", "A_min", "A_max"),
  all.x = TRUE,
  all.y = FALSE
)

Luminal_B_BCa_harmonized <- merge(
  TNFR2_sum_stats,
  Luminal_B_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
Luminal_B_BCa_harmonized = Luminal_B_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## ER_neg_BCa
ER_neg_BCa_outcome = ER_neg_BCa_outcome_data %>%
  filter(base_pair_location %in% TNFR2_sum_stats$Pos) %>%
  mutate(outcome = "ER_neg_BCa")

ER_neg_BCa_outcome = ER_neg_BCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(ER_neg_BCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                 "Effect_outcome", "StdErr_outcome",
                                 "A1", "A2", "outcome")

### Create sorted allele representations
ER_neg_BCa_outcome <- ER_neg_BCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

### Merge and harmonize
ER_neg_BCa_outcome <- merge(
  ER_neg_BCa_outcome,
  bim,
  by = c("Chr", "Pos", "A_min", "A_max"),
  all.x = TRUE,
  all.y = FALSE
)

ER_neg_BCa_harmonized <- merge(
  TNFR2_sum_stats,
  ER_neg_BCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
ER_neg_BCa_harmonized = ER_neg_BCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

## TNBCa
TNBCa_outcome = TNBCa_outcome_data %>%
  filter(base_pair_location %in% TNFR2_sum_stats$Pos) %>%
  mutate(outcome = "TNBCa",
         p_value = 10^(-neg_log_10_p_value))

TNBCa_outcome = TNBCa_outcome %>%
  select(c("chromosome","base_pair_location","p_value","beta",
           "standard_error","effect_allele","other_allele","outcome"))

colnames(TNBCa_outcome) = c("Chr", "Pos", "P_outcome", 
                                  "Effect_outcome", "StdErr_outcome",
                                  "A1", "A2", "outcome")

### Create sorted allele representations
TNBCa_outcome <- TNBCa_outcome %>%
  mutate(
    A_min = pmin(A1, A2),
    A_max = pmax(A1, A2)
  )

### Merge and harmonize
TNBCa_outcome <- merge(
  TNBCa_outcome,
  bim,
  by = c("Chr", "Pos", "A_min", "A_max"),
  all.x = TRUE,
  all.y = FALSE
)

TNBCa_harmonized <- merge(
  TNFR2_sum_stats,
  TNBCa_outcome,
  by = c("Chr", "Pos", "A_min", "A_max"),
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
TNBCa_harmonized = TNBCa_harmonized %>%
  mutate(A1 = toupper(A1.x),
         A2 = toupper(A2.x)) %>%
  select(-c(A1.x,A1.y,A2.x,A2.y))

# Join harmonized data
TNFR2_harmonized = rbind.data.frame(Overall_BCa_harmonized,
                                    HER2_pos_BCa_harmonized,
                                    Luminal_A_BCa_harmonized,
                              Luminal_B_BCa_harmonized,
                              ER_neg_BCa_harmonized,
                              TNBCa_harmonized) %>%
  filter(RSID %in% TNFR2_SNPs)

# Save data
write.table(
  CRP_harmonized,
  file = "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/BCa/CRP_harmonized_data.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  noBMI_CRP_harmonized,
  file = "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/BCa/CRP_noBMI_harmonized_data.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  TNFR2_harmonized,
  file = "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/BCa/TNFR2_harmonized_data.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
