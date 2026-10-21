# Setup
library(data.table)
library(R.utils)
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)
library(vcfR)

setDTthreads(percent = 100)

# Function to create MVMR exposure objects

# BMI GWAS summary statistics
BMI_sum_stats = vcfR::read.vcfR("/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/ukb-b-19953.vcf.gz",
                                verbose = FALSE)

data = data.table(CHR = BMI_sum_stats@fix[,1],
                  POS = BMI_sum_stats@fix[,2],
                  SNP = BMI_sum_stats@fix[,3],
                  REF = BMI_sum_stats@fix[,4],
                  ALT = BMI_sum_stats@fix[,5])

data$Effect <- sapply(strsplit(BMI_sum_stats@gt[,2], ":"), `[`, 1)
data$SE <- sapply(strsplit(BMI_sum_stats@gt[,2], ":"), `[`, 2)
data$log_P <- sapply(strsplit(BMI_sum_stats@gt[,2], ":"), `[`, 3)
data$EAF <- sapply(strsplit(BMI_sum_stats@gt[,2], ":"), `[`, 4)


data = data %>% mutate(P=10^(as.numeric(log_P)*-1)) 
data$Phenotype = "BMI"

data_for_IDs = data %>% select(CHR,POS,ALT,REF,SNP)


fwrite(data, file = "/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/BMI_summary_statistics.tsv",
       sep = " ",
       row.names = FALSE, quote = FALSE)

data = data %>% select("SNP","Effect","SE","ALT","REF","P","Phenotype","EAF")
colnames(data) = c("SNP","beta","se","effect_allele","other_allele","pval","Phenotype","eaf")


data = as.data.frame(data)

rm(BMI_sum_stats)

CRP_and_BMI_SNPs = fread("/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/CRP_and_BMI_snps.txt", header = FALSE)
BMI_SNPs = data %>% filter(SNP %in% CRP_and_BMI_SNPs$V1, pval < 5E-8) %>% pull(SNP)


# Load MD GWAS outcome data (no adjustment for BMI) 

## PMD
PMD_outcome_data = fread("/lindstroem/MD_GWAS/FilesforCameron/Total_NoBMI/PMD_NoBMI_Stderr_MetaAnalysis_Cleaned_Final.txt")
PMD_outcome_data$Allele1 = toupper(PMD_outcome_data$Allele1)
PMD_outcome_data$Allele2 = toupper(PMD_outcome_data$Allele2)

## DA 
DA_outcome_data = fread("/lindstroem/MD_GWAS/FilesforCameron/Total_NoBMI/DA_NoBMI_Stderr_MetaAnalysis_Cleaned_Final.txt.gz")
DA_outcome_data$Allele1 = toupper(DA_outcome_data$Allele1)
DA_outcome_data$Allele2 = toupper(DA_outcome_data$Allele2)

## NDA
NDA_outcome_data = fread("/lindstroem/MD_GWAS/FilesforCameron/Total_NoBMI/Rerun_20210325/NDA_NoBMI_Stderr_MetaAnalysis_Cleaned_Final.txt")
NDA_outcome_data$Allele1 = toupper(NDA_outcome_data$Allele1)
NDA_outcome_data$Allele2 = toupper(NDA_outcome_data$Allele2)


# IMD CIMs

## Generate list of filepaths and corresponding proteins

dir_path_filtered <- "/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/harmonized_fully_filtered/"
dir_path_unfiltered <- "/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/sum_stats/"


filepaths_filtered <- list.files(
  path = dir_path_filtered,
  pattern = "\\.tsv$",
  full.names = TRUE
)

filepaths_unfiltered <- list.files(
  path = dir_path_unfiltered,
  pattern = "\\.tbl.gz$",
  full.names = TRUE
)

### Keep only files with more than 1 line
filepaths_filtered <- filepaths_filtered[sapply(filepaths_filtered, function(fp) {
  length(readLines(fp, n = 2)) > 1
})]

proteins <- sub("_.*$", "", basename(filepaths_filtered))

filepaths_unfiltered <- filepaths_unfiltered[
  sub("_.*$", "", basename(filepaths_unfiltered)) %in% paste0(proteins,".tbl.gz")
]

## Harmonize with BMI data and then merge with MD outcome data

for (j in 1:length(filepaths_filtered)) {
  raw_sum_stats = fread(filepaths_unfiltered[j])
  sum_stats = raw_sum_stats %>%
    select(Chromosome,Position,Allele1,Allele2,Freq1,Effect,StdErr,N,"log(P)")
  rm(raw_sum_stats)
  
  sum_stats = sum_stats %>%
    mutate(
      Allele1 = toupper(Allele1),
      Allele2 = toupper(Allele2),
      P = 10^(`log(P)`),
      Phenotype = proteins[j]
    ) %>%
    select(-`log(P)`)
  
  colnames(sum_stats) = c("CHR","POS","ALT","REF","eaf","beta","se","N","pval","Phenotype")
  
  sum_stats$CHR = as.character(sum_stats$CHR)
  sum_stats$POS = as.character(sum_stats$POS)
  
  sum_stats <- merge(
    sum_stats,
    data_for_IDs %>% select(-c(ALT,REF)),
    by = c("CHR","POS")
  )
  
  filtered_CIM_SNPs <- fread(filepaths_filtered[j]) %>%
    pull(RSID) %>% unique()
  
  sum_stats = sum_stats %>% filter(SNP %in% c(BMI_SNPs,filtered_CIM_SNPs))
  sum_stats = sum_stats %>% select("SNP","beta","se","ALT","REF","pval","eaf","N","Phenotype")
  colnames(sum_stats) = c("SNP","beta","se","effect_allele","other_allele","pval","eaf","N","Phenotype")
  
  data_for_harmonization = data %>%
    filter(SNP %in% c(BMI_SNPs,filtered_CIM_SNPs))
  
  data_for_harmonization = data_for_harmonization %>% as.data.frame() %>%
    mutate(eaf = as.numeric(eaf)) %>%
    filter(eaf > 0.05)
  
  sum_stats = as.data.frame(sum_stats)
  
  filenames_exposure = list(data_for_harmonization,sum_stats)
  sep = " "
  phenotype_col = "Phenotype"
  snp_col = "SNP"
  beta_col = "beta"
  se_col = "se"
  eaf_col = "eaf"
  effect_allele_col = "effect_allele"
  other_allele_col = "other_allele"
  pval_col = "pval"
  samplesize_col = "N"
  
  flag <- if (inherits(filenames_exposure, "list")) {
    stopifnot(all(sapply(filenames_exposure, inherits, "data.frame")))
    "data.frame"
  } else "character"
  
  n <- length(filenames_exposure)
  replicate_args <- function(arg) if(length(arg)==1) rep(arg, n) else arg
  sep <- replicate_args(sep)
  phenotype_col <- replicate_args(phenotype_col)
  snp_col <- replicate_args(snp_col)
  beta_col <- replicate_args(beta_col)
  se_col <- replicate_args(se_col)
  eaf_col <- replicate_args(eaf_col)
  effect_allele_col <- replicate_args(effect_allele_col)
  other_allele_col <- replicate_args(other_allele_col)
  pval_col <- replicate_args(pval_col)
  samplesize_col <- replicate_args(samplesize_col)
  
  l_full <- list()
  
  for (i in seq_along(filenames_exposure)) {
    # Read or format exposure data
    l_full[[i]] <- if(flag=="character") {
      read_outcome_data(
        filenames_exposure[i],
        sep = sep[i],
        phenotype_col = phenotype_col[i],
        snp_col = snp_col[i],
        beta_col = beta_col[i],
        se_col = se_col[i],
        eaf_col = eaf_col[i],
        effect_allele_col = effect_allele_col[i],
        other_allele_col = other_allele_col[i],
        pval_col = pval_col[i],
        samplesize_col = samplesize_col[i],
      )
    } else {
      format_data(
        filenames_exposure[[i]],
        type = "outcome",
        phenotype_col = phenotype_col[i],
        snp_col = snp_col[i],
        beta_col = beta_col[i],
        se_col = se_col[i],
        eaf_col = eaf_col[i],
        effect_allele_col = effect_allele_col[i],
        other_allele_col = other_allele_col[i],
        pval_col = pval_col[i],
        samplesize_col = samplesize_col[i],
      )
    }
    
    if(l_full[[i]]$outcome[1]=="outcome") {
      l_full[[i]]$outcome <- paste0("exposure", i)
    }
    
  }
  
  # Combine all exposures
  exposure_dat <- dplyr::bind_rows(l_full)
  id_exposure <- unique(exposure_dat$id.outcome)
  
  d1 <- lapply(l_full, function(x) subset(x, SNP %in% exposure_dat$SNP)) %>% dplyr::bind_rows()

  d_CIM <- subset(d1, id.outcome == id_exposure[2])
  d_BMI <- subset(d1, id.outcome == id_exposure[1])
  d_CIM <- convert_outcome_to_exposure(subset(d_CIM, id.outcome == id_exposure[2]))
  
  # Harmonize against the first id
  d <- harmonise_data(d_CIM, d_BMI, action = 2)
  
  # Only keep SNPs that are present in all
  tab <- table(d$SNP)
  keepsnps <- names(tab)[tab == length(id_exposure)-1]
  d <- subset(d, SNP %in% keepsnps)
  
  # Reshape exposures
  dh1 <- subset(d, id.outcome==id.outcome[1],
                select=c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure,
                         eaf.exposure, beta.exposure, se.exposure, pval.exposure))
  dh2 <- subset(d,
                select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome,
                         eaf.outcome, beta.outcome, se.outcome, pval.outcome))
  names(dh2) <- gsub("outcome", "exposure", names(dh2))
  MVMR_data <- rbind(dh1, dh2)
  
  # PMD
  PMD_instrument_outcome_data = PMD_outcome_data %>% 
    filter(rsid %in% unique(MVMR_data$SNP)) %>% as.data.frame()
  PMD_instrument_outcome_data = PMD_instrument_outcome_data %>%
    select(c("chr","position","rsid","P-value","Effect",
             "StdErr","Allele2","Allele1"))
  colnames(PMD_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                                        "se","effect_allele","other_allele")
  PMD_instrument_outcome_data = format_data(PMD_instrument_outcome_data,
                                                        type = "outcome")
  PMD_instrument_outcome_data$id.outcome = "PMD"
  MVMR_PMD = mv_harmonise_data(MVMR_data, PMD_instrument_outcome_data,
                               harmonise_strictness  = 2)
  
  # DA
  DA_instrument_outcome_data = DA_outcome_data %>% 
    filter(rsid %in% unique(MVMR_data$SNP)) %>% as.data.frame()
  DA_instrument_outcome_data = DA_instrument_outcome_data %>%
    select(c("chr","position","rsid","P-value","Effect",
             "StdErr","Allele2","Allele1"))
  colnames(DA_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                            "se","effect_allele","other_allele")
  DA_instrument_outcome_data = format_data(DA_instrument_outcome_data,
                                            type = "outcome")
  DA_instrument_outcome_data$id.outcome = "DA"
  MVMR_DA = mv_harmonise_data(MVMR_data, DA_instrument_outcome_data,
                              harmonise_strictness  = 2)
  
  # NDA
  NDA_instrument_outcome_data = NDA_outcome_data %>% 
    filter(rsid %in% unique(MVMR_data$SNP)) %>% as.data.frame()
  NDA_instrument_outcome_data = NDA_instrument_outcome_data %>%
    select(c("chr","position","rsid","P-value","Effect",
             "StdErr","Allele2","Allele1"))
  colnames(NDA_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                            "se","effect_allele","other_allele")
  NDA_instrument_outcome_data = format_data(NDA_instrument_outcome_data,
                                            type = "outcome")
  NDA_instrument_outcome_data$id.outcome = "NDA"
  MVMR_NDA = mv_harmonise_data(MVMR_data, NDA_instrument_outcome_data,
                               harmonise_strictness  = 2)
  
  # Save
  saveRDS(MVMR_PMD,
          paste0("/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/IMD_CIMs/PMD_",proteins[j],"_MVMR_data.rds"))
  saveRDS(MVMR_DA,
          paste0("/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/IMD_CIMs/DA_",proteins[j],"_MVMR_data.rds"))
  saveRDS(MVMR_NDA,
          paste0("/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/IMD_CIMs/NDA_",proteins[j],"_MVMR_data.rds"))
}

# Open GWAS project CIMs

## Load and process data

### TNFR2

TNFR2 = vcfR::read.vcfR("/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/TNFR2_sum_stats.vcf.gz",
                       verbose = FALSE)

TNFR2_data = data.frame(CHR = TNFR2@fix[,1],
                       POS = TNFR2@fix[,2],
                       SNP = TNFR2@fix[,3],
                       REF = TNFR2@fix[,4],
                       ALT = TNFR2@fix[,5])

TNFR2_data$Effect <- sapply(strsplit(TNFR2@gt[,2], ":"), `[`, 1)
TNFR2_data$SE <- sapply(strsplit(TNFR2@gt[,2], ":"), `[`, 2)
TNFR2_data$log_P <- sapply(strsplit(TNFR2@gt[,2], ":"), `[`, 3)
TNFR2_data$EAF <- sapply(strsplit(TNFR2@gt[,2], ":"), `[`, 4)


TNFR2_data = TNFR2_data %>% mutate(P=10^(as.numeric(log_P)*-1)) 
TNFR2_data$Phenotype = "TNFR2"


TNFR2_data = TNFR2_data %>% select("SNP","Effect","SE","ALT","REF","P","Phenotype","EAF")
colnames(TNFR2_data) = c("SNP","beta","se","effect_allele","other_allele","pval","Phenotype","eaf")


TNFR2_data = as.data.frame(TNFR2_data)

rm(TNFR2)

filtered_CIM_SNPs = c("rs7655660", "rs972492")

data_for_harmonization = data %>%
  filter(SNP %in% c(BMI_SNPs,filtered_CIM_SNPs))

data_for_harmonization = data_for_harmonization %>%
  as.data.frame()

sum_stats = as.data.frame(TNFR2_data %>% filter(SNP %in% c(BMI_SNPs,filtered_CIM_SNPs)))

filenames_exposure = list(data_for_harmonization,sum_stats)
sep = " "
phenotype_col = "Phenotype"
snp_col = "SNP"
beta_col = "beta"
se_col = "se"
eaf_col = "eaf"
effect_allele_col = "effect_allele"
other_allele_col = "other_allele"
pval_col = "pval"
samplesize_col = "N"

flag <- if (inherits(filenames_exposure, "list")) {
  stopifnot(all(sapply(filenames_exposure, inherits, "data.frame")))
  "data.frame"
} else "character"

n <- length(filenames_exposure)
replicate_args <- function(arg) if(length(arg)==1) rep(arg, n) else arg
sep <- replicate_args(sep)
phenotype_col <- replicate_args(phenotype_col)
snp_col <- replicate_args(snp_col)
beta_col <- replicate_args(beta_col)
se_col <- replicate_args(se_col)
eaf_col <- replicate_args(eaf_col)
effect_allele_col <- replicate_args(effect_allele_col)
other_allele_col <- replicate_args(other_allele_col)
pval_col <- replicate_args(pval_col)
samplesize_col <- replicate_args(samplesize_col)

l_full <- list()

for (i in seq_along(filenames_exposure)) {
  # Read or format exposure data
  l_full[[i]] <- if(flag=="character") {
    read_outcome_data(
      filenames_exposure[i],
      sep = sep[i],
      phenotype_col = phenotype_col[i],
      snp_col = snp_col[i],
      beta_col = beta_col[i],
      se_col = se_col[i],
      eaf_col = eaf_col[i],
      effect_allele_col = effect_allele_col[i],
      other_allele_col = other_allele_col[i],
      pval_col = pval_col[i],
      samplesize_col = samplesize_col[i],
    )
  } else {
    format_data(
      filenames_exposure[[i]],
      type = "outcome",
      phenotype_col = phenotype_col[i],
      snp_col = snp_col[i],
      beta_col = beta_col[i],
      se_col = se_col[i],
      eaf_col = eaf_col[i],
      effect_allele_col = effect_allele_col[i],
      other_allele_col = other_allele_col[i],
      pval_col = pval_col[i],
      samplesize_col = samplesize_col[i],
    )
  }
  
  if(l_full[[i]]$outcome[1]=="outcome") {
    l_full[[i]]$outcome <- paste0("exposure", i)
  }
  
}

# Combine all exposures
exposure_dat <- dplyr::bind_rows(l_full)
id_exposure <- unique(exposure_dat$id.outcome)

d1 <- lapply(l_full, function(x) subset(x, SNP %in% exposure_dat$SNP)) %>% dplyr::bind_rows()

d_CIM <- subset(d1, id.outcome == id_exposure[2])
d_BMI <- subset(d1, id.outcome == id_exposure[1])
d_CIM <- convert_outcome_to_exposure(subset(d_CIM, id.outcome == id_exposure[2]))

# Harmonize against the first id
d <- harmonise_data(d_CIM, d_BMI, action = 2)

# Only keep SNPs that are present in all
tab <- table(d$SNP)
keepsnps <- names(tab)[tab == length(id_exposure)-1]
d <- subset(d, SNP %in% keepsnps)

# Reshape exposures
dh1 <- subset(d, id.outcome==id.outcome[1],
              select=c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure,
                       eaf.exposure, beta.exposure, se.exposure, pval.exposure))
dh2 <- subset(d,
              select=c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome,
                       eaf.outcome, beta.outcome, se.outcome, pval.outcome))
names(dh2) <- gsub("outcome", "exposure", names(dh2))
MVMR_data <- rbind(dh1, dh2)

# PMD
PMD_instrument_outcome_data = PMD_outcome_data %>% 
  filter(rsid %in% unique(MVMR_data$SNP)) %>% as.data.frame()
PMD_instrument_outcome_data = PMD_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(PMD_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                          "se","effect_allele","other_allele")
PMD_instrument_outcome_data = format_data(PMD_instrument_outcome_data,
                                          type = "outcome")
PMD_instrument_outcome_data$id.outcome = "PMD"
MVMR_PMD = mv_harmonise_data(MVMR_data, PMD_instrument_outcome_data,
                             harmonise_strictness = 2)

# DA
DA_instrument_outcome_data = DA_outcome_data %>% 
  filter(rsid %in% unique(MVMR_data$SNP)) %>% as.data.frame()
DA_instrument_outcome_data = DA_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(DA_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                         "se","effect_allele","other_allele")
DA_instrument_outcome_data = format_data(DA_instrument_outcome_data,
                                         type = "outcome")
DA_instrument_outcome_data$id.outcome = "DA"
MVMR_DA = mv_harmonise_data(MVMR_data, DA_instrument_outcome_data,
                            harmonise_strictness = 2)

# NDA
NDA_instrument_outcome_data = NDA_outcome_data %>% 
  filter(rsid %in% unique(MVMR_data$SNP)) %>% as.data.frame()
NDA_instrument_outcome_data = NDA_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(NDA_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                          "se","effect_allele","other_allele")
NDA_instrument_outcome_data = format_data(NDA_instrument_outcome_data,
                                          type = "outcome")
NDA_instrument_outcome_data$id.outcome = "NDA"
MVMR_NDA = mv_harmonise_data(MVMR_data, NDA_instrument_outcome_data,
                             harmonise_strictness = 2)

# Save
saveRDS(MVMR_PMD,
        "/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/OpenGWAS_project_CIMs/PMD_TNFR2_MVMR_data.rds")
saveRDS(MVMR_DA,
        "/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/OpenGWAS_project_CIMs/DA_TNFR2_MVMR_data.rds")
saveRDS(MVMR_NDA,
        "/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/OpenGWAS_project_CIMs/NDA_TNFR2_MVMR_data.rds")


