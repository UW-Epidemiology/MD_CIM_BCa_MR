# Setup
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)
setDTthreads(percent = 100)

options(ieugwasr_api = 'https://api.opengwas.io/')

jwt = "jwt" # JWT for IEU OpenGWAS Project API would go here

# Function to add European ancestry 1000G allele frequencies

add_eaf_1000G = function(data){
  af_data = afl2_rsid(data$SNP, opengwas_jwt = jwt)
  af_data <- af_data[order(as.numeric(af_data$chr), af_data$start), ]
  EUR_af = af_data$AF.EUR
  ref = af_data$ref
  
  # Check if the effect allele is the same as the 1000G reference allele
  if("exposure" %in% names(data)){
    diff = ifelse(data$effect_allele.exposure == ref, TRUE, FALSE)
    data$eaf.exposure = ifelse(diff, EUR_af, 1-EUR_af)
  } else {
    diff = ifelse(data$effect_allele.outcome == ref, TRUE, FALSE)
    data$eaf.outcome = ifelse(diff, EUR_af, 1-EUR_af)
    
  }
  
  
  return(data)
}

# Load instrument exposure data

## IL6
IL6_exposure_data = fread("/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/IL6_instrument_data.csv")
IL6_exposure_data = add_eaf_1000G(IL6_exposure_data)

### Filter to MAF > 0.05 and p < 5e-8
IL6_exposure_data <- IL6_exposure_data %>% filter(eaf.exposure > 0.05 & pval.exposure < 5e-8)

## TNF-a
#### No SNPs with sufficiently strong associations
TNF_a_exposure_data <- extract_instruments("ebi-a-GCST004426")


## CRP and BMI for MVMR
#### Using the largest ieugwasr European ancestry GWAS datasets for CRP and BMI
CRP_BMI_exposure_data = mv_extract_exposures(id_exposure = c("ukb-b-19953",
                                                                    "ebi-a-GCST90029070"),
                                                    harmonise_strictness = 2,
                                             opengwas_jwt = jwt,
                                             pop = "EUR")

## CRP for UVMR

### Without BMI adjustment 
CRP_exposure_data_no_BMI_adj = extract_instruments("ebi-a-GCST90029070",
                                        opengwas_jwt = jwt)

## IL8
IL8_exposure_data <- extract_instruments("ebi-a-GCST90011994",
                                         opengwas_jwt = jwt)

### Filter to MAF > 0.05 and p < 5e-8
IL8_exposure_data <- IL8_exposure_data %>% filter(eaf.exposure > 0.05 & pval.exposure < 5e-8)

## TGF_b
TGF_b_exposure_data <- extract_instruments("prot-a-2962")

#### No SNPs with sufficiently strong associations

## TNFR2
TNFR2_exposure_data <- extract_instruments("ebi-a-GCST90012026")

### Filter to MAF > 0.05 and p < 5e-8
TNFR2_exposure_data <- TNFR2_exposure_data %>% filter(eaf.exposure > 0.05 & pval.exposure < 5e-8)

## VEGF
VEGF_exposure_data <- extract_instruments("ebi-a-GCST90011995")

### Filter to MAF > 0.05 and p < 5e-8
VEGF_exposure_data <- VEGF_exposure_data %>% filter(eaf.exposure > 0.05 & pval.exposure < 5e-8)


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


# Filter and harmonize MD GWAS outcome data to instrument SNPs

## PMD and IL6

PMD_IL6_instrument_outcome_data = PMD_outcome_data %>% 
  filter(rsid %in% IL6_exposure_data$SNP) %>% as.data.frame()
PMD_IL6_instrument_outcome_data = PMD_IL6_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(PMD_IL6_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                             "se","effect_allele","other_allele")
PMD_IL6_instrument_outcome_data = format_data(PMD_IL6_instrument_outcome_data,
                                             type = "outcome")
PMD_IL6_instrument_outcome_data$id.outcome = "PMD"
PMD_IL6_instrument_outcome_data = add_eaf_1000G(PMD_IL6_instrument_outcome_data)
PMD_IL6_harmonized_data = harmonise_data(IL6_exposure_data,
                                        PMD_IL6_instrument_outcome_data)

## DA and IL6

DA_IL6_instrument_outcome_data = DA_outcome_data %>% 
  filter(rsid %in% IL6_exposure_data$SNP) %>% as.data.frame()
DA_IL6_instrument_outcome_data = DA_IL6_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(DA_IL6_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                             "se","effect_allele","other_allele")
DA_IL6_instrument_outcome_data = format_data(DA_IL6_instrument_outcome_data,
                                             type = "outcome")
DA_IL6_instrument_outcome_data$id.outcome = "DA"
DA_IL6_instrument_outcome_data = add_eaf_1000G(DA_IL6_instrument_outcome_data)
DA_IL6_harmonized_data = harmonise_data(IL6_exposure_data,
                                        DA_IL6_instrument_outcome_data)

## NDA and IL6

NDA_IL6_instrument_outcome_data = NDA_outcome_data %>% 
  filter(rsid %in% IL6_exposure_data$SNP) %>% as.data.frame()
NDA_IL6_instrument_outcome_data = NDA_IL6_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(NDA_IL6_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                             "se","effect_allele","other_allele")
NDA_IL6_instrument_outcome_data = format_data(NDA_IL6_instrument_outcome_data,
                                             type = "outcome")
NDA_IL6_instrument_outcome_data$id.outcome = "NDA"
NDA_IL6_instrument_outcome_data = add_eaf_1000G(NDA_IL6_instrument_outcome_data)
NDA_IL6_harmonized_data = harmonise_data(IL6_exposure_data,
                                        NDA_IL6_instrument_outcome_data)

## PMD, CRP, BMI
PMD_CRP_and_BMI_instrument_outcome_data = PMD_outcome_data %>% 
  filter(rsid %in% unique(CRP_BMI_leptin_exposure_data$SNP)) %>% as.data.frame()
PMD_CRP_and_BMI_instrument_outcome_data = PMD_CRP_and_BMI_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(PMD_CRP_and_BMI_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                              "se","effect_allele","other_allele")
PMD_CRP_and_BMI_instrument_outcome_data = format_data(PMD_CRP_and_BMI_instrument_outcome_data,
                                              type = "outcome")
PMD_CRP_and_BMI_instrument_outcome_data$id.outcome = "PMD"
PMD_CRP_and_BMI_instrument_outcome_data = add_eaf_1000G(PMD_CRP_and_BMI_instrument_outcome_data)
PMD_CRP_and_BMI_harmonized_data = mv_harmonise_data(CRP_BMI_leptin_exposure_data,
                                                    PMD_CRP_and_BMI_instrument_outcome_data)

## DA, CRP, BMI
DA_CRP_and_BMI_instrument_outcome_data = DA_outcome_data %>% 
  filter(rsid %in% unique(CRP_BMI_leptin_exposure_data$SNP)) %>% as.data.frame()
DA_CRP_and_BMI_instrument_outcome_data = DA_CRP_and_BMI_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(DA_CRP_and_BMI_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                                      "se","effect_allele","other_allele")
DA_CRP_and_BMI_instrument_outcome_data = format_data(DA_CRP_and_BMI_instrument_outcome_data,
                                                      type = "outcome")
DA_CRP_and_BMI_instrument_outcome_data$id.outcome = "DA"
DA_CRP_and_BMI_instrument_outcome_data = add_eaf_1000G(DA_CRP_and_BMI_instrument_outcome_data)
DA_CRP_and_BMI_harmonized_data = mv_harmonise_data(CRP_BMI_leptin_exposure_data,
                                                    DA_CRP_and_BMI_instrument_outcome_data)

## NDA, CRP, BMI
NDA_CRP_and_BMI_instrument_outcome_data = NDA_outcome_data %>% 
  filter(rsid %in% unique(CRP_BMI_leptin_exposure_data$SNP)) %>% as.data.frame()
NDA_CRP_and_BMI_instrument_outcome_data = NDA_CRP_and_BMI_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(NDA_CRP_and_BMI_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                                      "se","effect_allele","other_allele")
NDA_CRP_and_BMI_instrument_outcome_data = format_data(NDA_CRP_and_BMI_instrument_outcome_data,
                                                      type = "outcome")
NDA_CRP_and_BMI_instrument_outcome_data$id.outcome = "NDA"
NDA_CRP_and_BMI_instrument_outcome_data = add_eaf_1000G(NDA_CRP_and_BMI_instrument_outcome_data)
NDA_CRP_and_BMI_harmonized_data = mv_harmonise_data(CRP_BMI_leptin_exposure_data,
                                                    NDA_CRP_and_BMI_instrument_outcome_data)

## PMD and CRP no BMI adjustment

PMD_CRP_noBMIadj_instrument_outcome_data = PMD_outcome_data %>% 
  filter(rsid %in% CRP_exposure_data_no_BMI_adj$SNP) %>% as.data.frame()
PMD_CRP_noBMIadj_instrument_outcome_data = PMD_CRP_noBMIadj_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(PMD_CRP_noBMIadj_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                              "se","effect_allele","other_allele")
PMD_CRP_noBMIadj_instrument_outcome_data$eaf = NA

PMD_CRP_noBMIadj_instrument_outcome_data = format_data(PMD_CRP_noBMIadj_instrument_outcome_data,
                                              type = "outcome")
PMD_CRP_noBMIadj_instrument_outcome_data$id.outcome = "PMD"
PMD_CRP_noBMIadj_instrument_outcome_data = harmonise_data(CRP_exposure_data_no_BMI_adj,
                                                          PMD_CRP_noBMIadj_instrument_outcome_data)

## DA and CRP no BMI adjustment

DA_CRP_noBMIadj_instrument_outcome_data = DA_outcome_data %>% 
  filter(rsid %in% CRP_exposure_data_no_BMI_adj$SNP) %>% as.data.frame()
DA_CRP_noBMIadj_instrument_outcome_data = DA_CRP_noBMIadj_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(DA_CRP_noBMIadj_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                                       "se","effect_allele","other_allele")
DA_CRP_noBMIadj_instrument_outcome_data$eaf = NA

DA_CRP_noBMIadj_instrument_outcome_data = format_data(DA_CRP_noBMIadj_instrument_outcome_data,
                                                       type = "outcome")
DA_CRP_noBMIadj_instrument_outcome_data$id.outcome = "DA"
DA_CRP_noBMIadj_instrument_outcome_data = harmonise_data(CRP_exposure_data_no_BMI_adj,
                                                          DA_CRP_noBMIadj_instrument_outcome_data)

## NDA and CRP no BMI adjustment

NDA_CRP_noBMIadj_instrument_outcome_data = NDA_outcome_data %>% 
  filter(rsid %in% CRP_exposure_data_no_BMI_adj$SNP) %>% as.data.frame()
NDA_CRP_noBMIadj_instrument_outcome_data = NDA_CRP_noBMIadj_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(NDA_CRP_noBMIadj_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                                       "se","effect_allele","other_allele")
NDA_CRP_noBMIadj_instrument_outcome_data$eaf = NA

NDA_CRP_noBMIadj_instrument_outcome_data = format_data(NDA_CRP_noBMIadj_instrument_outcome_data,
                                                       type = "outcome")
NDA_CRP_noBMIadj_instrument_outcome_data$id.outcome = "NDA"
NDA_CRP_noBMIadj_instrument_outcome_data = harmonise_data(CRP_exposure_data_no_BMI_adj,
                                                          NDA_CRP_noBMIadj_instrument_outcome_data)

## PMD and IL8

PMD_IL8_instrument_outcome_data = PMD_outcome_data %>% 
  filter(rsid %in% IL8_exposure_data$SNP) %>% as.data.frame()
PMD_IL8_instrument_outcome_data = PMD_IL8_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(PMD_IL8_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                              "se","effect_allele","other_allele")
PMD_IL8_instrument_outcome_data = format_data(PMD_IL8_instrument_outcome_data,
                                              type = "outcome")
PMD_IL8_instrument_outcome_data$id.outcome = "PMD"
PMD_IL8_instrument_outcome_data = add_eaf_1000G(PMD_IL8_instrument_outcome_data)
PMD_IL8_harmonized_data = harmonise_data(IL8_exposure_data,
                                         PMD_IL8_instrument_outcome_data)

## DA and IL8

DA_IL8_instrument_outcome_data = DA_outcome_data %>% 
  filter(rsid %in% IL8_exposure_data$SNP) %>% as.data.frame()
DA_IL8_instrument_outcome_data = DA_IL8_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(DA_IL8_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                             "se","effect_allele","other_allele")
DA_IL8_instrument_outcome_data = format_data(DA_IL8_instrument_outcome_data,
                                             type = "outcome")
DA_IL8_instrument_outcome_data$id.outcome = "DA"
DA_IL8_instrument_outcome_data = add_eaf_1000G(DA_IL8_instrument_outcome_data)
DA_IL8_harmonized_data = harmonise_data(IL8_exposure_data,
                                        DA_IL8_instrument_outcome_data)

## NDA and IL8

NDA_IL8_instrument_outcome_data = NDA_outcome_data %>% 
  filter(rsid %in% IL8_exposure_data$SNP) %>% as.data.frame()
NDA_IL8_instrument_outcome_data = NDA_IL8_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(NDA_IL8_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                              "se","effect_allele","other_allele")
NDA_IL8_instrument_outcome_data = format_data(NDA_IL8_instrument_outcome_data,
                                              type = "outcome")
NDA_IL8_instrument_outcome_data$id.outcome = "NDA"
NDA_IL8_instrument_outcome_data = add_eaf_1000G(NDA_IL8_instrument_outcome_data)
NDA_IL8_harmonized_data = harmonise_data(IL8_exposure_data,
                                         NDA_IL8_instrument_outcome_data)
## PMD and TNFR2

PMD_TNFR2_instrument_outcome_data = PMD_outcome_data %>% 
  filter(rsid %in% TNFR2_exposure_data$SNP) %>% as.data.frame()
PMD_TNFR2_instrument_outcome_data = PMD_TNFR2_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(PMD_TNFR2_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                              "se","effect_allele","other_allele")
PMD_TNFR2_instrument_outcome_data = format_data(PMD_TNFR2_instrument_outcome_data,
                                              type = "outcome")
PMD_TNFR2_instrument_outcome_data$id.outcome = "PMD"
PMD_TNFR2_instrument_outcome_data = add_eaf_1000G(PMD_TNFR2_instrument_outcome_data)
PMD_TNFR2_harmonized_data = harmonise_data(TNFR2_exposure_data,
                                         PMD_TNFR2_instrument_outcome_data)

## DA and TNFR2

DA_TNFR2_instrument_outcome_data = DA_outcome_data %>% 
  filter(rsid %in% TNFR2_exposure_data$SNP) %>% as.data.frame()
DA_TNFR2_instrument_outcome_data = DA_TNFR2_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(DA_TNFR2_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                             "se","effect_allele","other_allele")
DA_TNFR2_instrument_outcome_data = format_data(DA_TNFR2_instrument_outcome_data,
                                             type = "outcome")
DA_TNFR2_instrument_outcome_data$id.outcome = "DA"
DA_TNFR2_instrument_outcome_data = add_eaf_1000G(DA_TNFR2_instrument_outcome_data)
DA_TNFR2_harmonized_data = harmonise_data(TNFR2_exposure_data,
                                        DA_TNFR2_instrument_outcome_data)

## NDA and TNFR2

NDA_TNFR2_instrument_outcome_data = NDA_outcome_data %>% 
  filter(rsid %in% TNFR2_exposure_data$SNP) %>% as.data.frame()
NDA_TNFR2_instrument_outcome_data = NDA_TNFR2_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(NDA_TNFR2_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                              "se","effect_allele","other_allele")
NDA_TNFR2_instrument_outcome_data = format_data(NDA_TNFR2_instrument_outcome_data,
                                              type = "outcome")
NDA_TNFR2_instrument_outcome_data$id.outcome = "NDA"
NDA_TNFR2_instrument_outcome_data = add_eaf_1000G(NDA_TNFR2_instrument_outcome_data)
NDA_TNFR2_harmonized_data = harmonise_data(TNFR2_exposure_data,
                                         NDA_TNFR2_instrument_outcome_data)

## PMD and VEGF

PMD_VEGF_instrument_outcome_data = PMD_outcome_data %>% 
  filter(rsid %in% VEGF_exposure_data$SNP) %>% as.data.frame()
PMD_VEGF_instrument_outcome_data = PMD_VEGF_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(PMD_VEGF_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                              "se","effect_allele","other_allele")
PMD_VEGF_instrument_outcome_data = format_data(PMD_VEGF_instrument_outcome_data,
                                              type = "outcome")
PMD_VEGF_instrument_outcome_data$id.outcome = "PMD"
PMD_VEGF_instrument_outcome_data = add_eaf_1000G(PMD_VEGF_instrument_outcome_data)
PMD_VEGF_harmonized_data = harmonise_data(VEGF_exposure_data,
                                         PMD_VEGF_instrument_outcome_data)

## DA and VEGF

DA_VEGF_instrument_outcome_data = DA_outcome_data %>% 
  filter(rsid %in% VEGF_exposure_data$SNP) %>% as.data.frame()
DA_VEGF_instrument_outcome_data = DA_VEGF_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(DA_VEGF_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                             "se","effect_allele","other_allele")
DA_VEGF_instrument_outcome_data = format_data(DA_VEGF_instrument_outcome_data,
                                             type = "outcome")
DA_VEGF_instrument_outcome_data$id.outcome = "DA"
DA_VEGF_instrument_outcome_data = add_eaf_1000G(DA_VEGF_instrument_outcome_data)
DA_VEGF_harmonized_data = harmonise_data(VEGF_exposure_data,
                                        DA_VEGF_instrument_outcome_data)

## NDA and VEGF

NDA_VEGF_instrument_outcome_data = NDA_outcome_data %>% 
  filter(rsid %in% VEGF_exposure_data$SNP) %>% as.data.frame()
NDA_VEGF_instrument_outcome_data = NDA_VEGF_instrument_outcome_data %>%
  select(c("chr","position","rsid","P-value","Effect",
           "StdErr","Allele2","Allele1"))
colnames(NDA_VEGF_instrument_outcome_data) = c("chr","pos","SNP","pval","beta",
                                              "se","effect_allele","other_allele")
NDA_VEGF_instrument_outcome_data = format_data(NDA_VEGF_instrument_outcome_data,
                                              type = "outcome")
NDA_VEGF_instrument_outcome_data$id.outcome = "NDA"
NDA_VEGF_instrument_outcome_data = add_eaf_1000G(NDA_VEGF_instrument_outcome_data)
NDA_VEGF_harmonized_data = harmonise_data(VEGF_exposure_data,
                                         NDA_VEGF_instrument_outcome_data)

# Write output files

## PMD and IL6
fwrite(PMD_IL6_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_IL6_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## DA and IL6
fwrite(DA_IL6_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/DA_IL6_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## NDA and IL6
fwrite(NDA_IL6_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/NDA_IL6_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## PMD, CRP, NDA, CRP without BMI adjustment and BMI
saveRDS(PMD_CRP_and_BMI_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_CRP_noBMIadj_BMIharmonized_data.rds")

## DA, NDA, CRP without BMI adjustment and BMI
saveRDS(DA_CRP_and_BMI_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/DA_CRP_noBMIadj_BMI_harmonized_data.rds")

## NDA, CRP without BMI adjustment and BMI
saveRDS(NDA_CRP_and_BMI_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/NDA_CRP_noBMIadj_BMI_harmonized_data.rds")

## PMD and CRP without BMI adjustment
saveRDS(PMD_CRP_noBMIadj_instrument_outcome_data,
        "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_CRP_noBMIadj_harmonized_data.rds")

## DA and CRP without BMI adjustment
saveRDS(DA_CRP_noBMIadj_instrument_outcome_data,
        "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/DA_CRP_noBMIadj_harmonized_data.rds")

## NDA and CRP without BMI adjustment
saveRDS(NDA_CRP_noBMIadj_instrument_outcome_data,
        "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/NDA_CRP_noBMIadj_harmonized_data.rds")

## PMD and IL8
fwrite(PMD_IL8_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_IL8_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## DA and IL8
fwrite(DA_IL8_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/DA_IL8_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## NDA and IL8
fwrite(NDA_IL8_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/NDA_IL8_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## PMD and TNFR2
fwrite(PMD_TNFR2_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_TNFR2_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## DA and TNFR2
fwrite(DA_TNFR2_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/DA_TNFR2_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## NDA and TNFR2
fwrite(NDA_TNFR2_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/NDA_TNFR2_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## PMD and VEGF
fwrite(PMD_VEGF_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_VEGF_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## DA and VEGF
fwrite(DA_VEGF_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/DA_VEGF_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)

## NDA and VEGF
fwrite(NDA_VEGF_harmonized_data,
       "/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/NDA_VEGF_harmonized_data.csv",
       row.names = FALSE, quote = FALSE)
