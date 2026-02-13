# Setup
library(tidyverse)
library(stringr)
library(data.table)

# Functions

rbind_fill <- function(..., fill = NA) {
  dfs <- list(...)
  
  # Get all column names present in any data.frame
  all_cols <- unique(unlist(lapply(dfs, names)))
  
  # Make sure each data.frame has all columns
  dfs_filled <- lapply(dfs, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0) {
      df[missing_cols] <- fill
    }
    # Reorder columns to be consistent
    df[all_cols]
  })
  
  # Bind them together
  do.call(rbind, dfs_filled)
}

read_VEP <- function(x) {
  VEP_text <- readLines(x)
  VEP_text <- VEP_text[!startsWith(VEP_text, "##")]
  df <- fread(text = VEP_text, sep = "\t", na.strings = "-", data.table = TRUE)
  if (startsWith(names(df)[1], "#")) {
    setnames(df, old = names(df)[1], new = sub("^#", "", names(df)[1]))
  }
  
  return(df)
}

# Load data 
IMD_data = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_GWAS_CIM_MR_results.tsv")
OpenGWAS_data = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_project_CIM_MR_results.tsv") 
OpenGWAS_data <- OpenGWAS_data %>%
  mutate(
    across(
      where(is.character),
      ~ str_replace_all(.x,
                        c("IL6" = "IL-6",
                          "IL8" = "IL-8",
                          "VEGF" = "VEGF_A"))
    ),
    Exposure = case_when(Exposure == "IL-6" ~ "IL-6*",
                         Exposure == "IL-8" ~ "IL-8*",
                         Exposure == "VEGF_A" ~ "VEGF_A*",
                         TRUE ~ Exposure) 
    
  )


OpenGWAS_data <- OpenGWAS_data %>%
  mutate(
    Source = case_when(Exposure == "CRP" ~ "PMC9033829",
                       Exposure == "CRP_noBMI" ~ "PMC9033829",
                       Exposure == "CRP_locus" ~ "PMC9033829",
                       Exposure == "BMI" ~ "MRC IEU UK Biobank GWAS pipeline version 2",
                       Exposure == "IL-6*" ~ "PMC8098112",
                       Exposure == "IL-8*" ~ "PMC7611474",
                       Exposure == "TNFR2" ~ "PMC7611474",
                       Exposure == "VEGF_A*" ~ "PMC7611474"),
    Instrument = SNP)

IMD_data <- IMD_data %>%
  mutate(
    Source = "PMC10457199",
    Exposure = sub("\\.", "-", Exposure),
    Instrument = sub("\\.", "-", SNP)
  )

data = rbind_fill(IMD_data,OpenGWAS_data)

CRP_BMI_MVMR = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_results.tsv")
PMD_qhet = read_rds("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_PMD_QHET.rds")
CI_bounds <- str_match(
  PMD_qhet$`95% CI`,
  "^(-?\\d*\\.?\\d+)-(-?\\d*\\.?\\d+)$"
)
PMD_qhet$CILower <- as.numeric(CI_bounds[, 2])
PMD_qhet$CIUpper <- as.numeric(CI_bounds[, 3])
PMD_qhet = PMD_qhet %>%
  mutate(Exposure = c("CRP","BMI"),
         Outcome = "PMD",
         Estimate = `Effect Estimates`,
         Method = "Q-statistic minimization") %>%
  select(-c(`Effect Estimates`,`95% CI`))

DA_qhet = read_rds("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_DA_QHET.rds")
CI_bounds <- str_match(
  DA_qhet$`95% CI`,
  "^(-?\\d*\\.?\\d+)-(-?\\d*\\.?\\d+)$"
)
DA_qhet$CILower <- as.numeric(CI_bounds[, 2])
DA_qhet$CIUpper <- as.numeric(CI_bounds[, 3])
DA_qhet = DA_qhet %>%
  mutate(Exposure = c("CRP","BMI"),
         Outcome = "DA",
         Estimate = `Effect Estimates`,
         Method = "Q-statistic minimization") %>%
  select(-c(`Effect Estimates`,`95% CI`))

NDA_qhet = read_rds("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_NDA_QHET.rds")
CI_bounds <- str_match(
  NDA_qhet$`95% CI`,
  "^(-?\\d*\\.?\\d+)-(-?\\d*\\.?\\d+)$"
)
NDA_qhet$CILower <- as.numeric(CI_bounds[, 2])
NDA_qhet$CIUpper <- as.numeric(CI_bounds[, 3])
NDA_qhet = NDA_qhet %>%
  mutate(Exposure = c("CRP","BMI"),
         Outcome = "NDA",
         Estimate = `Effect Estimates`,
         Method = "Q-statistic minimization") %>%
  select(-c(`Effect Estimates`,`95% CI`))

data_MVMR = rbind_fill(CRP_BMI_MVMR,PMD_qhet,DA_qhet,NDA_qhet) %>%
  mutate(Source = ifelse(Exposure == "CRP", "PMC9033829", "UMRC IEU UK Biobank GWAS pipeline version 2"),
         Instrument = SNP)

# Exposure instrument SNPs and F-statistics (Table S3)
MR_F_stats <- data %>%
  select(Instrument, Exposure, Source, Fstatistic, NSNPs) %>%
  filter(!is.na(Exposure), !is.na(Fstatistic)) %>%
  distinct()

MR_SNPs <- data %>%
  filter(Method == "MR Wald Ratio") %>%
  group_by(Exposure) %>%
  summarise(
    SNPs = paste(unique(SNP), collapse = ", "),
    .groups = "drop"
  )

MVMR_F_stats = data_MVMR %>%
  select(Instrument, Exposure, Source, Fstatistic, NSNPs) %>%
  filter(!is.na(Exposure), !is.na(Fstatistic)) %>%
  distinct()

CRP_and_BMI_VEP = read_VEP("/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/CRP_and_BMI_VEP.tsv")
CRP_and_BMI_common = CRP_and_BMI_VEP %>%
  filter(EUR_AF > 0.05) %>%
  select(Uploaded_variation) %>%
  pull()

MVMR_F_stats$SNPs = paste(CRP_and_BMI_common, collapse = ", ")

Table_S3 = merge(MR_F_stats, MR_SNPs, by = "Exposure") 

Table_S3 = rbind_fill(Table_S3, MVMR_F_stats) %>%
  arrange(Exposure)

openxlsx::write.xlsx(
  Table_S3,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S3.xlsx",
  rowNames = FALSE
)


# MR-Egger estimates (Table S4)
Table_S4 = data %>%
  filter(Method == "MR-Egger") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate",
           "StdError","CILower","CIUpper","Pvalue","Pvalue_MR_Egger_Int","I_squared")) 

temp = data_MVMR %>%
  filter(Method == "MVMR-Egger") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate",
           "StdError","CILower","CIUpper","Pvalue","Pvalue_MVMR_Egger_Int"))

Table_S4 = rbind_fill(Table_S4,temp) %>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S4,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S4.xlsx",
  rowNames = FALSE
)


# MR-PRESSO estimates (Table S5)
Table_S5 = data %>%
  filter(Method %in% c("MR-PRESSO distortion test","MR-PRESSO global test","MR-PRESSO outlier-corrected")) %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue","N_MR_PRESSO_outliers","MR_PRESSO_outliers"))

temp <- data_MVMR %>%
  filter(Method %in% c(
    "MR-PRESSO distortion test",
    "MR-PRESSO global test",
    "MR-PRESSO outlier-corrected",
    "MR-PRESSO raw"
  )) %>%
  select(
    Instrument, Exposure, Outcome, Method,
    Estimate, StdError, CILower, CIUpper,
    Pvalue, N_MR_PRESSO_outliers, MR_PRESSO_outliers
  ) %>%
  mutate(
    Method = paste0("MV", Method)
  )

Table_S5 = rbind_fill(Table_S5,temp) %>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S5,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S5.xlsx",
  rowNames = FALSE
)

# MR LASSO estimates (Table S6)
Table_S6 = data %>%
  filter(Method == "MR LASSO") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue","NSNPs","LASSO_N_valid","LASSO_valid_SNPs"))

temp = data_MVMR %>%
  filter(Method == "MVMR LASSO") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue","NSNPs","LASSO_N_valid","LASSO_valid_SNPs"))

Table_S6 = rbind.data.frame(Table_S6,temp) %>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S6,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S6.xlsx",
  rowNames = FALSE
)

# MR cML estimates (Table S7)
Table_S7 = data %>%
  filter(Method == "MR cML") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue","NSNPs","BIC_N_invalid"))

temp = data_MVMR %>%
  filter(Method == "MVMR cML") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue","NSNPs","cML_BIC_invalid"))

Table_S7 = rbind_fill(Table_S7,temp) %>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S7,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S7.xlsx",
  rowNames = FALSE
)


# MR IVW estimates (Table S8)
Table_S8 = data %>%
  filter(Method == "MR IVW") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue")) 

temp = data_MVMR %>%
  filter(Method == "MVMR IVW") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue"))

Table_S8 = rbind.data.frame(Table_S8,temp) %>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S8,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S8.xlsx",
  rowNames = FALSE
)


# MR Wald ratio estimates (Table S9)
Table_S9 = data %>%
  filter(Method == "MR Wald Ratio") %>%
  mutate(Method = "MR Wald ratio") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue")) %>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S9,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S9.xlsx",
  rowNames = FALSE
)

# MR median estimates (Table S10)
Table_S10 = data %>%
  filter(Method == "MR Median") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue")) %>%
  mutate(Method = "MR median")

temp = data_MVMR %>%
  filter(Method == "MVMR median") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue"))

Table_S10 = rbind.data.frame(Table_S10,temp) %>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S10,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S10.xlsx",
  rowNames = FALSE
)


# MVMR GMM estimates (Table S11)
Table_S11 = data_MVMR %>%
  filter(Method == "MVMR GMM") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","StdError","CILower","CIUpper","Pvalue"))%>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S11,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S11.xlsx",
  rowNames = FALSE
)

# MVMR Q-statistic minimization estimates (Table S12)
Table_S12 = data_MVMR %>%
  filter(Method == "Q-statistic minimization") %>%
  select(c("Instrument","Exposure","Outcome","Method","Estimate","CILower","CIUpper"))%>%
  arrange(Outcome,Exposure)

openxlsx::write.xlsx(
  Table_S12,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Tables/Table_S12.xlsx",
  rowNames = FALSE
)

