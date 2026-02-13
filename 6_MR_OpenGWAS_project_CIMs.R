library(data.table)
library(R.utils)
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)
library(MRPRESSO)

# Functions for MR analyses

## rbind.data.frame() where data.frames have different columns

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

## Read VEP annotation file

read_VEP <- function(x) {
  VEP_text <- readLines(x)
  VEP_text <- VEP_text[!startsWith(VEP_text, "##")]
  df <- fread(text = VEP_text, sep = "\t", na.strings = "-", data.table = TRUE)
  if (startsWith(names(df)[1], "#")) {
    setnames(df, old = names(df)[1], new = sub("^#", "", names(df)[1]))
  }
  
  return(df)
}

## MR IVW
run_MR_ivw <- function(MR_object, TSMR_object) {
  
  ### IVW (MendelianRandomization) 
  res_ivw <- MendelianRandomization::mr_ivw(MR_object, model = "random")
  
  ivw_df <- data.frame(
    Exposure    = MR_object@exposure[1],
    Outcome     = MR_object@outcome[1],
    SNP         = paste0(MR_object@exposure[1], " instrument") ,
    Estimate    = res_ivw@Estimate[1],
    StdError    = res_ivw@StdError[1],
    CILower     = res_ivw@CILower[1],
    CIUpper     = res_ivw@CIUpper[1],
    Pvalue      = res_ivw@Pvalue[1],
    Fstatistic  = res_ivw@Fstat[1],
    NSNPs       = res_ivw@SNPs[1],
    Method      = "MR IVW",
    stringsAsFactors = FALSE
  )
  
  ### Single SNP Wald Ratio MR (TwoSampleMR)
  res_single <- TwoSampleMR::mr_singlesnp(TSMR_object)
  
  single_df <- res_single %>%
    dplyr::transmute(
      Exposure    = exposure,
      Outcome     = outcome,
      SNP         = SNP,
      Estimate    = b,
      StdError    = se,
      CILower     = b - 1.96 * se,
      CIUpper     = b + 1.96 * se,
      Pvalue      = p,
      Fstatistic  = NA_real_,
      NSNPs       = 1,
      Method      = "MR Wald Ratio"
    ) %>%
    filter(!SNP %in% c("All - Inverse variance weighted",
                       "All - MR Egger") )
  
  ## ---- Combine IVW + Single SNP ----
  dplyr::bind_rows(ivw_df, single_df)
}

## MR Egger (if more than 2 variants; MendelianRandomization)
run_MR_Egger <- function(MR_object) {
  
  # Run MR-Egger
  res <- MendelianRandomization::mr_egger(MR_object)
  
  data.frame(
    Exposure     = MR_object@exposure[1],
    Outcome      = MR_object@outcome[1],
    SNP         = paste0(MR_object@exposure[1], " instrument") ,
    Estimate     = res@Estimate[1],
    StdError     = res@StdError.Est[1],
    CILower      = res@CILower.Est[1],
    CIUpper      = res@CIUpper.Est[1],
    Pvalue   = res@Pvalue.Est[1],  
    Pvalue_MR_Egger_Int   = res@Pvalue.Int[1],  
    I_squared   = res@I.sq[1],
    NSNPs        = res@SNPs[1],
    Method       = "MR-Egger",
    stringsAsFactors = FALSE
  )
  
}

## MR Median
run_MR_median <- function(MR_object) {
  
  ### Median (MendelianRandomization) 
  res_median <- MendelianRandomization::mr_median(MR_object)
  
  median_df <- data.frame(
    Exposure    = MR_object@exposure[1],
    Outcome     = MR_object@outcome[1],
    SNP         = paste0(MR_object@exposure[1], " instrument") ,
    Estimate    = res_median@Estimate[1],
    StdError    = res_median@StdError[1],
    CILower     = res_median@CILower[1],
    CIUpper     = res_median@CIUpper[1],
    Pvalue      = res_median@Pvalue[1],
    NSNPs       = res_median@SNPs[1],
    Method      = "MR Median",
    stringsAsFactors = FALSE
  )
}

## MR LASSO
run_MR_LASSO <- function(MR_object) {
  
  ### LASSO (MendelianRandomization) 
  res_LASSO <- MendelianRandomization::mr_lasso(MR_object)
  
  LASSO_df <- data.frame(
    Exposure    = MR_object@exposure[1],
    Outcome     = MR_object@outcome[1],
    SNP         = paste0(MR_object@exposure[1], " instrument") ,
    Estimate    = res_LASSO@Estimate[1],
    StdError    = res_LASSO@StdError[1],
    CILower     = res_LASSO@CILower[1],
    CIUpper     = res_LASSO@CIUpper[1],
    Pvalue      = res_LASSO@Pvalue[1],
    NSNPs       = res_LASSO@SNPs[1],
    LASSO_N_valid     = res_LASSO@Valid,
    LASSO_valid_SNPs  =  paste(res_LASSO@ValidSNPs, collapse = ","),
    Method      = "MR LASSO",
    stringsAsFactors = FALSE
  )
}

## MR cML
run_MR_cML <- function(MR_object, n_out, n_snp_valid) {
  
  res_cML <- MendelianRandomization::mr_cML(MR_object,
                                            n = n_out,
                                            K_vec = 0:n_snp_valid)
  
  cML_df <- data.frame(
    Exposure    = MR_object@exposure[1],
    Outcome     = MR_object@outcome[1],
    SNP         = paste0(MR_object@exposure[1], " instrument") ,
    Estimate    = res_cML@Estimate[1],
    StdError    = res_cML@StdError[1],
    CILower     = res_cML@CILower[1],
    CIUpper     = res_cML@CIUpper[1],
    Pvalue      = res_cML@Pvalue[1],
    NSNPs       = res_cML@SNPs[1],
    BIC_N_invalid = res_cML@BIC_invalid[1],
    Method      = "MR cML",
    stringsAsFactors = FALSE
  )
}

## MR-PRESSO
run_MR_PRESSO <- function(MR_object_TS) {
  
  res_MRPRESSO <- MRPRESSO::mr_presso(BetaOutcome = "beta.outcome",
                                 BetaExposure = "beta.exposure",
                                 SdOutcome = "se.outcome",
                                 SdExposure = "se.exposure",
                                 data = MR_object_TS,
                                 OUTLIERtest = TRUE,
                                 DISTORTIONtest = TRUE,
                                 SignifThreshold = 0.05,
                                 NbDistribution = 10000,
                                 seed = 40292)
  
  MR_PRESSO_df <- data.frame(
    Exposure    = rep(MR_object_TS$exposure,3),
    Outcome     = rep(MR_object_TS$outcome,3),
    SNP         = rep(paste0(MR_object_TS$exposure, " instrument"),3),
    Estimate    = c(res_MRPRESSO$`Main MR results`$`Causal Estimate`[2], NA, NA),
    StdError    = c(res_MRPRESSO$`Main MR results`$Sd[2],NA,NA),
    Pvalue      = c(res_MRPRESSO$`Main MR results`$`P-value`[2] %>% as.numeric() %>% round(3) %>% as.character(),
                    res_MRPRESSO$`MR-PRESSO results`$`Global Test`$Pvalue,
                    res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$Pvalue %>% as.numeric %>% round(3) %>% as.character()),
    NSNPs       = rep(MR_object_TS$SNP %>% length(),3),
    N_MR_PRESSO_outliers = rep(length(res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`), 3),
    MR_PRESSO_outliers = rep(paste(MR_object_TS$SNP[res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`], collapse = ","), 3) ,
    Method      = c("MR-PRESSO outlier-corrected",
                    "MR-PRESSO global test",
                    "MR-PRESSO distortion test"),
    stringsAsFactors = FALSE
  )
}

# Run univeriable MR analyses

dir_path <- "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/"

file_names <- c("PMD_IL6_harmonized_data.csv",
                "DA_IL6_harmonized_data.csv",
                "NDA_IL6_harmonized_data.csv",
                "PMD_IL8_harmonized_data.csv",
                "DA_IL8_harmonized_data.csv",
                "NDA_IL8_harmonized_data.csv",
                "PMD_VEGF_harmonized_data.csv",
                "DA_VEGF_harmonized_data.csv",
                "NDA_VEGF_harmonized_data.csv",
                "PMD_TNFR2_harmonized_data.csv",
                "DA_TNFR2_harmonized_data.csv",
                "NDA_TNFR2_harmonized_data.csv")

proteins <- c("IL6","IL8","VEGF","TNFR2","CRP")

MR_results = NULL

for (i in 1:length(proteins)) {
  if(proteins[i] != "CRP"){
    ## Exclude non-common SNPs
    VEP_data = read_VEP(paste0(dir_path,proteins[i],"_VEP.tsv")) %>%
      filter(EUR_AF > 0.05) %>%
      select(Uploaded_variation) %>%
      pull()
    
    ## Load and process harmonized data
    PMD = read.csv(paste0(dir_path,file_names[i + ((i-1)*2)])) %>%
      mutate(exposure = proteins[i], outcome = "PMD") %>%
      filter(SNP %in% VEP_data)
    
    DA = read.csv(paste0(dir_path,file_names[i + 1 + ((i-1)*2)])) %>%
      mutate(exposure = proteins[i], outcome = "DA") %>%
      filter(SNP %in% VEP_data)
    
    NDA = read.csv(paste0(dir_path,file_names[i + 2 + ((i-1)*2)])) %>%
      mutate(exposure = proteins[i], outcome = "NDA") %>%
      filter(SNP %in% VEP_data)
    
    
    
    PMD_MR_input = mr_input(bx = PMD$beta.exposure,
                            bxse = PMD$se.exposure,
                            by = PMD$beta.outcome,
                            byse = PMD$se.outcome,
                            exposure = PMD$exposure,
                            outcome = PMD$outcome,
                            snps = PMD$SNP)
    
    PMD_MR_input_TS = PMD
    
    DA_MR_input = mr_input(bx = DA$beta.exposure,
                           bxse = DA$se.exposure,
                           by = DA$beta.outcome,
                           byse = DA$se.outcome,
                           exposure = DA$exposure,
                           outcome = DA$outcome,
                           snps = DA$SNP)
    
    DA_MR_input_TS = DA
    
    NDA_MR_input = mr_input(bx = NDA$beta.exposure,
                            bxse = NDA$se.exposure,
                            by = NDA$beta.outcome,
                            byse = NDA$se.outcome,
                            exposure = NDA$exposure,
                            outcome = NDA$outcome,
                            snps = NDA$SNP)
    
    NDA_MR_input_TS = NDA
    
    ## Run IVW_MR analyses 
    if(nrow(PMD)< 3){
      PMD_IVW = run_MR_ivw(PMD_MR_input,PMD_MR_input_TS)
      DA_IVW = run_MR_ivw(DA_MR_input, DA_MR_input_TS)
      NDA_IVW = run_MR_ivw(NDA_MR_input, NDA_MR_input_TS)
      exposure_MR_results = rbind.data.frame(PMD_IVW,DA_IVW,NDA_IVW)
    }else{
      PMD_IVW = run_MR_ivw(PMD_MR_input,PMD_MR_input_TS)
      DA_IVW = run_MR_ivw(DA_MR_input, DA_MR_input_TS)
      NDA_IVW = run_MR_ivw(NDA_MR_input, NDA_MR_input_TS)
      
      PMD_Egger = run_MR_Egger(PMD_MR_input)
      DA_Egger = run_MR_Egger(DA_MR_input)
      NDA_Egger = run_MR_Egger(NDA_MR_input)
      
      PMD_Median = run_MR_median(PMD_MR_input)
      DA_Median = run_MR_median(DA_MR_input)
      NDA_Median = run_MR_median(NDA_MR_input)
      
      PMD_LASSO = run_MR_LASSO(PMD_MR_input)
      DA_LASSO = run_MR_LASSO(DA_MR_input)
      NDA_LASSO = run_MR_LASSO(NDA_MR_input)
      
      PMD_cML = run_MR_cML(PMD_MR_input, n_out = 17841, n_snp_valid = PMD_LASSO$LASSO_N_valid)
      DA_cML = run_MR_cML(DA_MR_input, n_out = 13965, n_snp_valid = DA_LASSO$LASSO_N_valid)
      NDA_cML = run_MR_cML(NDA_MR_input, n_out = 14036, n_snp_valid = NDA_LASSO$LASSO_N_valid)
      
      exposure_MR_results = rbind_fill(PMD_IVW, PMD_Egger, PMD_Median, PMD_LASSO, PMD_cML,
                                       DA_IVW, DA_Egger, DA_Median, DA_LASSO, DA_cML,
                                       NDA_IVW, NDA_Egger, NDA_Median, NDA_LASSO, NDA_cML
      )
      
    }
    
    MR_results = rbind_fill(MR_results,exposure_MR_results)
    

  }else{
    ## CRP instruments
    
    ### Exclude non-common SNPs
    VEP_data = read_VEP("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/CPR_VEP.tsv") %>%
      filter(EUR_AF > 0.05) %>%
      select(Uploaded_variation) %>%
      pull()
    
    ### Load and process harmonized data
    PMD = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_CRP_noBMIadj_harmonized_data.rds") %>%
      mutate(exposure = "CRP", outcome = "PMD") %>%
      filter(SNP %in% VEP_data)
    
    DA = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/DA_CRP_noBMIadj_harmonized_data.rds") %>%
      mutate(exposure = "CRP", outcome = "DA") %>%
      filter(SNP %in% VEP_data)
    
    NDA = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/NDA_CRP_noBMIadj_harmonized_data.rds") %>%
      mutate(exposure = "CRP", outcome = "NDA") %>%
      filter(SNP %in% VEP_data)
    
    PMD_MR_input = mr_input(bx = PMD$beta.exposure,
                            bxse = PMD$se.exposure,
                            by = PMD$beta.outcome,
                            byse = PMD$se.outcome,
                            exposure = PMD$exposure,
                            outcome = PMD$outcome,
                            snps = PMD$SNP)
    
    PMD_MR_input_TS = PMD
    
    DA_MR_input = mr_input(bx = DA$beta.exposure,
                           bxse = DA$se.exposure,
                           by = DA$beta.outcome,
                           byse = DA$se.outcome,
                           exposure = DA$exposure,
                           outcome = DA$outcome,
                           snps = DA$SNP)
    
    DA_MR_input_TS = DA
    
    NDA_MR_input = mr_input(bx = NDA$beta.exposure,
                            bxse = NDA$se.exposure,
                            by = NDA$beta.outcome,
                            byse = NDA$se.outcome,
                            exposure = NDA$exposure,
                            outcome = NDA$outcome,
                            snps = NDA$SNP)
    
    NDA_MR_input_TS = NDA
    
    ### Create subset CRP MR instruments
    
    #### Define CRP instrument with no BMI-associated SNPs 
    PMD_CRP_BMI = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_CRP_and_BMI_harmonized_data.rds")
    BMI_associated_SNPs = PMD_CRP_BMI$exposure_pval[,2] < 0.05
    BMI_CRP_SNPs = dimnames(PMD_CRP_BMI$exposure_pval)[[1]]
    BMI_associated_SNPs = BMI_CRP_SNPs[BMI_associated_SNPs]
    
    # Other BMI associated SNPs identified by checking LDtrait #
    other_BMI_associated_SNPs <- c("rs10498240","rs10498240","rs10797119","rs10871777",
                                   "rs112875651","rs1222209", "rs1260326", "rs12716972",
                                   "rs12929503","rs13066686","rs13294945","rs13416992",
                                   "rs1905505","rs2726033","rs28429148","rs288183",
                                   "rs34092660","rs2280406","rs197422","rs35695283",
                                   "rs3768321","rs3935032","rs429358","rs4655802",
                                   "rs630372","rs7794705","rs78428995","rs994596",
                                   "rs73137144","rs7103411","rs6992315")
    
    BMI_associated_SNPs = c(BMI_associated_SNPs,other_BMI_associated_SNPs)
    
    PMD_noBMI = PMD %>% 
      filter(!SNP %in% BMI_associated_SNPs)
    
    PMD_MR_input_noBMI = mr_input(bx = PMD_noBMI$beta.exposure,
                            bxse = PMD_noBMI$se.exposure,
                            by = PMD_noBMI$beta.outcome,
                            byse = PMD_noBMI$se.outcome,
                            exposure = "CRP_noBMI",
                            outcome = PMD_noBMI$outcome,
                            snps = PMD_noBMI$SNP)
    
    PMD_MR_input_TS_noBMI = PMD_noBMI %>%
      mutate(exposure = "CRP_noBMI")
    
    DA_noBMI = DA %>% 
      filter(!SNP %in% BMI_associated_SNPs)
    
    DA_MR_input_noBMI = mr_input(bx = DA_noBMI$beta.exposure,
                                  bxse = DA_noBMI$se.exposure,
                                  by = DA_noBMI$beta.outcome,
                                  byse = DA_noBMI$se.outcome,
                                  exposure = "CRP_noBMI",
                                  outcome = DA_noBMI$outcome,
                                  snps = DA_noBMI$SNP)
    
    DA_MR_input_TS_noBMI = DA_noBMI %>%
      mutate(exposure = "CRP_noBMI")
    
    NDA_noBMI = NDA %>% 
      filter(!SNP %in% BMI_associated_SNPs)
    
    NDA_MR_input_noBMI = mr_input(bx = NDA_noBMI$beta.exposure,
                                  bxse = NDA_noBMI$se.exposure,
                                  by = NDA_noBMI$beta.outcome,
                                  byse = NDA_noBMI$se.outcome,
                                  exposure = "CRP_noBMI",
                                  outcome = NDA_noBMI$outcome,
                                  snps = NDA_noBMI$SNP)
    
    NDA_MR_input_TS_noBMI = NDA_noBMI %>%
      mutate(exposure = "CRP_noBMI")
    
    
    #### Define CRP instrument restricted to the CRP locus 
    
    CRP_locus_SNPs = PMD %>% 
      filter(!SNP %in% BMI_associated_SNPs) %>%
      filter(chr.exposure == 1 & pos.exposure >= (159682079-500000) &  pos.exposure <= (159684379+500000)) %>%
      pull(SNP)
    
    PMD_CRP_locus = PMD %>% 
      filter(SNP %in% CRP_locus_SNPs)
    
    PMD_MR_input_CRP_locus = mr_input(bx = PMD_CRP_locus$beta.exposure,
                                  bxse = PMD_CRP_locus$se.exposure,
                                  by = PMD_CRP_locus$beta.outcome,
                                  byse = PMD_CRP_locus$se.outcome,
                                  exposure = "CRP_locus",
                                  outcome = PMD_CRP_locus$outcome,
                                  snps = PMD_CRP_locus$SNP)
    
    PMD_MR_input_TS_CRP_locus = PMD_CRP_locus %>%
      mutate(exposure = "CRP_locus")
    
    DA_CRP_locus = DA %>% 
      filter(SNP %in% CRP_locus_SNPs)
    
    DA_MR_input_CRP_locus = mr_input(bx = DA_CRP_locus$beta.exposure,
                                 bxse = DA_CRP_locus$se.exposure,
                                 by = DA_CRP_locus$beta.outcome,
                                 byse = DA_CRP_locus$se.outcome,
                                 exposure = "CRP_locus",
                                 outcome = DA_CRP_locus$outcome,
                                 snps = DA_CRP_locus$SNP)
    
    DA_MR_input_TS_CRP_locus = DA_CRP_locus %>%
      mutate(exposure = "CRP_locus")
    
    NDA_CRP_locus = NDA %>% 
      filter(SNP %in% CRP_locus_SNPs)
    
    NDA_MR_input_CRP_locus = mr_input(bx = NDA_CRP_locus$beta.exposure,
                                  bxse = NDA_CRP_locus$se.exposure,
                                  by = NDA_CRP_locus$beta.outcome,
                                  byse = NDA_CRP_locus$se.outcome,
                                  exposure = "CRP_locus",
                                  outcome = NDA_CRP_locus$outcome,
                                  snps = NDA_CRP_locus$SNP)
    
    NDA_MR_input_TS_CRP_locus = NDA_CRP_locus %>%
      mutate(exposure = "CRP_locus")
    
    # Run CRP MR analyses    
    ## Full set of CRP-associated SNPs
    
    PMD_IVW = run_MR_ivw(PMD_MR_input,PMD_MR_input_TS)
    DA_IVW = run_MR_ivw(DA_MR_input, DA_MR_input_TS)
    NDA_IVW = run_MR_ivw(NDA_MR_input, NDA_MR_input_TS)
    
    PMD_Egger = run_MR_Egger(PMD_MR_input)
    DA_Egger = run_MR_Egger(DA_MR_input)
    NDA_Egger = run_MR_Egger(NDA_MR_input)
    
    PMD_Median = run_MR_median(PMD_MR_input)
    DA_Median = run_MR_median(DA_MR_input)
    NDA_Median = run_MR_median(NDA_MR_input)
    
    PMD_LASSO = run_MR_LASSO(PMD_MR_input)
    DA_LASSO = run_MR_LASSO(DA_MR_input)
    NDA_LASSO = run_MR_LASSO(NDA_MR_input)
    
    PMD_cML = run_MR_cML(PMD_MR_input, n_out = 17841, n_snp_valid = PMD_LASSO$LASSO_N_valid)
    DA_cML = run_MR_cML(DA_MR_input, n_out = 13965, n_snp_valid = DA_LASSO$LASSO_N_valid)
    NDA_cML = run_MR_cML(NDA_MR_input, n_out = 14036, n_snp_valid = NDA_LASSO$LASSO_N_valid)
    
    PMD_MR_PRESSO = run_MR_PRESSO(PMD_MR_input_TS)
    DA_MR_PRESSO = run_MR_PRESSO(DA_MR_input_TS)
    NDA_MR_PRESSO = run_MR_PRESSO(NDA_MR_input_TS)
    
    full_CRP_MR_results = rbind_fill(PMD_IVW, PMD_Egger, PMD_Median, PMD_LASSO, PMD_cML, PMD_MR_PRESSO,
                                     DA_IVW, DA_Egger, DA_Median, DA_LASSO, DA_cML, DA_MR_PRESSO,
                                     NDA_IVW, NDA_Egger, NDA_Median, NDA_LASSO, NDA_cML, NDA_MR_PRESSO)
    
    ## CRP-associated SNPs excluding BMI-associated SNPs
    
    PMD_IVW = run_MR_ivw(PMD_MR_input_noBMI,PMD_MR_input_TS_noBMI)
    DA_IVW = run_MR_ivw(DA_MR_input_noBMI, DA_MR_input_TS_noBMI)
    NDA_IVW = run_MR_ivw(NDA_MR_input_noBMI, NDA_MR_input_TS_noBMI)
    
    PMD_Egger = run_MR_Egger(PMD_MR_input_noBMI)
    DA_Egger = run_MR_Egger(DA_MR_input_noBMI)
    NDA_Egger = run_MR_Egger(NDA_MR_input_noBMI)
    
    PMD_Median = run_MR_median(PMD_MR_input_noBMI)
    DA_Median = run_MR_median(DA_MR_input_noBMI)
    NDA_Median = run_MR_median(NDA_MR_input_noBMI)
    
    PMD_LASSO = run_MR_LASSO(PMD_MR_input_noBMI)
    DA_LASSO = run_MR_LASSO(DA_MR_input_noBMI)
    NDA_LASSO = run_MR_LASSO(NDA_MR_input_noBMI)
    
    PMD_cML = run_MR_cML(PMD_MR_input_noBMI, n_out = 17841, n_snp_valid = PMD_LASSO$LASSO_N_valid)
    DA_cML = run_MR_cML(DA_MR_input_noBMI, n_out = 13965, n_snp_valid = DA_LASSO$LASSO_N_valid)
    NDA_cML = run_MR_cML(NDA_MR_input_noBMI, n_out = 14036, n_snp_valid = NDA_LASSO$LASSO_N_valid)
    
    PMD_MR_PRESSO = run_MR_PRESSO(PMD_MR_input_TS_noBMI)
    DA_MR_PRESSO = run_MR_PRESSO(DA_MR_input_TS_noBMI)
    NDA_MR_PRESSO = run_MR_PRESSO(NDA_MR_input_TS_noBMI)
    
    noBMI_CRP_MR_results = rbind_fill(PMD_IVW, PMD_Egger, PMD_Median, PMD_LASSO, PMD_cML, PMD_MR_PRESSO,
                                     DA_IVW, DA_Egger, DA_Median, DA_LASSO, DA_cML, DA_MR_PRESSO,
                                     NDA_IVW, NDA_Egger, NDA_Median, NDA_LASSO, NDA_cML, NDA_MR_PRESSO)
    
    ## CRP-locus SNP 
    
    PMD_IVW = run_MR_ivw(PMD_MR_input_CRP_locus,PMD_MR_input_TS_CRP_locus)
    DA_IVW = run_MR_ivw(DA_MR_input_CRP_locus, DA_MR_input_TS_CRP_locus)
    NDA_IVW = run_MR_ivw(NDA_MR_input_CRP_locus, NDA_MR_input_TS_CRP_locus)
    

   CRP_locus_MR_results = rbind_fill(PMD_IVW,
                                      DA_IVW,
                                      NDA_IVW)
    
    
  }
  
}

# Merge

MR_results = rbind_fill(MR_results,full_CRP_MR_results,noBMI_CRP_MR_results,CRP_locus_MR_results)


# Remove extra rows
MR_results = MR_results[-1,]

MR_results <- MR_results[!duplicated(MR_results), ]

# Save results

fwrite(x = MR_results, file = "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_project_CIM_MR_results.tsv",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
)

xlsx::write.xlsx(MR_results, "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_project_CIM_MR_results.xlsx", row.names = FALSE,
                 showNA = FALSE)




