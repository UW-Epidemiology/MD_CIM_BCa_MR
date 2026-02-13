library(data.table)
library(R.utils)
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)
library(MRPRESSO)
library(MVMR)

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

# Load instrument data

PMD_CRP_BMI = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/PMD_CRP_and_BMI_harmonized_data.rds")
DA_CRP_BMI = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/DA_CRP_and_BMI_harmonized_data.rds")
NDA_CRP_BMI = readRDS("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/NDA_CRP_and_BMI_harmonized_data.rds")

# Exclude SNPs with low MAF

CRP_and_BMI_VEP = read_VEP("/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/instrument_data/CRP_and_BMI_VEP.tsv")
CRP_and_BMI_common = CRP_and_BMI_VEP %>%
  filter(EUR_AF > 0.05) %>%
  select(Uploaded_variation) %>%
  pull()

filter_vec = dimnames(PMD_CRP_BMI$exposure_beta)[1] %>% unlist() %in% CRP_and_BMI_common

PMD_CRP_BMI[1:6] <- lapply(seq_along(PMD_CRP_BMI[1:6]), function(i) {
  x <- PMD_CRP_BMI[[i]]
  if (is.matrix(x) || is.data.frame(x)) {
    filtered_x <- x[filter_vec, , drop = FALSE]
    if (!is.null(rownames(x))) {
      rownames(filtered_x) <- rownames(x)[filter_vec]
    }
    return(filtered_x)
  } else if (is.vector(x)) {
    return(x[filter_vec])  
  } else {
    return(x)  
  }
})

DA_CRP_BMI[1:6] <- lapply(seq_along(DA_CRP_BMI[1:6]), function(i) {
  x <- DA_CRP_BMI[[i]]
  if (is.matrix(x) || is.data.frame(x)) {
    filtered_x <- x[filter_vec, , drop = FALSE]
    if (!is.null(rownames(x))) {
      rownames(filtered_x) <- rownames(x)[filter_vec]
    }
    return(filtered_x)
  } else if (is.vector(x)) {
    return(x[filter_vec])  
  } else {
    return(x)  
  }
})

NDA_CRP_BMI[1:6] <- lapply(seq_along(NDA_CRP_BMI[1:6]), function(i) {
  x <- NDA_CRP_BMI[[i]]
  if (is.matrix(x) || is.data.frame(x)) {
    filtered_x <- x[filter_vec, , drop = FALSE]
    if (!is.null(rownames(x))) {
      rownames(filtered_x) <- rownames(x)[filter_vec]
    }
    return(filtered_x)
  } else if (is.vector(x)) {
    return(x[filter_vec])  
  } else {
    return(x)  
  }
})

## Prepare MVMR format dataset for Mendelian Randomization package

exp_sample_sizes = c(575531,461460)
out_sample_sizes = c(17841,13965,14036)

PMD_MVMR_input = mr_mvinput(bx = PMD_CRP_BMI$exposure_beta,
                            bxse = PMD_CRP_BMI$exposure_se,
                            by = PMD_CRP_BMI$outcome_beta,
                            byse = PMD_CRP_BMI$outcome_se,
                            exposure = c("CRP","BMI"),
                            outcome = "PMD",
                            snps = dimnames(PMD_CRP_BMI$exposure_se)[[1]])

DA_MVMR_input = mr_mvinput(bx = DA_CRP_BMI$exposure_beta,
                           bxse = DA_CRP_BMI$exposure_se,
                           by = DA_CRP_BMI$outcome_beta,
                           byse = DA_CRP_BMI$outcome_se,
                           exposure = c("CRP","BMI"),
                           outcome = "DA",
                           snps = dimnames(DA_CRP_BMI$exposure_se)[[1]])

NDA_MVMR_input = mr_mvinput(bx = NDA_CRP_BMI$exposure_beta,
                            bxse = NDA_CRP_BMI$exposure_se,
                            by = NDA_CRP_BMI$outcome_beta,
                            byse = NDA_CRP_BMI$outcome_se,
                            exposure = c("CRP","BMI"),
                            outcome = "NDA",
                            snps = dimnames(NDA_CRP_BMI$exposure_se)[[1]])

## Prepare MVMR format for the MVMR package
### Correlation of 0.43 is from overall correlation between CRP and BMI in UK Biobank

CRP_BMI_covar = phenocov_mvmr(matrix(c(1,0.43,0.43,1), nrow = 2, ncol = 2),
                              PMD_CRP_BMI$exposure_se)

PMD_MVMR_sensitivity_input = format_mvmr(BXGs = PMD_CRP_BMI$exposure_beta,
                                         seBXGs = PMD_CRP_BMI$exposure_se,
                                         BYG = PMD_CRP_BMI$outcome_beta,
                                         seBYG = PMD_CRP_BMI$outcome_se,
                                         RSID = dimnames(PMD_CRP_BMI$exposure_beta)[1])

DA_MVMR_sensitivity_input = format_mvmr(BXGs = DA_CRP_BMI$exposure_beta,
                                        seBXGs = DA_CRP_BMI$exposure_se,
                                        BYG = DA_CRP_BMI$outcome_beta,
                                        seBYG = DA_CRP_BMI$outcome_se,
                                        RSID = dimnames(DA_CRP_BMI$exposure_beta)[1])

NDA_MVMR_sensitivity_input = format_mvmr(BXGs = NDA_CRP_BMI$exposure_beta,
                                         seBXGs = NDA_CRP_BMI$exposure_se,
                                         BYG = NDA_CRP_BMI$outcome_beta,
                                         seBYG = NDA_CRP_BMI$outcome_se,
                                         RSID = dimnames(NDA_CRP_BMI$exposure_beta)[1])


## MendelianRandomization MVMR IVW
PMD_IVW <- mr_mvivw(PMD_MVMR_input, robust = TRUE,
                    nx = exp_sample_sizes,
                    correl.x = matrix(c(1,0.43,0.43,1), nrow = 2, ncol = 2))
DA_IVW <- mr_mvivw(DA_MVMR_input, robust = TRUE,
                    nx = exp_sample_sizes,
                    correl.x = matrix(c(1,0.43,0.43,1), nrow = 2, ncol = 2))
NDA_IVW <- mr_mvivw(NDA_MVMR_input, robust = TRUE,
                    nx = exp_sample_sizes,
                    correl.x = matrix(c(1,0.43,0.43,1), nrow = 2, ncol = 2))

IVW_data = data.frame(Exposure = c(PMD_IVW@Exposure, DA_IVW@Exposure, NDA_IVW@Exposure),
                      Outcome = c("PMD","PMD","DA","DA","NDA","NDA"),
                      SNP = rep("CRP_BMI_instrument",6),
                      Estimate = c(PMD_IVW@Estimate, DA_IVW@Estimate, NDA_IVW@Estimate),
                      StdError = c(PMD_IVW@StdError, DA_IVW@StdError, NDA_IVW@StdError),
                      CILower = c(PMD_IVW@CILower, DA_IVW@CILower, NDA_IVW@CILower),
                      CIUpper = c(PMD_IVW@CIUpper, DA_IVW@CIUpper, NDA_IVW@CIUpper),
                      Pvalue = c(PMD_IVW@Pvalue, DA_IVW@Pvalue, NDA_IVW@Pvalue),
                      Fstatistic = c(PMD_IVW@CondFstat, DA_IVW@CondFstat, NDA_IVW@CondFstat),
                      NSNPs = c(PMD_IVW@SNPs, DA_IVW@SNPs, NDA_IVW@SNPs),
                      Method = rep("MVMR IVW",6),
                      Hetero_statistic = c(PMD_IVW@Heter.Stat[1], PMD_IVW@Heter.Stat[1], DA_IVW@Heter.Stat[1], DA_IVW@Heter.Stat[1], NDA_IVW@Heter.Stat[1], NDA_IVW@Heter.Stat[1]),
                      Hetero_Pvalue = c(PMD_IVW@Heter.Stat[2], PMD_IVW@Heter.Stat[2], DA_IVW@Heter.Stat[2], DA_IVW@Heter.Stat[2], NDA_IVW@Heter.Stat[2], NDA_IVW@Heter.Stat[2]))

## MendelianRandomization MVMR MR-Egger
PMD_Egger <- mr_mvegger(PMD_MVMR_input)
DA_Egger <- mr_mvegger(DA_MVMR_input)
NDA_Egger <- mr_mvegger(NDA_MVMR_input)

Egger_data = data.frame(Exposure = c(PMD_Egger@Exposure, DA_Egger@Exposure, NDA_Egger@Exposure),
                      Outcome = c("PMD","PMD","DA","DA","NDA","NDA"),
                      SNP = rep("CRP_BMI_instrument",6),
                      Estimate = c(PMD_Egger@Estimate, DA_Egger@Estimate, NDA_Egger@Estimate),
                      StdError = c(PMD_Egger@StdError.Est, DA_Egger@StdError.Est, NDA_Egger@StdError.Est),
                      CILower = c(PMD_Egger@CILower.Est, DA_Egger@CILower.Est, NDA_Egger@CILower.Est),
                      CIUpper = c(PMD_Egger@CIUpper.Est, DA_Egger@CIUpper.Est, NDA_Egger@CIUpper.Est),
                      Pvalue = c(PMD_Egger@Pvalue.Est, DA_Egger@Pvalue.Est, NDA_Egger@Pvalue.Est),
                      NSNPs = c(PMD_Egger@SNPs, DA_Egger@SNPs, NDA_Egger@SNPs),
                      Method = rep("MVMR-Egger",6),
                      Pvalue_MVMR_Egger_Int = c(PMD_Egger@Pvalue.Int, DA_Egger@Pvalue.Int, NDA_Egger@Pvalue.Int),
                      Hetero_statistic = c(PMD_Egger@Heter.Stat[1], PMD_Egger@Heter.Stat[1], DA_Egger@Heter.Stat[1], DA_Egger@Heter.Stat[1], NDA_Egger@Heter.Stat[1], NDA_Egger@Heter.Stat[1]),
                      Hetero_Pvalue = c(PMD_Egger@Heter.Stat[2], PMD_Egger@Heter.Stat[2], DA_Egger@Heter.Stat[2], DA_Egger@Heter.Stat[2], NDA_Egger@Heter.Stat[2], NDA_Egger@Heter.Stat[2]))


## MendelianRandomization MVMR MR-Egger
PMD_Egger <- mr_mvegger(PMD_MVMR_input)
DA_Egger <- mr_mvegger(DA_MVMR_input)
NDA_Egger <- mr_mvegger(NDA_MVMR_input)

Egger_data = data.frame(Exposure = c(PMD_Egger@Exposure, DA_Egger@Exposure, NDA_Egger@Exposure),
                        Outcome = c("PMD","PMD","DA","DA","NDA","NDA"),
                        SNP = rep("CRP_BMI_instrument",6),
                        Estimate = c(PMD_Egger@Estimate, DA_Egger@Estimate, NDA_Egger@Estimate),
                        StdError = c(PMD_Egger@StdError.Est, DA_Egger@StdError.Est, NDA_Egger@StdError.Est),
                        CILower = c(PMD_Egger@CILower.Est, DA_Egger@CILower.Est, NDA_Egger@CILower.Est),
                        CIUpper = c(PMD_Egger@CIUpper.Est, DA_Egger@CIUpper.Est, NDA_Egger@CIUpper.Est),
                        Pvalue = c(PMD_Egger@Pvalue.Est, DA_Egger@Pvalue.Est, NDA_Egger@Pvalue.Est),
                        NSNPs = c(PMD_Egger@SNPs, DA_Egger@SNPs, NDA_Egger@SNPs),
                        Method = rep("MVMR-Egger",6),
                        Pvalue_MVMR_Egger_Int = c(PMD_Egger@Pvalue.Int, DA_Egger@Pvalue.Int, NDA_Egger@Pvalue.Int),
                        Hetero_statistic = c(PMD_Egger@Heter.Stat[1], PMD_Egger@Heter.Stat[1], DA_Egger@Heter.Stat[1], DA_Egger@Heter.Stat[1], NDA_Egger@Heter.Stat[1], NDA_Egger@Heter.Stat[1]),
                        Hetero_Pvalue = c(PMD_Egger@Heter.Stat[2], PMD_Egger@Heter.Stat[2], DA_Egger@Heter.Stat[2], DA_Egger@Heter.Stat[2], NDA_Egger@Heter.Stat[2], NDA_Egger@Heter.Stat[2]))


## MendelianRandomization MVMR weighted median
PMD_median <- mr_mvmedian(PMD_MVMR_input)
DA_median <- mr_mvmedian(DA_MVMR_input)
NDA_median <- mr_mvmedian(PMD_MVMR_input)

Median_data = data.frame(Exposure = c(PMD_median@Exposure, DA_median@Exposure, NDA_median@Exposure),
                        Outcome = c("PMD","PMD","DA","DA","NDA","NDA"),
                        SNP = rep("CRP_BMI_instrument",6),
                        Estimate = c(PMD_median@Estimate, DA_median@Estimate, NDA_median@Estimate),
                        StdError = c(PMD_median@StdError, DA_median@StdError, NDA_median@StdError),
                        CILower = c(PMD_median@CILower, DA_median@CILower, NDA_median@CILower),
                        CIUpper = c(PMD_median@CIUpper, DA_median@CIUpper, NDA_median@CIUpper),
                        Pvalue = c(PMD_median@Pvalue, DA_median@Pvalue, NDA_median@Pvalue),
                        NSNPs = c(PMD_median@SNPs, DA_median@SNPs, NDA_median@SNPs),
                        Method = rep("MVMR median",6))

## MendelianRandomization MVMR LASSO

PMD_LASSO <- mr_mvlasso(PMD_MVMR_input)
DA_LASSO <- mr_mvlasso(DA_MVMR_input)
NDA_LASSO <- mr_mvlasso(NDA_MVMR_input)

LASSO_data = data.frame(Exposure = c(PMD_LASSO@Exposure, DA_LASSO@Exposure, NDA_LASSO@Exposure),
                         Outcome = c("PMD","PMD","DA","DA","NDA","NDA"),
                         SNP = rep("CRP_BMI_instrument",6),
                         Estimate = c(PMD_LASSO@Estimate, DA_LASSO@Estimate, NDA_LASSO@Estimate),
                         StdError = c(PMD_LASSO@StdError, DA_LASSO@StdError, NDA_LASSO@StdError),
                         CILower = c(PMD_LASSO@CILower, DA_LASSO@CILower, NDA_LASSO@CILower),
                         CIUpper = c(PMD_LASSO@CIUpper, DA_LASSO@CIUpper, NDA_LASSO@CIUpper),
                         Pvalue = c(PMD_LASSO@Pvalue, DA_LASSO@Pvalue, NDA_LASSO@Pvalue),
                         NSNPs = c(PMD_LASSO@SNPs, DA_LASSO@SNPs, NDA_LASSO@SNPs),
                         Method = rep("MVMR LASSO",6),
                         LASSO_N_valid = c(PMD_LASSO@Valid, PMD_LASSO@Valid, DA_LASSO@Valid, DA_LASSO@Valid, NDA_LASSO@Valid, NDA_LASSO@Valid),
                         LASSO_valid_SNPs = c(paste(PMD_LASSO@ValidSNPs, collapse = ","),
                                              paste(PMD_LASSO@ValidSNPs, collapse = ","),
                                              paste(DA_LASSO@ValidSNPs, collapse = ","),
                                              paste(DA_LASSO@ValidSNPs, collapse = ","),
                                              paste(NDA_LASSO@ValidSNPs, collapse = ","),
                                              paste(NDA_LASSO@ValidSNPs, collapse = ",")))

## MendelianRandomization MVMR cML
PMD_cML <- mr_mvcML(PMD_MVMR_input, n = out_sample_sizes[1],
                    K_vec = 0:(PMD_LASSO@SNPs-PMD_LASSO@Valid)*2,
                    seed = 32509)
DA_cML <- mr_mvcML(DA_MVMR_input, n = out_sample_sizes[2],
                    K_vec = 0:(DA_LASSO@SNPs-DA_LASSO@Valid)*2,
                   seed = 32509)
NDA_cML <- mr_mvcML(NDA_MVMR_input, n = out_sample_sizes[3],
                    K_vec = 0:(NDA_LASSO@SNPs-NDA_LASSO@Valid)*2,
                    seed = 32509)


cML_data = data.frame(Exposure = c(PMD_cML@Exposure, DA_cML@Exposure, NDA_cML@Exposure),
                        Outcome = c("PMD","PMD","DA","DA","NDA","NDA"),
                        SNP = rep("CRP_BMI_instrument",6),
                        Estimate = c(PMD_cML@Estimate, DA_cML@Estimate, NDA_cML@Estimate),
                        StdError = c(PMD_cML@StdError, DA_cML@StdError, NDA_cML@StdError),
                        CILower = c(PMD_cML@CILower, DA_cML@CILower, NDA_cML@CILower),
                        CIUpper = c(PMD_cML@CIUpper, DA_cML@CIUpper, NDA_cML@CIUpper),
                        Pvalue = c(PMD_cML@Pvalue, DA_cML@Pvalue, NDA_cML@Pvalue),
                        NSNPs = c(PMD_cML@SNPs, DA_cML@SNPs, NDA_cML@SNPs),
                        Method = rep("MVMR cML",6),
                        cML_BIC_invalid = c(PMD_cML@BIC_invalid[1],
                                            PMD_cML@BIC_invalid[2],
                                            DA_cML@BIC_invalid[1],
                                            DA_cML@BIC_invalid[2],
                                            NDA_cML@BIC_invalid[1],
                                            NDA_cML@BIC_invalid[2]))

## MendelianRandomization MVMR GMM 

PMD_GMM <- mr_mvgmm(PMD_MVMR_input, nx = exp_sample_sizes,
                    ny = out_sample_sizes[1],
                    cor.x = matrix(c(1,0.43,0.43,1),
                                   nrow = 2, ncol = 2))
DA_GMM <- mr_mvgmm(DA_MVMR_input, nx = exp_sample_sizes,
                    ny = out_sample_sizes[2],
                    cor.x = matrix(c(1,0.43,0.43,1),
                                   nrow = 2, ncol = 2))
NDA_GMM <- mr_mvgmm(DA_MVMR_input, nx = exp_sample_sizes,
                   ny = out_sample_sizes[3],
                   cor.x = matrix(c(1,0.43,0.43,1),
                                  nrow = 2, ncol = 2))

GMM_data = data.frame(Exposure = c(PMD_GMM@Exposure, DA_GMM@Exposure, NDA_GMM@Exposure),
                      Outcome = c("PMD","PMD","DA","DA","NDA","NDA"),
                      SNP = rep("CRP_BMI_instrument",6),
                      Estimate = c(PMD_GMM@Estimate, DA_GMM@Estimate, NDA_GMM@Estimate),
                      StdError = c(PMD_GMM@StdError, DA_GMM@StdError, NDA_GMM@StdError),
                      CILower = c(PMD_GMM@CILower, DA_GMM@CILower, NDA_GMM@CILower),
                      CIUpper = c(PMD_GMM@CIUpper, DA_GMM@CIUpper, NDA_GMM@CIUpper),
                      Pvalue = c(PMD_GMM@Pvalue, DA_GMM@Pvalue, NDA_GMM@Pvalue),
                      Method = rep("MVMR GMM",6))


## MR-PRESSO MVMR
### Prepare data.frames for MVMR

PMD = data.frame(outcome.beta = PMD_CRP_BMI$outcome_beta,
                 exposure1.beta = PMD_CRP_BMI$exposure_beta[,1],
                 exposure2.beta = PMD_CRP_BMI$exposure_beta[,2],
                 outcome.se = PMD_CRP_BMI$outcome_se,
                 exposure1.se = PMD_CRP_BMI$exposure_se[,1],
                 exposure2.se = PMD_CRP_BMI$exposure_se[,2])

DA = data.frame(outcome.beta = DA_CRP_BMI$outcome_beta,
                 exposure1.beta = DA_CRP_BMI$exposure_beta[,1],
                 exposure2.beta = DA_CRP_BMI$exposure_beta[,2],
                 outcome.se = DA_CRP_BMI$outcome_se,
                 exposure1.se = DA_CRP_BMI$exposure_se[,1],
                 exposure2.se = DA_CRP_BMI$exposure_se[,2])

NDA = data.frame(outcome.beta = NDA_CRP_BMI$outcome_beta,
                 exposure1.beta = NDA_CRP_BMI$exposure_beta[,1],
                 exposure2.beta = NDA_CRP_BMI$exposure_beta[,2],
                 outcome.se = NDA_CRP_BMI$outcome_se,
                 exposure1.se = NDA_CRP_BMI$exposure_se[,1],
                 exposure2.se = NDA_CRP_BMI$exposure_se[,2])


PMD_MRPRESSO <- MRPRESSO::mr_presso(BetaOutcome = "outcome.beta",
                                    BetaExposure = c("exposure1.beta","exposure2.beta"),
                                    SdOutcome = "outcome.se",
                                    SdExposure = c("exposure1.se","exposure2.se"),
                                    data = PMD,
                                    OUTLIERtest = TRUE,
                                    DISTORTIONtest = TRUE,
                                    SignifThreshold = 0.05,
                                    NbDistribution = 10000,
                                    seed = 40292)

DA_MRPRESSO <- MRPRESSO::mr_presso(BetaOutcome = "outcome.beta",
                                    BetaExposure = c("exposure1.beta","exposure2.beta"),
                                    SdOutcome = "outcome.se",
                                    SdExposure = c("exposure1.se","exposure2.se"),
                                    data = DA,
                                    OUTLIERtest = TRUE,
                                    DISTORTIONtest = TRUE,
                                    SignifThreshold = 0.05,
                                    NbDistribution = 10000,
                                    seed = 40292)

NDA_MRPRESSO <- MRPRESSO::mr_presso(BetaOutcome = "outcome.beta",
                                    BetaExposure = c("exposure1.beta","exposure2.beta"),
                                    SdOutcome = "outcome.se",
                                    SdExposure = c("exposure1.se","exposure2.se"),
                                    data = NDA,
                                    OUTLIERtest = TRUE,
                                    DISTORTIONtest = TRUE,
                                    SignifThreshold = 0.05,
                                    NbDistribution = 10000,
                                    seed = 40292)

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    formatC(as.numeric(p), format = "e", digits = 3)
  )
}


make_MRPRESSO_df <- function(MR_object, res_MRPRESSO) {
  
  data.frame(
    Exposure    = rep(MR_object@exposure, 3),
    Outcome     = rep(MR_object@outcome, 6),
    SNP         = rep("CRP_BMI_instrument",6),
    
    Estimate = c(
      res_MRPRESSO$`Main MR results`$`Causal Estimate`[1],
      res_MRPRESSO$`Main MR results`$`Causal Estimate`[2],
      res_MRPRESSO$`Main MR results`$`Causal Estimate`[3],
      res_MRPRESSO$`Main MR results`$`Causal Estimate`[4],
      NA,
      NA
    ),
    
    StdError = c(
      res_MRPRESSO$`Main MR results`$Sd[1],
      res_MRPRESSO$`Main MR results`$`Sd`[2],
      res_MRPRESSO$`Main MR results`$`Sd`[3],
      res_MRPRESSO$`Main MR results`$`Sd`[4],
      NA,
      NA
    ),
    
    Pvalue = c(
      fmt_p(res_MRPRESSO$`Main MR results`$`P-value`[1]),
      fmt_p(res_MRPRESSO$`Main MR results`$`P-value`[2]),
      fmt_p(res_MRPRESSO$`Main MR results`$`P-value`[3]),
      fmt_p(res_MRPRESSO$`Main MR results`$`P-value`[4]),
      fmt_p(res_MRPRESSO$`MR-PRESSO results`$`Global Test`$Pvalue),
      fmt_p(res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$Pvalue)
    ),
    
    NSNPs = rep(length(MR_object@snps), 6),
    
    N_MR_PRESSO_outliers = rep(
      length(res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`)-1,
      6
    ),
    
    MR_PRESSO_outliers = rep(
      paste(
        MR_object@snps[
          res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
        ],
        collapse = ","
      ),
      6
    ),
    
    Method = c(
      "MR-PRESSO raw",
      "MR-PRESSO raw",
      "MR-PRESSO outlier-corrected",
      "MR-PRESSO outlier-corrected",
      "MR-PRESSO global test",
      "MR-PRESSO distortion test"
    ),
    
    stringsAsFactors = FALSE
  ) %>%
    mutate(across(where(is.character), ~ na_if(.x, "NA")))
}

MR_PRESSO_data <- rbind(
  make_MRPRESSO_df(PMD_MVMR_input, PMD_MRPRESSO),
  make_MRPRESSO_df(DA_MVMR_input,  DA_MRPRESSO),
  make_MRPRESSO_df(NDA_MVMR_input, NDA_MRPRESSO)
)

# Merge

MR_results = rbind_fill(IVW_data,
                        Egger_data,
                        Median_data,
                        LASSO_data,
                        cML_data,
                        GMM_data,
                        MR_PRESSO_data)

# Save results

fwrite(x = MR_results, file = "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_results.tsv",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
)

xlsx::write.xlsx(MR_results, "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_results.xlsx", row.names = FALSE,
                 showNA = FALSE)



## MVMR package Q-statistic minimization 

PMD_QHET <- qhet_mvmr(r_input = PMD_MVMR_sensitivity_input,
                      pcor = matrix(c(1,0.43,0.43,1), nrow = 2, ncol = 2),
                      CI=TRUE, iterations = 1000)

DA_QHET <- qhet_mvmr(r_input = DA_MVMR_sensitivity_input,
                      pcor = matrix(c(1,0.43,0.43,1), nrow = 2, ncol = 2),
                      CI=TRUE, iterations = 1000)

NDA_QHET <- qhet_mvmr(r_input = NDA_MVMR_sensitivity_input,
                     pcor = matrix(c(1,0.43,0.43,1), nrow = 2, ncol = 2),
                     CI=TRUE, iterations = 1000)

saveRDS(PMD_QHET, "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_PMD_QHET.rds")
saveRDS(DA_QHET, "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_DA_QHET.rds")
saveRDS(NDA_QHET, "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_NDA_QHET.rds")


##

