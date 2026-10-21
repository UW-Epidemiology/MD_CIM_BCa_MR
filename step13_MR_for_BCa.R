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
  
  ### cML (MendelianRandomization) 
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
  
  res_MRPRESSO <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      data = MR_object_TS,
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      SignifThreshold = 0.05,
      NbDistribution = 10000,
      seed = 40292
    ),
    error = function(e) {
      message("MR-PRESSO failed: ", e$message)
      return(NULL)
    }
  )
  
  # If MR-PRESSO failed, fill with NAs
  if (is.null(res_MRPRESSO)) {
    ## Extract single values safely
    exposure_val <- unique(MR_object_TS$exposure)[1]
    outcome_val  <- unique(MR_object_TS$outcome)[1]
    nsnps_val    <- length(unique(MR_object_TS$SNP))
    
    ## Build dataframe with fixed-length columns
    MR_PRESSO_df <- data.frame(
      Exposure    = rep(exposure_val, 3),
      Outcome     = rep(outcome_val, 3),
      SNP         = rep(paste0(exposure_val, " instrument"), 3),
      Estimate    = rep(NA_real_, 3),
      StdError    = rep(NA_real_, 3),
      Pvalue      = rep(NA_character_, 3),
      NSNPs       = rep(nsnps_val, 3),
      N_MR_PRESSO_outliers = rep(NA_integer_, 3),
      MR_PRESSO_outliers   = rep(NA_character_, 3),
      Method      = c(
        "MR-PRESSO outlier-corrected",
        "MR-PRESSO global test",
        "MR-PRESSO distortion test"
      ),
      stringsAsFactors = FALSE
    )
    
    ## Keep only unique rows
    MR_PRESSO_df <- unique(MR_PRESSO_df)
  } else {
    ## Safely extract MR-PRESSO results
    estimate <- if (!is.null(res_MRPRESSO$`Main MR results`$`Causal Estimate`)) {
      res_MRPRESSO$`Main MR results`$`Causal Estimate`[2]
    } else NA
    
    stderr <- if (!is.null(res_MRPRESSO$`Main MR results`$Sd)) {
      res_MRPRESSO$`Main MR results`$Sd[2]
    } else NA
    
    p_main <- if (!is.null(res_MRPRESSO$`Main MR results`$`P-value`)) {
      res_MRPRESSO$`Main MR results`$`P-value`[2] %>% 
        as.numeric() %>% round(3) %>% as.character()
    } else NA
    
    p_global <- if (!is.null(res_MRPRESSO$`MR-PRESSO results`$`Global Test`$Pvalue)) {
      res_MRPRESSO$`MR-PRESSO results`$`Global Test`$Pvalue
    } else NA
    
    p_dist <- if (!is.null(res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$Pvalue)) {
      res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$Pvalue %>% 
        as.numeric() %>% round(3) %>% as.character()
    } else NA
    
    outlier_idx <- res_MRPRESSO$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`
    
    n_outliers <- if (!is.null(outlier_idx)) length(outlier_idx) else NA
    
    outlier_snps <- if (!is.null(outlier_idx) && length(outlier_idx) > 0) {
      paste(MR_object_TS$SNP[outlier_idx], collapse = ",")
    } else NA
    
    
    ## Build dataframe with fixed-length columns
    MR_PRESSO_df <- data.frame(
      Exposure    = rep(unique(MR_object_TS$exposure)[1], 3),
      Outcome     = rep(unique(MR_object_TS$outcome)[1], 3),
      SNP         = rep(paste0(unique(MR_object_TS$exposure)[1], " instrument"), 3),
      Estimate    = c(estimate, NA, NA),
      StdError    = c(stderr, NA, NA),
      Pvalue      = c(p_main, p_global, p_dist),
      NSNPs       = rep(length(unique(MR_object_TS$SNP)), 3),
      N_MR_PRESSO_outliers = rep(n_outliers, 3),
      MR_PRESSO_outliers   = rep(outlier_snps, 3),
      Method      = c(
        "MR-PRESSO outlier-corrected",
        "MR-PRESSO global test",
        "MR-PRESSO distortion test"
      ),
      stringsAsFactors = FALSE
    )
    
    ## Keep only unique rows
    MR_PRESSO_df <- unique(MR_PRESSO_df)
  }
  
  return(MR_PRESSO_df)
}


# Run MR analyses

dir_path <- "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/BCa/harmonized_fully_filtered/"

filepaths <- list.files(
  path = dir_path,
  pattern = "\\.tsv$",
  full.names = TRUE
)

proteins <- sub("_.*$", "", basename(filepaths))


# File for all  results 
MR_results = NULL

for (i in 1:length(proteins)) {
  
  ## Load and process harmonized data
  data = fread(filepaths[i])
  
  if(nrow(data)< 2){
    next
  } 
  
  Overall_BCa = data %>% 
    filter(outcome == "Overall_BCa")
  
  HER2_pos_BCa = data %>% 
    filter(outcome == "HER2_pos_BCa")
  
  Luminal_A_BCa = data %>% 
    filter(outcome == "Luminal_A_BCa")
  
  Luminal_B_BCa = data %>% 
    filter(outcome == "Luminal_B_BCa")
  
  ER_neg_BCa = data %>% 
    filter(outcome == "ER_neg_BCa")
  
  TNBCa = data %>% 
    filter(outcome == "TNBCa")
  
  Overall_BCa_MR_input = mr_input(bx = Overall_BCa$Effect_exposure,
                          bxse = Overall_BCa$StdErr_exposure,
                          by = Overall_BCa$Effect_outcome,
                          byse = Overall_BCa$StdErr_outcome,
                          exposure = Overall_BCa$exposure,
                          outcome = Overall_BCa$outcome,
                          snps = Overall_BCa$RSID)
  
  Overall_BCa_MR_input_TS = data.frame(SNP = Overall_BCa$RSID,
                            effect_allele.exposure = Overall_BCa$A1,
                            other_allele.exposure = Overall_BCa$A2,
                            effect_allele.outcome = Overall_BCa$A1,
                            other_allele.outcome = Overall_BCa$A2,
                            beta.exposure = Overall_BCa$Effect_exposure,
                            beta.outcome = Overall_BCa$Effect_outcome,
                            eaf.exposure = rep(NA, nrow(Overall_BCa)),
                            eaf.outcome = rep(NA, nrow(Overall_BCa)),
                            remove = rep(FALSE, nrow(Overall_BCa)),
                            palindromic = rep(FALSE, nrow(Overall_BCa)),
                            ambigouous = rep(NA, nrow(Overall_BCa)),
                            id.outcome = Overall_BCa$outcome,
                            chr.outcome = Overall_BCa$Chr,
                            pos.outcome = Overall_BCa$Pos,
                            pval.outcome = Overall_BCa$P_outcome,
                            se.outcome = Overall_BCa$StdErr_outcome,
                            outcome = rep("Overall_BCa", nrow(Overall_BCa)),
                            mr_keep.outcome = rep(TRUE, nrow(Overall_BCa)),
                            pval_origin.outcome = rep("reported", nrow(Overall_BCa)),
                            id.exposure = rep("Overall_BCa", nrow(Overall_BCa)),
                            chr.exposure = Overall_BCa$Chr,
                            pos.exposure = Overall_BCa$Pos,
                            pval.exposure = Overall_BCa$P_exposure,
                            se.exposure = Overall_BCa$StdErr_exposure,
                            exposure = Overall_BCa$exposure,
                            mr_keep.exposure = rep(TRUE, nrow(Overall_BCa)),
                            pval_origin.exposure = rep("reported", nrow(Overall_BCa)),
                            data_source.exposure = rep("IMD_GWAS", nrow(Overall_BCa)),
                            action = rep(2, nrow(Overall_BCa)),
                            SNP_index = rep(1, nrow(Overall_BCa)),
                            mr_keep = rep(TRUE, nrow(Overall_BCa)),
                            samplesize.outcome = rep(247173, nrow(Overall_BCa)))
  
  HER2_pos_BCa_MR_input = mr_input(bx = HER2_pos_BCa$Effect_exposure,
                                  bxse = HER2_pos_BCa$StdErr_exposure,
                                  by = HER2_pos_BCa$Effect_outcome,
                                  byse = HER2_pos_BCa$StdErr_outcome,
                                  exposure = HER2_pos_BCa$exposure,
                                  outcome = HER2_pos_BCa$outcome,
                                  snps = HER2_pos_BCa$RSID)
  
  HER2_pos_BCa_MR_input_TS = data.frame(SNP = HER2_pos_BCa$RSID,
                                       effect_allele.exposure = HER2_pos_BCa$A1,
                                       other_allele.exposure = HER2_pos_BCa$A2,
                                       effect_allele.outcome = HER2_pos_BCa$A1,
                                       other_allele.outcome = HER2_pos_BCa$A2,
                                       beta.exposure = HER2_pos_BCa$Effect_exposure,
                                       beta.outcome = HER2_pos_BCa$Effect_outcome,
                                       eaf.exposure = rep(NA, nrow(HER2_pos_BCa)),
                                       eaf.outcome = rep(NA, nrow(HER2_pos_BCa)),
                                       remove = rep(FALSE, nrow(HER2_pos_BCa)),
                                       palindromic = rep(FALSE, nrow(HER2_pos_BCa)),
                                       ambigouous = rep(NA, nrow(HER2_pos_BCa)),
                                       id.outcome = HER2_pos_BCa$outcome,
                                       chr.outcome = HER2_pos_BCa$Chr,
                                       pos.outcome = HER2_pos_BCa$Pos,
                                       pval.outcome = HER2_pos_BCa$P_outcome,
                                       se.outcome = HER2_pos_BCa$StdErr_outcome,
                                       outcome = rep("HER2_pos_BCa", nrow(HER2_pos_BCa)),
                                       mr_keep.outcome = rep(TRUE, nrow(HER2_pos_BCa)),
                                       pval_origin.outcome = rep("reported", nrow(HER2_pos_BCa)),
                                       id.exposure = rep("HER2_pos_BCa", nrow(HER2_pos_BCa)),
                                       chr.exposure = HER2_pos_BCa$Chr,
                                       pos.exposure = HER2_pos_BCa$Pos,
                                       pval.exposure = HER2_pos_BCa$P_exposure,
                                       se.exposure = HER2_pos_BCa$StdErr_exposure,
                                       exposure = HER2_pos_BCa$exposure,
                                       mr_keep.exposure = rep(TRUE, nrow(HER2_pos_BCa)),
                                       pval_origin.exposure = rep("reported", nrow(HER2_pos_BCa)),
                                       data_source.exposure = rep("IMD_GWAS", nrow(HER2_pos_BCa)),
                                       action = rep(2, nrow(HER2_pos_BCa)),
                                       SNP_index = rep(1, nrow(HER2_pos_BCa)),
                                       mr_keep = rep(TRUE, nrow(HER2_pos_BCa)),
                                       samplesize.outcome = rep(97877, nrow(HER2_pos_BCa)))
  
  Luminal_A_BCa_MR_input = mr_input(bx = Luminal_A_BCa$Effect_exposure,
                                  bxse = Luminal_A_BCa$StdErr_exposure,
                                  by = Luminal_A_BCa$Effect_outcome,
                                  byse = Luminal_A_BCa$StdErr_outcome,
                                  exposure = Luminal_A_BCa$exposure,
                                  outcome = Luminal_A_BCa$outcome,
                                  snps = Luminal_A_BCa$RSID)
  
  Luminal_A_BCa_MR_input_TS = data.frame(SNP = Luminal_A_BCa$RSID,
                                       effect_allele.exposure = Luminal_A_BCa$A1,
                                       other_allele.exposure = Luminal_A_BCa$A2,
                                       effect_allele.outcome = Luminal_A_BCa$A1,
                                       other_allele.outcome = Luminal_A_BCa$A2,
                                       beta.exposure = Luminal_A_BCa$Effect_exposure,
                                       beta.outcome = Luminal_A_BCa$Effect_outcome,
                                       eaf.exposure = rep(NA, nrow(Luminal_A_BCa)),
                                       eaf.outcome = rep(NA, nrow(Luminal_A_BCa)),
                                       remove = rep(FALSE, nrow(Luminal_A_BCa)),
                                       palindromic = rep(FALSE, nrow(Luminal_A_BCa)),
                                       ambigouous = rep(NA, nrow(Luminal_A_BCa)),
                                       id.outcome = Luminal_A_BCa$outcome,
                                       chr.outcome = Luminal_A_BCa$Chr,
                                       pos.outcome = Luminal_A_BCa$Pos,
                                       pval.outcome = Luminal_A_BCa$P_outcome,
                                       se.outcome = Luminal_A_BCa$StdErr_outcome,
                                       outcome = rep("Luminal_A_BCa", nrow(Luminal_A_BCa)),
                                       mr_keep.outcome = rep(TRUE, nrow(Luminal_A_BCa)),
                                       pval_origin.outcome = rep("reported", nrow(Luminal_A_BCa)),
                                       id.exposure = rep("Luminal_A_BCa", nrow(Luminal_A_BCa)),
                                       chr.exposure = Luminal_A_BCa$Chr,
                                       pos.exposure = Luminal_A_BCa$Pos,
                                       pval.exposure = Luminal_A_BCa$P_exposure,
                                       se.exposure = Luminal_A_BCa$StdErr_exposure,
                                       exposure = Luminal_A_BCa$exposure,
                                       mr_keep.exposure = rep(TRUE, nrow(Luminal_A_BCa)),
                                       pval_origin.exposure = rep("reported", nrow(Luminal_A_BCa)),
                                       data_source.exposure = rep("IMD_GWAS", nrow(Luminal_A_BCa)),
                                       action = rep(2, nrow(Luminal_A_BCa)),
                                       SNP_index = rep(1, nrow(Luminal_A_BCa)),
                                       mr_keep = rep(TRUE, nrow(Luminal_A_BCa)),
                                       samplesize.outcome = rep(148877, nrow(Luminal_A_BCa)))
                            
  Luminal_B_BCa_MR_input = mr_input(bx = Luminal_B_BCa$Effect_exposure,
                         bxse = Luminal_B_BCa$StdErr_exposure,
                         by = Luminal_B_BCa$Effect_outcome,
                         byse = Luminal_B_BCa$StdErr_outcome,
                         exposure = Luminal_B_BCa$exposure,
                         outcome = Luminal_B_BCa$outcome,
                         snps = Luminal_B_BCa$RSID)
  
  Luminal_B_BCa_MR_input_TS = data.frame(SNP = Luminal_B_BCa$RSID,
                               effect_allele.exposure = Luminal_B_BCa$A1,
                               other_allele.exposure = Luminal_B_BCa$A2,
                               effect_allele.outcome = Luminal_B_BCa$A1,
                               other_allele.outcome = Luminal_B_BCa$A2,
                               beta.exposure = Luminal_B_BCa$Effect_exposure,
                               beta.outcome = Luminal_B_BCa$Effect_outcome,
                               eaf.exposure = rep(NA, nrow(Luminal_B_BCa)),
                               eaf.outcome = rep(NA, nrow(Luminal_B_BCa)),
                               remove = rep(FALSE, nrow(Luminal_B_BCa)),
                               palindromic = rep(FALSE, nrow(Luminal_B_BCa)),
                               ambigouous = rep(NA, nrow(Luminal_B_BCa)),
                               id.outcome = Luminal_B_BCa$outcome,
                               chr.outcome = Luminal_B_BCa$Chr,
                               pos.outcome = Luminal_B_BCa$Pos,
                               pval.outcome = Luminal_B_BCa$P_outcome,
                               se.outcome = Luminal_B_BCa$StdErr_outcome,
                               outcome = rep("Luminal_B_BCa", nrow(Luminal_B_BCa)),
                               mr_keep.outcome = rep(TRUE, nrow(Luminal_B_BCa)),
                               pval_origin.outcome = rep("reported", nrow(Luminal_B_BCa)),
                               id.exposure = rep("Luminal_B_BCa", nrow(Luminal_B_BCa)),
                               chr.exposure = Luminal_B_BCa$Chr,
                               pos.exposure = Luminal_B_BCa$Pos,
                               pval.exposure = Luminal_B_BCa$P_exposure,
                               se.exposure = Luminal_B_BCa$StdErr_exposure,
                               exposure = Luminal_B_BCa$exposure,
                               mr_keep.exposure = rep(TRUE, nrow(Luminal_B_BCa)),
                               pval_origin.exposure = rep("reported", nrow(Luminal_B_BCa)),
                               data_source.exposure = rep("IMD_GWAS", nrow(Luminal_B_BCa)),
                               action = rep(2, nrow(Luminal_B_BCa)),
                               SNP_index = rep(1, nrow(Luminal_B_BCa)),
                               mr_keep = rep(TRUE, nrow(Luminal_B_BCa)),
                               samplesize.outcome = rep(105277, nrow(Luminal_B_BCa)))
  
  ER_neg_BCa_MR_input = mr_input(bx = ER_neg_BCa$Effect_exposure,
                          bxse = ER_neg_BCa$StdErr_exposure,
                          by = ER_neg_BCa$Effect_outcome,
                          byse = ER_neg_BCa$StdErr_outcome,
                          exposure = ER_neg_BCa$exposure,
                          outcome = ER_neg_BCa$outcome,
                          snps = ER_neg_BCa$RSID)
  
  ER_neg_BCa_MR_input_TS = data.frame(SNP = ER_neg_BCa$RSID,
                               effect_allele.exposure = ER_neg_BCa$A1,
                               other_allele.exposure = ER_neg_BCa$A2,
                               effect_allele.outcome = ER_neg_BCa$A1,
                               other_allele.outcome = ER_neg_BCa$A2,
                               beta.exposure = ER_neg_BCa$Effect_exposure,
                               beta.outcome = ER_neg_BCa$Effect_outcome,
                               eaf.exposure = rep(NA, nrow(ER_neg_BCa)),
                               eaf.outcome = rep(NA, nrow(ER_neg_BCa)),
                               remove = rep(FALSE, nrow(ER_neg_BCa)),
                               palindromic = rep(FALSE, nrow(ER_neg_BCa)),
                               ambigouous = rep(NA, nrow(ER_neg_BCa)),
                               id.outcome = ER_neg_BCa$outcome,
                               chr.outcome = ER_neg_BCa$Chr,
                               pos.outcome = ER_neg_BCa$Pos,
                               pval.outcome = ER_neg_BCa$P_outcome,
                               se.outcome = ER_neg_BCa$StdErr_outcome,
                               outcome = rep("ER_neg_BCa", nrow(ER_neg_BCa)),
                               mr_keep.outcome = rep(TRUE, nrow(ER_neg_BCa)),
                               pval_origin.outcome = rep("reported", nrow(ER_neg_BCa)),
                               id.exposure = rep("ER_neg_BCa", nrow(ER_neg_BCa)),
                               chr.exposure = ER_neg_BCa$Chr,
                               pos.exposure = ER_neg_BCa$Pos,
                               pval.exposure = ER_neg_BCa$P_exposure,
                               se.exposure = ER_neg_BCa$StdErr_exposure,
                               exposure = ER_neg_BCa$exposure,
                               mr_keep.exposure = rep(TRUE, nrow(ER_neg_BCa)),
                               pval_origin.exposure = rep("reported", nrow(ER_neg_BCa)),
                               data_source.exposure = rep("IMD_GWAS", nrow(ER_neg_BCa)),
                               action = rep(2, nrow(ER_neg_BCa)),
                               SNP_index = rep(1, nrow(ER_neg_BCa)),
                               mr_keep = rep(TRUE, nrow(ER_neg_BCa)),
                               samplesize.outcome = rep(72261, nrow(ER_neg_BCa)))
  
  TNBCa_MR_input = mr_input(bx = TNBCa$Effect_exposure,
                                  bxse = TNBCa$StdErr_exposure,
                                  by = TNBCa$Effect_outcome,
                                  byse = TNBCa$StdErr_outcome,
                                  exposure = TNBCa$exposure,
                                  outcome = TNBCa$outcome,
                                  snps = TNBCa$RSID)
  
  TNBCa_MR_input_TS = data.frame(SNP = TNBCa$RSID,
                                       effect_allele.exposure = TNBCa$A1,
                                       other_allele.exposure = TNBCa$A2,
                                       effect_allele.outcome = TNBCa$A1,
                                       other_allele.outcome = TNBCa$A2,
                                       beta.exposure = TNBCa$Effect_exposure,
                                       beta.outcome = TNBCa$Effect_outcome,
                                       eaf.exposure = rep(NA, nrow(TNBCa)),
                                       eaf.outcome = rep(NA, nrow(TNBCa)),
                                       remove = rep(FALSE, nrow(TNBCa)),
                                       palindromic = rep(FALSE, nrow(TNBCa)),
                                       ambigouous = rep(NA, nrow(TNBCa)),
                                       id.outcome = TNBCa$outcome,
                                       chr.outcome = TNBCa$Chr,
                                       pos.outcome = TNBCa$Pos,
                                       pval.outcome = TNBCa$P_outcome,
                                       se.outcome = TNBCa$StdErr_outcome,
                                       outcome = rep("TNBCa", nrow(TNBCa)),
                                       mr_keep.outcome = rep(TRUE, nrow(TNBCa)),
                                       pval_origin.outcome = rep("reported", nrow(TNBCa)),
                                       id.exposure = rep("TNBCa", nrow(TNBCa)),
                                       chr.exposure = TNBCa$Chr,
                                       pos.exposure = TNBCa$Pos,
                                       pval.exposure = TNBCa$P_exposure,
                                       se.exposure = TNBCa$StdErr_exposure,
                                       exposure = TNBCa$exposure,
                                       mr_keep.exposure = rep(TRUE, nrow(TNBCa)),
                                       pval_origin.exposure = rep("reported", nrow(TNBCa)),
                                       data_source.exposure = rep("IMD_GWAS", nrow(TNBCa)),
                                       action = rep(2, nrow(TNBCa)),
                                       SNP_index = rep(1, nrow(TNBCa)),
                                       mr_keep = rep(TRUE, nrow(TNBCa)),
                                       samplesize.outcome = rep(106377, nrow(TNBCa)))
  
  
  ## Run IVW_MR analyses 
  if(nrow(Overall_BCa)< 3){
    Overall_BCa_IVW = run_MR_ivw(Overall_BCa_MR_input,Overall_BCa_MR_input_TS)
    HER2_pos_BCa_IVW = run_MR_ivw(HER2_pos_BCa_MR_input,HER2_pos_BCa_MR_input_TS)
    Luminal_A_BCa_IVW = run_MR_ivw(Luminal_A_BCa_MR_input, Luminal_A_BCa_MR_input_TS)
    Luminal_B_BCa_IVW = run_MR_ivw(Luminal_B_BCa_MR_input, Luminal_B_BCa_MR_input_TS)
    ER_neg_BCa_IVW = run_MR_ivw(ER_neg_BCa_MR_input, ER_neg_BCa_MR_input_TS)
    TNBCa_IVW = run_MR_ivw(TNBCa_MR_input,TNBCa_MR_input_TS)
    exposure_MR_results = rbind.data.frame(Overall_BCa_IVW,HER2_pos_BCa_IVW,
                                           Luminal_A_BCa_IVW,Luminal_B_BCa_IVW,
                                           ER_neg_BCa_IVW,TNBCa_IVW)
  }else{
    safe_run <- function(expr) {
      tryCatch(
        expr,
        error = function(e) NA
      )
    }
    
    # IVW
    Overall_BCa_IVW <- safe_run(run_MR_ivw(Overall_BCa_MR_input, Overall_BCa_MR_input_TS))
    HER2_pos_BCa_IVW <- safe_run(run_MR_ivw(HER2_pos_BCa_MR_input, HER2_pos_BCa_MR_input_TS))
    Luminal_A_BCa_IVW <- safe_run(run_MR_ivw(Luminal_A_BCa_MR_input, Luminal_A_BCa_MR_input_TS))
    Luminal_B_BCa_IVW <- safe_run(run_MR_ivw(Luminal_B_BCa_MR_input, Luminal_B_BCa_MR_input_TS))
    ER_neg_BCa_IVW <- safe_run(run_MR_ivw(ER_neg_BCa_MR_input, ER_neg_BCa_MR_input_TS))
    TNBCa_IVW <- safe_run(run_MR_ivw(TNBCa_MR_input, TNBCa_MR_input_TS))
    
    # MR-Egger
    Overall_BCa_Egger <- safe_run(run_MR_Egger(Overall_BCa_MR_input))
    HER2_pos_BCa_Egger <- safe_run(run_MR_Egger(HER2_pos_BCa_MR_input))
    Luminal_A_BCa_Egger <- safe_run(run_MR_Egger(Luminal_A_BCa_MR_input))
    Luminal_B_BCa_Egger <- safe_run(run_MR_Egger(Luminal_B_BCa_MR_input))
    ER_neg_BCa_Egger <- safe_run(run_MR_Egger(ER_neg_BCa_MR_input))
    TNBCa_Egger <- safe_run(run_MR_Egger(TNBCa_MR_input))
    
    # MR Median
    Overall_BCa_Median <- safe_run(run_MR_median(Overall_BCa_MR_input))
    HER2_pos_BCa_Median <- safe_run(run_MR_median(HER2_pos_BCa_MR_input))
    Luminal_A_BCa_Median <- safe_run(run_MR_median(Luminal_A_BCa_MR_input))
    Luminal_B_BCa_Median <- safe_run(run_MR_median(Luminal_B_BCa_MR_input))
    ER_neg_BCa_Median <- safe_run(run_MR_median(ER_neg_BCa_MR_input))
    TNBCa_Median <- safe_run(run_MR_median(TNBCa_MR_input))
    
    # LASSO
    Overall_BCa_LASSO <- safe_run(run_MR_LASSO(Overall_BCa_MR_input))
    HER2_pos_BCa_LASSO <- safe_run(run_MR_LASSO(HER2_pos_BCa_MR_input))
    Luminal_A_BCa_LASSO <- safe_run(run_MR_LASSO(Luminal_A_BCa_MR_input))
    Luminal_B_BCa_LASSO <- safe_run(run_MR_LASSO(Luminal_B_BCa_MR_input))
    ER_neg_BCa_LASSO <- safe_run(run_MR_LASSO(ER_neg_BCa_MR_input))
    TNBCa_LASSO <- safe_run(run_MR_LASSO(TNBCa_MR_input))
    
    # cML
    Overall_BCa_cML <- safe_run(
      run_MR_cML(
        Overall_BCa_MR_input,
        n_out = 247173,
        n_snp_valid = ifelse(is.list(Overall_BCa_LASSO), Overall_BCa_LASSO$LASSO_N_valid, NA)
      )
    )
    
    HER2_pos_BCa_cML <- safe_run(
      run_MR_cML(
        HER2_pos_BCa_MR_input,
        n_out = 97877,
        n_snp_valid = ifelse(is.list(HER2_pos_BCa_LASSO), HER2_pos_BCa_LASSO$LASSO_N_valid, NA)
      )
    )
    
    Luminal_A_BCa_cML <- safe_run(
      run_MR_cML(
        Luminal_A_BCa_MR_input,
        n_out = 148877,
        n_snp_valid = ifelse(is.list(Luminal_A_BCa_LASSO), Luminal_A_BCa_LASSO$LASSO_N_valid, NA)
      )
    )
    
    Luminal_B_BCa_cML <- safe_run(
      run_MR_cML(
        Luminal_B_BCa_MR_input,
        n_out = 105277,
        n_snp_valid = ifelse(is.list(Luminal_B_BCa_LASSO), Luminal_B_BCa_LASSO$LASSO_N_valid, NA)
      )
    )
    
    ER_neg_BCa_cML <- safe_run(
      run_MR_cML(
        ER_neg_BCa_MR_input,
        n_out = 72261,
        n_snp_valid = ifelse(is.list(ER_neg_BCa_LASSO), ER_neg_BCa_LASSO$LASSO_N_valid, NA)
      )
    )
    
    TNBCa_cML <- safe_run(
      run_MR_cML(
        TNBCa_MR_input,
        n_out = 106377,
        n_snp_valid = ifelse(is.list(TNBCa_LASSO), TNBCa_LASSO$LASSO_N_valid, NA)
      )
    )
    
    # Collect all results
    obj_names <- c(
      "Overall_BCa_IVW", "Overall_BCa_Egger", "Overall_BCa_Median", "Overall_BCa_LASSO", "Overall_BCa_cML",
      "HER2_pos_BCa_IVW", "HER2_pos_BCa_Egger", "HER2_pos_BCa_Median", "HER2_pos_BCa_LASSO", "HER2_pos_BCa_cML",
      "Luminal_A_BCa_IVW", "Luminal_A_BCa_Egger", "Luminal_A_BCa_Median", "Luminal_A_BCa_LASSO", "Luminal_A_BCa_cML",
      "Luminal_B_BCa_IVW", "Luminal_B_BCa_Egger", "Luminal_B_BCa_Median", "Luminal_B_BCa_LASSO", "Luminal_B_BCa_cML",
      "ER_neg_BCa_IVW", "ER_neg_BCa_Egger", "ER_neg_BCa_Median", "ER_neg_BCa_LASSO", "ER_neg_BCa_cML",
      "TNBCa_IVW", "TNBCa_Egger", "TNBCa_Median", "TNBCa_LASSO", "TNBCa_cML"
    )
    
    existing_objs <- mget(obj_names[sapply(obj_names, exists)], inherits = TRUE)
    
    exposure_MR_results <- do.call(rbind_fill, existing_objs)
    
  }
  
 MR_results = rbind_fill(MR_results,exposure_MR_results)
  
}

MR_results = MR_results[-1,]

# Save results

fwrite(x = MR_results, file = "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_GWAS_CIM_BCa_MR_results.tsv",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
)

xlsx::write.xlsx(MR_results, "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_GWAS_CIM_BCa_MR_results.xlsx", row.names = FALSE,
                 showNA = FALSE)

# Run univaraible MR analyses

dir_path <- "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/OpenGWAS_project_CIMs/BCa/"

filepaths <- list.files(
  path = dir_path,
  pattern = "\\.tsv$",
  full.names = TRUE
)

proteins <- c("CRP","CRP_noBMI","TNFR2")

MR_results = NULL

for (i in 1:length(proteins)) {
    ## Load and process harmonized data
    data = fread(filepaths[i])
    
    if(nrow(data)< 2){
      next
    } 
    
    Overall_BCa = data %>% 
      filter(outcome == "Overall_BCa")
    
    HER2_pos_BCa = data %>% 
      filter(outcome == "HER2_pos_BCa")
    
    Luminal_A_BCa = data %>% 
      filter(outcome == "Luminal_A_BCa")
    
    Luminal_B_BCa = data %>% 
      filter(outcome == "Luminal_B_BCa")
    
    ER_neg_BCa = data %>% 
      filter(outcome == "ER_neg_BCa")
    
    TNBCa = data %>% 
      filter(outcome == "TNBCa")
    
    Overall_BCa_MR_input = mr_input(bx = Overall_BCa$Effect_exposure,
                                    bxse = Overall_BCa$StdErr_exposure,
                                    by = Overall_BCa$Effect_outcome,
                                    byse = Overall_BCa$StdErr_outcome,
                                    exposure = Overall_BCa$exposure,
                                    outcome = Overall_BCa$outcome,
                                    snps = Overall_BCa$RSID)
    
    Overall_BCa_MR_input_TS = data.frame(SNP = Overall_BCa$RSID,
                                         effect_allele.exposure = Overall_BCa$A1,
                                         other_allele.exposure = Overall_BCa$A2,
                                         effect_allele.outcome = Overall_BCa$A1,
                                         other_allele.outcome = Overall_BCa$A2,
                                         beta.exposure = Overall_BCa$Effect_exposure,
                                         beta.outcome = Overall_BCa$Effect_outcome,
                                         eaf.exposure = rep(NA, nrow(Overall_BCa)),
                                         eaf.outcome = rep(NA, nrow(Overall_BCa)),
                                         remove = rep(FALSE, nrow(Overall_BCa)),
                                         palindromic = rep(FALSE, nrow(Overall_BCa)),
                                         ambigouous = rep(NA, nrow(Overall_BCa)),
                                         id.outcome = Overall_BCa$outcome,
                                         chr.outcome = Overall_BCa$Chr,
                                         pos.outcome = Overall_BCa$Pos,
                                         pval.outcome = Overall_BCa$P_outcome,
                                         se.outcome = Overall_BCa$StdErr_outcome,
                                         outcome = rep("Overall_BCa", nrow(Overall_BCa)),
                                         mr_keep.outcome = rep(TRUE, nrow(Overall_BCa)),
                                         pval_origin.outcome = rep("reported", nrow(Overall_BCa)),
                                         id.exposure = rep("Overall_BCa", nrow(Overall_BCa)),
                                         chr.exposure = Overall_BCa$Chr,
                                         pos.exposure = Overall_BCa$Pos,
                                         pval.exposure = Overall_BCa$P_exposure,
                                         se.exposure = Overall_BCa$StdErr_exposure,
                                         exposure = Overall_BCa$exposure,
                                         mr_keep.exposure = rep(TRUE, nrow(Overall_BCa)),
                                         pval_origin.exposure = rep("reported", nrow(Overall_BCa)),
                                         data_source.exposure = rep("IMD_GWAS", nrow(Overall_BCa)),
                                         action = rep(2, nrow(Overall_BCa)),
                                         SNP_index = rep(1, nrow(Overall_BCa)),
                                         mr_keep = rep(TRUE, nrow(Overall_BCa)),
                                         samplesize.outcome = rep(247173, nrow(Overall_BCa)))
    
    HER2_pos_BCa_MR_input = mr_input(bx = HER2_pos_BCa$Effect_exposure,
                                     bxse = HER2_pos_BCa$StdErr_exposure,
                                     by = HER2_pos_BCa$Effect_outcome,
                                     byse = HER2_pos_BCa$StdErr_outcome,
                                     exposure = HER2_pos_BCa$exposure,
                                     outcome = HER2_pos_BCa$outcome,
                                     snps = HER2_pos_BCa$RSID)
    
    HER2_pos_BCa_MR_input_TS = data.frame(SNP = HER2_pos_BCa$RSID,
                                          effect_allele.exposure = HER2_pos_BCa$A1,
                                          other_allele.exposure = HER2_pos_BCa$A2,
                                          effect_allele.outcome = HER2_pos_BCa$A1,
                                          other_allele.outcome = HER2_pos_BCa$A2,
                                          beta.exposure = HER2_pos_BCa$Effect_exposure,
                                          beta.outcome = HER2_pos_BCa$Effect_outcome,
                                          eaf.exposure = rep(NA, nrow(HER2_pos_BCa)),
                                          eaf.outcome = rep(NA, nrow(HER2_pos_BCa)),
                                          remove = rep(FALSE, nrow(HER2_pos_BCa)),
                                          palindromic = rep(FALSE, nrow(HER2_pos_BCa)),
                                          ambigouous = rep(NA, nrow(HER2_pos_BCa)),
                                          id.outcome = HER2_pos_BCa$outcome,
                                          chr.outcome = HER2_pos_BCa$Chr,
                                          pos.outcome = HER2_pos_BCa$Pos,
                                          pval.outcome = HER2_pos_BCa$P_outcome,
                                          se.outcome = HER2_pos_BCa$StdErr_outcome,
                                          outcome = rep("HER2_pos_BCa", nrow(HER2_pos_BCa)),
                                          mr_keep.outcome = rep(TRUE, nrow(HER2_pos_BCa)),
                                          pval_origin.outcome = rep("reported", nrow(HER2_pos_BCa)),
                                          id.exposure = rep("HER2_pos_BCa", nrow(HER2_pos_BCa)),
                                          chr.exposure = HER2_pos_BCa$Chr,
                                          pos.exposure = HER2_pos_BCa$Pos,
                                          pval.exposure = HER2_pos_BCa$P_exposure,
                                          se.exposure = HER2_pos_BCa$StdErr_exposure,
                                          exposure = HER2_pos_BCa$exposure,
                                          mr_keep.exposure = rep(TRUE, nrow(HER2_pos_BCa)),
                                          pval_origin.exposure = rep("reported", nrow(HER2_pos_BCa)),
                                          data_source.exposure = rep("IMD_GWAS", nrow(HER2_pos_BCa)),
                                          action = rep(2, nrow(HER2_pos_BCa)),
                                          SNP_index = rep(1, nrow(HER2_pos_BCa)),
                                          mr_keep = rep(TRUE, nrow(HER2_pos_BCa)),
                                          samplesize.outcome = rep(97877, nrow(HER2_pos_BCa)))
    
    Luminal_A_BCa_MR_input = mr_input(bx = Luminal_A_BCa$Effect_exposure,
                                      bxse = Luminal_A_BCa$StdErr_exposure,
                                      by = Luminal_A_BCa$Effect_outcome,
                                      byse = Luminal_A_BCa$StdErr_outcome,
                                      exposure = Luminal_A_BCa$exposure,
                                      outcome = Luminal_A_BCa$outcome,
                                      snps = Luminal_A_BCa$RSID)
    
    Luminal_A_BCa_MR_input_TS = data.frame(SNP = Luminal_A_BCa$RSID,
                                           effect_allele.exposure = Luminal_A_BCa$A1,
                                           other_allele.exposure = Luminal_A_BCa$A2,
                                           effect_allele.outcome = Luminal_A_BCa$A1,
                                           other_allele.outcome = Luminal_A_BCa$A2,
                                           beta.exposure = Luminal_A_BCa$Effect_exposure,
                                           beta.outcome = Luminal_A_BCa$Effect_outcome,
                                           eaf.exposure = rep(NA, nrow(Luminal_A_BCa)),
                                           eaf.outcome = rep(NA, nrow(Luminal_A_BCa)),
                                           remove = rep(FALSE, nrow(Luminal_A_BCa)),
                                           palindromic = rep(FALSE, nrow(Luminal_A_BCa)),
                                           ambigouous = rep(NA, nrow(Luminal_A_BCa)),
                                           id.outcome = Luminal_A_BCa$outcome,
                                           chr.outcome = Luminal_A_BCa$Chr,
                                           pos.outcome = Luminal_A_BCa$Pos,
                                           pval.outcome = Luminal_A_BCa$P_outcome,
                                           se.outcome = Luminal_A_BCa$StdErr_outcome,
                                           outcome = rep("Luminal_A_BCa", nrow(Luminal_A_BCa)),
                                           mr_keep.outcome = rep(TRUE, nrow(Luminal_A_BCa)),
                                           pval_origin.outcome = rep("reported", nrow(Luminal_A_BCa)),
                                           id.exposure = rep("Luminal_A_BCa", nrow(Luminal_A_BCa)),
                                           chr.exposure = Luminal_A_BCa$Chr,
                                           pos.exposure = Luminal_A_BCa$Pos,
                                           pval.exposure = Luminal_A_BCa$P_exposure,
                                           se.exposure = Luminal_A_BCa$StdErr_exposure,
                                           exposure = Luminal_A_BCa$exposure,
                                           mr_keep.exposure = rep(TRUE, nrow(Luminal_A_BCa)),
                                           pval_origin.exposure = rep("reported", nrow(Luminal_A_BCa)),
                                           data_source.exposure = rep("IMD_GWAS", nrow(Luminal_A_BCa)),
                                           action = rep(2, nrow(Luminal_A_BCa)),
                                           SNP_index = rep(1, nrow(Luminal_A_BCa)),
                                           mr_keep = rep(TRUE, nrow(Luminal_A_BCa)),
                                           samplesize.outcome = rep(148877, nrow(Luminal_A_BCa)))
    
    Luminal_B_BCa_MR_input = mr_input(bx = Luminal_B_BCa$Effect_exposure,
                                      bxse = Luminal_B_BCa$StdErr_exposure,
                                      by = Luminal_B_BCa$Effect_outcome,
                                      byse = Luminal_B_BCa$StdErr_outcome,
                                      exposure = Luminal_B_BCa$exposure,
                                      outcome = Luminal_B_BCa$outcome,
                                      snps = Luminal_B_BCa$RSID)
    
    Luminal_B_BCa_MR_input_TS = data.frame(SNP = Luminal_B_BCa$RSID,
                                           effect_allele.exposure = Luminal_B_BCa$A1,
                                           other_allele.exposure = Luminal_B_BCa$A2,
                                           effect_allele.outcome = Luminal_B_BCa$A1,
                                           other_allele.outcome = Luminal_B_BCa$A2,
                                           beta.exposure = Luminal_B_BCa$Effect_exposure,
                                           beta.outcome = Luminal_B_BCa$Effect_outcome,
                                           eaf.exposure = rep(NA, nrow(Luminal_B_BCa)),
                                           eaf.outcome = rep(NA, nrow(Luminal_B_BCa)),
                                           remove = rep(FALSE, nrow(Luminal_B_BCa)),
                                           palindromic = rep(FALSE, nrow(Luminal_B_BCa)),
                                           ambigouous = rep(NA, nrow(Luminal_B_BCa)),
                                           id.outcome = Luminal_B_BCa$outcome,
                                           chr.outcome = Luminal_B_BCa$Chr,
                                           pos.outcome = Luminal_B_BCa$Pos,
                                           pval.outcome = Luminal_B_BCa$P_outcome,
                                           se.outcome = Luminal_B_BCa$StdErr_outcome,
                                           outcome = rep("Luminal_B_BCa", nrow(Luminal_B_BCa)),
                                           mr_keep.outcome = rep(TRUE, nrow(Luminal_B_BCa)),
                                           pval_origin.outcome = rep("reported", nrow(Luminal_B_BCa)),
                                           id.exposure = rep("Luminal_B_BCa", nrow(Luminal_B_BCa)),
                                           chr.exposure = Luminal_B_BCa$Chr,
                                           pos.exposure = Luminal_B_BCa$Pos,
                                           pval.exposure = Luminal_B_BCa$P_exposure,
                                           se.exposure = Luminal_B_BCa$StdErr_exposure,
                                           exposure = Luminal_B_BCa$exposure,
                                           mr_keep.exposure = rep(TRUE, nrow(Luminal_B_BCa)),
                                           pval_origin.exposure = rep("reported", nrow(Luminal_B_BCa)),
                                           data_source.exposure = rep("IMD_GWAS", nrow(Luminal_B_BCa)),
                                           action = rep(2, nrow(Luminal_B_BCa)),
                                           SNP_index = rep(1, nrow(Luminal_B_BCa)),
                                           mr_keep = rep(TRUE, nrow(Luminal_B_BCa)),
                                           samplesize.outcome = rep(105277, nrow(Luminal_B_BCa)))
    
    ER_neg_BCa_MR_input = mr_input(bx = ER_neg_BCa$Effect_exposure,
                                   bxse = ER_neg_BCa$StdErr_exposure,
                                   by = ER_neg_BCa$Effect_outcome,
                                   byse = ER_neg_BCa$StdErr_outcome,
                                   exposure = ER_neg_BCa$exposure,
                                   outcome = ER_neg_BCa$outcome,
                                   snps = ER_neg_BCa$RSID)
    
    ER_neg_BCa_MR_input_TS = data.frame(SNP = ER_neg_BCa$RSID,
                                        effect_allele.exposure = ER_neg_BCa$A1,
                                        other_allele.exposure = ER_neg_BCa$A2,
                                        effect_allele.outcome = ER_neg_BCa$A1,
                                        other_allele.outcome = ER_neg_BCa$A2,
                                        beta.exposure = ER_neg_BCa$Effect_exposure,
                                        beta.outcome = ER_neg_BCa$Effect_outcome,
                                        eaf.exposure = rep(NA, nrow(ER_neg_BCa)),
                                        eaf.outcome = rep(NA, nrow(ER_neg_BCa)),
                                        remove = rep(FALSE, nrow(ER_neg_BCa)),
                                        palindromic = rep(FALSE, nrow(ER_neg_BCa)),
                                        ambigouous = rep(NA, nrow(ER_neg_BCa)),
                                        id.outcome = ER_neg_BCa$outcome,
                                        chr.outcome = ER_neg_BCa$Chr,
                                        pos.outcome = ER_neg_BCa$Pos,
                                        pval.outcome = ER_neg_BCa$P_outcome,
                                        se.outcome = ER_neg_BCa$StdErr_outcome,
                                        outcome = rep("ER_neg_BCa", nrow(ER_neg_BCa)),
                                        mr_keep.outcome = rep(TRUE, nrow(ER_neg_BCa)),
                                        pval_origin.outcome = rep("reported", nrow(ER_neg_BCa)),
                                        id.exposure = rep("ER_neg_BCa", nrow(ER_neg_BCa)),
                                        chr.exposure = ER_neg_BCa$Chr,
                                        pos.exposure = ER_neg_BCa$Pos,
                                        pval.exposure = ER_neg_BCa$P_exposure,
                                        se.exposure = ER_neg_BCa$StdErr_exposure,
                                        exposure = ER_neg_BCa$exposure,
                                        mr_keep.exposure = rep(TRUE, nrow(ER_neg_BCa)),
                                        pval_origin.exposure = rep("reported", nrow(ER_neg_BCa)),
                                        data_source.exposure = rep("IMD_GWAS", nrow(ER_neg_BCa)),
                                        action = rep(2, nrow(ER_neg_BCa)),
                                        SNP_index = rep(1, nrow(ER_neg_BCa)),
                                        mr_keep = rep(TRUE, nrow(ER_neg_BCa)),
                                        samplesize.outcome = rep(72261, nrow(ER_neg_BCa)))
    
    TNBCa_MR_input = mr_input(bx = TNBCa$Effect_exposure,
                              bxse = TNBCa$StdErr_exposure,
                              by = TNBCa$Effect_outcome,
                              byse = TNBCa$StdErr_outcome,
                              exposure = TNBCa$exposure,
                              outcome = TNBCa$outcome,
                              snps = TNBCa$RSID)
    
    TNBCa_MR_input_TS = data.frame(SNP = TNBCa$RSID,
                                   effect_allele.exposure = TNBCa$A1,
                                   other_allele.exposure = TNBCa$A2,
                                   effect_allele.outcome = TNBCa$A1,
                                   other_allele.outcome = TNBCa$A2,
                                   beta.exposure = TNBCa$Effect_exposure,
                                   beta.outcome = TNBCa$Effect_outcome,
                                   eaf.exposure = rep(NA, nrow(TNBCa)),
                                   eaf.outcome = rep(NA, nrow(TNBCa)),
                                   remove = rep(FALSE, nrow(TNBCa)),
                                   palindromic = rep(FALSE, nrow(TNBCa)),
                                   ambigouous = rep(NA, nrow(TNBCa)),
                                   id.outcome = TNBCa$outcome,
                                   chr.outcome = TNBCa$Chr,
                                   pos.outcome = TNBCa$Pos,
                                   pval.outcome = TNBCa$P_outcome,
                                   se.outcome = TNBCa$StdErr_outcome,
                                   outcome = rep("TNBCa", nrow(TNBCa)),
                                   mr_keep.outcome = rep(TRUE, nrow(TNBCa)),
                                   pval_origin.outcome = rep("reported", nrow(TNBCa)),
                                   id.exposure = rep("TNBCa", nrow(TNBCa)),
                                   chr.exposure = TNBCa$Chr,
                                   pos.exposure = TNBCa$Pos,
                                   pval.exposure = TNBCa$P_exposure,
                                   se.exposure = TNBCa$StdErr_exposure,
                                   exposure = TNBCa$exposure,
                                   mr_keep.exposure = rep(TRUE, nrow(TNBCa)),
                                   pval_origin.exposure = rep("reported", nrow(TNBCa)),
                                   data_source.exposure = rep("IMD_GWAS", nrow(TNBCa)),
                                   action = rep(2, nrow(TNBCa)),
                                   SNP_index = rep(1, nrow(TNBCa)),
                                   mr_keep = rep(TRUE, nrow(TNBCa)),
                                   samplesize.outcome = rep(106377, nrow(TNBCa)))
    
    ## Run IVW_MR analyses 
    if(nrow(Overall_BCa)< 3){
      Overall_BCa_IVW = run_MR_ivw(Overall_BCa_MR_input,Overall_BCa_MR_input_TS)
      HER2_pos_BCa_IVW = run_MR_ivw(HER2_pos_BCa_MR_input,HER2_pos_BCa_MR_input_TS)
      Luminal_A_BCa_IVW = run_MR_ivw(Luminal_A_BCa_MR_input, Luminal_A_BCa_MR_input_TS)
      Luminal_B_BCa_IVW = run_MR_ivw(Luminal_B_BCa_MR_input, Luminal_B_BCa_MR_input_TS)
      ER_neg_BCa_IVW = run_MR_ivw(ER_neg_BCa_MR_input, ER_neg_BCa_MR_input_TS)
      TNBCa_IVW = run_MR_ivw(TNBCa_MR_input,TNBCa_MR_input_TS)
      exposure_MR_results = rbind.data.frame(Overall_BCa_IVW,HER2_pos_BCa_IVW,
                                             Luminal_A_BCa_IVW,Luminal_B_BCa_IVW,
                                             ER_neg_BCa_IVW,TNBCa_IVW)
    }else{
      safe_run <- function(expr) {
        tryCatch(
          expr,
          error = function(e) NA
        )
      }
      
      Overall_BCa_IVW <- safe_run(run_MR_ivw(Overall_BCa_MR_input, Overall_BCa_MR_input_TS))
      HER2_pos_BCa_IVW <- safe_run(run_MR_ivw(HER2_pos_BCa_MR_input, HER2_pos_MR_input_TS))
      Luminal_A_BCa_IVW <- safe_run(run_MR_ivw(Luminal_A_BCa_MR_input, Luminal_A_BCa_MR_input_TS))
      Luminal_B_BCa_IVW <- safe_run(run_MR_ivw(Luminal_B_BCa_MR_input, Luminal_B_BCa_MR_input_TS))
      ER_neg_BCa_IVW <- safe_run(run_MR_ivw(ER_neg_BCa_MR_input, ER_neg_BCa_MR_input_TS))
      TNBCa_IVW <- safe_run(run_MR_ivw(TNBCa_IVW_MR_input, TNBCa_IVW_MR_input_TS))
      
      Overall_BCa_Egger <- safe_run(run_MR_Egger(Overall_BCa_MR_input))
      HER2_pos_BCa_Egger <- safe_run(run_MR_Egger(HER2_pos_BCa_MR_input))
      Luminal_A_BCa_Egger <- safe_run(run_MR_Egger(Luminal_A_BCa_MR_input))
      Luminal_B_BCa_Egger <- safe_run(run_MR_Egger(Luminal_B_BCa_MR_input))
      ER_neg_BCa_IVW_Egger <- safe_run(run_MR_Egger(ER_neg_BCa_IVW_MR_input))
      TNBCa_Egger <- safe_run(run_MR_Egger(TNBCa_MR_input))
      
      Overall_BCa_Median <- safe_run(run_MR_median(Overall_BCa_MR_input))
      HER2_pos_BCa_Median <- safe_run(run_MR_median(HER2_pos_BCa_MR_input))
      Luminal_A_BCa_Median <- safe_run(run_MR_median(Luminal_A_BCa_MR_input))
      Luminal_B_BCa_Median <- safe_run(run_MR_median(Luminal_B_BCa_MR_input))
      ER_neg_BCa_Median <- safe_run(run_MR_median(ER_neg_BCa_MR_input))
      TNBCa_Median <- safe_run(run_MR_median(TNBCa_MR_input))
      
      Overall_BCa_LASSO <- safe_run(run_MR_LASSO(Overall_BCa_MR_input))
      HER2_pos_BCa_LASSO <- safe_run(run_MR_LASSO(HER2_pos_BCa_MR_input))
      Luminal_A_BCa_LASSO <- safe_run(run_MR_LASSO(Luminal_A_BCa_MR_input))
      Luminal_B_BCa_LASSO <- safe_run(run_MR_LASSO(Luminal_B_BCa_MR_input))
      ER_neg_BCa_LASSO <- safe_run(run_MR_LASSO(ER_neg_BCa_MR_input))
      TNBCa_LASSO <- safe_run(run_MR_LASSO(TNBCa_MR_input))
      
      Overall_BCa_cML = run_MR_cML(Overall_BCa_MR_input, n_out = 247173, n_snp_valid = Overall_BCa_LASSO$LASSO_N_valid)
      HER2_pos_BCa_cML = run_MR_cML(HER2_pos_BCa_MR_input, n_out = 97877, n_snp_valid = HER2_pos_BCa_LASSO$LASSO_N_valid)
      Luminal_A_BCa_cML = run_MR_cML(Luminal_A_BCa_MR_input, n_out = 148877, n_snp_valid = Luminal_A_BCa_LASSO$LASSO_N_valid)
      Luminal_B_BCa_cML = run_MR_cML(Luminal_B_BCa_MR_input, n_out = 105277, n_snp_valid = Luminal_B_BCa_LASSO$LASSO_N_valid)
      ER_neg_BCa_cML = run_MR_cML(ER_neg_BCa_MR_input, n_out = 72261, n_snp_valid = ER_neg_BCa_LASSO$LASSO_N_valid)
      TNBCa_cML = run_MR_cML(TNBCa_MR_input, n_out = 106377, n_snp_valid = TNBCa_LASSO$LASSO_N_valid)
      
      Overall_BCa_PRESSO = run_MR_PRESSO(Overall_BCa_MR_input_TS)
      HER2_pos_BCa_PRESSO = run_MR_PRESSO(HER2_pos_BCa_MR_input_TS)
      Luminal_A_BCa_MR_PRESSO = run_MR_PRESSO(Luminal_A_BCa_MR_input_TS)
      Luminal_B_BCa_MR_PRESSO = run_MR_PRESSO(Luminal_B_BCa_MR_input_TS)
      ER_neg_BCa_PRESSO = run_MR_PRESSO(ER_neg_BCa_MR_input_TS)
      TNBCa_PRESSO = run_MR_PRESSO(TNBCa_MR_input_TS)
      
      
    }
    
    obj_names <- c(
      "Overall_BCa_IVW", "Overall_BCa_Egger", "Overall_BCa_Median", "Overall_BCa_LASSO", "Overall_BCa_cML","Overall_BCa_PRESSO",
      "HER2_pos_BCa_IVW", "HER2_pos_BCa_Egger", "HER2_pos_BCa_Median", "HER2_pos_BCa_LASSO", "HER2_pos_BCa_cML","HER2_pos_BCa_PRESSO",
      "Luminal_A_BCa_IVW", "Luminal_A_BCa_Egger", "Luminal_A_BCa_Median", "Luminal_A_BCa_LASSO", "Luminal_A_BCa_cML","Luminal_A_BCa_MR_PRESSO",
      "Luminal_B_BCa_IVW", "Luminal_B_BCa_Egger", "Luminal_B_BCa_Median", "Luminal_B_BCa_LASSO", "Luminal_B_BCa_cML","Luminal_B_BCa_MR_PRESSO",
      "ER_neg_BCa_IVW", "ER_neg_BCa_Egger", "ER_neg_BCa_Median", "ER_neg_BCa_LASSO", "ER_neg_BCa_cML","ER_neg_BCa_PRESSO",
      "TNBCa_IVW", "TNBCa_Egger", "TNBCa_Median", "TNBCa_LASSO", "TNBCa_cML","TNBCa_PRESSO"
    )
    
    existing_objs <- mget(obj_names[sapply(obj_names, exists)], inherits = TRUE)
    
    exposure_MR_results <- do.call(rbind_fill, existing_objs)
    MR_results = rbind_fill(MR_results,exposure_MR_results)

}

MR_results = MR_results[-1,]


# Save results

fwrite(x = MR_results, file = "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_BCa_MR_results.tsv",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
)

xlsx::write.xlsx(MR_results, "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_BCa_MR_results.xlsx", row.names = FALSE,
                 showNA = FALSE)
