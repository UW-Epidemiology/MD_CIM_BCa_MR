library(data.table)
library(R.utils)
library(tidyverse)
library(MendelianRandomization)
library(TwoSampleMR)

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


# Run MR analyses

dir_path <- "/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/harmonized_fully_filtered/"

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
  
  PMD = data %>% 
    filter(outcome == "PMD")
  
  DA = data %>% 
    filter(outcome == "DA")
  
  NDA = data %>% 
    filter(outcome == "NDA")
  
  PMD_MR_input = mr_input(bx = PMD$Effect_exposure,
                          bxse = PMD$StdErr_exposure,
                          by = PMD$Effect_outcome,
                          byse = PMD$StdErr_outcome,
                          exposure = PMD$exposure,
                          outcome = PMD$outcome,
                          snps = PMD$RSID)
  
  PMD_MR_input_TS = data.frame(SNP = PMD$RSID,
                            effect_allele.exposure = PMD$A1,
                            other_allele.exposure = PMD$A2,
                            effect_allele.outcome = PMD$A1,
                            other_allele.outcome = PMD$A2,
                            beta.exposure = PMD$Effect_exposure,
                            beta.outcome = PMD$Effect_outcome,
                            eaf.exposure = rep(NA, nrow(PMD)),
                            eaf.outcome = rep(NA, nrow(PMD)),
                            remove = rep(FALSE, nrow(PMD)),
                            palindromic = rep(FALSE, nrow(PMD)),
                            ambigouous = rep(NA, nrow(PMD)),
                            id.outcome = PMD$outcome,
                            chr.outcome = PMD$Chr,
                            pos.outcome = PMD$Pos,
                            pval.outcome = PMD$P_outcome,
                            se.outcome = PMD$StdErr_outcome,
                            outcome = rep("PMD", nrow(PMD)),
                            mr_keep.outcome = rep(TRUE, nrow(PMD)),
                            pval_origin.outcome = rep("reported", nrow(PMD)),
                            id.exposure = rep("PMD", nrow(PMD)),
                            chr.exposure = PMD$Chr,
                            pos.exposure = PMD$Pos,
                            pval.exposure = PMD$P_exposure,
                            se.exposure = PMD$StdErr_exposure,
                            exposure = PMD$exposure,
                            mr_keep.exposure = rep(TRUE, nrow(PMD)),
                            pval_origin.exposure = rep("reported", nrow(PMD)),
                            data_source.exposure = rep("IMD_GWAS", nrow(PMD)),
                            action = rep(2, nrow(PMD)),
                            SNP_index = rep(1, nrow(PMD)),
                            mr_keep = rep(TRUE, nrow(PMD)),
                            samplesize.outcome = rep(17841, nrow(PMD)))
                            
  DA_MR_input = mr_input(bx = DA$Effect_exposure,
                         bxse = DA$StdErr_exposure,
                         by = DA$Effect_outcome,
                         byse = DA$StdErr_outcome,
                         exposure = DA$exposure,
                         outcome = DA$outcome,
                         snps = DA$RSID)
  
  DA_MR_input_TS = data.frame(SNP = DA$RSID,
                               effect_allele.exposure = DA$A1,
                               other_allele.exposure = DA$A2,
                               effect_allele.outcome = DA$A1,
                               other_allele.outcome = DA$A2,
                               beta.exposure = DA$Effect_exposure,
                               beta.outcome = DA$Effect_outcome,
                               eaf.exposure = rep(NA, nrow(DA)),
                               eaf.outcome = rep(NA, nrow(DA)),
                               remove = rep(FALSE, nrow(DA)),
                               palindromic = rep(FALSE, nrow(DA)),
                               ambigouous = rep(NA, nrow(DA)),
                               id.outcome = DA$outcome,
                               chr.outcome = DA$Chr,
                               pos.outcome = DA$Pos,
                               pval.outcome = DA$P_outcome,
                               se.outcome = DA$StdErr_outcome,
                               outcome = rep("DA", nrow(DA)),
                               mr_keep.outcome = rep(TRUE, nrow(DA)),
                               pval_origin.outcome = rep("reported", nrow(DA)),
                               id.exposure = rep("DA", nrow(DA)),
                               chr.exposure = DA$Chr,
                               pos.exposure = DA$Pos,
                               pval.exposure = DA$P_exposure,
                               se.exposure = DA$StdErr_exposure,
                               exposure = DA$exposure,
                               mr_keep.exposure = rep(TRUE, nrow(DA)),
                               pval_origin.exposure = rep("reported", nrow(DA)),
                               data_source.exposure = rep("IMD_GWAS", nrow(DA)),
                               action = rep(2, nrow(DA)),
                               SNP_index = rep(1, nrow(DA)),
                               mr_keep = rep(TRUE, nrow(DA)),
                               samplesize.outcome = rep(13965, nrow(DA)))
  
  NDA_MR_input = mr_input(bx = NDA$Effect_exposure,
                          bxse = NDA$StdErr_exposure,
                          by = NDA$Effect_outcome,
                          byse = NDA$StdErr_outcome,
                          exposure = NDA$exposure,
                          outcome = NDA$outcome,
                          snps = NDA$RSID)
  
  NDA_MR_input_TS = data.frame(SNP = NDA$RSID,
                               effect_allele.exposure = NDA$A1,
                               other_allele.exposure = NDA$A2,
                               effect_allele.outcome = NDA$A1,
                               other_allele.outcome = NDA$A2,
                               beta.exposure = NDA$Effect_exposure,
                               beta.outcome = NDA$Effect_outcome,
                               eaf.exposure = rep(NA, nrow(NDA)),
                               eaf.outcome = rep(NA, nrow(NDA)),
                               remove = rep(FALSE, nrow(NDA)),
                               palindromic = rep(FALSE, nrow(NDA)),
                               ambigouous = rep(NA, nrow(NDA)),
                               id.outcome = NDA$outcome,
                               chr.outcome = NDA$Chr,
                               pos.outcome = NDA$Pos,
                               pval.outcome = NDA$P_outcome,
                               se.outcome = NDA$StdErr_outcome,
                               outcome = rep("NDA", nrow(NDA)),
                               mr_keep.outcome = rep(TRUE, nrow(NDA)),
                               pval_origin.outcome = rep("reported", nrow(NDA)),
                               id.exposure = rep("NDA", nrow(NDA)),
                               chr.exposure = NDA$Chr,
                               pos.exposure = NDA$Pos,
                               pval.exposure = NDA$P_exposure,
                               se.exposure = NDA$StdErr_exposure,
                               exposure = NDA$exposure,
                               mr_keep.exposure = rep(TRUE, nrow(NDA)),
                               pval_origin.exposure = rep("reported", nrow(NDA)),
                               data_source.exposure = rep("IMD_GWAS", nrow(NDA)),
                               action = rep(2, nrow(NDA)),
                               SNP_index = rep(1, nrow(NDA)),
                               mr_keep = rep(TRUE, nrow(NDA)),
                               samplesize.outcome = rep(14036, nrow(NDA)))
  
  ## Run IVW_MR analyses 
  if(nrow(PMD)< 3){
    PMD_IVW = run_MR_ivw(PMD_MR_input,PMD_MR_input_TS)
    DA_IVW = run_MR_ivw(DA_MR_input, DA_MR_input_TS)
    NDA_IVW = run_MR_ivw(NDA_MR_input, NDA_MR_input_TS)
    exposure_MR_results = rbind.data.frame(PMD_IVW,DA_IVW,NDA_IVW)
  }else{
    safe_run <- function(expr) {
      tryCatch(
        expr,
        error = function(e) NA
      )
    }
    
    PMD_IVW <- safe_run(run_MR_ivw(PMD_MR_input, PMD_MR_input_TS))
    DA_IVW <- safe_run(run_MR_ivw(DA_MR_input, DA_MR_input_TS))
    NDA_IVW <- safe_run(run_MR_ivw(NDA_MR_input, NDA_MR_input_TS))
    
    PMD_Egger <- safe_run(run_MR_Egger(PMD_MR_input))
    DA_Egger <- safe_run(run_MR_Egger(DA_MR_input))
    NDA_Egger <- safe_run(run_MR_Egger(NDA_MR_input))
    
    PMD_Median <- safe_run(run_MR_median(PMD_MR_input))
    DA_Median <- safe_run(run_MR_median(DA_MR_input))
    NDA_Median <- safe_run(run_MR_median(NDA_MR_input))
    
    PMD_LASSO <- safe_run(run_MR_LASSO(PMD_MR_input))
    DA_LASSO <- safe_run(run_MR_LASSO(DA_MR_input))
    NDA_LASSO <- safe_run(run_MR_LASSO(NDA_MR_input))
    
    PMD_cML <- safe_run(
      run_MR_cML(
        PMD_MR_input,
        n_out = 17841,
        n_snp_valid = ifelse(is.list(PMD_LASSO), PMD_LASSO$LASSO_N_valid, NA)
      )
    )
    
    DA_cML <- safe_run(
      run_MR_cML(
        DA_MR_input,
        n_out = 13965,
        n_snp_valid = ifelse(is.list(DA_LASSO), DA_LASSO$LASSO_N_valid, NA)
      )
    )
    
    NDA_cML <- safe_run(
      run_MR_cML(
        NDA_MR_input,
        n_out = 14036,
        n_snp_valid = ifelse(is.list(NDA_LASSO), NDA_LASSO$LASSO_N_valid, NA)
      )
    )
    
    obj_names <- c(
      "PMD_IVW", "PMD_Egger", "PMD_Median", "PMD_LASSO", "PMD_cML",
      "DA_IVW", "DA_Egger", "DA_Median", "DA_LASSO", "DA_cML",
      "NDA_IVW", "NDA_Egger", "NDA_Median", "NDA_LASSO", "NDA_cML"
    )
    
    existing_objs <- mget(obj_names[sapply(obj_names, exists)], inherits = TRUE)
    
    exposure_MR_results <- do.call(rbind_fill, existing_objs)
    
  }
  
 MR_results = rbind_fill(MR_results,exposure_MR_results)
  
}

MR_results = MR_results[-1,]

# Save results

fwrite(x = MR_results, file = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_GWAS_CIM_MR_results.tsv",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
)

xlsx::write.xlsx(MR_results, "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_GWAS_CIM_MR_results.xlsx", row.names = FALSE,
                 showNA = FALSE)
