# Define function
run_MVMR_pipeline <- function(
    file_paths,
    exposure_names,
    outcome_names,
    exp_sample_sizes,
    out_sample_sizes,
    output_prefix
) {
  
  library(data.table)
  library(dplyr)
  library(MendelianRandomization)
  library(MVMR)
  library(MRPRESSO)
  
  ### Helper functions

  rbind_fill <- function(..., fill = NA) {
    dfs <- list(...)
    all_cols <- unique(unlist(lapply(dfs, names)))
    dfs_filled <- lapply(dfs, function(df) {
      missing_cols <- setdiff(all_cols, names(df))
      if (length(missing_cols) > 0) df[missing_cols] <- fill
      df[all_cols]
    })
    do.call(rbind, dfs_filled)
  }
  
  ### Load and instruments

  data_list <- lapply(file_paths, readRDS)
  
  ### Build MVMR inputs

  MVMR_inputs <- lapply(seq_along(data_list), function(i) {
    dat <- data_list[[i]]
    
    # Align exposure with names and create MVMR input
    
    exp_order = colnames(dat$exposure_beta)
    exp_map = dat$expname
    
    exposure = c(exp_map %>%
                   filter(id.exposure == exp_order[1]) %>%
                   pull(exposure),
                 exp_map %>%
                   filter(id.exposure == exp_order[2]) %>%
                   pull(exposure))
    
    
    mr_mvinput(
      bx = dat$exposure_beta,
      bxse = dat$exposure_se,
      by = dat$outcome_beta,
      byse = dat$outcome_se,
      exposure = exposure,
      outcome = outcome_names[i],
      snps = rownames(dat$exposure_beta)
    )
    
  })
  
  ### Run MVMR Methods

  results_list <- list()
  
  for (i in seq_along(MVMR_inputs)) {
    
    input <- MVMR_inputs[[i]]
    outcome <- outcome_names[i]
    snp_label <- paste0(exposure_names[i], "_BMI_instrument")
    
    ### IVW

    ### Compute conditional F-statistics using MVMR package
    Fstat <- tryCatch({
      
      mvmr_formatted <- MVMR::format_mvmr(
        BXGs = input@betaX,
        BYG = input@betaY,
        seBXGs = input@betaXse,
        seBYG = input@betaYse,
        RSID = rownames(input@betaX)
      )
      
      strength_res <- MVMR::strength_mvmr(mvmr_formatted)
      
      # Extract conditional F-statistics
      data.frame(Fstatistic = c(strength_res$exposure1,strength_res$exposure2)) %>%
        arrange(Fstatistic) %>%
      mutate(Exposure = c(exposure_names[i],"BMI"))

    }, error = function(e) {
      rep(NA, ncol(input@betaX))  # one per exposure
    })
    
    ### IVW
    ivw_df <- tryCatch({
      
      ivw <- mr_mvivw(
        input,
        robust = TRUE,
        nx = exp_sample_sizes[c(2*i - 1, 2*i)]
      )
      
      data.frame(
        Exposure = ivw@Exposure,
        Outcome = outcome,
        SNP =  snp_label,
        Estimate = ivw@Estimate,
        StdError = ivw@StdError,
        CILower = ivw@CILower,
        CIUpper = ivw@CIUpper,
        Pvalue = ivw@Pvalue,
        NSNPs = ivw@SNPs,
        Method = "MVMR IVW",
        Hetero_statistic = ivw@Heter.Stat[1],
        Hetero_Pvalue = ivw@Heter.Stat[2]
      )
      
    }, error = function(e) {
      data.frame(
        Exposure = NA, Outcome = outcome, SNP = paste0(exposure_names[i],"_BMI_instrument"),
        Estimate = NA, StdError = NA, CILower = NA, CIUpper = NA,
        Pvalue = NA, Fstatistic = NA, NSNPs = NA,
        Method = "MVMR IVW",
        Hetero_statistic = NA, Hetero_Pvalue = NA
      )
    })
    
    ivw_df = merge(ivw_df, Fstat, by= "Exposure")
    
    ### Egger
    egger <- mr_mvegger(input)
    
    egger_df <- data.frame(
      Exposure = egger@Exposure,
      Outcome = outcome,
      SNP = snp_label,
      Estimate = egger@Estimate,
      StdError = egger@StdError.Est,
      CILower = egger@CILower.Est,
      CIUpper = egger@CIUpper.Est,
      Pvalue = egger@Pvalue.Est,
      NSNPs = egger@SNPs,
      Method = "MVMR Egger",
      Egger_intercept_p = egger@Pvalue.Int,
      Hetero_statistic = egger@Heter.Stat[1],
      Hetero_Pvalue = egger@Heter.Stat[2]
    )
    
    ### Median
    median <- mr_mvmedian(input)
    
    median_df <- data.frame(
      Exposure = median@Exposure,
      Outcome = outcome,
      SNP = snp_label,
      Estimate = median@Estimate,
      StdError = median@StdError,
      CILower = median@CILower,
      CIUpper = median@CIUpper,
      Pvalue = median@Pvalue,
      NSNPs = median@SNPs,
      Method = "MVMR median"
    )
    
    ### LASSO
    lasso <- mr_mvlasso(input)
    
    lasso_df <- data.frame(
      Exposure = lasso@Exposure,
      Outcome = outcome,
      SNP = snp_label,
      Estimate = lasso@Estimate,
      StdError = lasso@StdError,
      CILower = lasso@CILower,
      CIUpper = lasso@CIUpper,
      Pvalue = lasso@Pvalue,
      NSNPs = lasso@SNPs,
      Method = "MVMR LASSO",
      Valid_SNPs = paste(lasso@ValidSNPs, collapse = ",")
    )
    
    ### cML
    cml <- mr_mvcML(
      input,
      n = out_sample_sizes[i],
      K_vec = 0:(lasso@SNPs - lasso@Valid) * 2
    )
    
    cml_df <- data.frame(
      Exposure = cml@Exposure,
      Outcome = outcome,
      SNP = snp_label,
      Estimate = cml@Estimate,
      StdError = cml@StdError,
      CILower = cml@CILower,
      CIUpper = cml@CIUpper,
      Pvalue = cml@Pvalue,
      NSNPs = cml@SNPs,
      Method = "MVMR cML"
    )
    
    ### GMM
    gmm <- mr_mvgmm(input, nx = exp_sample_sizes[c(2*i - 1,2*i)], ny = out_sample_sizes[i])
    
    gmm_df <- data.frame(
      Exposure = gmm@Exposure,
      Outcome = outcome,
      SNP = snp_label,
      Estimate = gmm@Estimate,
      StdError = gmm@StdError,
      CILower = gmm@CILower,
      CIUpper = gmm@CIUpper,
      Pvalue = gmm@Pvalue,
      Method = "MVMR GMM"
    )
    
    results_list[[i]] <- rbind_fill(ivw_df, egger_df, median_df, lasso_df, cml_df, gmm_df)
  }
  
  ### Combine and save

  final_results <- do.call(rbind, results_list)
  
  fwrite(final_results, paste0(output_prefix, ".tsv"), sep = "\t")
  }

# Define input

file_paths <- c(
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/OpenGWAS_project_CIMs/PMD_TNFR2_MVMR_data.rds",
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/OpenGWAS_project_CIMs/DA_TNFR2_MVMR_data.rds",
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/OpenGWAS_project_CIMs/NDA_TNFR2_MVMR_data.rds"
)

outcome_names <- c("PMD", "DA", "NDA")

exposure_names <- c(rep("TNFR2",3))

exp_sample_sizes <- c(rep(c(21758, 461460),3))

out_sample_sizes <- c(rep(c(17841, 13965, 14036),3))

output_prefix = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_project_CIM_BMI_MVMR_results"

# Run

run_MVMR_pipeline(file_paths = file_paths,
                  outcome_names = outcome_names,
                  out_sample_sizes = out_sample_sizes,
                  exposure_names = exposure_names,
                  exp_sample_sizes = exp_sample_sizes,
                  output_prefix = output_prefix)
