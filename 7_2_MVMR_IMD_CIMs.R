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
    
    ## Match exposure names to the order of columns in input dataset
    exposure_names <- dat$expname$exposure[
      match(colnames(dat$exposure_beta), dat$expname$id.exposure)
    ]
    
    mr_mvinput(
      bx = dat$exposure_beta,
      bxse = dat$exposure_se,
      by = dat$outcome_beta,
      byse = dat$outcome_se,
      exposure = exposure_names,
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
    
    exposure = exposure_names[i]
    print(paste0("Running MVMR for ", exposure ," and ", outcome))
    
    
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
    egger_df <- tryCatch({
      egger <- mr_mvegger(input)
      
      data.frame(
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
    }, error = function(e) {
      data.frame(
        Exposure = NA, Outcome = outcome, SNP = snp_label,
        Estimate = NA, StdError = NA, CILower = NA, CIUpper = NA,
        Pvalue = NA, NSNPs = NA,
        Method = "MVMR Egger",
        Egger_intercept_p = NA,
        Hetero_statistic = NA, Hetero_Pvalue = NA
      )
    })
    
    ### Median
    median_df <- tryCatch({
      median <- mr_mvmedian(input)
      
      data.frame(
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
    }, error = function(e) {
      data.frame(
        Exposure = NA, Outcome = outcome, SNP = snp_label,
        Estimate = NA, StdError = NA, CILower = NA, CIUpper = NA,
        Pvalue = NA, NSNPs = NA,
        Method = "MVMR median"
      )
    })
    
    ### LASSO
    lasso_df <- tryCatch({
      lasso <- mr_mvlasso(input)
      
      data.frame(
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
    }, error = function(e) {
      data.frame(
        Exposure = NA, Outcome = outcome, SNP = snp_label,
        Estimate = NA, StdError = NA, CILower = NA, CIUpper = NA,
        Pvalue = NA, NSNPs = NA,
        Method = "MVMR LASSO",
        Valid_SNPs = NA
      )
    })
    
    ### cML
    cml_df <- tryCatch({
      lasso <- mr_mvlasso(input)  # needed for K_vec
      cml <- mr_mvcML(
        input,
        n = out_sample_sizes[i],
        K_vec = 0:(lasso@SNPs - lasso@Valid) * 2
      )
      
      data.frame(
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
    }, error = function(e) {
      data.frame(
        Exposure = NA, Outcome = outcome, SNP = snp_label,
        Estimate = NA, StdError = NA, CILower = NA, CIUpper = NA,
        Pvalue = NA, NSNPs = NA,
        Method = "MVMR cML"
      )
    })
    
    ### GMM
    gmm <- tryCatch({
      
      mr_mvgmm(
        input,
        nx = exp_sample_sizes[c(1:2)],
        ny = out_sample_sizes[i]
      )
      
    }, error = function(e) {
      
      cat("\n========== MR GMM ERROR ==========\n", file = stderr())
      cat("Iteration:", i, "\n", file = stderr())
      cat("Outcome:", outcome, "\n", file = stderr())
      cat("SNP label:", snp_label, "\n", file = stderr())
      cat("Error message:", conditionMessage(e), "\n", file = stderr())
      
      cat("\nSample sizes:\n", file = stderr())
      cat("nx:", paste(exp_sample_sizes[c(2*i - 1, 2*i)], collapse = ", "), "\n", file = stderr())
      cat("ny:", out_sample_sizes[i], "\n", file = stderr())
      
      cat("\nInput structure:\n", file = stderr())
      cat("Class:", class(input), "\n", file = stderr())
      cat("Slots:", paste(slotNames(input), collapse = ", "), "\n", file = stderr())
      
      for (slot in slotNames(input)) {
        obj <- tryCatch(slot(input, slot), error = function(x) NULL)
        
        if (!is.null(obj)) {
          cat("\n--- Slot:", slot, "---\n", file = stderr())
          cat("Class:", class(obj), "\n", file = stderr())
          
          if (is.matrix(obj) || is.data.frame(obj)) {
            cat("Dim:", paste(dim(obj), collapse = " x "), "\n", file = stderr())
            cat("NA:", sum(is.na(obj)), "\n", file = stderr())
            cat("NaN:", sum(is.nan(as.matrix(obj))), "\n", file = stderr())
            cat("Inf:", sum(is.infinite(as.matrix(obj))), "\n", file = stderr())
          } else if (is.numeric(obj)) {
            cat("Length:", length(obj), "\n", file = stderr())
            cat("NA:", sum(is.na(obj)), "\n", file = stderr())
            cat("NaN:", sum(is.nan(obj)), "\n", file = stderr())
            cat("Inf:", sum(is.infinite(obj)), "\n", file = stderr())
            cat("Range:", paste(range(obj, na.rm = TRUE), collapse = " to "), "\n", file = stderr())
          }
        }
      }
      
      cat("================================\n\n", file = stderr())
      
      return(NULL)
    })
    
    if (!is.null(gmm)) {
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
    } else {
      gmm_df <- NULL
    }    
    
    results_list[[i]] <- rbind_fill(
      ivw_df, egger_df, median_df, lasso_df, cml_df, gmm_df
    )
  }
  
  ### Combine and save
  
  final_results <- do.call(rbind, results_list)
  
  fwrite(final_results, paste0(output_prefix, ".tsv"), sep = "\t")
}

# Define input

## Detect file paths

dir_path <- "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/additional_MVMR_summary_statistics/IMD_CIMs/"


filepaths_PMD <- list.files(
  path = dir_path,
  pattern = "^PMD.*\\.rds$",
  full.names = TRUE
)

filepaths_DA <- list.files(
  path = dir_path,
  pattern = "^DA.*\\.rds$",
  full.names = TRUE
)

filepaths_NDA <- list.files(
  path = dir_path,
  pattern = "^NDA.*\\.rds$",
  full.names = TRUE
)

file_paths = c(filepaths_PMD,filepaths_DA,filepaths_NDA)

proteins <- sub("^PMD_(.*)_MVMR.*$", "\\1", basename(filepaths_PMD))

outcome_names <- c(rep("PMD",length(filepaths_PMD)),
                   rep("DA",length(filepaths_PMD)),
                   rep("NDA",length(filepaths_PMD)))

exposure_names <- c(rep(proteins,3))

exp_sample_sizes <- c(rep(c(14824, 461460),length(filepaths_PMD)))

out_sample_sizes <- c(rep(17841,length(filepaths_PMD)),rep(13965,length(filepaths_PMD)), rep(14036,length(filepaths_PMD)))
output_prefix = "/projects/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_CIM_BMI_MVMR_results"

# Run
run_MVMR_pipeline(file_paths = file_paths,
                  outcome_names = outcome_names,
                  out_sample_sizes = out_sample_sizes,
                  exposure_names = exposure_names,
                  exp_sample_sizes = exp_sample_sizes,
                  output_prefix = output_prefix)
