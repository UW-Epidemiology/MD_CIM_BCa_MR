# Setup
library(tidyverse)
library(patchwork)
library(stringr)
library(data.table)

# Function 

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

######################## Univariable MR #######################################

# Load data 
IMD_data = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_GWAS_CIM_MR_results.tsv")
OpenGWAS_data = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_project_CIM_MR_results.tsv")

OpenGWAS_data <- OpenGWAS_data %>%
  mutate(
    Exposure = paste0(Exposure, "*")
  )

data = rbind_fill(IMD_data,OpenGWAS_data)

# Filter to IVW results and relevant columns
data_IVW = data %>%
  filter(Method == "MR IVW") %>%
  select(c(Exposure,Outcome,Estimate,CILower,CIUpper,Fstatistic,NSNPs))

# Modify  plotting

data_IVW = data_IVW %>%
  mutate(Exposure = sub("\\.", "-", Exposure))

data_PMD <- data_IVW %>%
  filter(Outcome == "PMD") %>%
  mutate(
    Exposure_label = paste0(
        Exposure,
        "  (", NSNPs ,")" 
      )
    ) %>%
  arrange(Estimate) %>%
  mutate(
    Exposure_label = factor(Exposure_label, levels = Exposure_label)
  )

data_DA <- data_IVW %>%
  filter(Outcome == "DA") %>%
  mutate(
    Exposure_label = paste0(
      Exposure,
      "  (", NSNPs ,")" 
    )
  ) %>%
  arrange(Estimate) %>%
  mutate(
    Exposure_label = factor(Exposure_label, levels = Exposure_label)
  )


data_NDA <- data_IVW %>%
  filter(Outcome == "NDA") %>%
  mutate(
    Exposure_label = paste0(
      Exposure,
      "  (", NSNPs ,")" 
    )
  ) %>%
  arrange(Estimate) %>%
  mutate(
    Exposure_label = factor(Exposure_label, levels = Exposure_label)
  )


# Plots

## PMD 

plot_PMD = ggplot(data_PMD, aes(x = Estimate, y = Exposure_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2
  ) +
  geom_point(size = 2) +
  labs(
    x = "MR IVW estimate for PMD (sqrt(%))",
    y = "MR instrument (N SNPs)"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


## DA

plot_DA = ggplot(data_DA, aes(x = Estimate, y = Exposure_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2
  ) +
  geom_point(size = 2) +
  labs(
    x = expression("MR IVW estimate for DA (sqrt(cm"^2*"))"),
    y = "MR instrument (N SNPs)"
  )+
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

## NDA
plot_NDA = ggplot(data_NDA, aes(x = Estimate, y = Exposure_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2
  ) +
  geom_point(size = 2) +
  labs(
    x = expression("MR IVW estimate for NDA (sqrt(cm"^2*"))"),
    y = "MR instrument (N SNPs)"
  )+
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

## Combined plot

combined = plot_PMD + plot_DA + plot_NDA +
  plot_layout(ncol = 3)

ggsave(combined, filename = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Figures/Figure_2.jpg",
       width = 12, height = 8, dpi = 300)


# Sensitivity analyses for significant CIMs #

# Filter to IVW results and relevant columns
data_IVW = data %>%
  filter(Method == "MR IVW") %>%
  select(c(Exposure,Outcome,Estimate,CILower,CIUpper,Fstatistic,NSNPs,Pvalue,Method)) 

# Prepare for plotting
sig_CIMs = data_IVW %>% 
  filter(Pvalue < 0.05) %>%
  pull(Exposure) %>%
  unique

data_sig = data %>%
  filter(Exposure %in% sig_CIMs, !Method %in% c("MR Wald Ratio",
                                                "MR-PRESSO distortion test",
                                                "MR-PRESSO global test",
                                                "MR-PRESSO outlier-corrected") ) %>%
  mutate(Exposure = sub("\\.", "-", Exposure)) %>%
  mutate(Method = case_when(Method == "MR Median" ~ "MR median",
                            TRUE ~ Method))

data_PMD <- data_sig %>%
  filter(Outcome == "PMD")

data_DA <- data_sig %>%
  filter(Outcome == "DA")


data_NDA <- data_sig %>%
  filter(Outcome == "NDA")

# Plots

pd <- position_dodge(width = 0.6)

## PMD
data_PMD <- data_PMD %>%
  mutate(
    Exposure = fct_reorder(
      Exposure,
      ifelse(Method == "MR IVW", Estimate, NA_real_),
      .fun = function(x) mean(x, na.rm = TRUE)
    )
  )

plot_PMD = ggplot(data_PMD,
                  aes(x = Estimate,
                      y = Exposure,
                      color = Method)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2,
    position = pd
  ) +
  geom_point(
    size = 2,
    position = pd
  ) +
  labs(
    x = "MR estimate for PMD (sqrt(%))",
    y = "MR instrument",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )


## DA
data_DA <- data_DA %>%
  mutate(
    Exposure = fct_reorder(
      Exposure,
      ifelse(Method == "MR IVW", Estimate, NA_real_),
      .fun = function(x) mean(x, na.rm = TRUE)
    )
  )

plot_DA = ggplot(data_DA,
                 aes(x = Estimate,
                     y = Exposure,
                     color = Method)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2,
    position = pd
  ) +
  geom_point(
    size = 2,
    position = pd
  ) +
  labs(
    x = "MR estimate for DA (sqrt(%))",
    y = "MR instrument",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )


## NDA
data_NDA <- data_NDA %>%
  mutate(
    Exposure = fct_reorder(
      Exposure,
      ifelse(Method == "MR IVW", Estimate, NA_real_),
      .fun = function(x) mean(x, na.rm = TRUE)
    )
  )

plot_NDA = ggplot(data_NDA,
                  aes(x = Estimate,
                      y = Exposure,
                      color = Method)) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2,
    position = pd
  ) +
  geom_point(
    size = 2,
    position = pd
  ) +
  labs(
    x = "MR estimate for NDA (sqrt(%))",
    y = "MR instrument",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## Combined

combined =
  plot_PMD + plot_DA + plot_NDA +
  plot_layout(ncol = 3, guides = "collect") &
  theme(
    legend.position = "top",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

combined

ggsave(combined, filename = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Figures/Figure_3.jpg",
       width = 12, height = 8, dpi = 300)




######################## Multivariable MR #######################################
CRP_BMI_MVMR = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_results.tsv")
PMD_qhet = read_rds("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_PMD_QHET.rds")
CI_bounds <- str_match(
  PMD_qhet$`95% CI`,
  "^(-?\\d*\\.?\\d+)-(-?\\d*\\.?\\d+)$"
)
PMD_qhet$CI_lower <- as.numeric(CI_bounds[, 2])
PMD_qhet$CI_upper <- as.numeric(CI_bounds[, 3])

DA_qhet = read_rds("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_DA_QHET.rds")
CI_bounds <- str_match(
  DA_qhet$`95% CI`,
  "^(-?\\d*\\.?\\d+)-(-?\\d*\\.?\\d+)$"
)
DA_qhet$CI_lower <- as.numeric(CI_bounds[, 2])
DA_qhet$CI_upper <- as.numeric(CI_bounds[, 3])

NDA_qhet = read_rds("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_NDA_QHET.rds")
CI_bounds <- str_match(
  NDA_qhet$`95% CI`,
  "^(-?\\d*\\.?\\d+)-(-?\\d*\\.?\\d+)$"
)
NDA_qhet$CI_lower <- as.numeric(CI_bounds[, 2])
NDA_qhet$CI_upper <- as.numeric(CI_bounds[, 3])


qhet = data.frame(Exposure = c("CRP","BMI","CRP","BMI","CRP","BMI"),
                  Outcome = c("PMD","PMD","DA","DA","NDA","NDA"),
                  SNP = rep("CRP_BMI_instrument",6),
                  Estimate = c(PMD_qhet$`Effect Estimates`[1], PMD_qhet$`Effect Estimates`[2],
                               DA_qhet$`Effect Estimates`[1], DA_qhet$`Effect Estimates`[2],
                               NDA_qhet$`Effect Estimates`[1], NDA_qhet$`Effect Estimates`[2]),
                  CILower = c(PMD_qhet$CI_lower[1], PMD_qhet$CI_lower[2],
                              DA_qhet$CI_lower[1], DA_qhet$CI_lower[2],
                              NDA_qhet$CI_lower[1], NDA_qhet$CI_lower[2]),
                  CIUpper = c(PMD_qhet$CI_upper[1], PMD_qhet$CI_upper[2],
                              DA_qhet$CI_upper[1], DA_qhet$CI_upper[2],
                              NDA_qhet$CI_upper[1], NDA_qhet$CI_upper[2]),
                  Method = rep("Q_statistic_mimimization",6))

data = rbind_fill(CRP_BMI_MVMR,qhet)

fwrite(x = data, file = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_results_updated.tsv",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
)

openxlsx::write.xlsx(
  data,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_results_updated.xlsx",
  rowNames = FALSE
)

# Plots 
data_IVW = data %>%
  filter(Method == "MVMR IVW") %>%
  select(c(Exposure,Outcome,Estimate,CILower,CIUpper,Fstatistic,NSNPs))

# Modify  plotting

data_IVW = data_IVW %>%
  mutate(Exposure = sub("\\.", "-", Exposure))

data_PMD <- data_IVW %>%
  filter(Outcome == "PMD") %>%
  arrange(Estimate) 

data_DA <- data_IVW %>%
  filter(Outcome == "DA") %>%
  arrange(Estimate)


data_NDA <- data_IVW %>%
  filter(Outcome == "NDA") %>%
  arrange(Estimate) 

# Plots

## PMD 

plot_PMD = ggplot(data_PMD, aes(x = Estimate, y = Exposure)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2
  ) +
  geom_point(size = 2) +
  labs(
    x = "IVW MVMR estimate for PMD (sqrt(%))",
    y = "Exposure"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

## DA

plot_DA = ggplot(data_DA, aes(x = Estimate, y = Exposure)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2
  ) +
  geom_point(size = 2) +
  labs(
    x = expression("IVW MVMR estimate for DA (sqrt(cm"^2*"))"),
    y = "Exposure"
  )+
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

## NDA
plot_NDA = ggplot(data_NDA, aes(x = Estimate, y = Exposure)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(
    aes(xmin = CILower, xmax = CIUpper),
    height = 0.2
  ) +
  geom_point(size = 2) +
  labs(
    x = expression("MVMR IVW estimate for NDA (sqrt(cm"^2*"))"),
    y = "Exposure"
  )+
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

## Combined plot

combined = plot_PMD + plot_DA + plot_NDA +
  plot_layout(nrow = 3)

ggsave(combined, filename = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Figures/Figure_4.jpg",
       width = 12, height = 8, dpi = 300)
