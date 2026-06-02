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
  mutate(Exposure = sub("\\.", "-", Exposure),
         Exposure = sub("VEGF-A", "VEGF", Exposure))

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

# Ensure Exposure ordering is alphabetical
data_PMD$Exposure <- factor(data_PMD$Exposure,
                            levels = rev(sort(unique(data_PMD$Exposure))))

data_DA$Exposure  <- factor(data_DA$Exposure,
                            levels = rev(sort(unique(data_DA$Exposure))))

data_NDA$Exposure <- factor(data_NDA$Exposure,
                            levels = rev(sort(unique(data_NDA$Exposure))))

# Get all unique methods across datasets
all_methods <- unique(c(data_PMD$Method,
                        data_DA$Method,
                        data_NDA$Method))

# Create default brewer palette
default_cols <- scales::brewer_pal(palette = "Set1")(length(all_methods))
names(default_cols) <- all_methods

# Force MR IVW to black
default_cols["MR IVW"] <- "black"

# Add function to check if error 95% CI for MR estimate crosses 0
add_cross_flag <- function(df) {
  df$crosses0 <- with(df, CILower <= 0 & CIUpper >= 0)
  df
}

# Run function 

data_PMD <- add_cross_flag(data_PMD)
data_DA  <- add_cross_flag(data_DA)
data_NDA <- add_cross_flag(data_NDA)


# Plots

pd <- position_dodge(width = 0.6)
shape_vals <- c(`FALSE` = 16, `TRUE` = 1)

## PMD
plot_PMD = ggplot(data_PMD,
                  aes(x = Estimate,
                      y = Exposure,
                      color = Method,
                      shape = crosses0)) +
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
    stroke = 1,
    position = pd
  ) +
  scale_color_manual(values = default_cols) +
  scale_shape_manual(values = shape_vals, guide = "none") +
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
plot_DA = ggplot(data_DA,
                 aes(x = Estimate,
                     y = Exposure,
                     color = Method,
                     shape = crosses0)) +
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
    stroke = 1,
    position = pd
  ) +
  scale_color_manual(values = default_cols) +
  scale_shape_manual(values = shape_vals, guide = "none") +
  labs(
    x = expression("MR estimate for DA (sqrt(cm"^2*"))"),
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
plot_NDA = ggplot(data_NDA,
                  aes(x = Estimate,
                      y = Exposure,
                      color = Method,
                      shape = crosses0)) +
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
    stroke = 1,
    position = pd
  ) +
  scale_color_manual(values = default_cols) +
  scale_shape_manual(values = shape_vals, guide = "none") +
  labs(
    x = expression("MR estimate for NDA (sqrt(cm"^2*"))"),
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

ggsave(combined, filename = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Figures/Figure_2.jpg",
       width = 12, height = 8, dpi = 300)

######################## Multivariable MR #######################################
CRP_BMI_MVMR = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/CRP_BMI_MVMR_results.tsv")

CRP_BMI_MVMR <- CRP_BMI_MVMR %>%
  mutate(
    Exposure = paste0(Exposure, "*")
  )

OpenGWAS_BMI_MVMR = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_project_CIM_BMI_MVMR_results.tsv")

OpenGWAS_BMI_MVMR <- OpenGWAS_BMI_MVMR %>%
  mutate(
    Exposure = paste0(Exposure, "*"),
    Method = gsub("MVMR Egger", "MVMR-Egger", Method)
  )

IMD_CIM_BMI_MVMR = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_CIM_BMI_MVMR_results.tsv") %>%
  mutate(
    Method = gsub("MVMR Egger", "MVMR-Egger", Method)
  )

MVMR = rbind_fill(CRP_BMI_MVMR,OpenGWAS_BMI_MVMR,IMD_CIM_BMI_MVMR)  

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

qhet = data.frame(Exposure = c("CRP*","BMI","CRP*","BMI","CRP*","BMI"),
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

data = rbind_fill(MVMR,qhet)

data = data %>%
  mutate(Exposure = sub("VEGF.A", "VEGF", Exposure),
         SNP = sub("VEGF.A_BMI_instrument","VEGF_BMI_instrument",SNP))


fwrite(x = data, file = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/BMI_MVMR_results_updated.tsv",
       sep = "\t",
       quote = FALSE,
       row.names = FALSE,
)

openxlsx::write.xlsx(
  data,
  "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/BMI_MVMR_results_updated.xlsx",
  rowNames = FALSE
)


data = data %>%
  filter(!Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected",
                        "MR-PRESSO raw"), 
         Exposure != "BMI*",
         Exposure != "BMI") %>%
  mutate(Exposure = sub("\\.", "-", Exposure)) %>%
  mutate(Method = case_when(Method == "MVMR Median" ~ "MVMR median",
                            Method == "Q_statistic_mimimization" ~ "MVMR Q-stat min",
                            TRUE ~ Method))

data_IVW = data %>% filter(Method == "MVMR IVW")

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
    x = "MVMR IVW estimate for PMD (sqrt(%))",
    y = "CIM exposure in MVMR (N SNPs in instrument)"
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
    x = expression("MVMR IVW estimate for DA (sqrt(cm"^2*"))"),
    y = "CIM exposure in MVMR (N SNPs in instrument)"
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
    x = expression("MVMR IVW estimate for NDA (sqrt(cm"^2*"))"),
    y = "CIM exposure in MVMR (N SNPs in instrument)"
  )+
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

## Combined plot

combined = plot_PMD + plot_DA + plot_NDA +
  plot_layout(ncol = 3)

#ggsave(combined, filename = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Figures/Figure_4.jpg",
       #width = 12, height = 8, dpi = 300)

# Sensitivity analyses for significant CIMs in MVMR analyses #

# Prepare for plotting
sig_CIMs = data_IVW %>% 
  filter(Pvalue < 0.05) %>%
  pull(Exposure) %>%
  unique()

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

# Ensure Exposure ordering is alphabetical
data_PMD$Exposure <- factor(data_PMD$Exposure,
                            levels = rev(sort(unique(data_PMD$Exposure))))

data_DA$Exposure  <- factor(data_DA$Exposure,
                            levels = rev(sort(unique(data_DA$Exposure))))

data_NDA$Exposure <- factor(data_NDA$Exposure,
                            levels = rev(sort(unique(data_NDA$Exposure))))

# Get all unique methods across datasets
all_methods <- unique(c(data_PMD$Method,
                        data_DA$Method,
                        data_NDA$Method))

# Create default brewer palette
default_cols <- scales::brewer_pal(palette = "Paired")(length(all_methods))
names(default_cols) <- all_methods

# Force MVMR IVW to black
default_cols["MVMR IVW"] <- "black"

# Add function to check if error 95% CI for MR estimate crosses 0
add_cross_flag <- function(df) {
  df$crosses0 <- with(df, CILower <= 0 & CIUpper >= 0)
  df
}

# Run function 

data_PMD <- add_cross_flag(data_PMD)
data_DA  <- add_cross_flag(data_DA)
data_NDA <- add_cross_flag(data_NDA)


# Plots

pd <- position_dodge(width = 0.6)
shape_vals <- c(`FALSE` = 16, `TRUE` = 1)

## PMD
plot_PMD = ggplot(data_PMD,
                  aes(x = Estimate,
                      y = Exposure,
                      color = Method,
                      shape = crosses0)) +
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
    stroke = 1,
    position = pd
  ) +
  scale_color_manual(values = default_cols) +
  scale_shape_manual(values = shape_vals, guide = "none") +
  labs(
    x = "MVMR estimate for PMD (sqrt(%))",
    y = "CIM exposure in MVMR",
    color = "MVMR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## DA
plot_DA = ggplot(data_DA,
                 aes(x = Estimate,
                     y = Exposure,
                     color = Method,
                     shape = crosses0)) +
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
    stroke = 1,
    position = pd
  ) +
  scale_color_manual(values = default_cols) +
  scale_shape_manual(values = shape_vals, guide = "none") +
  labs(
    x = expression("MVMR estimate for DA (sqrt(cm"^2*"))"),
    y = "CIM exposure in MVMR",
    color = "MVMR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## NDA
plot_NDA = ggplot(data_NDA,
                  aes(x = Estimate,
                      y = Exposure,
                      color = Method,
                      shape = crosses0)) +
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
    stroke = 1,
    position = pd
  ) +
  scale_color_manual(values = default_cols) +
  scale_shape_manual(values = shape_vals, guide = "none") +
  labs(
    x = expression("MVMR estimate for NDA (sqrt(cm"^2*"))"),
    y = "CIM exposure in MVMR",
    color = "MVMR method"
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

ggsave(combined, filename = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Figures/Figure_3.jpg",
       width = 12, height = 8, dpi = 300)
