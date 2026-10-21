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

######################## MD plots #######################################

# Load data 
MD_data = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_GWAS_CIM_MR_results.tsv")
MD_OpenGWAS_data = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_project_CIM_MR_results.tsv")
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

MVMR_data = rbind_fill(MVMR,qhet)

MVMR_data = MVMR_data %>%
  mutate(Exposure = sub("VEGF.A", "VEGF", Exposure),
         SNP = sub("VEGF.A_BMI_instrument","VEGF_BMI_instrument",SNP))

MD_data = rbind_fill(MD_data,MD_OpenGWAS_data,MVMR_data)


# Sensitivity analyses for significant CIMs #

# Filter to IVW results and relevant columns
data_IVW_MD = MD_data %>%
  filter(Method %in% c("MR IVW", "MVMR IVW")) %>%
  select(c(Exposure,Outcome,Estimate,CILower,CIUpper,Fstatistic,NSNPs,Pvalue,Method)) 

# Prepare for plotting
sig_CIMs_MD = data_IVW_MD %>% 
  filter(Pvalue < 0.05) %>%
  pull(Exposure) %>%
  unique()


MD_data = MD_data %>%
  filter(Exposure %in% sig_CIMs_MD, !Method %in% c("MR Wald Ratio",
                                                   "MR-PRESSO distortion test",
                                                   "MR-PRESSO global test",
                                                   "MR-PRESSO outlier-corrected",
                                                   "MR-PRESSO raw") ) %>%
  mutate(Exposure = sub("\\.", "-", Exposure)) %>%
  mutate(Exposure = sub("\\*", "", Exposure)) %>%
  mutate(Method = case_when(Method == "MR Median" ~ "MR median",
                            TRUE ~ Method))

MD_data = MD_data %>%
  mutate(method = case_when(Method %in% c("MR cML", "MVMR cML") ~ "cML",
         Method %in% c("MR IVW", "MVMR IVW") ~ "IVW",
         Method %in% c("MR LASSO", "MVMR LASSO") ~ "LASSO",
         Method %in% c("MR median", "MVMR median") ~ "median",
         Method %in% c("MR-Egger", "MVMR-Egger") ~ "Egger",
         Method == "MVMR GMM"  ~ "GMM",
         Method == "Q_statistic_minimization"  ~ "Q-statistic minimization"),
         MR_type = ifelse(Method %in% c("MR IVW", "MR cML", "MR LASSO",
                                        "MR median", "MR-Egger"), "MR", "MVMR"))

# Plot data
data_DA <- MD_data %>%
  filter(Outcome == "DA",
         Exposure %in% c(sig_CIMs_MD),
         Exposure != "BMI",
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected"))

data_NDA <- MD_data %>%
  filter(Outcome == "NDA",
         Exposure %in% c(sig_CIMs_MD),
         Exposure != "BMI",
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected"))


data_PMD <- MD_data %>%
  filter(Outcome == "PMD",
         Exposure %in% c(sig_CIMs_MD),
         Exposure != "BMI",
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected"))

# Alphabetize exposures
data_DA$Exposure <- factor(
  data_DA$Exposure,
  levels = rev(sort(unique(data_DA$Exposure)))
)

data_NDA$Exposure <- factor(
  data_NDA$Exposure,
  levels = rev(sort(unique(data_NDA$Exposure)))
)

data_PMD$Exposure <- factor(
  data_PMD$Exposure,
  levels = rev(sort(unique(data_PMD$Exposure)))
)


# Identify whether the 95% CI crosses zero
data_DA <- data_DA %>%
  mutate(CI_crosses_zero = CILower <= 0 & CIUpper >= 0)

data_NDA <- data_NDA %>%
  mutate(CI_crosses_zero = CILower <= 0 & CIUpper >= 0)

data_PMD <- data_PMD %>%
  mutate(CI_crosses_zero = CILower <= 0 & CIUpper >= 0)


# Alphabetize exposures
data_DA$Exposure <- factor(
  data_DA$Exposure,
  levels = rev(sort(unique(data_DA$Exposure)))
)

data_NDA$Exposure <- factor(
  data_NDA$Exposure,
  levels = rev(sort(unique(data_NDA$Exposure)))
)

data_PMD$Exposure <- factor(
  data_PMD$Exposure,
  levels = rev(sort(unique(data_PMD$Exposure)))
)


# Identify whether the 95% CI crosses zero
add_shape_group <- function(data) {
  data %>%
    mutate(
      CI_crosses_zero = CILower <= 0 & CIUpper >= 0,
      shape_group = interaction(
        MR_type,
        CI_crosses_zero,
        sep = "_",
        drop = TRUE
      )
    )
}

data_DA  <- add_shape_group(data_DA)
data_NDA <- add_shape_group(data_NDA)
data_PMD <- add_shape_group(data_PMD)


# Get all unique methods across datasets
all_methods <- unique(c(
  as.character(data_DA$method),
  as.character(data_NDA$method),
  as.character(data_PMD$method)
))

# Assign colors
default_cols <- scales::brewer_pal(palette = "Dark2")(
  max(length(all_methods), 3)
)
default_cols <- default_cols[seq_along(all_methods)]
names(default_cols) <- all_methods

# Force IVW to black
default_cols["IVW"] <- "black"


# Define filled and open shapes for each MR type
all_mr_types <- unique(c(
  as.character(data_DA$MR_type),
  as.character(data_NDA$MR_type),
  as.character(data_PMD$MR_type)
))

filled_shapes <- c(16, 17, 15, 18, 8, 3, 4)
open_shapes   <- c(1, 2, 0, 5, 6, 7, 4)

shape_values <- unlist(lapply(seq_along(all_mr_types), function(i) {
  c(
    filled_shapes[i],
    open_shapes[i]
  )
}))


names(shape_values) <- unlist(lapply(all_mr_types, function(x) {
  c(
    paste0(x, "_FALSE"),
    paste0(x, "_TRUE")
  )
}))


pd <- position_dodge(width = 0.6)


# Plotting function
make_forest_plot <- function(data, outcome_label) {
  
  ggplot(
    data,
    aes(
      x = Estimate,
      y = Exposure,
      color = method,
      shape = shape_group
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
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
    scale_shape_manual(
      values = shape_values,
      breaks = unlist(lapply(all_mr_types, function(x) {
        c(paste0(x, "_FALSE"), paste0(x, "_TRUE"))
      })),
      labels = unlist(lapply(all_mr_types, function(x) {
        c(
          paste0(x, " (p<0.05)"),
          paste0(x)
        )
      })),
      name = "MR type (significance)"
    ) +
    labs(
      x = bquote(beta ~ " estimate for " ~ .(outcome_label)),
      y = "Exposure",
      color = "MR method"
    ) +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}


# Generate plots
plot_DA  <- make_forest_plot(data_DA, "DA")
plot_NDA <- make_forest_plot(data_NDA, "NDA")
plot_PMD <- make_forest_plot(data_PMD, "PMD")

## Combined plot
combined =
  plot_DA + plot_NDA + plot_PMD +
  plot_layout(ncol = 3, nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.tag = element_text(size = 12, face = "bold")
  )

ggsave(combined, filename = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Figures/Figure_3_new.jpg",
       width = 12, height = 16, dpi = 500)

####################### BCa plots #####################################

BCa_data = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/IMD_GWAS_CIM_BCa_MR_results.tsv")
BCa_OpenGWAS_data = read_tsv("/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/OpenGWAS_BCa_MR_results.tsv")


data = rbind_fill(MD_data,MD_OpenGWAS_data,BCa_data,BCa_OpenGWAS_data)

# Restrict to instruments with 3 or more SNPs

data = data %>% filter(NSNPs >2)

# Sensitivity analyses for significant CIMs #

# Filter to IVW results and relevant columns
data_IVW = data %>%
  filter(Method == "MR IVW") %>%
  select(c(Exposure,Outcome,Estimate,CILower,CIUpper,Fstatistic,NSNPs,Pvalue,Method)) 

# Prepare for plotting
sig_CIMs_BCa = data_IVW %>% 
  filter(Pvalue < 0.05, !Outcome %in% c("DA","NDA","PMD")) %>%
  pull(Exposure) %>%
  unique

data_sig_BCa = data %>%
  filter(Exposure %in% sig_CIMs_BCa, !Method %in% c("MR Wald Ratio",
                                                "MR-PRESSO distortion test",
                                                "MR-PRESSO global test",
                                                "MR-PRESSO outlier-corrected") ) %>%
  mutate(Exposure = sub("\\.", "-", Exposure)) %>%
  mutate(Method = case_when(Method == "MR Median" ~ "MR median",
                            TRUE ~ Method))

# Modify  plotting

data = data %>%
  mutate(Exposure = sub("\\.", "-", Exposure),
         Exposure = sub("VEGF-A", "VEGF", Exposure))

data_Overall <- data %>%
  filter(Outcome == "Overall_BCa",
         Exposure %in% data_sig_BCa$Exposure,
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected"))

data_HER2_pos <- data %>%
  filter(Outcome == "HER2_pos_BCa",
         Exposure %in% data_sig_BCa$Exposure,
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected"))


data_Luminal_A <- data %>%
  filter(Outcome == "Luminal_A_BCa",
         Exposure %in% data_sig_BCa$Exposure,
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected"))


data_Luminal_B <- data %>%
  filter(Outcome == "Luminal_B_BCa",
         Exposure %in% data_sig_BCa$Exposure,
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected"))


data_ER_neg <- data %>%
  filter(Outcome == "ER_neg_BCa",
         Exposure %in% data_sig_BCa$Exposure,
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected")) 


data_TN <- data %>%
  filter(Outcome == "TNBCa",
         Exposure %in% data_sig_BCa$Exposure,
         !Method %in% c("MR Wald Ratio",
                        "MR-PRESSO distortion test",
                        "MR-PRESSO global test",
                        "MR-PRESSO outlier-corrected"))

# Ensure exposure ordering is alphabetical

data_Overall$Exposure <- factor(data_Overall$Exposure,
                                levels = rev(sort(unique(data_Overall$Exposure))))

data_HER2_pos$Exposure <- factor(data_HER2_pos$Exposure,
                                 levels = rev(sort(unique(data_HER2_pos$Exposure))))

data_Luminal_A$Exposure  <- factor(data_Luminal_A$Exposure,
                                   levels = rev(sort(unique(data_Luminal_A$Exposure))))

data_Luminal_B$Exposure  <- factor(data_Luminal_B$Exposure,
                                   levels = rev(sort(unique(data_Luminal_B$Exposure))))

data_ER_neg$Exposure <- factor(data_ER_neg$Exposure,
                               levels = rev(sort(unique(data_ER_neg$Exposure))))

data_TN$Exposure <- factor(data_TN$Exposure,
                           levels = rev(sort(unique(data_TN$Exposure))))

# Get all unique methods across datasets
all_methods <- unique(c(data_Overall$Method,
                        data_HER2_pos$Method,
                        data_Luminal_A$Method,
                        data_Luminal_B$Method,
                        data_ER_neg$Method,
                        data_TN$Method))

# Create default brewer palette
default_cols <- scales::brewer_pal(palette = "Set1")(length(all_methods))
names(default_cols) <- all_methods

# Force MR IVW to black
default_cols["MR IVW"] <- "black"

# Run cross flag function 
add_cross_flag <- function(df) {
  df$crosses0 <- with(df, CILower <= 0 & CIUpper >= 0)
  df
}

data_Overall <- add_cross_flag(data_Overall)
data_HER2_pos <- add_cross_flag(data_HER2_pos)
data_Luminal_A  <- add_cross_flag(data_Luminal_A)
data_Luminal_B  <- add_cross_flag(data_Luminal_B)
data_ER_neg <- add_cross_flag(data_ER_neg)
data_TN <- add_cross_flag(data_TN)


# Plots

pd <- position_dodge(width = 0.6)
shape_vals <- c(`FALSE` = 16, `TRUE` = 1)

## Overall BCa
plot_Overall = ggplot(data_Overall,
                      aes(x = exp(Estimate),
                          y = Exposure,
                          color = Method,
                          shape = crosses0)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = exp(CILower), xmax = exp(CIUpper)),
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
    x = "OR estimate for overall BCa",
    y = "Exposure",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## ER_neg BCa
plot_ER_neg = ggplot(data_ER_neg,
                      aes(x = exp(Estimate),
                          y = Exposure,
                          color = Method,
                          shape = crosses0)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = exp(CILower), xmax = exp(CIUpper)),
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
    x = "OR estimate for ER- BCa",
    y = "Exposure",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )


## HER2 positive BCa
plot_HER2_pos = ggplot(data_HER2_pos,
                      aes(x = exp(Estimate),
                          y = Exposure,
                          color = Method,
                          shape = crosses0)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = exp(CILower), xmax = exp(CIUpper)),
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
    x = "OR estimate for HER2+ BCa",
    y = "Exposure",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## Luminal A
plot_Luminal_A = ggplot(data_Luminal_A,
                      aes(x = exp(Estimate),
                          y = Exposure,
                          color = Method,
                          shape = crosses0)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = exp(CILower), xmax = exp(CIUpper)),
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
    x = "OR estimate for LumA BCa",
    y = "Exposure",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## Luminal B
plot_Luminal_B = ggplot(data_Luminal_B,
                      aes(x = exp(Estimate),
                          y = Exposure,
                          color = Method,
                          shape = crosses0)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = exp(CILower), xmax = exp(CIUpper)),
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
    x = "OR estimate for LumB BCa",
    y = "Exposure",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## TN BCa
plot_TN = ggplot(data_TN,
                      aes(x = exp(Estimate),
                          y = Exposure,
                          color = Method,
                          shape = crosses0)) +
  geom_vline(xintercept = 1,
             linetype = "dashed",
             color = "grey50") +
  geom_errorbarh(
    aes(xmin = exp(CILower), xmax = exp(CIUpper)),
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
    x = "OR estimate for TN BCa",
    y = "Exposure",
    color = "MR method"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

## Combined plot
combined =
  plot_Overall  + plot_ER_neg + plot_HER2_pos +
  plot_Luminal_A + plot_Luminal_B + plot_TN +
  plot_layout(ncol = 2, nrow = 3, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    plot.tag = element_text(size = 16, face = "bold")
  )

ggsave(combined, filename = "/lindstroem/MD_MR_Inflammatory_Biomarkers/MR_results/Figures/Figure_2_new.jpg",
       width = 12, height = 16, dpi = 500)

