# Setup
library(data.table)
library(R.utils)
library(tidyverse)
setDTthreads(percent = 100)


# Function to extract vectors of RSIDs from clumping files

extract_rsids <- function(filepaths) {
  
  rs_ids <- character()
  
  for (fp in filepaths) {
    
    lines <- readLines(fp, warn = FALSE)
    
    for (line in lines) {
      
      fields <- strsplit(line, "[\t ]+")[[1]]
      
      idx <- grep("^rs", fields)
      if (length(idx) > 0) {
        rs_ids <- c(rs_ids, fields[min(idx)])
      }
    }
  }
  
  unique(rs_ids)
}


# Generate list of filepaths and corresponding proteins

dir_path <- "/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/BCa/harmonized_p_filtered/"

filepaths <- list.files(
  path = dir_path,
  pattern = "\\.tsv$",
  full.names = TRUE
)

proteins <- sub("_.*$", "", basename(filepaths))

# Load harmonized dataset, filter to MAF-filtered and LD pruned SNPs in the final instruments, and save

for (i in 1:length(proteins)) {
  ## Load harmonized datasets
  
  harmonized_data <- fread(filepaths[i])
  
  ## Apply LD pruning and MAF filters
  
  LD_prune_path <- "/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/BCa/harmonized_fully_filtered/post_clumping_SNPs/"
  
  LD_filepaths <- list.files(
    path = LD_prune_path,
    pattern = paste0("\\Q", proteins[i], "\\E.*\\.clumps$"),
    full.names = TRUE
  )
  
  clump_rsids = extract_rsids(LD_filepaths)
  
  filtered_data = harmonized_data %>%
    filter(RSID %in% clump_rsids)
  
  ## Save
  fwrite(x = filtered_data, file = paste0("/lindstroem/MD_MR_Inflammatory_Biomarkers/IMD_GWAS_CIMs/BCa/harmonized_fully_filtered/",proteins[i],"_filtered_harmonized_data.tsv"),
         sep = "\t",
         quote = FALSE,
         row.names = FALSE,
  )
  
}

