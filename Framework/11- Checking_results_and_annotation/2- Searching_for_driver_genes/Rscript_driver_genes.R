### Filtering common genes in two variables from two different files which tier 1 COSMIC Cancer driver genes are in our results

# Set the working directory to the folder containing the relevant files
setwd("F:/Higor/Cancer_metabolism_analysis_06/11- Checking_results_and_annotation/2- Searching_for_driver_genes")

# Load required libraries
library(rio)      
library(dplyr)    

# Import the COSMIC Cancer driver genes dataset
df01_driver_genes_COSMIC <- import("Census_allTue Jun 11 16_27_10 2024.tsv")

# Import the results dataset containing potential metabolic cell death genes
df02_results <- import("/Higor/Cancer_metabolism_analysis_06/11- Checking_results_and_annotation/1- Searching_for_metabolic_cell_death_gene/Results.rds")

# Transform the 'Target' column in df02_results to standardize gene names
df02_results$Target <- df02_results$Target %>%
  sub("^hsa-", "", .) %>%       # Remove "hsa-" prefix if present
  sub("miR", "MIR", .) %>%      # Replace "miR" with "MIR" for consistency
  gsub("-", "", .) %>%          # Remove all hyphens from gene names
  sub("(5p|3p)$", "", .)        # Remove "-5p" or "-3p" suffix if present

# Filter the COSMIC dataset to keep only genes that are present in the transformed 'Target' column of df02_results
df03_filtered_driver_genes <- df01_driver_genes_COSMIC %>% 
  filter(`Gene Symbol` %in% df02_results$Target)

# Export the filtered dataset to a TSV (tab-separated values) file for further analysis
export(df03_filtered_driver_genes, "Driver_gene_target.tsv")
