##### RCD DATABASE - Searching for Metabolic cell death target

# Load required libraries
library(rio)      
library(dplyr)    

# Set the working directory to the folder where the data files are located
setwd("C:/Users/quiqu/OneDrive/Área de Trabalho/RCD_database")

# Import cell death pathway datasets from CSV files
df1_Alkaliptosis <- import("Alkaliptosis.csv")
df2_Apoptosis <- import("Apoptosis.csv")
df3_Autophagy_dependent_cell_death <- import("Autophagy_dependent_cell_death.csv")
df4_Cuproptosis <- import("Cuproptosis.csv")
df5_Disulfidptosis <- import("Disulfidptosis.csv")
df6_Entotic_cell_death <- import("Entotic_cell_death.csv")
df7_Ferroptosis <- import("Ferroptosis.csv")
df8_Immunogenic_cell_death <- import("Immunogenic_cell_death.csv")
df9_Lysosome_dependent_cell_death <- import("Lysosome_dependent_cell_death.csv")
df10_MPT_driven_necrosis <- import("MPT_driven_necrosis.csv")
df11_Necroptosis <- import("Necroptosis.csv")
df12_NETotic_cell_death <- import("NETotic_cell_death.csv")
df13_Oxeiptosis <- import("Oxeiptosis.csv")
df14_Parthanatos <- import("Parthanatos.csv")
df15_Pyroptosis <- import("Pyroptosis.csv")

# Combine all individual cell death pathway datasets into a single dataframe
df16_all_RCD <- rbind(df1_Alkaliptosis, df2_Apoptosis, df3_Autophagy_dependent_cell_death, df4_Cuproptosis,
                      df5_Disulfidptosis, df6_Entotic_cell_death, df7_Ferroptosis, df8_Immunogenic_cell_death,
                      df9_Lysosome_dependent_cell_death, df10_MPT_driven_necrosis, df11_Necroptosis, df12_NETotic_cell_death,
                      df13_Oxeiptosis, df14_Parthanatos, df15_Pyroptosis)

# Import results dataset containing the target genes for filtering
df17_metabolic_target <- import("Results.rds")

# Transform Target column
df17_metabolic_target$Target <- df17_metabolic_target$Target %>%
  sub("^hsa-", "", .) %>%       # Remove "hsa-" if present
  sub("miR", "MIR", .) %>%      # Replace "miR" with "MIR"
  gsub("-", "", .) %>%          # Remove all hyphens
  sub("(5p|3p)$", "", .)      # Remove only "-5p" or "-3p" at the end

# Filter df1_all_RCD to keep only rows where the 'gene' column matches the 'Target' column in df2_results
df18_filtered_RCD <- df16_all_RCD %>% 
  filter(gene %in% df17_metabolic_target$Target)

# Export the filtered dataset to a TSV (tab-separated values) file for further analysis
export(df18_filtered_RCD, "Metabolic_cell_death_target.tsv")

