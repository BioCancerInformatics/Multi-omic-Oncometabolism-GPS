#### RSCRIPT TO BUILD SIGNATURES ACCORDING TO METABOLISM AND PATHWAYS
#### HIGOR ALMEIDA

library(dplyr)
library(tidyr)
library(rio)

################ gene, lncRNA, miRNA, and protein signatures ##################

# 1. Construction of coding gene signatures 
# ----------------------------------------------------------------
setwd("E:/Higor/Cancer_metabolism_project/11- Checking_results_and_annotation/3- Checking_results_and_annotation/")

df001 <- import("df063_gene_protein_mirna_lncrna_results.rds")

# Filter rows where Correlation_rho is not between -0.1 and 0.1 
df001 <- df001 %>%
  filter(Correlation_rho <= -0.1 | Correlation_rho >= 0.1)

# 2. Construction of gene signatures based on results and metabolic pathways
# ----------------------------------------------------------------

# Function to format the Target signature according to requirements
format_gene_signature <- function(Target) {
  formatted_Target <- sapply(Target, function(gene) {
    if (grepl("-", gene)) {
      paste0("", gene, "")
    } else {
      gene
    }
  })
  signature <- paste(formatted_Target, collapse = " + ")
  paste0("(", signature, ")")
}

# Function to extract common interactions between group elements
find_common_interactions <- function(Interactions) {
  interaction_list <- strsplit(Interactions, " / ")  # Split components
  
  # If there is only one row in the group, return the original value directly
  if (length(interaction_list) == 1) {
    return(Interactions[1])
  }
  
  # Check if all interactions are identical
  if (all(sapply(interaction_list, function(x) identical(x, interaction_list[[1]])))) {
    return(Interactions[1])  # Keep the original set if they are identical
  }
  
  # Identify intersection among elements
  common <- Reduce(intersect, interaction_list)
  
  # If there is an intersection, return the common values; otherwise, return "No common interactions"
  if (length(common) > 0) {
    return(paste(common, collapse = " / "))  # Return only the common elements
  } else {
    return("No common interactions")  # Return explicit message instead of an empty string
  }
}

# Update of the main aggregating function
aggregate_Target <- function(df001) {
  # Internal function to group data
  aggregate_and_format <- function(df) {
    df %>%
      group_by(Molecular_class, Metabolism, Pathways, Metabolic_cell_death, Cancer_types, Tumor_vs_normal, Omic_layer, Phenotypic_layer, Cox_OS_type, Cox_DSS_type, Cox_DFI_type, Cox_PFI_type,
               OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group, 
               Immune_classification) %>%
      summarise(
        Target = format_gene_signature(Target),
        Interactions = find_common_interactions(Interactions),  # Correct processing of Interactions
        .groups = 'drop'
      ) 
  }
  
  # Separate data with positive and negative correlations
  positive <- df001 %>% filter(Correlation_rho > 0)
  negative <- df001 %>% filter(Correlation_rho < 0)
  
  # Aggregate and format for positive and negative correlations
  positive_aggregated <- aggregate_and_format(positive)
  negative_aggregated <- aggregate_and_format(negative)
  
  # Add a column to identify the correlation type
  positive_aggregated <- positive_aggregated %>%
    mutate(Correlation_type = "positive")
  negative_aggregated <- negative_aggregated %>%
    mutate(Correlation_type = "negative")
  
  # Combine both aggregated dataframes
  aggregated_data <- bind_rows(positive_aggregated, negative_aggregated)
  
  return(aggregated_data)
}

# Apply the function to the dataframe
df002_signatures <- aggregate_Target(df001)

# Add a column that counts the number of Targets in each signature
df002_signatures <- df002_signatures %>%
  mutate(Members = sapply(strsplit(Target, " \\+ "), length))

# Sort by 'Cancer' and 'Members' in descending order within each cancer type
df002_signatures <- df002_signatures %>%
  arrange(Cancer_types, desc(Members))

# 3. Reorder and rename columns
# ----------------------------------------------------------------
New_column_order_01 <- c("Members", "Target", "Molecular_class", "Interactions", "Metabolism", "Pathways", "Metabolic_cell_death", "Cancer_types", "Tumor_vs_normal", 
                         "Omic_layer", "Phenotypic_layer", "Correlation_type", "Cox_OS_type", "Cox_DSS_type", "Cox_DFI_type", "Cox_PFI_type",
                         "OS_worst_prognosis_group", "DSS_worst_prognosis_group", "DFI_worst_prognosis_group", "PFI_worst_prognosis_group", 
                         "Immune_classification")

df002_signatures <- df002_signatures[New_column_order_01]

# Rename columns
df002_signatures <- df002_signatures %>%
  rename(
    Multiomics_Signature = Target,
    Common_interaction = Interactions
  )

# 10. Export constructed signatures
# ----------------------------------------------------------------
setwd("E:/Higor/Cancer_metabolism_project/13- Building_metabolic_signatures/")

export(df002_signatures, "Gene_lncRNA_miRNA_Protein_signatures.tsv")

saveRDS(df002_signatures, "Gene_lncRNA_miRNA_Protein_signatures.rds")


##################### transcript signatures #######################

# 1. Construction of coding gene signatures 
# ----------------------------------------------------------------
setwd("E:/Higor/Cancer_metabolism_project/11- Checking_results_and_annotation/3- Checking_results_and_annotation/")

df003_transcript <- import("df065_transcript_results.rds")

# Filter rows where Correlation_rho is not between -0.1 and 0.1 
df003_transcript <- df003_transcript %>%
  filter(Correlation_rho <= -0.1 | Correlation_rho >= 0.1)

# 2. Construction of gene signatures based on results and metabolic pathways
# ----------------------------------------------------------------

# Function to format the Target signature according to requirements
format_gene_signature <- function(Target) {
  formatted_Target <- sapply(Target, function(gene) {
    if (grepl("-", gene)) {
      paste0("", gene, "")
    } else {
      gene
    }
  })
  signature <- paste(formatted_Target, collapse = " + ")
  paste0("(", signature, ")")
}

# Function to extract common interactions between group elements
find_common_interactions <- function(Interactions) {
  interaction_list <- strsplit(Interactions, " / ")  # Split components
  
  # If there is only one row in the group, return the original value directly
  if (length(interaction_list) == 1) {
    return(Interactions[1])
  }
  
  # Check if all interactions are identical
  if (all(sapply(interaction_list, function(x) identical(x, interaction_list[[1]])))) {
    return(Interactions[1])  # Keep the original set if they are identical
  }
  
  # Identify intersection among elements
  common <- Reduce(intersect, interaction_list)
  
  # If there is an intersection, return the common values; otherwise, return "No common interactions"
  if (length(common) > 0) {
    return(paste(common, collapse = " / "))  # Return only the common elements
  } else {
    return("No common interactions")  # Return explicit message instead of an empty string
  }
}

# Update of the main aggregating function
aggregate_Target <- function(df003_transcript) {
  # Internal function to group data
  aggregate_and_format <- function(df) {
    df %>%
      group_by(Molecular_class, Metabolism, Pathways, Metabolic_cell_death, Cancer_types, Tumor_vs_normal, Omic_layer, Phenotypic_layer, Cox_OS_type, Cox_DSS_type, Cox_DFI_type, Cox_PFI_type,
               OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group, Immune_classification) %>%
      summarise(
        Target = format_gene_signature(Target),
        TranscriptGene = format_gene_signature(TranscriptGene),
        Interactions = find_common_interactions(Interactions),  # Correct processing of Interactions
        .groups = 'drop'
      ) 
  }
  
  # Separate data with positive and negative correlations
  positive <- df003_transcript %>% filter(Correlation_rho > 0)
  negative <- df003_transcript %>% filter(Correlation_rho < 0)
  
  # Aggregate and format for positive and negative correlations
  positive_aggregated <- aggregate_and_format(positive)
  negative_aggregated <- aggregate_and_format(negative)
  
  # Add a column to identify the correlation type
  positive_aggregated <- positive_aggregated %>%
    mutate(Correlation_type = "positive")
  negative_aggregated <- negative_aggregated %>%
    mutate(Correlation_type = "negative")
  
  # Combine both aggregated dataframes
  aggregated_data <- bind_rows(positive_aggregated, negative_aggregated)
  
  return(aggregated_data)
}

# Apply the function to the dataframe
df004_transcript_signatures <- aggregate_Target(df003_transcript)

# Add a column that counts the number of Targets in each signature
df004_transcript_signatures <- df004_transcript_signatures %>%
  mutate(Members = sapply(strsplit(Target, " \\+ "), length))

# Sort by 'Cancer' and 'Members' in descending order within each cancer type
df004_transcript_signatures <- df004_transcript_signatures %>%
  arrange(Cancer_types, desc(Members))

# 3. Reorder and rename columns
# ----------------------------------------------------------------
New_column_order_02 <- c("Members", "Target", "TranscriptGene", "Molecular_class", "Interactions", "Metabolism", "Pathways", "Metabolic_cell_death", "Cancer_types", "Tumor_vs_normal", 
                         "Omic_layer", "Phenotypic_layer", "Correlation_type", "Cox_OS_type", "Cox_DSS_type", "Cox_DFI_type", "Cox_PFI_type",
                         "OS_worst_prognosis_group", "DSS_worst_prognosis_group", "DFI_worst_prognosis_group", "PFI_worst_prognosis_group", 
                         "Immune_classification")

df004_transcript_signatures <- df004_transcript_signatures[New_column_order_02]

# Rename columns
df004_transcript_signatures <- df004_transcript_signatures %>%
  rename(
    Multiomics_Signature = Target,
    Common_interaction = Interactions
  )

# 10. Export constructed signatures
# ----------------------------------------------------------------
setwd("E:/Higor/Cancer_metabolism_project/13- Building_metabolic_signatures/")

export(df004_transcript_signatures, "Transcript_signatures.tsv")

saveRDS(df004_transcript_signatures, "Transcript_signatures.rds")


################################################################################

# 1. Criar pasta principal para armazenar assinaturas ômicas
# ----------------------------------------------------------------
setwd("E:/Higor/Cancer_metabolism_project/13- Building_metabolic_signatures/")
dir.create("Signatures_by_omic_layer")
setwd("E:/Higor/Cancer_metabolism_project/13- Building_metabolic_signatures/Signatures_by_omic_layer/")

# 3. Filtrar assinaturas específicas para cada camada ômica e salvar
# ----------------------------------------------------------------
Gene_expression_signature <- df002_signatures %>% filter(Omic_layer == "Gene expression")
saveRDS(Gene_expression_signature, "Gene_expression_signature.rds")

Methylation_signature <- df002_signatures %>% filter(Omic_layer == "Methylation")
saveRDS(Methylation_signature, "Methylation_signature.rds")

miRNA_signature <- df002_signatures %>% filter(Omic_layer == "miRNA expression")
saveRDS(miRNA_signature, "miRNA_signature.rds")

Protein_signature <- df002_signatures %>% filter(Omic_layer == "Protein expression")
saveRDS(Protein_signature, "Protein_signature.rds")

Mutation_signature <- df002_signatures %>% filter(Omic_layer == "Mutation")
saveRDS(Mutation_signature, "Mutation_signature.rds")

CNV_signature <- df002_signatures %>% filter(Omic_layer == "CNV")
saveRDS(CNV_signature, "CNV_signature.rds")

Transcript_signature <- df004_transcript_signatures %>% filter(Omic_layer == "Transcript expression")
saveRDS(Transcript_signature, "Transcript_signature.rds")

