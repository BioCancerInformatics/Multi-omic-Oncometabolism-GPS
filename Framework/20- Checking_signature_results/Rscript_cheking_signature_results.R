##### Checking results and annotation

library(dplyr)
library(tidyr)
library(rio)

setwd("E:/Oncometabolism_GPS/19- Tumor_immune_and_microenvironment_signature_classification/")

##### CODING GENE, PROTEIN, lncRNA, miRNA SIGNATURE RESULTS

# List of gene results to be imported
filepaths_gene <- c(
  "Gene/CNV/Immune_classification.tsv",
  "Gene/Methylation/Immune_classification.tsv",
  "Gene/mRNA/Immune_classification.tsv",
  "Gene/Mutation/Immune_classification.tsv",
  "miRNA/Immune_classification.tsv",
  "Protein/Immune_classification.tsv"
)

# Import and combine all files 
df001_results <- filepaths_gene %>%
  lapply(rio::import) %>% # Importa cada arquivo da lista
  bind_rows()             # Combina os arquivos em um único data.frame

# Remove specific columns
df002_results_filtered <- df001_results %>%
  select(-log10_correlation, -Expression_p.sgnif, -log10_p.adj, -immune_cells, -cor, -padj, -Infiltrate_profile)

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
df002_results_filtered <- df002_results_filtered %>%
  distinct()

# Replace empty strings with "no data"
df002_results_filtered <- df002_results_filtered %>%
  mutate(across(
    c(Expression, Expression_p.adj, p.value_Cox_OS, Type_Cox_OS, p.value_Cox_DSS, Type_Cox_DSS, p.value_Cox_DFI, Type_Cox_DFI, 
      p.value_Cox_PFI, Type_Cox_PFI, OS_log_rank_chisq, OS_p_val, OS_worst_prognosis_group, DSS_log_rank_chisq, DSS_p_val, DSS_worst_prognosis_group, 
      DFI_log_rank_chisq, DFI_p_val, DFI_worst_prognosis_group, PFI_log_rank_chisq, PFI_p_val, PFI_worst_prognosis_group, 
      Classification, Microenvironment_score_details, Immune_classification, Immune_score_details),
    ~ ifelse(is.na(.) | . == "" | . == "NA", "No data", as.character(.))
  ))


# Change column names
df003_renamed <- df002_results_filtered %>%
  rename(
    Multiomics_Signature = genes,
    Cancer_types = cancer_types,
    Tumor_vs_normal = Expression,
    Tumor_vs_normal_p.adj = Expression_p.adj,
    Omic_layer = Genotypic_var,
    Phenotypic_layer = Phenotypic_var,
    Correlation_rho = correlation,
    Correlation_p.adj = Correlation_p.adj,
    Cox_OS_type = Type_Cox_OS,
    Cox_OS_p.value = p.value_Cox_OS,
    Cox_DSS_type = Type_Cox_DSS,
    Cox_DSS_p.value = p.value_Cox_DSS,
    Cox_DFI_type = Type_Cox_DFI,
    Cox_DFI_p.value = p.value_Cox_DFI,
    Cox_PFI_type = Type_Cox_PFI,
    Cox_PFI_p.value = p.value_Cox_PFI,
    OS_log_rank_chisq = OS_log_rank_chisq,
    OS_p.value = OS_p_val,
    OS_worst_prognosis_group = OS_worst_prognosis_group,
    DSS_log_rank_chisq = DSS_log_rank_chisq,
    DSS_p.value = DSS_p_val,
    DSS_worst_prognosis_group = DSS_worst_prognosis_group,
    DFI_log_rank_chisq = DFI_log_rank_chisq,
    DFI_p.value = DFI_p_val,
    DFI_worst_prognosis_group = DFI_worst_prognosis_group,
    PFI_log_rank_chisq = PFI_log_rank_chisq,
    PFI_p.value = PFI_p_val,
    PFI_worst_prognosis_group = PFI_worst_prognosis_group,
    Microenvironment_classification = Classification,
    Microenvironment_score_details = Microenvironment_score_details,
    Immune_classification = Immune_classification,
    Immune_score_details = Immune_score_details
  )


# Change names
df003_renamed <- df003_renamed %>%
  mutate(Phenotypic_layer = case_when(
    Phenotypic_layer == "MSI" ~ "Microsatellite instability",
    Phenotypic_layer == "TMB" ~ "Tumor mutational burden",
    TRUE ~ Phenotypic_layer # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df003_renamed <- df003_renamed %>%
  mutate(Immune_classification = case_when(
    Immune_classification == "Frio" ~ "Cold",
    Immune_classification == "Quente" ~ "Hot",
    Immune_classification == "Variável" ~ "Variable",
    TRUE ~ Immune_classification # Mantém os valores originais caso não sejam alterados
  ))

# Import signature information
df004_signature_information <- import("/Oncometabolism_GPS/13- Building_metabolic_signatures/Gene_lncRNA_miRNA_Protein_signatures.rds")

# Filtering columns df_003
df005_signature_information_filtered <- df004_signature_information %>% 
  select(Members, Multiomics_Signature, Molecular_class, Common_interaction, Metabolism, Pathways, Metabolic_cell_death, 
         Cancer_types, Omic_layer, Phenotypic_layer)


# Função para processar assinaturas corretamente
process_gene_name <- function(gene) {
  # Remover crases se presentes
  gene <- gsub("`", "", gene)
  
  # Verificar se há múltiplos membros na assinatura
  if (grepl("\\+", gene)) {
    # Remover parênteses existentes antes de processar
    gene <- gsub("^\\(|\\)$", "", gene)
    
    # Dividir os elementos da assinatura
    members <- unlist(strsplit(gene, " \\+ "))
    
    # Adicionar crase nos membros que contêm "-"
    members <- sapply(members, function(x) {
      if (grepl("-", x)) {
        return(paste0("`", x, "`"))  # Adiciona crases
      } else {
        return(x)
      }
    })
    
    # Reconstruir a assinatura formatada com os membros processados
    formatted_gene <- paste(members, collapse = " + ")
    
    # Manter um único par de parênteses
    return(paste0("(", formatted_gene, ")"))
  } else {
    # Se for uma assinatura única, remover parênteses
    return(gsub("\\(|\\)", "", gene))
  }
}

# Substituir os valores na coluna Multiomics_Signature usando a função process_gene_name()
df005_signature_information_filtered$Multiomics_Signature_Original <- df005_signature_information_filtered$Multiomics_Signature  # Armazena original para comparação
df005_signature_information_filtered$Multiomics_Signature <- sapply(df005_signature_information_filtered$Multiomics_Signature, process_gene_name)


# Define the key columns for merging
key_cols <- c("Multiomics_Signature", "Cancer_types", "Omic_layer", "Phenotypic_layer")

# Perform the merge
df006_merged_information <- merge(df005_signature_information_filtered, df003_renamed, by = key_cols, all.x = TRUE)

# Filtering columns df_006
df007_merged_information_filtered <- df006_merged_information %>% 
  select(-Multiomics_Signature_Original, -OS_log_rank_chisq, -DSS_log_rank_chisq, -DFI_log_rank_chisq, -PFI_log_rank_chisq)

# Reorder the columns
New_column_order_01 <- c("Members", "Multiomics_Signature", "Molecular_class", "Common_interaction", "Metabolism", "Pathways", "Metabolic_cell_death", "Cancer_types", 
                         "Tumor_vs_normal",  "Tumor_vs_normal_p.adj", "Omic_layer", "Phenotypic_layer", "Correlation_rho", "Correlation_p.adj", "Cox_OS_type", "Cox_OS_p.value", 
                         "Cox_DSS_type", "Cox_DSS_p.value", "Cox_DFI_type", "Cox_DFI_p.value", "Cox_PFI_type", "Cox_PFI_p.value", "OS_worst_prognosis_group", "OS_p.value", 
                         "DSS_worst_prognosis_group", "DSS_p.value", "DFI_worst_prognosis_group", "DFI_p.value", "PFI_worst_prognosis_group", "PFI_p.value",
                         "Microenvironment_classification", "Microenvironment_score_details", "Immune_classification", "Immune_score_details")

df007_merged_information_filtered <- df007_merged_information_filtered[New_column_order_01]

df008_final <- df007_merged_information_filtered %>%
  mutate(across(c(Tumor_vs_normal, Tumor_vs_normal_p.adj), as.character))

df008_final[df008_final$Omic_layer %in% c("Methylation", "CNV", "Mutation", "Protein expression"), 
                                                 c("Tumor_vs_normal", "Tumor_vs_normal_p.adj")] <- "No data"
# # Verifying results 
# 
# # First, keep only the columns of interest
# df004_sub <- df004_signature_information[, c("Cox_OS_type", "Cox_DSS_type", "Cox_DFI_type", "Cox_PFI_type",
#                                              "OS_worst_prognosis_group", "DSS_worst_prognosis_group", "DFI_worst_prognosis_group", "PFI_worst_prognosis_group")]
# df008_sub <- df008_final[, c("Cox_OS_type", "Cox_DSS_type", "Cox_DFI_type", "Cox_PFI_type",
#                              "OS_worst_prognosis_group", "DSS_worst_prognosis_group", "DFI_worst_prognosis_group", "PFI_worst_prognosis_group")]
# 
# # Make sure the data frames are the same size and in the same order (optional: depends on your case)
# df004_sub <- df004_sub %>% arrange(Cox_OS_type, Cox_DSS_type, Cox_DFI_type, Cox_PFI_type,
#                                    OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group)
# 
# df008_sub <- df008_sub %>% arrange(Cox_OS_type, Cox_DSS_type, Cox_DFI_type, Cox_PFI_type,
#                                    OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group)
# 
# # Compare row by row
# comparison <- df004_sub == df008_sub
# 
# # Summarize the differences
# differences_per_column <- colSums(!comparison, na.rm = TRUE)
# total_differences <- sum(!comparison, na.rm = TRUE)
# 
# # Output
# print("Differences per column:")
# print(differences_per_column)
# 
# print(paste("Total differences:", total_differences))


setwd("/Oncometabolism_GPS/20- Checking_signature_results/")

export(df008_final, "df_gene_protein_mirna_lncrna_signatures.tsv")

saveRDS(df008_final, "df_gene_protein_mirna_lncrna_signatures.rds")
             
##### CODING GENE TRANSCRIPT AND lncRNA TRANSCRIPT SIGNATURE RESULTS

setwd("E:/Oncometabolism_GPS/19- Tumor_immune_and_microenvironment_signature_classification/")

# Lista de arquivos a serem importados
filepaths_transcript_gene <- c(
  "Transcript/Immune_classification.tsv"
)

# Import and combine all files 
df009_results <- filepaths_transcript_gene %>%
  lapply(rio::import) %>% # Importa cada arquivo da lista
  bind_rows()             # Combina os arquivos em um único data.frame


# Remove specific columns
df010_results_filtered <- df009_results %>%
  select(-log10_correlation, -Expression_p.sgnif, -log10_p.adj, -immune_cells, -cor, -padj, -Infiltrate_profile)

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
df010_results_filtered <- df010_results_filtered %>%
  distinct()

# Replace empty strings with "no data"
df010_results_filtered <- df010_results_filtered %>%
  mutate(across(
    c(Expression, Expression_p.adj, p.value_Cox_OS, Type_Cox_OS, p.value_Cox_DSS, Type_Cox_DSS, p.value_Cox_DFI, Type_Cox_DFI, 
      p.value_Cox_PFI, Type_Cox_PFI, OS_log_rank_chisq, OS_p_val, OS_worst_prognosis_group, DSS_log_rank_chisq, DSS_p_val, DSS_worst_prognosis_group, 
      DFI_log_rank_chisq, DFI_p_val, DFI_worst_prognosis_group, PFI_log_rank_chisq, PFI_p_val, PFI_worst_prognosis_group, 
      Classification, Microenvironment_score_details, Immune_classification, Immune_score_details),
    ~ ifelse(is.na(.) | . == "" | . == "NA", "No data", as.character(.))
  ))

# Change column names
df011_renamed <- df010_results_filtered %>%
  rename(
    Multiomics_Signature = genes,
    Cancer_types = cancer_types,
    Tumor_vs_normal = Expression,
    Tumor_vs_normal_p.adj = Expression_p.adj,
    Omic_layer = Genotypic_var,
    Phenotypic_layer = Phenotypic_var,
    Correlation_rho = correlation,
    Correlation_p.adj = Correlation_p.adj,
    Cox_OS_type = Type_Cox_OS,
    Cox_OS_p.value = p.value_Cox_OS,
    Cox_DSS_type = Type_Cox_DSS,
    Cox_DSS_p.value = p.value_Cox_DSS,
    Cox_DFI_type = Type_Cox_DFI,
    Cox_DFI_p.value = p.value_Cox_DFI,
    Cox_PFI_type = Type_Cox_PFI,
    Cox_PFI_p.value = p.value_Cox_PFI,
    OS_log_rank_chisq = OS_log_rank_chisq,
    OS_p.value = OS_p_val,
    OS_worst_prognosis_group = OS_worst_prognosis_group,
    DSS_log_rank_chisq = DSS_log_rank_chisq,
    DSS_p.value = DSS_p_val,
    DSS_worst_prognosis_group = DSS_worst_prognosis_group,
    DFI_log_rank_chisq = DFI_log_rank_chisq,
    DFI_p.value = DFI_p_val,
    DFI_worst_prognosis_group = DFI_worst_prognosis_group,
    PFI_log_rank_chisq = PFI_log_rank_chisq,
    PFI_p.value = PFI_p_val,
    PFI_worst_prognosis_group = PFI_worst_prognosis_group,
    Microenvironment_classification = Classification,
    Microenvironment_score_details = Microenvironment_score_details,
    Immune_classification = Immune_classification,
    Immune_score_details = Immune_score_details
  )


# Change names
df011_renamed <- df011_renamed %>%
  mutate(Phenotypic_layer = case_when(
    Phenotypic_layer == "MSI" ~ "Microsatellite instability",
    Phenotypic_layer == "TMB" ~ "Tumor mutational burden",
    TRUE ~ Phenotypic_layer # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df011_renamed <- df011_renamed %>%
  mutate(Immune_classification = case_when(
    Immune_classification == "Frio" ~ "Cold",
    Immune_classification == "Quente" ~ "Hot",
    Immune_classification == "Variável" ~ "Variable",
    TRUE ~ Immune_classification # Mantém os valores originais caso não sejam alterados
  ))

# Import signature information
df012_signature_information <- import("/Oncometabolism_GPS/13- Building_metabolic_signatures/Transcript_signatures.rds")

# Filtering columns df_003
df013_signature_information_filtered <- df012_signature_information %>% 
  select(Members, Multiomics_Signature, TranscriptGene, Molecular_class, Common_interaction, Metabolism, Pathways, Metabolic_cell_death, 
         Cancer_types, Omic_layer, Phenotypic_layer)


# Função para processar assinaturas corretamente
process_gene_name <- function(gene) {
  # Remover crases se presentes
  gene <- gsub("`", "", gene)
  
  # Verificar se há múltiplos membros na assinatura
  if (grepl("\\+", gene)) {
    # Remover parênteses existentes antes de processar
    gene <- gsub("^\\(|\\)$", "", gene)
    
    # Dividir os elementos da assinatura
    members <- unlist(strsplit(gene, " \\+ "))
    
    # Adicionar crase nos membros que contêm "-"
    members <- sapply(members, function(x) {
      if (grepl("-", x)) {
        return(paste0("`", x, "`"))  # Adiciona crases
      } else {
        return(x)
      }
    })
    
    # Reconstruir a assinatura formatada com os membros processados
    formatted_gene <- paste(members, collapse = " + ")
    
    # Manter um único par de parênteses
    return(paste0("(", formatted_gene, ")"))
  } else {
    # Se for uma assinatura única, remover parênteses
    return(gsub("\\(|\\)", "", gene))
  }
}

# Substituir os valores na coluna Multiomics_Signature usando a função process_gene_name()
df013_signature_information_filtered$Multiomics_Signature_Original <- df013_signature_information_filtered$Multiomics_Signature  # Armazena original para comparação
df013_signature_information_filtered$Multiomics_Signature <- sapply(df013_signature_information_filtered$Multiomics_Signature, process_gene_name)

# Define the key columns for merging
key_cols <- c("Multiomics_Signature", "Cancer_types", "Omic_layer", "Phenotypic_layer")

# Perform the merge
df014_merged_information <- merge(df013_signature_information_filtered, df011_renamed, by = key_cols, all.x = TRUE)

# Filtering columns df_006
df015_merged_information_filtered <- df014_merged_information %>% 
  select(-Multiomics_Signature_Original, -OS_log_rank_chisq, -DSS_log_rank_chisq, -DFI_log_rank_chisq, -PFI_log_rank_chisq)

# Reorder the columns
New_column_order_01 <- c("Members", "Multiomics_Signature", "TranscriptGene", "Molecular_class", "Common_interaction", "Metabolism", "Pathways", "Metabolic_cell_death", "Cancer_types", 
                         "Tumor_vs_normal",  "Tumor_vs_normal_p.adj", "Omic_layer", "Phenotypic_layer", "Correlation_rho", "Correlation_p.adj", "Cox_OS_type", "Cox_OS_p.value", 
                         "Cox_DSS_type", "Cox_DSS_p.value", "Cox_DFI_type", "Cox_DFI_p.value", "Cox_PFI_type", "Cox_PFI_p.value", "OS_worst_prognosis_group", "OS_p.value", 
                         "DSS_worst_prognosis_group", "DSS_p.value", "DFI_worst_prognosis_group", "DFI_p.value", "PFI_worst_prognosis_group", "PFI_p.value",
                         "Microenvironment_classification", "Microenvironment_score_details", "Immune_classification", "Immune_score_details")

df015_merged_information_filtered <- df015_merged_information_filtered[New_column_order_01]

df016_final <- df015_merged_information_filtered %>%
  mutate(across(c(Tumor_vs_normal, Tumor_vs_normal_p.adj), as.character))

setwd("/Oncometabolism_GPS/20- Checking_signature_results/")

export(df016_final, "df_transcript_signatures.tsv")

saveRDS(df016_final, "df_transcript_signatures.rds")

# Binding all signature results
df017_transcript_filtered <- df016_final %>% 
  select(-TranscriptGene)

df018_all_signatures <- rbind(df017_transcript_filtered, df008_final)

saveRDS(df018_all_signatures, "df_all_signatures.rds")
