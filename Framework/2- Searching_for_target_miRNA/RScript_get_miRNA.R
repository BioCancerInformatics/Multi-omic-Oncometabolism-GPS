#### RSCRIPT TO SEARCH AND GET miRNA RELATED TO THE METABOLIC GENE (mirecords, mirtarbase, and tarbase)
#### HIGOR ALMEIDA, M.sC
#### LAST VERSION - 13/08/2024

# Carregamento dos pacotes necessários
library(BiocManager)
library(multiMiR)
library(rio)
library(dplyr)

setwd("D:/Cancer_metabolism_analysis_06/2- Searching_for_target_miRNA")

# # Importando conjunto de genes
# data <- import("/Cancer_metabolism_analysis_06/1- Searching_for_target_coding_genes/Target_genes.tsv")
# 
# # Extraindo a coluna de genes como um vetor de caracteres
# genes <- as.character(data$Genes)
# 
# # Utilizando o multiMiR para prever os miRNAs alvos dos mRNAs
# mirna_targets <- get_multimir(
#   org = 'hsa', 
#   target = genes,
#   table = "validated",
#   summary = FALSE
# )
# 
# # Extraindo os dados
# mirna_targets_data <- mirna_targets@data
# 
# # Salvando os resultados em um arquivo tsv
# export(mirna_targets_data, "mirna_targets_data.tsv")
# 
# # Filtrando as colunas relevantes para mostrar a relação entre miRNA e mRNA
# Target_mirna <- mirna_targets_data[, c("database", "mature_mirna_id", "target_symbol", "pubmed_id")]
# 
# # Renomear as colunas 
# Target_mirna <- Target_mirna %>%
#   rename(
#     Database = database,
#     Mature_mirna_id = mature_mirna_id,
#     Target_genes = target_symbol,
#     Pubmed_id = pubmed_id,
#   )
# 
# # Ornanizando os Target genes em ordem alfabética
# All_target_mirna <- Target_mirna[order(Target_mirna$Target_genes), ]
# 
# # Salvando os resultados em um arquivo tsv
# export(Target_mirna, "All_target_mirna.tsv")
# 
# # Remove rows where Pubmed_id is NA or an empty string
# Target_mirna <- Target_mirna[!(is.na(Target_mirna$Pubmed_id) | Target_mirna$Pubmed_id == ""), ]
# 
# # To remove rows that are duplicated across three columns
# Target_mirna <- Target_mirna %>% 
#   distinct(Database, Mature_mirna_id, Target_genes, Pubmed_id, .keep_all = TRUE)
# 
# # To combine the different Database values into a single row for each unique combination of Mature_mirna_id and Target_gene, while keeping the Pubmed_id column
# Target_mirna <- Target_mirna %>%
#   group_by(Mature_mirna_id, Target_genes) %>%
#   summarise(
#     Database = paste(unique(Database), collapse = " / "),
#     Pubmed_id = paste(unique(Pubmed_id), collapse = " / "),
#     .groups = 'drop'
#   )
# 
# # Salvando os resultados em um arquivo tsv
# export(Target_mirna, "Target_mirna_pubmed_id.tsv")
# 
# # Removendo miRNA duplicados 
# Target_mirna <- Target_mirna %>% distinct(Mature_mirna_id, .keep_all = TRUE)
# 
# # Keep only the column named "Mature_mirna_id"
# Target_mirna <- Target_mirna %>% select(Mature_mirna_id)
# 
# # Salvando os resultados em um arquivo tsv
# export(Target_mirna, "Target_mirna.tsv")

######### New Approach

# Importando conjunto de genes
df001_genes <- import("/Cancer_metabolism_analysis_06/1- Searching_for_target_coding_genes/Target_genes.tsv")

# Importando conjunto de miRNAs
df002_mirna <- import("/Cancer_metabolism_analysis_06/2- Searching_for_target_miRNA/hsa_MTI.csv")

# Identificar as linhas onde miRNA corresponde a gene
df003_matching <- df002_mirna[df002_mirna$`Target Gene` %in% df001_genes$Genes, ]

# ✅ Remove rows where 'PMID' is NA
df004_clean <- df003_matching %>%
  # Remove rows with NA in the PMID column
  filter(!is.na(`References (PMID)`)) %>%
  # Keep only rows where specie is "hsa"
  filter(`Species (miRNA)` == "hsa")

# Exporting results
export(df004_clean, "miRNA-related_genes.tsv")

# Removing duplicated 
df005_target_lncrna <- df004_clean %>% distinct(miRNA, .keep_all = TRUE)

# Exporting results
export(df005_target_lncrna, "New_Target_mirna.tsv")
