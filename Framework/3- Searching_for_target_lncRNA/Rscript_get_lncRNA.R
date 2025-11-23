# RSCRIPT TO CHECK WHICH LNCRNA AVAILABLE IN THE EVlncRNAs 3.0 DATABASE INTERACT WITH METABOLIC GENES AND MATURE MIRNA
# AUTHOR: HIGOR ALMEIDA - PhD WORK
# LAST VERSION: 06/09/2024

# Loading library
library("dplyr")
library("rio")

# Setting directory 
setwd("F:/Higor/Cancer_metabolism_analysis_06/3- Searching for target lncRNA")

# Importing data
lncrna <- import("function_human.xls")
miRNA <- import("/Higor/Cancer_metabolism_analysis_06/2- Searching for target mature miRNA/Target_mirna.tsv")
gene <- import("/Higor/Cancer_metabolism_analysis_06/1- Searching for target coding genes/Target_genes.tsv")

# # Remove o prefixo "hsa-" de miRNA$Mature_mirna_id
# miRNA$Mature_mirna_id <- sub("^hsa-", "", miRNA$Mature_mirna_id)
# 
# # Identifica as linhas onde lncrna$`Interaction target` corresponde a Target_genes ou Mature_mirna_id
# matching_rows <- lncrna[lncrna$`Interaction target` %in% c(gene$Genes, miRNA$Mature_mirna_id), ]

# Criar uma coluna Matches com os itens correspondentes
lncrna$Matches <- ifelse(
  lncrna$`Interaction target` %in% c(gene$Genes, miRNA$Mature_mirna_id),
  lncrna$`Interaction target`, 
  NA
)

# Filtrar as linhas onde houve correspondência
matching_rows <- lncrna[!is.na(lncrna$Matches), ]


# Função para verificar correspondência e retornar os itens encontrados
find_matches <- function(row, genes, mirnas) {
  # Remover espaços extras ao redor de "|"
  row <- gsub("\\s*\\|\\s*", "|", row)
  
  # Encontrar genes correspondentes
  matching_genes <- genes[sapply(genes, function(g) grepl(paste0("\\b", g, "\\b"), row, ignore.case = TRUE))]
  
  # Encontrar miRNAs correspondentes
  matching_mirnas <- mirnas[sapply(mirnas, function(m) grepl(paste0("\\b", m, "\\b"), row, ignore.case = TRUE))]
  
  # Combinar resultados em uma string
  paste(c(matching_genes, matching_mirnas), collapse = ", ")
}

# Adicionar uma coluna indicando correspondências encontradas
lncrna$Matches <- sapply(lncrna$`Detailed Pathway`, find_matches, gene$Genes, miRNA$Mature_mirna_id)

# Filtrar as linhas onde há correspondência
matching_rows2 <- lncrna[lncrna$Matches != "", ]


# Combina os dois data frames
combined_rows <- rbind(matching_rows, matching_rows2)

# Remove as linhas duplicadas
unique_rows <- combined_rows[!duplicated(combined_rows), ]

# Exporting results
export(unique_rows, "lncRNA-related_metabolic_genes.tsv")

# Keeping only the column named "Lncrna name"
Target_lncrna <- matching_rows %>% select(`LncRNA name`)

# Renaming the column name
Target_lncrna <- Target_lncrna %>% rename(Target_lncrna = `LncRNA name`)

# Removing duplicated 
Target_lncrna <- Target_lncrna %>% distinct(Target_lncrna, .keep_all = TRUE)

# Exporting results
export(Target_lncrna, "Target_lncrna.tsv")

# -------------------------

# RSCRIPT TO CHECK WHICH LNCRNA AVAILABLE IN THE EVlncRNAs 3.0 DATABASE INTERACT WITH METABOLIC GENES AND MATURE MIRNA
# AUTHOR: HIGOR ALMEIDA - PhD WORK
# LAST VERSION: 06/09/2024
  
# Loading library
library("dplyr")
library("rio")

# Setting directory 
setwd("F:/Higor/Cancer_metabolism_analysis_06/3- Searching for target lncRNA")

# Importing data
lncrna <- import("function_human.xls")
# miRNA <- import("/Higor/Cancer_metabolism_analysis_06/2- Searching for target mature miRNA/Target_mirna.tsv")
gene <- import("/Higor/Cancer_metabolism_analysis_06/1- Searching for target coding genes/Target_genes.tsv")

# # Remove o prefixo "hsa-" de miRNA$Mature_mirna_id
# miRNA$Mature_mirna_id <- sub("^hsa-", "", miRNA$Mature_mirna_id)

# Identificar as linhas onde lncrna$`Interaction target` corresponde a gene$Genes
matching_rows_genes <- lncrna[lncrna$`Interaction target` %in% gene$Genes, ]

# Exporting results
export(matching_rows_genes, "lncRNA-related_genes.tsv")

# # Identificar as linhas onde lncrna$`Interaction target` corresponde a miRNA$Mature_mirna_id
# matching_rows_mirna <- lncrna[lncrna$`Interaction target` %in% miRNA$Mature_mirna_id, ]
# 
# # Exporting results
# export(matching_rows_mirna, "lncRNA-related_mirna.tsv")
# 
# # Combina os dois data frames
# combined_matching <- rbind(matching_rows_genes, matching_rows_mirna)
# 
# # Exporting results
# export(combined_matching, "lncRNA-related_genes_mirna.tsv")

# Keeping only the column named "Lncrna name"
Target_lncrna <- matching_rows_genes %>% select(`LncRNA name`)

# Renaming the column name
Target_lncrna <- Target_lncrna %>% rename(Target_lncrna = `LncRNA name`)

# Removing duplicated 
Target_lncrna <- Target_lncrna %>% distinct(Target_lncrna, .keep_all = TRUE)

# Exporting results
export(Target_lncrna, "New_Target_lncrna.tsv")
