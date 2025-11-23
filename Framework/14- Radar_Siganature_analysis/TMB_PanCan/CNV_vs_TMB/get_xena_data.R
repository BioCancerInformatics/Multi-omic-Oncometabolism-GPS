### UCSCXenaShiny Pancan genotype-phenotype correlations
### Authors: Juan Carlo S Silva, Emanuell Rodrigues de Souza, Leonardo Henrique da Silva, Higor A C Nogueira, Ana Beatriz Garcia, Enrique Medina-Acosta
### Updated: 01/08/2024

### GET_XENA_DATA - CNV vs TMB

# Load required libraries
library(UCSCXenaShiny)
library(UCSCXenaTools)
library(readr)
library(lubridate)
library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(fmsb)
library(rio)

# Print start time
start_time <- now()
cat("Processing started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

setwd("E:/Higor/Cancer_metabolism_project/14- Radar_Siganature_analysis/TMB_PanCan/CNV_vs_TMB/")

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

# Carregar os dados do gene a partir do arquivo RDS
genes <- tryCatch({
  readRDS("/Higor/Cancer_metabolism_project//13- Building_metabolic_signatures/Signatures_by_omic_layer/CNV_signature.rds")
}, error = function(e) {
  cat("Error reading gene data: ", e$message, "\n")
  return(NULL)  # Retorna NULL em caso de falha
})


# Filter rows where Phenotypic_layer contains "Tumor mutational burden"
genes <- genes[grep("Tumor mutational burden", genes$Phenotypic_layer, ignore.case = TRUE), ]



# Verificar se os dados foram carregados corretamente
if (is.null(genes)) {
  cat("Gene data could not be loaded. Check the error message above for details.\n")
} else {
  # Remover genes duplicados
  genes <- genes[!duplicated(genes$Multiomics_Signature), ]
  
  # Ordenar os símbolos gênicos alfabeticamente
  genes$Multiomics_Signature <- sort(genes$Multiomics_Signature)
  
  # Substituir os valores na coluna Multiomics_Signature usando a função process_gene_name()
  genes$Multiomics_Signature_Original <- genes$Multiomics_Signature  # Armazena original para comparação
  genes$Multiomics_Signature <- sapply(genes$Multiomics_Signature, process_gene_name)
  
  # Exibir comparação entre antes e depois
  cat("\nPré-visualização das modificações:\n")
  resultado <- data.frame(
    "Original" = genes$Multiomics_Signature_Original,
    "Processado" = genes$Multiomics_Signature,
    stringsAsFactors = FALSE
  )
  print(resultado)
  
  # Visualizar o DataFrame atualizado
  View(genes)
  
  # Subset the genes column
  current_genes <- unique(genes$Multiomics_Signature)
  
}

# Criar um vetor para armazenar genes com NULL (caso necessário)
genes_with_null <- c()

# Source custom functions
source('source_functions.R')
rds_xena_tmb_name_output <- "int/Target_genes.rds"

# Set the number of genes for each subset
N_FOR_SPLIT <- 500

# Function to handle data retrieval for each gene
get_gene_data <- function(gene) {
  data <- try(get_tmb_data(gene)$data, silent = TRUE)
  if (inherits(data, "try-error")) {
    genes_with_null <<- c(genes_with_null, gene)
    return(NULL)
  } else {
    return(data)
  }
}
# Check if the number of genes is less than or equal to N_FOR_SPLIT
if (length(current_genes) <= N_FOR_SPLIT) {
  # Data with <= 500 genes
  current_tmb_genes <- get_tmb_data(current_genes)
  gene_data <- lapply(current_tmb_genes, function(x) x$data)
  
  saveRDS(gene_data, rds_xena_tmb_name_output)
  
  print("Bingo! I ran Xena data for all genes!")
} else {
  # Data with > 500 genes
  subsets <- split(seq_along(current_genes), 
                   rep(1:ceiling(length(current_genes)/N_FOR_SPLIT),
                       each = N_FOR_SPLIT, 
                       length.out = length(current_genes)))
  full_list <- list()
  
  for (sub in subsets) {
    current_subset <- current_genes[sub]
    current_tmb_genes <- get_tmb_data(current_subset)
    gene_data <- lapply(current_tmb_genes, function(x) x$data)
    
    full_list <- c(full_list, gene_data)
    
    saveRDS(full_list, rds_xena_tmb_name_output)
  }
  
  print("Bingo! I succeeded running Xena data for the gene set!")
}

# Log genes with NULL data into a text file
if (length(genes_with_null) > 0) {
  # Create a data frame with the genes_with_null vector
  genes_null_df <- data.frame(Genes = genes_with_null)
  
  # Write the data frame to a TSV file
  write.table(genes_null_df, "int/genes_with_null.tsv", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  cat("Genes with NULL data logged in 'genes_with_null.tsv'\n")
}

# Print end time
end_time <- now()
cat("Processing ended at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

# Calculate processing time
processing_time <- interval(start_time, end_time)
cat("Processing time:", as.duration(processing_time), "\n")

# Print highest peak of RAM memory used
peak_memory <- round(max(memory.size()) / 1024, 2)
cat("Highest peak of RAM memory used:", peak_memory, "MB", "\n")

# Create output text file
output_file <- "processing_stats.txt"

# Open the file in append mode
file_conn <- file(output_file, "a")

# Write the processing statistics to the file
writeLines(paste("Start Time:", format(start_time, "%Y-%m-%d %H:%M:%S")), file_conn)
writeLines(paste("End Time:", format(end_time, "%Y-%m-%d %H:%M:%S")), file_conn)
writeLines(paste("Processing Time:", as.duration(processing_time)), file_conn)
writeLines(paste("Peak RAM Memory Used:", peak_memory, "MB"), file_conn)

# Close the file connection
close(file_conn)

# Save the current get_xena_workspace to a file named 'get_xena_workspace.RData' in the current working directory
save.image(file = "get_xena_workspace.RData")

# Alternatively, if you want to confirm the file path where the get_xena_workspace is saved, you can use:
current_working_directory <- getwd()  # Get the current working directory
get_xena_workspace_path <- file.path(current_working_directory, "get_xena_workspace.RData")  # Build the file path
save.image(file = get_xena_workspace_path)  # Save the get_xena_workspace

cat("get_xena_workspace saved to: ", get_xena_workspace_path, "\n")  # Print confirmation message with file path

# Print end time
end_time <- now()
cat("Processing ended at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

# Calculate processing time
processing_time <- interval(start_time, end_time)
cat("Processing time:", as.duration(processing_time), "\n")

# Save the session information to a file
get_xena_session_info_path <- file.path(current_working_directory, "get_xena_session_info.txt")
sink(get_xena_session_info_path)  # Redirect output to file
sessionInfo()  # Print session information
sink()  # Stop redirecting output

cat("Session information saved to: ", get_xena_session_info_path, "\n")

