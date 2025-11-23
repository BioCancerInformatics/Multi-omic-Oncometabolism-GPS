## UCSCXenaShiny Pancan omics-phenotype correlations
## Authors: Juan Carlo S Silva, Emanuell Rodrigues de Souza, Leonardo Henrique da Silva, Higor Almeida C Nogueira, Ana Beatriz Garcia, Enrique Medina-Acosta
## Updated: 03/12/2024

### MAIN CODE - CNV vs MSI

# Load required libraries
library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(fmsb)
library(UCSCXenaShiny)

setwd("E:/Higor/Cancer_metabolism_project/14- Radar_Siganature_analysis/MSI_PanCan/CNV_vs_MSI/")

start_time <- now()
cat("Processing started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")

source('source_functions.R')
setTimeLimit(cpu=Inf)

######### GET GENE-PANCANCER CORRELATIONS ########
msi_gene_data = read_rds("int/Target_genes.rds")

####
# Filter out genes with null data
empty_gene_names <- names(msi_gene_data)[sapply(msi_gene_data, function(x) all(is.na(x)))]

# Check if any empty gene names were found
if (length(empty_gene_names) > 0) {
  cat("Empty gene names found:", empty_gene_names, "\n")
  # Save the empty gene names to a .txt file in the working directory
  writeLines(as.character(empty_gene_names), "int/empty_gene_names.txt")
} else {
  cat("No empty gene names found.\n")
}

####
empty_gene_names = names(msi_gene_data[which(sapply(msi_gene_data, is.null))])
if (length(empty_gene_names) > 0){
  msi_gene_data = msi_gene_data[-which(sapply(msi_gene_data, is.null))]
}

gene_sig_list = lapply(msi_gene_data, pivot_wider_manager)
gene_sig_list = gene_sig_list[lengths(gene_sig_list) > 0]

gene_sig_table = do.call(bind_rows, gene_sig_list)
gene_sig_table = as.data.frame(gene_sig_table)

row.names(gene_sig_table) = names(gene_sig_list)

rows_with_na = row.names(gene_sig_table[
  which(apply(gene_sig_table, 1, function(row) all(is.na(row)))),])

###### Save
write.table(empty_gene_names,
            file = "int/null_genes_pancan_msi_radar.tsv",
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

# Add the row names as a new column at the beginning of the data frame
gene_sig_table <- data.frame(Genes = rownames(gene_sig_table), gene_sig_table)

# Reset the row names to NULL since they are now part of the data
rownames(gene_sig_table) <- NULL

# Write the modified gene_sig_table to a file
write.table(gene_sig_table,
            file = "results/sig_genes_corr_pancan_msi.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)  # No longer using row.names, as "Genes" is now a regular column

##################################### SET NETWORK ----
network_table <- fread('results/sig_genes_corr_pancan_msi.tsv')
colnames(network_table)[1] <- 'genes'

cancers <- colnames(network_table[, 2:ncol(network_table)])
first_cancer <- cancers[1]
last_cancer <- cancers[length(cancers)]

##### Set up edges ----
network_table <- network_table %>% 
  mutate(across(all_of(first_cancer):all_of(last_cancer), as.numeric)) %>%
  pivot_longer(cols = all_of(first_cancer):all_of(last_cancer), 
               names_to = 'cancer_types',
               values_to = 'correlation')

# INÍCIO DA PARTE ADICIONADA #
network_table = network_table[!is.na(network_table$correlation),]

# Add p.value and padj columns to network_table
network_table$Correlation_p.adj <- NA

# Loop through each gene and add p.adj values from msi_gene_data
for (gene in unique(network_table$genes)) {
  if (gene %in% names(msi_gene_data)) {
    # Extract all p.adj values for the current gene
    Correlation_p.adj <- msi_gene_data[[gene]][["padj"]]
    
    # Get the corresponding rows in network_table for the current gene
    gene_rows <- network_table$genes == gene
    
    # Assign p.adj values to corresponding rows in network_table
    network_table[gene_rows, "Correlation_p.adj"] <- Correlation_p.adj
  }
}

# Verificar se a coluna "Correlation_p.adj" existe na tabela
if("Correlation_p.adj" %in% colnames(network_table)) {
  # Adicionar uma nova coluna para armazenar os valores transformados pelo log10
  network_table$log10_correlation <- -log10(network_table$Correlation_p.adj)
  
  # Imprimir uma mensagem de conclusão
  print("coluna 'log10(p.adj)' adicionada com sucesso!")
}

#### DIFFERENTIAL EXPRESSION ANALYSIS USING "Visualize Gene TPM in Single Cancer Type (TCGA vs GTEX) ####

# Inicialize as colunas de expressão como NA por padrão
network_table$Expression <- rep(NA, nrow(network_table))
network_table$Expression_p.sgnif <- rep(NA, nrow(network_table))
network_table$Expression_p.adj <- rep(NA, nrow(network_table))

# Loop para analisar cada gene de acordo com o tipo de câncer
for (i in 1:nrow(network_table)) {
  cancer_type <- network_table$cancer_types[i]
  genes <- network_table$genes[i]
  
  # Use tryCatch para lidar com possíveis erros ao chamar vis_toil_TvsN_cancer
  tryCatch({
    # Analisar a expressão diferencial usando a função vis_toil_TvsN_cancer
    Results_expression <- vis_toil_TvsN_cancer(
      Gene = genes,
      Mode = "Violinplot",
      data_type = "mRNA",
      Show.P.value = TRUE,
      Show.P.label = FALSE,
      Method = "wilcox.test",
      values = c("#DF2020", "#DDDF21"),
      TCGA.only = FALSE,
      Cancer = cancer_type
    )
    
    # Inicializar a expressão como "Unchanged"
    network_table$Expression[i] <- "Unchanged"
    
    # Verificar se há dados de expressão retornados
    if (!is.null(Results_expression[["data"]])) {
      # Extrair dados de expressão
      tpm_data <- Results_expression[["data"]][["value"]]
      type2_data <- Results_expression[["data"]][["type2"]]
      
      # Verificar se há dados suficientes para calcular as medianas
      if (sum(type2_data == "tumor", na.rm = TRUE) > 0 && sum(type2_data == "normal", na.rm = TRUE) > 0) {
        # Calcular as medianas dos grupos tumor e normal
        tumor_median <- median(tpm_data[type2_data == "tumor"], na.rm = TRUE)
        normal_median <- median(tpm_data[type2_data == "normal"], na.rm = TRUE)
        
        # Determinar se é Overexpression ou Underexpression
        if (tumor_median > normal_median) {
          network_table$Expression[i] <- "Overexpression"
        } else {
          network_table$Expression[i] <- "Underexpression"
        }
      }
    }
    
    # Verificar se há p-values retornados
    if (length(Results_expression$layers[[3]]$data) > 0) {
      # Extrair os p-values
      p_values <- Results_expression$layers[[3]]$data
      
      # Adicionar os p-values às colunas correspondentes
      network_table$Expression_p.sgnif[i] <- p_values$p.signif
      network_table$Expression_p.adj[i] <- p_values$p.adj
      
      # Classificar como "Unchanged" se o valor de p for maior ou igual a 0.05
      if (p_values$p.adj >= 0.05) {
        network_table$Expression[i] <- "Unchanged"
      }
    }
  }, error = function(e) {
    # Captura e imprime o erro, mas continua a execução do loop
    print(paste("Erro ao analisar câncer:", cancer_type))
    print(e)
  })
}

# Verificar se a coluna "Expression_p.adj" existe na tabela
if("Expression_p.adj" %in% colnames(network_table)) {
  # Adicionar uma nova coluna para armazenar os valores transformados pelo log10
  network_table$log10_p.adj <- -log10(network_table$Expression_p.adj)
  
  # Imprimir uma mensagem de conclusão
  print("coluna 'log10(p.adj)' adicionada com sucesso!")
}

# Remover linhas duplicadas, se houver
network_table <- unique(network_table)

# FINAL DA PARTE ADICIONADA #

# Lista dos valores de corte de correlação
correlation_cutoffs <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)

# Initialize an empty vector to store status messages
status_messages <- c()

# Loop through the correlation cutoffs
for (cutoff in correlation_cutoffs) {
  init_message <- paste("Processing for correlation cutoff:", cutoff)
  
  # Append the initial message to status_messages
  status_messages <- c(status_messages, init_message)
  
  # Print the initial message to the console
  cat(init_message, "\n")
  
  cat("Processing for correlation cutoff:", cutoff, "\n")
  
  # Atualize a variável de corte de correlação no código
  correlation_cutoff = cutoff
  
  # Filtre com base no valor de corte de correlação
  network_table_filtered <- subset(network_table, correlation >= correlation_cutoff | correlation <= -correlation_cutoff)
  
  # Verifique se há pelo menos um resultado antes de salvar
  if (nrow(network_table_filtered) > 0) {
    
    tryCatch({
      ##### Set up nodes ----
      nodes = get_gene_cancer_node_table(network_table_filtered)
      
      ###### Get cancer types for each gene
      cancers_by_genes = aggregate(cancer_types ~ genes,
                                   data = network_table_filtered,
                                   FUN = function(x) paste(unique(x),
                                                           collapse = "/"))
      nodes = left_join(nodes, cancers_by_genes, by=join_by(nodes == genes))
      
      # Function to check signature sign directions
      sign_verification <- function(network_table) {
        has_positive <- any(network_table$correlation > 0)
        has_negative <- any(network_table$correlation < 0)
        
        if (has_positive & has_negative) {
          return(c('positive', 'negative'))
        } else if (has_positive) {
          return('positive')
        } else if (has_negative) {
          return('negative')
        } else {
          return(character(0))  # Return an empty vector if there are no correlations.
        }
      }
      
      # Salve os resultados com o nome de arquivo apropriado
      network_filename <- paste0('results/gene_cancer_correlation_network_', cutoff, '.tsv')
      nodes_filename <- paste0('results/gene_cancer_nodes_network_', cutoff, '.tsv')
      
      write.table(network_table_filtered,
                  network_filename,
                  sep = '\t',
                  quote = FALSE,
                  row.names = FALSE)
      
      write.table(nodes,
                  nodes_filename,
                  sep = '\t',
                  quote = FALSE,
                  row.names = FALSE)
      
      completion_message <- paste("Processing for correlation cutoff:", cutoff, "completed.")
      cat(completion_message, "\n")
      
      # Append the completion message to status_messages
      status_messages <- c(status_messages, completion_message)
      
    }, error = function(e) {
      no_result_message <- paste("No results for correlation cutoff:", cutoff)
      cat(no_result_message, "\n")
      
      # Append the no result message to status_messages
      status_messages <- c(status_messages, no_result_message)
    })
    
  } else {
    no_result_message <- paste("No results for correlation cutoff:", cutoff)
    cat(no_result_message, "\n")
    
    # Append the no result message to status_messages
    status_messages <- c(status_messages, no_result_message)
  }
}

  # Save the status messages
  write.table(status_messages, 
              "results/cutoff_status.tsv", 
              sep = "\t", 
              quote = FALSE, 
              col.names = FALSE, 
              row.names = FALSE)

  # Define the path to the 'results' subdirectory
  results_directory <- file.path(getwd(), "results")
  
  # Locate the files with the highest numerical suffix in the 'results' directory
  files <- list.files(path = results_directory,
                      pattern = "gene_cancer_correlation_network_\\d+\\.\\d+\\.tsv",
                      full.names = TRUE)
  
  # Process to find the highest numerical suffix
  if (length(files) > 0) {
    suffixes <- gsub(".*gene_cancer_correlation_network_", "", files)
    suffixes <- gsub("\\.tsv", "", suffixes)
    numeric_suffixes <- as.numeric(suffixes)
    
    if (length(numeric_suffixes) > 0) {
      highest_suffix <- max(numeric_suffixes, na.rm = TRUE)
      highest_serial_file <- files[which(numeric_suffixes == highest_suffix)]
      
      # Read the file with the highest serial number
      highest_data <- read.delim(highest_serial_file, stringsAsFactors = FALSE)
      highest_suffix_found <- highest_suffix
    }
  }
  
  # Optional: Write the highest suffix found to a text file and print it
  highest_suffix_found_path <- file.path(results_directory, "highest_suffix_found.txt")
  write.table(highest_suffix_found,  # assuming 'highest_suffix_found' is your data
              highest_suffix_found_path,
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE)
  
  # Optional: Print the highest suffix found
  cat("The highest suffix found is: ", highest_suffix_found, "\n")
  
end_time <- now()
cat("Processing ended at:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")

processing_time <- interval(start_time, end_time)
cat("Processing time:", as.duration(processing_time), "\n")

peak_memory <- round(max(memory.size()) / 1024, 2)
cat("Highest peak of RAM memory used:", peak_memory, "MB", "\n")

output_file <- "processing_stats.txt"

file_conn <- file(output_file, "a")

writeLines(paste("Start Time:", format(start_time, "%Y-%m-%d %H:%M:%S")), file_conn)
writeLines(paste("End Time:", format(end_time, "%Y-%m-%d %H:%M:%S")), file_conn)
writeLines(paste("Processing Time:", as.duration(processing_time)), file_conn)
writeLines(paste("Peak RAM Memory Used:", peak_memory, "MB"), file_conn)

# Close the file connection
close(file_conn)

# Save the current main_workspace to a file named 'main_workspace.RData' in the current working directory
save.image(file = "main_workspace.RData")

# Alternatively, if you want to confirm the file path where the main_workspace is saved, you can use:
current_working_directory <- getwd()  # Get the current working directory
main_workspace_path <- file.path(current_working_directory, "main_workspace.RData")  # Build the file path
save.image(file = main_workspace_path)  # Save the main_workspace

cat("main_workspace saved to: ", main_workspace_path, "\n")  # Print confirmation message with file path

# Save the session information to a file
main_session_info_path <- file.path(current_working_directory, "main_session_info.txt")
sink(main_session_info_path)  # Redirect output to file
sessionInfo()  # Print session information
sink()  # Stop redirecting output

cat("Session information saved to: ", main_session_info_path, "\n")

