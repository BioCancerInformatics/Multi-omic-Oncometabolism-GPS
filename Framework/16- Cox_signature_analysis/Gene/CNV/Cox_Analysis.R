#### RSCRIPT COX ANALYSIS 
#### HIGOR ALMEIDA, EMANNUEL, ENRIQUE MEDINA-ACOSTA 
#### LAST VERSION - 12/08/2024

# Diretório
setwd("E:/Higor/Cancer_metabolism_project/16- Cox_signature_analysis/Gene/CNV/")

# Carregar os pacotes necessários
library(dplyr)
library(rio)
library(UCSCXenaShiny)

# Dataframe
genes <- import("/Higor/Cancer_metabolism_project/15- FIltering_signature_results/mRNA_miRNA_lncRNA_Protein_signatures/All_correlation_signature_results.tsv")
genes <- genes[order(genes$genes), ]

# Filter the genes to include only rows where genotypic is "CNV"
genes <- genes %>% filter(Genotypic_var == "CNV")

# Função para realizar análise de Cox para uma métrica específica
perform_cox_analysis <- function(genes, measure) {
  results <- data.frame()  # Cria um data frame vazio para armazenar os resultados
  
  # Usar unique para evitar repetições desnecessárias
  unique_genes <- unique(genes$genes)
  
  for (gene in unique_genes) {
    cox_result <- vis_unicox_tree(
      Gene = gene,
      measure = measure,
      data_type = "cnv",
      use_optimal_cutoff = FALSE,
      opt_pancan = .opt_pancan
    )
    
    # Extrai os valores de p.value, Type e cancer
    p_value <- cox_result[["data"]][["p.value"]]
    Type <- cox_result[["data"]][["Type"]]
    cancer <- cox_result[["data"]][["cancer"]]
    
    # Cria um data frame temporário com os resultados atuais
    temp_df <- data.frame(
      Gene = gene,
      p.value = p_value,
      Type = Type,
      cancer = cancer
    )
    
    # Adiciona o data frame temporário aos resultados principais
    results <- rbind(results, temp_df)
  }
  
  # Renomeia as colunas em results para evitar conflitos ao mesclar
  colnames(results) <- c("Gene", paste0("p.value_Cox_", measure), paste0("Type_Cox_", measure), "cancer")
  
  # Mescla os data frames genes e results pelo Gene e cancer_types
  merged_genes <- merge(genes, results[, c("Gene", paste0("p.value_Cox_", measure), paste0("Type_Cox_", measure), "cancer")], 
                        by.x = c("genes", "cancer_types"),
                        by.y = c("Gene", "cancer"),
                        all.x = TRUE)
  
  return(merged_genes)
}

# Aplicar a função para cada métrica e armazenar os resultados no data frame genes
measures <- c("OS", "DSS", "DFI", "PFI")
genes_merged <- genes

# Lista para armazenar os horários de conclusão
end_times <- list()

for (measure in measures) {
  cat("Analisando métrica:", measure, "\n")  # Indicar a métrica atual
  genes_merged <- perform_cox_analysis(genes_merged, measure)
  end_times[[measure]] <- Sys.time()  # Armazenar o horário de conclusão
}

# Salva o resultado final em um arquivo tSV
export(genes_merged, "Cox_results.tsv")

# Imprime os horários de conclusão
print(end_times)

