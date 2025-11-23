#### RSCRIPT 
#### HIGOR ALMEIDA, M.sC
#### LAST VERSION - 12/08/2024

library(dplyr)
library(tidyr)
library(rio)

setwd("E:/Higor/Cancer_metabolism_project/14- Radar_Siganature_analysis/")


##### BUSCANDO OS RESULTADOS DE TODAS AS CORRELAÇÕES EM DIFERENTES PASTAS E CRIANDO UMA TABELA GERAL DE RESULTADOS

# Função para importar e adicionar colunas de identificação
import_and_label <- function(filepath, genotypic_var, phenotypic_var) {
  import(filepath) %>%
    mutate(Genotypic_var = genotypic_var, Phenotypic_var = phenotypic_var)
}

# Lista de arquivos e suas respectivas colunas de identificação
file_info <- list(
  list("MSI_PanCan/CNV_vs_MSI/results/gene_cancer_correlation_network_0.1.tsv", "CNV", "Microsatellite instability"),
  list("MSI_PanCan/Methylation_vs_MSI/results/gene_cancer_correlation_network_0.1.tsv", "Methylation", "Microsatellite instability"),
  list("MSI_PanCan/mRNA_vs_MSI/results/gene_cancer_correlation_network_0.1.tsv", "Gene expression", "Microsatellite instability"),
  list("MSI_PanCan/Mutation_vs_MSI/results/gene_cancer_msi_correlation_network_0.1.tsv", "Mutation", "Microsatellite instability"),
  list("MSI_PanCan/Protein_vs_MSI/results/gene_cancer_msi_correlation_network_0.1.tsv", "Protein expression", "Microsatellite instability"),
  list("MSI_PanCan/miRNA_vs_MSI/results/gene_cancer_correlation_network_0.1.tsv", "miRNA expression", "Microsatellite instability"),
  list("Stemness_PanCan/CNV_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "CNV", "Stemness"),
  list("Stemness_PanCan/Methylation_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "Methylation", "Stemness"),
  list("Stemness_PanCan/mRNA_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "Gene expression", "Stemness"),
  list("Stemness_PanCan/Mutation_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "Mutation", "Stemness"),
  list("Stemness_PanCan/Protein_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "Protein expression", "Stemness"),
  list("Stemness_PanCan/miRNA_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "miRNA expression", "Stemness"),
  list("TMB_PanCan/CNV_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "CNV", "Tumor mutational burden"),
  list("TMB_PanCan/Methylation_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "Methylation", "Tumor mutational burden"),
  list("TMB_PanCan/mRNA_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "Gene expression", "Tumor mutational burden"),
  list("TMB_PanCan/Mutation_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "Mutation", "Tumor mutational burden"),
  list("TMB_PanCan/Protein_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "Protein expression", "Tumor mutational burden"),
  list("TMB_PanCan/miRNA_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "miRNA expression", "Tumor mutational burden")
)

# Importar e combinar todos os dados
All_results <- file_info %>%
  lapply(function(info) {
    filepath <- info[[1]]
    genotypic_var <- info[[2]]
    phenotypic_var <- info[[3]]
    import_and_label(filepath, genotypic_var, phenotypic_var)
  }) %>%
  bind_rows()

genes <- import("/Higor/Cancer_metabolism_project/13- Building_metabolic_signatures/Gene_lncRNA_miRNA_Protein_signatures.tsv")

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
genes$Multiomics_Signature_Original <- genes$Multiomics_Signature  # Armazena original para comparação
genes$Multiomics_Signature <- sapply(genes$Multiomics_Signature, process_gene_name)

# Remove duplicates based on the specified columns
Signatures <- genes %>%
  distinct(Multiomics_Signature, Cancer_types, Omic_layer, Phenotypic_layer, .keep_all = TRUE)

# Rename columns in Signatures to match All_results for join
Signatures_renamed <- Signatures %>%
  rename(
    genes = Multiomics_Signature,
    cancer_types = Cancer_types,
    Genotypic_var = Omic_layer,
    Phenotypic_var = Phenotypic_layer
  )

# Perform the semi-join to keep only matching rows
Filtered_results <- All_results %>%
  semi_join(Signatures_renamed, by = c("genes", "cancer_types", "Genotypic_var", "Phenotypic_var"))

setwd("E:/Higor/Cancer_metabolism_project/15- FIltering_signature_results/mRNA_miRNA_lncRNA_Protein_signatures/")

export(Filtered_results, "All_correlation_signature_results.tsv")







# 
# 
# ## ============================================================
# ## Comparar combinações (Gene, Genotypic_var, Phenotypic_var)
# ## entre dois data frames
# ## ============================================================
# 
# # Seleciona apenas as 3 colunas de interesse em cada df
# df1_sel <- Filtered_results %>%
#   dplyr::select(genes, Genotypic_var, Phenotypic_var) %>%
#   dplyr::distinct()
# 
# df2_sel <- Signatures_renamed %>%
#   dplyr::select(genes, Genotypic_var, Phenotypic_var) %>%
#   dplyr::distinct()
# 
# # Combinações presentes só em df1
# only_in_df1 <- dplyr::anti_join(df1_sel, df2_sel,
#                                 by = c("genes", "Genotypic_var", "Phenotypic_var")) %>%
#   dplyr::mutate(Source = "Filtered_results")
# 
# # Combinações presentes só em df2
# only_in_df2 <- dplyr::anti_join(df2_sel, df1_sel,
#                                 by = c("genes", "Genotypic_var", "Phenotypic_var")) %>%
#   dplyr::mutate(Source = "Signatures_renamed")
# 
# # Junta os dois conjuntos em uma tabela final
# diff_combinations <- dplyr::bind_rows(only_in_df1, only_in_df2)
# 
# print(diff_combinations)
# 
# ## Se quiser salvar
# rio::export(diff_combinations, "diff_combinations.tsv")
# 
# 
# 
