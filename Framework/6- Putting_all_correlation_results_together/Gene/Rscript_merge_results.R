#### RSCRIPT TO BUILD A FINAL TABLE WITH THE RESULTS OF THE CORRELATIONS, CATEGORIZE THE METABOLISM AND PATHWAY OF EACH GENE AND BUILD SIGNATURES ACCORDING TO METABOLISM OR PATHWAYS
#### HIGOR ALMEIDA, M.sC
#### LAST VERSION - 12/08/2024

library(dplyr)
library(tidyr)
library(rio)

setwd("F:/Higor/Cancer_metabolism_analysis_06/5- Radar_Analysis/")


##### BUSCANDO OS RESULTADOS DE TODAS AS CORRELAÇÕES EM DIFERENTES PASTAS E CRIANDO UMA TABELA GERAL DE RESULTADOS

# Função para importar e adicionar colunas de identificação
import_and_label <- function(filepath, genotypic_var, phenotypic_var) {
  import(filepath) %>%
    mutate(Genotypic_var = genotypic_var, Phenotypic_var = phenotypic_var)
}

# Lista de arquivos e suas respectivas colunas de identificação
file_info <- list(
  list("MSI_PanCan/CNV_vs_MSI/results/gene_cancer_correlation_network_0.1.tsv", "CNV", "MSI"),
  list("MSI_PanCan/Methylation_vs_MSI/results/gene_cancer_correlation_network_0.1.tsv", "Methylation", "MSI"),
  list("MSI_PanCan/mRNA_vs_MSI/results/gene_cancer_correlation_network_0.1.tsv", "Gene expression", "MSI"),
  list("MSI_PanCan/Mutation_vs_MSI/results/gene_cancer_msi_correlation_network_0.1.tsv", "Mutation", "MSI"),
  list("MSI_PanCan/Protein_vs_MSI/results/gene_cancer_msi_correlation_network_0.1.tsv", "Protein expression", "MSI"),
  list("Stemness_PanCan/CNV_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "CNV", "Stemness"),
  list("Stemness_PanCan/Methylation_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "Methylation", "Stemness"),
  list("Stemness_PanCan/mRNA_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "Gene expression", "Stemness"),
  list("Stemness_PanCan/Mutation_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "Mutation", "Stemness"),
  list("Stemness_PanCan/Protein_vs_Stemness/results/gene_cancer_correlation_network_0.1.tsv", "Protein expression", "Stemness"),
  list("TMB_PanCan/CNV_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "CNV", "TMB"),
  list("TMB_PanCan/Methylation_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "Methylation", "TMB"),
  list("TMB_PanCan/mRNA_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "Gene expression", "TMB"),
  list("TMB_PanCan/Mutation_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "Mutation", "TMB"),
  list("TMB_PanCan/Protein_vs_TMB/results/gene_cancer_correlation_network_0.1.tsv", "Protein expression", "TMB")
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

setwd("F:/Higor/Cancer_metabolism_analysis_06/6- Merge all correlation results/Gene/")

export(All_results, "All_results_gene.tsv")
