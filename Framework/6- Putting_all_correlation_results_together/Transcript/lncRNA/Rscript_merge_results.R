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
  list("MSI_PanCan/Transcript_vs_MSI/lncRNA/results/gene_cancer_correlation_network_0.1.tsv", "Transcript expression", "MSI"),
  list("TMB_PanCan/Transcript_vs_TMB/lncRNA/results/gene_cancer_correlation_network_0.1.tsv", "Transcript expression", "TMB"),
  list("Stemness_PanCan/Transcript_vs_Stemness/lncRNA/results/gene_cancer_correlation_network_0.1.tsv", "Transcript expression", "Stemness")
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

setwd("F:/Higor/Cancer_metabolism_analysis_06/6- Merge all correlation results/Transcript/lncRNA/")

export(All_results, "All_results_transcript.tsv")
