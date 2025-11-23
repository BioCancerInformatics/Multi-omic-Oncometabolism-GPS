## UCSCXenaShiny TIL Immune infiltrates Pan_Cancer correlates
## Authors: Higor Almeida
## Last update: 20/12/2024

# Carregando as bibliotecas necessárias
library(UCSCXenaShiny)
library(UCSCXenaTools)
library(dplyr)
library(tidyr)
library(purrr)
library(rio)

# Definindo o diretório
setwd("F:/Higor/Cancer_metabolism_analysis_06/9- Tumor_cell_infiltration_analysis/lncRNA/CNV/")

# Importando o arquivo
Table <- import("/Higor/Cancer_metabolism_analysis_06/8- Survival_analysis/lncRNA/CNV/Updated_Cox_results_with_survival_analysis.tsv")

# Remover as duplicatas de genes em 'Target', mantendo apenas uma linha por gene
Target <- Table %>%
  distinct(genes, .keep_all = TRUE)

# Inicializando um vetor para acompanhar genes que falharam devido à falta de dados
null_genes <- character(0)

# Inicializando uma lista para armazenar os resultados das correlações
til_results_list <- list()
null_genes <- character(0)

# Loop para analisar as correlações dos genes TIL
for (i in seq_along(Target$genes)) {
  gene <- Target$genes[i]  # Usando o gene diretamente sem processamento
  
  # Usando tryCatch apenas na chamada à função vis_gene_TIL_cor
  tryCatch({
    p <- vis_gene_TIL_cor(
      Gene = gene,
      cor_method = "spearman",
      data_type = "cnv",
      sig = c("B cell memory_CIBERSORT",
              "B cell naive_CIBERSORT",
              "B cell plasma_CIBERSORT",
              "Cancer associated fibroblast_XCELL",
              "Class-switched memory B cell_XCELL",
              "Common lymphoid progenitor_XCELL",
              "Endothelial cell_XCELL",
              "Eosinophil_CIBERSORT",
              "Granulocyte-monocyte progenitor_XCELL",
              "Hematopoietic stem cell_XCELL",
              "Macrophage M0_CIBERSORT",
              "Macrophage M1_CIBERSORT",
              "Macrophage M2_CIBERSORT",
              "Mast cell activated_CIBERSORT",
              "Monocyte_CIBERSORT",
              "Myeloid dendritic cell activated_CIBERSORT",
              "Myeloid dendritic cell resting_CIBERSORT",
              "Neutrophil_CIBERSORT",
              "NK cell activated_CIBERSORT",
              "NK cell resting_CIBERSORT",
              "T cell CD4+ memory activated_CIBERSORT",
              "T cell CD4+ memory resting_CIBERSORT",
              "T cell CD4+ naive_CIBERSORT",
              "T cell CD4+ Th1_XCELL",
              "T cell CD4+ Th2_XCELL",
              "T cell CD8+_CIBERSORT",
              "T cell follicular helper_CIBERSORT",
              "T cell gamma delta_CIBERSORT",
              "T cell regulatory (Tregs)_CIBERSORT"),
      Plot = "TRUE"
    )
    
    # Verificando se p[["data"]] existe e é um dataframe
    if (!is.null(p) && inherits(p[["data"]], "data.frame")) {
      # Armazenando os resultados temporariamente em uma lista
      til_results_list[[length(til_results_list) + 1]] <- p[["data"]]
    } else {
      warning("No data for gene: ", gene)
      null_genes <- append(null_genes, gene)
    }
  }, error = function(e) {
    cat("Error processing gene: ", gene, " with error message: ", e$message, "\n")
    null_genes <- append(null_genes, gene)
  })
}

# Verificando se há genes que falharam e escrevendo-os em um arquivo
if (length(null_genes) > 0) {
  null_genes_file_path <- file.path(getwd(), "null_genes_that_failed.txt")
  write.table(null_genes, file = null_genes_file_path, col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("List of genes with no data or errors written to ", null_genes_file_path, "\n")
}

# Combining the results into a single data frame
if (length(til_results_list) > 0) {
  til_results <- do.call(rbind, til_results_list)  # Combining the list into a data frame
} else {
  til_results <- data.frame()  # Se não houver resultados, cria um data frame vazio
}

# Ajuste de p-values por gene e tipo de câncer, levando em consideração os 30 tipos de células imunes
til_results <- til_results %>%
  group_by(gene, cancer) %>%  # Agrupar por gene e tipo de câncer
  mutate(padj = p.adjust(p.value, method = "holm")) %>%  # Ajuste de p-values dentro de cada grupo (gene e câncer)
  ungroup()  # Desagrupar após o ajuste

# Filtrando os dados para incluir apenas linhas onde o FDR seja menor ou igual a 0.05
til_filtered <- til_results %>%
  filter(padj <= 0.05)

# Pegar os dados das colunas immune_cells, cor, p.value do data frame til_filtered e unir ao data frame Target
data <- merge(
  Table, 
  til_filtered %>% select(cancer, gene, immune_cells, cor, padj), 
  by.x = c("cancer_types", "genes"), 
  by.y = c("cancer", "gene"), 
  all.x = TRUE
)

# Exportando os resultados gerais para um arquivo
export(til_filtered, "TIL_results.tsv")

# Exportando os resultados filtrados para um arquivo
export(data, "TIL_filtered_results.tsv")

cat("Analysis complete")
