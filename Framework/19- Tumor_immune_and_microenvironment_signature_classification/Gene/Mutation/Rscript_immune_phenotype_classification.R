#### RSCRIPT TO CLASSIFY THE TUMOR IMMUNE PHENOTYPE ADAPTATIVELY

# Carregar pacotes necessários
library(dplyr)
library(tidyr)
library(purrr)
library(rio)

setwd("E:/Oncometabolism_GPS/19- Tumor_immune_and_microenvironment_signature_classification/Gene/Mutation/")

# Importar dados
data <- import("Microenvironment_classification.tsv")

# Definir os critérios de "hot" e "cold"
hot_criteria <- c("T cell CD8+_CIBERSORT", "NK cell activated_CIBERSORT", 
                  "Macrophage M1_CIBERSORT", "T cell CD4+ memory activated_CIBERSORT", 
                  "Myeloid dendritic cell activated_CIBERSORT")
cold_criteria <- c("Macrophage M2_CIBERSORT", "T cell regulatory (Tregs)_CIBERSORT")

# Função para classificar tumores com base nos critérios
tumor_classification <- function(cor_data) {
  # Inicializar contadores
  hot_score <- 0
  cold_score <- 0
  hot_direct <- FALSE
  cold_direct <- FALSE
  score_details <- ""
  
  # Ajustar valores de 'Expression' vazios
  cor_data$Expression <- ifelse(trimws(cor_data$Expression) == "", "Unchanged", cor_data$Expression)
  
  # Verificar e ajustar valores de 'Genotypic_var'
  cor_data$Expression <- ifelse(!cor_data$Genotypic_var %in% c("Gene expression", "miRNA expression", "Transcript expression"),
                                "Unchanged", cor_data$Expression)
  
  # Analisar cada célula imune
  for (i in 1:nrow(cor_data)) {
    row <- cor_data[i, ]
    
    cell_type <- row$immune_cells
    rho <- row$cor
    expression_type <- row$Expression
    
    # Determinar tipo de correlação
    correlation_type <- ifelse(
      (expression_type %in% c("Overexpression", "Unchanged") && rho > 0) ||
        (expression_type == "Underexpression" && rho < 0),
      "Direta",
      "Inversa"
    )
    
    # Determinar o tipo de célula
    cell_classification <- ifelse(cell_type %in% hot_criteria, "Imunoefetora", 
                                  ifelse(cell_type %in% cold_criteria, "Imunossupressora", NA))
    
    # Atualizar contadores e detalhes
    if (!is.na(cell_classification)) {
      if (cell_classification == "Imunoefetora" && correlation_type == "Direta") {
        hot_score <- hot_score + 1
        hot_direct <- TRUE
      } else if (cell_classification == "Imunossupressora" && correlation_type == "Direta") {
        cold_score <- cold_score + 1
        cold_direct <- TRUE
      } else if (cell_classification == "Imunoefetora" && correlation_type == "Inversa") {
        cold_score <- cold_score + 1
      } else if (cell_classification == "Imunossupressora" && correlation_type == "Inversa") {
        hot_score <- hot_score + 1
      }
      
      # Adicionar detalhe formatado à string
      score_details <- paste0(score_details, 
                              sprintf("immune_cells: %s, rho: %.6f, correlation_type: %s, cell_type: %s\n", 
                                      cell_type, rho, correlation_type, cell_classification))
    }
  }
  
  # Classificar baseado nos critérios
  if (hot_direct && !cold_direct && cold_score == 0) {
    classification <- "Quente"
  } else if (!hot_direct && cold_direct && hot_score == 0) {
    classification <- "Frio"
  } else if (hot_direct && cold_direct) {
    classification <- "Variável"
  } else {
    classification <- "Frio" # Caso padrão para ausência de correlações válidas
    score_details <- "No significant results for immune infiltrates"
  }
  
  return(list(classification = classification, score_details = score_details))
}

# Agrupar dados por gene e câncer e aplicar a classificação
result <- data %>%
  group_by(genes, cancer_types) %>%
  nest() %>%
  mutate(classification_result = map(data, ~tumor_classification(.x))) %>%
  mutate(classification = map_chr(classification_result, "classification"),
         score_details = map_chr(classification_result, "score_details")) %>%
  select(-classification_result) %>%
  unnest(cols = c(data)) %>%
  mutate(classification = classification[match(paste(genes, cancer_types), paste(genes, cancer_types))],
         score_details = score_details[match(paste(genes, cancer_types), paste(genes, cancer_types))])

# Renomear as colunas conforme solicitado
result <- result %>%
  rename(
    Immune_classification = classification,
    Immune_score_details = score_details
  )

# Exportar dados
export(result, "Immune_classification.tsv")

