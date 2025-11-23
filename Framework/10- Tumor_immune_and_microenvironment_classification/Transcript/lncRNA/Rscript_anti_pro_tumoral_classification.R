# Carregar as bibliotecas necessárias
library(tidyverse)
library(rio)

setwd("E:/Higor/Cancer_metabolism_project/10- Tumor_immune_and_microenvironment_classification/Transcript/lncRNA/")

# Carregar os dados
data <- import("/Higor/Cancer_metabolism_project/9- Tumor_cell_infiltration_analysis/Transcript/lncRNA/TIL_filtered_results_v1.tsv")
data <- data %>%
  filter(cor >= 0.1 | cor <= -0.1)

# Verificar a distribuição das expressões
expression_dist <- data %>% 
  group_by(Expression) %>% 
  summarize(Count = n())

# Definir a tabela de classificação das células imunes
cell_classification <- tibble(
  Cell_Type = c("B cell memory_CIBERSORT", "B cell naive_CIBERSORT", "B cell plasma_CIBERSORT",
                "Cancer associated fibroblast_XCELL", "Class-switched memory B cell_XCELL",
                "Common lymphoid progenitor_XCELL", "Endothelial cell_XCELL", "Eosinophil_CIBERSORT",
                "Granulocyte-monocyte progenitor_XCELL", "Hematopoietic stem cell_XCELL", 
                "Macrophage M0_CIBERSORT", "Macrophage M1_CIBERSORT", "Macrophage M2_CIBERSORT", 
                "Mast cell activated_CIBERSORT", "Monocyte_CIBERSORT", "Myeloid dendritic cell activated_CIBERSORT",
                "Myeloid dendritic cell resting_CIBERSORT", "Neutrophil_CIBERSORT", 
                "NK cell activated_CIBERSORT", "NK cell resting_CIBERSORT", 
                "T cell CD4+ memory activated_CIBERSORT", "T cell CD4+ memory resting_CIBERSORT",
                "T cell CD4+ naive_CIBERSORT", "T cell CD4+ Th1_XCELL", "T cell CD4+ Th2_XCELL", 
                "T cell CD8+_CIBERSORT", "T cell follicular helper_CIBERSORT", 
                "T cell gamma delta_CIBERSORT", "T cell regulatory (Tregs)_CIBERSORT"),
  Classification = c("dual", "dual", "dual", "pro-tumoral", "dual", "dual",
                     "pro-tumoral", "dual", "pro-tumoral", "pro-tumoral", "dual", "anti-tumoral",
                     "pro-tumoral", "dual", "dual", "anti-tumoral", "dual", "dual", 
                     "anti-tumoral", "dual", "anti-tumoral", "dual", "dual", "anti-tumoral", 
                     "pro-tumoral", "anti-tumoral", "dual", "dual", "pro-tumoral")
)

# Combinar dados de correlação com a classificação das células imunes
merged_data <- data %>% 
  left_join(cell_classification, by = c("immune_cells" = "Cell_Type"))

# Função para classificar assinaturas com base nas correlações e retornar detalhes da pontuação
classify_signature <- function(data) {
  if (nrow(data) == 0 || all(is.na(data$cor))) {
    return(list(Classification = "NS", Score_Details = "No data"))
  }
  
  anti_tumoral_score <- 0
  pro_tumoral_score <- 0
  dual_score <- 0
  
  anti_tumoral_cells <- c()
  pro_tumoral_cells <- c()
  dual_cells <- c()
  
  for (i in 1:nrow(data)) {
    row <- data[i, ]
    
    # Check if Genotypic_var is valid
    if (!row$Genotypic_var %in% c("Gene expression", "miRNA expression", "Transcript expression")) {
      row$Expression <- "Unchanged"
    }
    
    # Tratar valores ausentes ou strings vazias na coluna Expression
    if (is.na(row$Expression) || row$Expression == "") {
      row$Expression <- "Unchanged"
    }
    
    # Skip if any critical column is missing
    if (is.na(row$cor) || is.na(row$Classification)) {
      next
    }
    
    cell_info <- paste(row$immune_cells, "(", row$cor, ")", sep = "")
    correlation_value <- row$cor
    
    # Determine direct vs. inverse correlation
    is_direct <- (row$Expression == "Overexpression" && correlation_value > 0) ||
      (row$Expression == "Underexpression" && correlation_value < 0) ||
      (row$Expression == "Unchanged" && correlation_value > 0)
    is_inverse <- (row$Expression == "Overexpression" && correlation_value < 0) ||
      (row$Expression == "Underexpression" && correlation_value > 0) ||
      (row$Expression == "Unchanged" && correlation_value < 0)
    
    # Update scores based on direct or inverse correlation
    if (row$Classification == "anti-tumoral") {
      if (is_direct) {
        anti_tumoral_score <- anti_tumoral_score + abs(correlation_value)
      } else if (is_inverse) {
        anti_tumoral_score <- anti_tumoral_score - abs(correlation_value)
      }
      anti_tumoral_cells <- c(anti_tumoral_cells, cell_info)
    } else if (row$Classification == "pro-tumoral") {
      if (is_direct) {
        pro_tumoral_score <- pro_tumoral_score + abs(correlation_value)
      } else if (is_inverse) {
        pro_tumoral_score <- pro_tumoral_score - abs(correlation_value)
      }
      pro_tumoral_cells <- c(pro_tumoral_cells, cell_info)
    } else if (row$Classification == "dual") {
      if (is_direct) {
        dual_score <- dual_score + abs(correlation_value)
      } else if (is_inverse) {
        dual_score <- dual_score - abs(correlation_value)
      }
      dual_cells <- c(dual_cells, cell_info)
    }
  }
  
  # Determine final classification
  classification <- if (anti_tumoral_score > pro_tumoral_score & anti_tumoral_score > dual_score) {
    "anti-tumoral"
  } else if (pro_tumoral_score > anti_tumoral_score & pro_tumoral_score > dual_score) {
    "pro-tumoral"
  } else {
    "dual"
  }
  
  # Prepare detailed output
  score_details <- paste("Anti-tumoral:", anti_tumoral_score, 
                         "(", paste(anti_tumoral_cells, collapse = ", "), ")",
                         "pro-tumoral:", pro_tumoral_score, 
                         "(", paste(pro_tumoral_cells, collapse = ", "), ")",
                         "dual:", dual_score, 
                         "(", paste(dual_cells, collapse = ", "), ")")
  
  return(list(Classification = classification, Score_Details = score_details))
}


# Aplicar a função de classificação aos dados agrupados por assinatura e tipo de câncer
classification_result <- merged_data %>%
  group_by(genes, cancer_types) %>%
  mutate(Result = list(classify_signature(pick(everything())))) %>%
  unnest_wider(Result, names_sep = "_result") %>% # Add a suffix to avoid name conflicts
  ungroup()

# Renomear as colunas conforme solicitado
classification_result <- classification_result %>%
  rename(
    Infiltrate_profile = Classification,
    Classification = Result_resultClassification,
    Microenvironment_score_details = Result_resultScore_Details
  )

# Exportar dados
export(classification_result, "Microenvironment_classification.tsv")

