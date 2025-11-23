# =========================================================
# Load Required Packages (deduplicated and organized)
# =========================================================
suppressPackageStartupMessages({
  
  # Core Shiny
  library(shiny)
  library(shinycssloaders)
  
  # Data handling and utilities
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(purrr)
  library(readr)
  library(readxl)
  library(glue)
  library(memoise)
  
  # Só mantenha qs se realmente ainda for usar em outro lugar
  # library(qs)
  
  # Plotting and visualization
  library(ggplot2)
  library(ggsci)
  library(cowplot)
  library(grid)
  library(gridtext)
  
  # Network & Graph visualization
  library(igraph)
  library(tidygraph)
  library(ggraph)
  
  # Biological / Survival Analysis
  library(survival)
  library(survminer)
  library(fmsb)
  library(ggradar)
  
  # Shiny Data Tables
  library(DT)
  library(htmltools)
  
  # UCSC Xena data utilities
  library(UCSCXenaShiny)
  library(UCSCXenaTools)
  
  # Development tools
  library(devtools)
})

# Pacotes que não precisam estar carregados em produção:
# - rsconnect: usado apenas para deploy via RStudio
# - devtools: apenas desenvolvimento
# Retire-os para não carregar coisas desnecessárias no worker
# library(rsconnect)
# library(devtools)

# ---------------------- Configurações do Shiny ----------------------
options(shiny.maxRequestSize = 500 * 1024^2)  # até 500 MB upload

# Em produção é melhor sanitizar erros (não expor stack trace bruto)
options(shiny.sanitize.errors = TRUE)


# ## -------------------------
# ## 1) Definir caminhos
# ## -------------------------
# arquivos <- c(
#   "data/Dataset_S1.tsv",
#   "data/Dataset_S2.tsv",
#   "data/Dataset_S3.tsv",
#   "data/Dataset_S4.tsv",
#   "data/Dataset_S5.tsv",
#   "data/Search_your_target.rds"
# )
# 
# ## -------------------------
# ## 2) Função de leitura
# ## -------------------------
# ler_arquivo <- function(x){
#   if (grepl("\\.tsv$", x)) {
#     read_tsv(x, show_col_types = FALSE)
#   } else if (grepl("\\.rds$", x)) {
#     readRDS(x)
#   } else {
#     stop(paste("Formato não suportado:", x))
#   }
# }
# 
# ## -------------------------
# ## 3) Ler todos os arquivos  
# ## -------------------------
# all_data <- lapply(arquivos, ler_arquivo)
# 
# ## -------------------------
# ## 4) Nomear a lista
# ## -------------------------
# names(all_data) <- file_path_sans_ext(basename(arquivos))
# 
# ## -------------------------
# ## 5) Salvar
# ## -------------------------
# saveRDS(all_data, file = "All_data.rds")


# ---------------------- Carregamento dos Dados ----------------------
# Estratégia: usar .rds em vez de .qs (evitando o qread pesado)
# e encapsular em uma função memoizada.

load_all_data <- memoise(function() {
  file_path <- file.path("data", "All_data.rds")
  if (!file.exists(file_path)) {
    stop("Arquivo 'data/All_data.rds' não encontrado. Verifique o deploy.")
  }
  message("Carregando All_data.rds...")
  readRDS(file_path)
})

# Atenção:
# Não estamos chamando load_all_data() aqui.
# Ele será chamado no início do server(), o que facilita debugging
# e evita tentar carregar em situações onde o app nem precisa subir.

# TCGA types — provavelmente leve, pode ficar aqui
tcga_types <- readxl::read_excel("data/TCGA_Cancer_types.xlsx")

# ---------------------- Importação de Módulos ----------------------
source("modules/module_user_manual.R")
# source("modules/module_atlas.R")
source("modules/module_search_your_target.R")
source("modules/module_nomenclature.R")
source("modules/module_all_signatures.R")
source("modules/module_signature.R")
# source("modules/module_meaningful_interaction.R")
# source("modules/module_regulatory_circuitry.R")
source("modules/module_interaction_network.R")
# source("modules/module_custom_signatures.R")
source("modules/module_correlation_analysis.R")
source("modules/module_tumor_normal_analysis.R")
source("modules/module_cox_analysis.R")
source("modules/module_survival_analysis.R")
source("modules/module_infiltrates_analysis.R")
source("modules/module_web_resources.R")
# source("modules/module_article.R")
source("modules/module_developers.R")
