# =========================================================
# server.R
# =========================================================

# (Opcional) Contador de visitas: arquivo local
counter_file <- "counter.txt"

update_counter <- function() {
  # Em shinyapps.io, o sistema de arquivos é efêmero,
  # mas é gravável durante a sessão. Ainda assim,
  # envolvemos em tryCatch para evitar que erros de I/O
  # impeçam o app de subir.
  tryCatch({
    if (!file.exists(counter_file)) {
      writeLines("0", counter_file)
    }
    count <- as.integer(readLines(counter_file))
    if (is.na(count)) count <- 0
    count <- count + 1
    writeLines(as.character(count), counter_file)
    count
  }, error = function(e) {
    message("Falha ao atualizar o contador: ", conditionMessage(e))
    NA_integer_
  })
}

server <- function(input, output, session) {
  session$allowReconnect(TRUE)  # ok, pode manter
  
  # ---------------------- Carregar all_data ----------------------
  # Aqui usamos a função memoizada definida no global.R
  all_data <- NULL
  all_data <- tryCatch(
    {
      load_all_data()
    },
    error = function(e) {
      # Essa mensagem aparecerá no log do shinyapps.io
      message("Erro ao carregar all_data.rds: ", conditionMessage(e))
      NULL
    }
  )
  
  if (is.null(all_data)) {
    # Se chegou aqui, provavelmente o arquivo está faltando
    # ou estourando limite. Você pode exibir uma mensagem amigável
    # na interface se quiser.
    showModal(
      modalDialog(
        title = "Erro de inicialização",
        "Não foi possível carregar os dados necessários (all_data.rds). 
        Verifique se o arquivo foi incluído corretamente no deploy 
        e se o tamanho não excede os limites do shinyapps.io.",
        easyClose = TRUE
      )
    )
    # Podemos retornar cedo, evitando erros em cascata:
    return(invisible(NULL))
  }
  
  # ---------------------- Contador de visitas ----------------------
  visitor_count <- update_counter()
  
  output$visitor_count <- renderText({
    if (is.na(visitor_count)) {
      "Total: (contador indisponível no servidor)"
    } else {
      paste("Total:", visitor_count)
    }
  })
  
  # ---------------------- Reativos principais ----------------------
  # Aqui assumo que all_data é uma lista com esses nomes.
  # Se o nome de algum elemento for diferente, é só ajustar.
  
  search_your_target      <- reactiveVal(all_data$Search_your_target)
  Target                  <- reactiveVal(all_data$Dataset_S3)
  
  # cnv_signature           <- reactiveVal(all_data$cnv_signature)
  # gene_signature          <- reactiveVal(all_data$gene_signature)
  # methylation_signature   <- reactiveVal(all_data$methylation_signature)
  # mutation_signature      <- reactiveVal(all_data$mutation_signature)
  # mirna_signature         <- reactiveVal(all_data$mirna_signature)
  # protein_signature       <- reactiveVal(all_data$protein_signature)
  # transcript_signature    <- reactiveVal(all_data$transcript_signature)
  
  meaningful_interaction  <- reactiveVal(all_data$Dataset_S4)
  regulatory_circuitry    <- reactiveVal(all_data$Dataset_S5)
  # atlas                  <- reactiveVal(all_data$Dataset_S1)
  
  # ---------------------- Chamada dos módulos ----------------------
  # Atlas (se/quando for reativado)
  # mod_atlas_server("atlas", atlas_reactive = atlas, tcga_types = tcga_types)
  
  mod_search_your_target_server("search_your_target", search_your_target, Target, tcga_types)
  
  mod_all_signatures_server("all_signatures", Target)
  mod_signature_server("cnv_signature",         Target)
  mod_signature_server("gene_signature",        Target)
  mod_signature_server("methylation_signature", Target)
  mod_signature_server("mutation_signature",    Target)
  mod_signature_server("mirna_signature",       Target)
  mod_signature_server("protein_signature",     Target)
  mod_signature_server("transcript_signature",  Target)
  
  # mod_meaningful_interaction_server("Mean_int", meaningful_interaction)
  # mod_regulatory_circuitry_server("reg_circ",   regulatory_circuitry)
  
  mod_interaction_network_server(
    id                 = "interaction_net",
    meaningful_dataset = meaningful_interaction,   # Dataset_S3
    regulatory_dataset = regulatory_circuitry,     # Dataset_S4
    tcga_types         = tcga_types
  )
  
  mod_correlation_analysis_server("correlation_analysis", Target)
  mod_tumor_normal_analysis_server("tumor_normal_analysis", Target)
  mod_cox_analysis_server("cox_analysis", Target)
  mod_survival_analysis_server("survival_analysis", Target)
  mod_infiltrates_analysis_server("infiltrates_analysis", Target)
  
  # mod_web_resources_server("Web_resources")
  # mod_article_server("article")
  mod_developers_server("Developers")
  
  # ---------------------- Limpeza em fim de sessão ----------------------
  # NÃO usar rm(list = ls()) aqui: isso mexe no ambiente do worker inteiro.
  # Use apenas o que for realmente necessário (ex.: desconectar banco, fechar conexões, etc.).
  session$onSessionEnded(function() {
    message("Session ended.")
    # Se você criar objetos específicos da sessão (ex. conexões DB),
    # limpe-os aqui explicitamente.
    gc()
  })
}

