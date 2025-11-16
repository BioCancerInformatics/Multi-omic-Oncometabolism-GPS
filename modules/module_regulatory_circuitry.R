# =========================================================
# Module: Regulatory Circuitry Plotter (spec-compliant)
# File: modules/module_regulatory_circuitry.R
# =========================================================

# UI --------------------------------------------------------------------------
mod_regulatory_circuitry_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Regulatory Circuitry"),
    fluidRow(
      column(
        width = 3,
        wellPanel(
          textInput(ns("nomenclature"), "Nomenclature", placeholder = "Type/paste a Nomenclature..."),
          actionButton(ns("plot_btn"), "Plot circuitry", class = "btn-primary"),
          tags$hr(),
          checkboxGroupInput(
            ns("layers_on"),
            "Edge groups to include:",
            choices = c(
              "Biological processes"  = "bio",
              "Tumor phenotype"       = "phen",
              "Immune profile"        = "immune",
              "Prognosis association" = "prog"
            ),
            selected = c("bio","phen","immune","prog")
          ),
          selectInput(
            ns("layout"), "Layout",
            choices = c("circle", "stress", "kk", "fr"), selected = "circle"
          ),
          sliderInput(ns("text_size"), "Label size", min = 2, max = 6, value = 3, step = 0.5),
          sliderInput(ns("node_size"), "Node size", min = 4, max = 12, value = 7, step = 1),
          tags$hr(),
          downloadButton(ns("download_png"), "Download PNG")
        )
      ),
      column(
        width = 9,
        wellPanel(
          uiOutput(ns("case_label")),
          shinycssloaders::withSpinner(plotOutput(ns("plot"), height = 600), type = 4)
        ),
        verbatimTextOutput(ns("diag"))
      )
    )
  )
}


# SERVER ----------------------------------------------------------------------
mod_regulatory_circuitry_server <- function(id, df_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ---------------------- Helpers -------------------------------------------
    `%||%` <- function(a, b) if (!is.null(a)) a else b
    
    split_regulators <- function(x) {
      if (is.null(x) || all(is.na(x))) return(character(0))
      raw <- as.character(x[1])
      parts <- unlist(strsplit(raw, "\\s*[+;\\|/]\\s*"))
      unique(parts[nzchar(parts)])
    }
    
    # NÃO dividir Pathways_signature (rótulo único, preserva vírgulas)
    as_single_label <- function(x, fallback = NA_character_) {
      z <- as.character(x[1]); z <- stringr::str_squish(z)
      if (!nzchar(z)) fallback else z
    }
    
    norm_cd <- function(v) {
      v <- trimws(tolower(as.character(v)))
      v <- ifelse(v %in% c("convergent","divergent","ns"), v, NA_character_)
      factor(v, levels = c("convergent","divergent","ns"))
    }
    
    mk_edge <- function(from, to, type, concordance = NA_character_) {
      if (is.na(from) || is.na(to) || from == "" || to == "") return(NULL)
      if (!is.na(concordance) && identical(tolower(concordance), "ns")) return(NULL)
      tibble::tibble(from = from, to = to, type = type, concordance = concordance)
    }
    
    phenotype_from_rho <- function(phenotype_label, rho) {
      if (is.null(rho) || is.na(rho)) return(list(label = NA_character_, dir = NA_character_))
      rho_num <- suppressWarnings(as.numeric(rho))
      if (is.na(rho_num)) return(list(label = NA_character_, dir = NA_character_))
      plab <- tolower(trimws(as.character(phenotype_label)))
      dir <- ifelse(rho_num > 0, "pos", ifelse(rho_num < 0, "neg", NA_character_))
      
      if (stringr::str_detect(plab, "stemness")) {
        lab <- if (rho_num > 0) "Higher stemness" else if (rho_num < 0) "Lower stemness" else NA_character_
        return(list(label = lab, dir = dir))
      }
      if (stringr::str_detect(plab, "tumor mutational burden|\\bTMB\\b")) {
        lab <- if (rho_num > 0) "Higher TMB" else if (rho_num < 0) "Lower TMB" else NA_character_
        return(list(label = lab, dir = dir))
      }
      if (stringr::str_detect(plab, "microsatellite instability|\\bMSI\\b")) {
        lab <- if (rho_num > 0) "Higher MSI" else if (rho_num < 0) "Lower MSI" else NA_character_
        return(list(label = lab, dir = dir))
      }
      list(label = NA_character_, dir = dir)
    }
    
    pick_immune_node <- function(label) {
      lab <- tolower(trimws(as.character(label)))
      if (lab %in% c("hot"))      return("Hot immune profile")
      if (lab %in% c("cold"))     return("Cold immune profile")
      if (lab %in% c("variable")) return("Variable immune profile")
      NA_character_
    }
    
    map_cox_to_node <- function(v) {
      if (is.null(v) || all(is.na(v))) return(NA_character_)
      s <- tolower(trimws(as.character(v[1])))
      if (!nzchar(s)) return(NA_character_)
      if (grepl("protect|low", s))  return("Favorable prognosis")
      if (grepl("risk|high", s))    return("Worse prognosis")
      if (s %in% c("protective","low")) return("Favorable prognosis")
      if (s %in% c("risky","high"))     return("Worse prognosis")
      NA_character_
    }
    
    prog_node_from_keywords <- function(text) {
      if (is.null(text) || all(is.na(text))) return(NA_character_)
      tx <- tolower(paste(na.omit(as.character(text)), collapse = " | "))
      if (tx == "") return(NA_character_)
      if (stringr::str_detect(tx, "\\brisk\\b|\\bhigh\\b")) return("Worse prognosis")
      if (stringr::str_detect(tx, "\\bprotect\\b|\\blow\\b")) return("Favorable prognosis")
      NA_character_
    }
    
    collect_signature_text <- function(row) {
      cols <- grep("signature", names(row), ignore.case = TRUE, value = TRUE)
      cols <- union(cols, c("Signatures"))
      unlist(row[, intersect(cols, names(row)), drop = TRUE])
    }
    collect_regulator_text <- function(row) {
      cols <- grep("interaction", names(row), ignore.case = TRUE, value = TRUE)
      cols <- union(cols, c("Meaningful_interaction"))
      unlist(row[, intersect(cols, names(row)), drop = TRUE])
    }
    
    any_present <- function(x) any(!is.na(x) & nzchar(as.character(x)))
    
    # ---------------------- Core builder ---------------------------------------
    build_graph <- function(x, layers_on = c("bio","phen","immune","prog"),
                            layout_name = "circle", node_size = 7, text_size = 3) {
      
      validate(need(all(c("Nomenclature","Meaningful_interaction") %in% names(x)),
                    "Required columns are missing."))
      
      NODE_SIG <- as_single_label(x$Nomenclature)
      
      # Reguladores
      regulators <- split_regulators(x$Meaningful_interaction)
      if (!length(regulators)) regulators <- "(no regulator)"
      
      # Biological processes (Pathways + RCD) — rótulos únicos (sem split por vírgula)
      pathways_sig <- as_single_label(x$Pathways_signature, fallback = "Metabolic pathway")
      mcd_sig_raw  <- as_single_label(x$Metabolic_cell_death_signature, fallback = NA_character_)
      mcd_sig <- if (!is.na(mcd_sig_raw) && tolower(mcd_sig_raw) != "unrelated") mcd_sig_raw else NA_character_
      bio_nodes <- unique(na.omit(c(pathways_sig, mcd_sig)))
      
      # Fenótipo por ρ
      PH_sig <- phenotype_from_rho(x$Phenotypic_layer_signature %||% NA, x$rho_signature %||% NA)
      PH_ent <- phenotype_from_rho(x$Phenotypic_layer_signature %||% NA, x$rho_interaction %||% NA)
      phen_conc <- if (is.na(PH_sig$dir) || is.na(PH_ent$dir)) NA_character_
      else if (PH_sig$dir == PH_ent$dir) "convergent" else "divergent"
      
      # Imune
      IMM_sig_node <- pick_immune_node(x$Immune_classification_signature %||% NA)
      IMM_ent_node <- pick_immune_node(x$Immune_classification_interaction %||% NA)
      immune_cd    <- norm_cd(x$Immune_concordance %||% NA)
      
      # Prognóstico (Cox)
      cox_sig <- c(
        OS  = map_cox_to_node(x$Cox_OS_type_signature  %||% x$OS_worst_prognosis_group_signature  %||% NA),
        DSS = map_cox_to_node(x$Cox_DSS_type_signature %||% x$DSS_worst_prognosis_group_signature %||% NA),
        PFI = map_cox_to_node(x$Cox_PFI_type_signature %||% x$PFI_worst_prognosis_group_signature %||% NA),
        DFI = map_cox_to_node(x$Cox_DFI_type_signature %||% x$DFI_worst_prognosis_group_signature %||% NA)
      )
      cox_int <- c(
        OS  = map_cox_to_node(x$Cox_OS_type_interaction  %||% x$OS_worst_prognosis_group_interaction  %||% NA),
        DSS = map_cox_to_node(x$Cox_DSS_type_interaction %||% x$DSS_worst_prognosis_group_interaction %||% NA),
        PFI = map_cox_to_node(x$Cox_PFI_type_interaction %||% x$PFI_worst_prognosis_group_interaction %||% NA),
        DFI = map_cox_to_node(x$Cox_DFI_type_interaction %||% x$DFI_worst_prognosis_group_interaction %||% NA)
      )
      
      # Fallback por keywords se Cox ausente em ambos os lados
      if (!any_present(cox_sig) && !any_present(cox_int)) {
        sig_kw <- prog_node_from_keywords(collect_signature_text(x))
        int_kw <- prog_node_from_keywords(collect_regulator_text(x))
        cox_sig[] <- sig_kw
        cox_int[] <- int_kw
      }
      
      # ---------------------- Edges --------------------------------------------
      edges <- list()
      
      # 1) Regulation
      edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, NODE_SIG, "Regulation", NA_character_)))
      
      # 2) Biological processes — regra fixa: mesma via => "convergent"
      if ("bio" %in% layers_on && length(bio_nodes)) {
        # assinatura -> bio (concordance registrada como 'divergent', porém cor será fixa por "sig_assoc")
        edges <- append(edges, lapply(bio_nodes, function(bn) mk_edge(NODE_SIG, bn, "Association", "convergent")))
        # reguladores -> bio (divergent)
        for (bn in bio_nodes) {
          edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, bn, "Association", "convergent")))
        }
      }
      
      # 3) Fenótipo (ρ)
      if ("phen" %in% layers_on) {
        if (!is.na(PH_sig$label)) edges <- append(edges, list(mk_edge(NODE_SIG, PH_sig$label, "Association", phen_conc)))
        if (!is.na(PH_ent$label)) edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, PH_ent$label, "Association", phen_conc)))
      }
      
      # 4) Imune (somente se Immune_concordance != NS)
      if ("immune" %in% layers_on && !is.na(immune_cd) && !identical(as.character(immune_cd), "ns")) {
        if (!is.na(IMM_sig_node)) edges <- append(edges, list(mk_edge(NODE_SIG, IMM_sig_node, "Association", as.character(immune_cd))))
        if (!is.na(IMM_ent_node)) edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, IMM_ent_node, "Association", as.character(immune_cd))))
      }
      
      # 5) Prognóstico (Cox) por endpoint
      if ("prog" %in% layers_on) {
        for (ep in names(cox_sig)) {
          sig_node <- cox_sig[[ep]]
          int_node <- cox_int[[ep]]
          if (!is.na(sig_node)) edges <- append(edges, list(mk_edge(NODE_SIG, sig_node, "Association",
                                                                    if (!is.na(int_node) && int_node == sig_node) "convergent"
                                                                    else if (!is.na(int_node) && int_node != sig_node) "divergent" else NA_character_)))
          if (!is.na(int_node)) edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, int_node, "Association",
                                                                                               if (!is.na(sig_node) && sig_node == int_node) "convergent"
                                                                                               else if (!is.na(sig_node) && sig_node != int_node) "divergent" else NA_character_)))
        }
      }
      
      edges <- dplyr::bind_rows(purrr::compact(edges))
      validate(need(nrow(edges) > 0, "No edges to plot for this Nomenclature and current layer selection."))
      edges <- edges %>% dplyr::distinct(from, to, type, concordance, .keep_all = TRUE)
      
      # --------- Mapeamento visual: assinatura fixa, interação por concordance
      edges_for_plot <- edges %>%
        dplyr::mutate(
          conc_plot = dplyr::case_when(
            .data$type == "Regulation"                 ~ "regulation",
            .data$type == "Association" & .data$from == NODE_SIG ~ "sig_assoc",    # cor fixa p/ arestas da assinatura
            !is.na(.data$concordance)                  ~ as.character(.data$concordance),  # interação: verde/vermelho
            TRUE                                       ~ "association"                      # neutra (se sobrar NA)
          )
        )
      
      # ---------------------- Nodes & Layers -----------------------------------
      layer_of <- function(name) {
        if (name %in% regulators)                                 return("Interaction")
        if (name %in% c(NODE_SIG))                                return("Signature")
        if (name %in% bio_nodes)                                  return("Biological processes")
        if (grepl("immune profile", tolower(name)))               return("Immune phenotype")
        if (grepl("stemness|tmb|msi", tolower(name)))             return("Tumor phenotype")
        if (name %in% c("Favorable prognosis","Worse prognosis")) return("Prognosis")
        "Biological processes"
      }
      
      nodes <- tibble::tibble(name = unique(c(edges$from, edges$to))) %>%
        dplyr::mutate(layer = vapply(name, layer_of, character(1)))
      
      # ---------------------- Aesthetics ---------------------------------------
      layer_colors <- c(
        "Interaction"          = "#E67E22",
        "Signature"            = "#2ECC71",
        "Biological processes" = "#5DADE2",
        "Immune phenotype"     = "gold",
        "Tumor phenotype"      = "#8E44AD",
        "Prognosis"            = "#C0392B"
      )
      
      edge_colours <- c(
        "sig_assoc"  = "#34495E",  # cor fixa p/ associações da assinatura
        "convergent" = "#27AE60",
        "divergent"  = "#C0392B",
        "association"= "#7f8c8d",  # neutra (quase nunca usada agora)
        "regulation" = "#E67E22"
      )
      
      etype_lty <- c("Regulation" = "dashed", "Association" = "solid")
      etype_w   <- c("Regulation" = 1.2,      "Association" = 1.2)
      
      g <- tidygraph::tbl_graph(nodes = nodes, edges = edges_for_plot, directed = TRUE)
      
      set.seed(123)
      p <- ggraph::ggraph(g, layout = layout_name) +
        ggraph::geom_edge_arc(
          ggplot2::aes(edge_colour = conc_plot, edge_linetype = type, edge_width = type),
          arrow      = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
          strength   = 0.8,
          edge_alpha = 0.9
        ) +
        ggraph::geom_node_point(ggplot2::aes(fill = layer), size = node_size, shape = 21, color = "white") +
        ggraph::geom_node_text(ggplot2::aes(label = name), vjust = -1.5, size = text_size, fontface = "bold") +
        ggplot2::scale_fill_manual(values = layer_colors, name = "Dimension") +
        ggraph::scale_edge_colour_manual(
          values = edge_colours,
          name   = "Interpretation",
          breaks = c("sig_assoc","convergent","divergent")  # oculta "regulation" e "association" da legenda
        ) +
        ggraph::scale_edge_width_manual(values = etype_w, guide = "none") +
        ggraph::scale_edge_linetype_manual(values = etype_lty, name = "Relation type") +
        ggplot2::theme_void() +
        ggplot2::coord_equal(clip = "off") +
        ggplot2::theme(
          legend.position   = "left",
          legend.box.margin = ggplot2::margin(r = 30, l = 30, t = 20, b = 20),
          plot.title        = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1.5),
          plot.title.position = "panel",
          plot.margin       = ggplot2::margin(t = 20, r = 40, b = 20, l = 55)
        )
      
      case_label <- x$Signatures %||% x$Nomenclature %||% "Selected circuitry"
      list(plot = p, label = case_label, nodes = nodes, edges = edges_for_plot)
    }
    
    # ---------------------- Seleção por Nomenclature ---------------------------
    picked_row <- eventReactive(input$plot_btn, {
      df <- df_reactive(); req(df, nrow(df) > 0)
      validate(need("Nomenclature" %in% names(df), "'Nomenclature' column is missing in the dataset."))
      
      nm <- trimws(input$nomenclature)
      validate(need(nzchar(nm), "Please type a Nomenclature."))
      
      hits <- dplyr::filter(df, .data$Nomenclature == nm)
      validate(need(nrow(hits) > 0, sprintf("No row found for Nomenclature = '%s'.", nm)))
      if (nrow(hits) > 1) showNotification("Multiple rows match this Nomenclature — using the first.", type = "warning", duration = 6)
      hits[1, , drop = FALSE]
    }, ignoreInit = TRUE)
    
    built <- reactive({
      x <- picked_row(); req(x)
      tryCatch({
        build_graph(
          x = x,
          layers_on  = input$layers_on,
          layout_name = input$layout,
          node_size   = input$node_size,
          text_size   = input$text_size
        )
      }, error = function(e) list(error = TRUE, message = conditionMessage(e)))
    })
    
    output$case_label <- renderUI({
      b <- built(); req(b)
      if (isTRUE(b$error)) return(tags$p(tags$b("Error:"), b$message, style = "color:#b30000;"))
      tags$h4(glue::glue("Regulatory Circuitry — {b$label}"))
    })
    
    output$plot <- renderPlot({
      b <- built(); req(b); validate(need(!isTRUE(b$error), b$message)); b$plot
    })
    
    output$diag <- renderPrint({
      b <- built(); req(b)
      if (isTRUE(b$error)) return(b$message)
      list(
        nodes = dplyr::as_tibble(b$nodes),
        edges = dplyr::as_tibble(b$edges) %>% head(25)
      )
    })

    # ------------------------ Downloads -----------------------------------------
    fname_safe <- reactive({
      b <- built(); req(b)
      base <- if (isTRUE(b$error)) "Regulatory_Circuitry" else paste0("Regulatory_Circuitry_", gsub("[^A-Za-z0-9]+","_", b$label))
      paste0(base)
    })
    
    output$download_png <- downloadHandler(
      filename = function() paste0(fname_safe(), ".png"),
      content = function(file) {
        b <- built(); validate(need(!isTRUE(b$error), b$message))
        ggplot2::ggsave(file, plot = b$plot, width = 16, height = 9, units = "in", dpi = 300, bg = "white")
      }
    )
  })
}

