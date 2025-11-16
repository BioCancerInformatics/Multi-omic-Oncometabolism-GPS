# =========================================================
# Module: Meaningful Interactions (by Nomenclature)
# File: modules/module_meaningful_interaction.R
# =========================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(ggplot2); library(igraph); library(tidygraph); library(ggraph)
  library(htmltools); library(DT); library(readr); library(glue); library(grid)
})

# =========================== UI ==============================================
mod_meaningful_interaction_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Meaningful Interactions — Signature-centric"),
    fluidRow(
      column(
        width = 3,
        wellPanel(
          selectizeInput(
            ns("nomen_input"),
            "Signature Nomenclature:",
            choices = NULL, multiple = FALSE,
            options = list(placeholder = "Type to search…")
          ),
          actionButton(ns("run"), "Load", class = "btn-primary"),
          tags$hr(),
          downloadButton(ns("dl_row"),  "Download row (.csv)"),
          downloadButton(ns("dl_plot_png"), "Plot (.png)"),
          downloadButton(ns("dl_plot_pdf"), "Plot (.pdf)"),
          downloadButton(ns("dl_summary"),  "Summary (.txt)")
        )
      ),
      column(
        width = 9,
        wellPanel(
          tags$details(
            open = NA,
            tags$summary(tags$strong("Integrative summary")),
            uiOutput(ns("summaryText"))
          )
        ),
        wellPanel(
          h4("Signature regulatory circuitry"),
          plotOutput(ns("netplot"), height = "620px")
        ),
        wellPanel(
          h4("Selected row"),
          DTOutput(ns("row_table"))
        )
      )
    )
  )
}

# ======================= Helpers (safe, column-aware) =========================
# Split meaningful interactions list: separators + ; / |
.mi_split <- function(x) {
  if (is.null(x) || all(is.na(x))) return(character(0))
  x <- as.character(x)[1]
  if (!nzchar(x) || tolower(x) == "no meaningful interaction") return(character(0))
  x %>%
    str_replace_all("[(){}\\[\\]]", "") %>%
    str_split("\\s*[+;\\|/]\\s*") %>% pluck(1) %>%
    discard(~ .x == "" | is.na(.x)) %>% unique()
}

# Immune node mapping
.pick_immune <- function(label) {
  lab <- tolower(trimws(as.character(label)))
  if (!nzchar(lab) || is.na(lab)) return(NA_character_)
  if (lab %in% c("hot"))      return("Hot immune profile")
  if (lab %in% c("cold"))     return("Cold immune profile")
  if (lab %in% c("variable")) return("Variable immune profile")
  NA_character_
}

# Phenotype label from phenotype column + rho sign
.phenotype_from_rho <- function(phenotype_label, rho) {
  if (is.null(phenotype_label) || is.na(phenotype_label)) return(list(label = NA_character_, dir = NA_character_))
  rho_num <- suppressWarnings(as.numeric(rho))
  if (is.na(rho_num)) return(list(label = NA_character_, dir = NA_character_))
  dir <- ifelse(rho_num > 0, "pos", ifelse(rho_num < 0, "neg", NA_character_))
  plab <- tolower(trimws(as.character(phenotype_label)))
  if (!nzchar(plab)) return(list(label = NA_character_, dir = NA_character_))
  
  make_lab <- function(base, pos_txt, neg_txt) {
    if (rho_num > 0) glue("{pos_txt} {base}") else glue("{neg_txt} {base}")
  }
  
  if (str_detect(plab, "stemness"))
    return(list(label = make_lab("stemness", "Higher", "Lower"), dir = dir))
  if (str_detect(plab, "tumor mutational burden") || str_detect(plab, "\\bTMB\\b"))
    return(list(label = make_lab("TMB", "Higher", "Lower"), dir = dir))
  if (str_detect(plab, "microsatellite instability") || str_detect(plab, "\\bMSI\\b"))
    return(list(label = make_lab("MSI", "Higher", "Lower"), dir = dir))
  
  # default: do not fabricate unknown phenotypes
  list(label = NA_character_, dir = dir)
}

# Conform concordance flag to factor with expected levels
.norm_cd <- function(v) {
  v <- trimws(tolower(as.character(v)))
  v <- ifelse(v %in% c("convergent","divergent"), v, NA_character_)
  factor(v, levels = c("convergent","divergent"))
}

# Convert risk/protective keywords to prognosis nodes
.prog_node_from_type <- function(type_str) {
  s <- tolower(trimws(as.character(type_str)))
  if (!nzchar(s)) return(NA_character_)
  if (s %in% c("protective","low","lower","favorable")) return("Favorable prognosis")
  if (s %in% c("risky","risk","high","higher","worse"))  return("Worse prognosis")
  NA_character_
}

# Prefer explicit *_concordance_* columns; otherwise compute from rhos (for phenotype)
.pick_concordance_safe <- function(pref, fallback = NA_character_) {
  cd <- .norm_cd(pref)
  if (!is.na(cd)) as.character(cd) else fallback
}

# Deduplicate edges
.dedupe_edges <- function(df) {
  df %>% distinct(from, to, type, conc_plot, .keep_all = TRUE)
}

# =========================== Plot constructor =================================
.build_network_from_row <- function(x) {
  stopifnot(nrow(x) == 1)
  
  # --- Canonical nodes ---
  sig_node <- as.character(x$Nomenclature)
  
  # Regulators
  regulators <- .mi_split(x$Meaningful_interaction)
  if (!length(regulators)) regulators <- "(no regulator)"
  
  # Biological processes nodes
  pathway_lab <- as.character(x$Pathways)
  mcd_lab <- as.character(x$Metabolic_cell_death)
  mcd_valid <- nzchar(mcd_lab) && !is.na(mcd_lab) && tolower(mcd_lab) != "unrelated"
  
  # Phenotype nodes (signature & interaction)
  ph_sig <- .phenotype_from_rho(x$Phenotypic_layer_signature, x$rho_signature)
  ph_int <- .phenotype_from_rho(x$Phenotypic_layer_interaction, x$rho_interaction)
  
  phen_cd <- .pick_concordance_safe(x$Phenotypic_concordance,
                                    if (!is.na(ph_sig$dir) && !is.na(ph_int$dir) && nzchar(ph_sig$dir) && nzchar(ph_int$dir)) {
                                      if (identical(ph_sig$dir, ph_int$dir)) "convergent" else "divergent"
                                    } else NA_character_)
  
  # Immune nodes + concordance
  imm_sig <- .pick_immune(x$Immune_classification_signature)
  imm_int <- .pick_immune(x$Immune_classification_interaction)
  imm_cd  <- .norm_cd(x$Immune_concordance) %>% as.character()
  
  # Prognosis per endpoint: prefer Cox_*_type_* with Cox_concordance_*
  # If all Cox types are NA/NS, fall back to *_worst_prognosis_group_* with Survival_concordance_*
  cox_types_sig <- c(OS = x$Cox_OS_type_signature, DSS = x$Cox_DSS_type_signature,
                     DFI = x$Cox_DFI_type_signature, PFI = x$Cox_PFI_type_signature)
  cox_types_int <- c(OS = x$Cox_OS_type_interaction, DSS = x$Cox_DSS_type_interaction,
                     DFI = x$Cox_DFI_type_interaction, PFI = x$Cox_PFI_type_interaction)
  cox_cd <- c(OS = x$Cox_concordance_OS, DSS = x$Cox_concordance_DSS,
              DFI = x$Cox_concordance_DFI, PFI = x$Cox_concordance_PFI)
  
  has_cox <- any(tolower(na.omit(cox_types_sig)) != "ns") || any(tolower(na.omit(cox_types_int)) != "ns")
  
  surv_grp_sig <- c(OS = x$OS_worst_prognosis_group_signature, DSS = x$DSS_worst_prognosis_group_signature,
                    DFI = x$DFI_worst_prognosis_group_signature, PFI = x$PFI_worst_prognosis_group_signature)
  surv_grp_int <- c(OS = x$OS_worst_prognosis_group_interaction, DSS = x$DSS_worst_prognosis_group_interaction,
                    DFI = x$DFI_worst_prognosis_group_interaction, PFI = x$PFI_worst_prognosis_group_interaction)
  surv_cd <- c(OS = x$Survival_concordance_OS, DSS = x$Survival_concordance_DSS,
               DFI = x$Survival_concordance_DFI, PFI = x$Survival_concordance_PFI)
  
  # ---------- Build edges ----------
  edges <- list()
  
  # (1) Regulation (regulator -> signature), dashed orange
  edges <- append(edges, lapply(regulators, function(rg)
    tibble(from = rg, to = sig_node, type = "Regulation", conc_plot = "regulation")
  ))
  
  # (2) Biological processes association (signature and regulators to pathways/MCD)
  if (nzchar(pathway_lab) && !is.na(pathway_lab)) {
    edges <- append(edges, list(tibble(from = sig_node, to = pathway_lab, type = "Association", conc_plot = "association")))
    edges <- append(edges, lapply(regulators, function(rg)
      tibble(from = rg, to = pathway_lab, type = "Association", conc_plot = "association")))
  }
  if (mcd_valid) {
    edges <- append(edges, list(tibble(from = sig_node, to = mcd_lab, type = "Association", conc_plot = "association")))
    edges <- append(edges, lapply(regulators, function(rg)
      tibble(from = rg, to = mcd_lab, type = "Association", conc_plot = "association")))
  }
  
  # (3) Phenotype (from rho) with concordance by sign equality or provided Phenotypic_concordance
  if (!is.na(ph_sig$label)) edges <- append(edges, list(
    tibble(from = sig_node, to = ph_sig$label, type = "Association", conc_plot = phen_cd)
  ))
  if (!is.na(ph_int$label)) {
    edges <- append(edges, lapply(regulators, function(rg)
      tibble(from = rg, to = ph_int$label, type = "Association", conc_plot = phen_cd)))
  }
  
  # (4) Immune (only when Immune_concordance ∈ {convergent, divergent})
  if (!is.na(imm_cd) && imm_cd %in% c("convergent","divergent")) {
    if (!is.na(imm_sig)) edges <- append(edges, list(
      tibble(from = sig_node, to = imm_sig, type = "Association", conc_plot = imm_cd)
    ))
    if (!is.na(imm_int)) {
      edges <- append(edges, lapply(regulators, function(rg)
        tibble(from = rg, to = imm_int, type = "Association", conc_plot = imm_cd)))
    }
  }
  
  # (5) Prognosis per endpoint
  add_prog_edges <- function(sig_type, int_type, cd, label) {
    # Map types -> nodes
    n_sig <- .prog_node_from_type(sig_type)
    n_int <- .prog_node_from_type(int_type)
    cd2   <- .pick_concordance_safe(cd, NA_character_)
    out <- list()
    if (!is.na(n_sig)) out <- append(out, list(tibble(from = sig_node, to = n_sig, type = "Association", conc_plot = cd2)))
    if (!is.na(n_int)) out <- append(out, lapply(regulators, function(rg)
      tibble(from = rg, to = n_int, type = "Association", conc_plot = cd2)))
    out
  }
  
  if (has_cox) {
    edges <- append(edges, add_prog_edges(cox_types_sig["OS"],  cox_types_int["OS"],  cox_cd["OS"],  "OS"))
    edges <- append(edges, add_prog_edges(cox_types_sig["DSS"], cox_types_int["DSS"], cox_cd["DSS"], "DSS"))
    edges <- append(edges, add_prog_edges(cox_types_sig["DFI"], cox_types_int["DFI"], cox_cd["DFI"], "DFI"))
    edges <- append(edges, add_prog_edges(cox_types_sig["PFI"], cox_types_int["PFI"], cox_cd["PFI"], "PFI"))
  } else {
    edges <- append(edges, add_prog_edges(surv_grp_sig["OS"],  surv_grp_int["OS"],  surv_cd["OS"],  "OS"))
    edges <- append(edges, add_prog_edges(surv_grp_sig["DSS"], surv_grp_int["DSS"], surv_cd["DSS"], "DSS"))
    edges <- append(edges, add_prog_edges(surv_grp_sig["DFI"], surv_grp_int["DFI"], surv_cd["DFI"], "DFI"))
    edges <- append(edges, add_prog_edges(surv_grp_sig["PFI"], surv_grp_int["PFI"], surv_cd["PFI"], "PFI"))
  }
  
  edges <- bind_rows(compact(edges)) %>%
    mutate(
      type = factor(type, levels = c("Regulation","Association")),
      conc_plot = factor(conc_plot, levels = c("convergent","divergent","association","regulation"))
    ) %>%
    .dedupe_edges()
  
  # ---------- Nodes with layers ----------
  layer_of <- function(name) {
    if (name %in% regulators)                                return("Interaction")
    if (identical(name, sig_node))                            return("Signature")
    if (grepl("immune profile", tolower(name)))               return("Immune phenotype")
    if (grepl("^higher|^lower", tolower(name)))               return("Tumor phenotype")
    if (name %in% c("Favorable prognosis","Worse prognosis")) return("Prognosis")
    # default: Biological processes (pathways + MCD)
    "Biological processes"
  }
  
  nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
    mutate(layer = vapply(name, layer_of, character(1)))
  
  # ---------- Aesthetics ----------
  layer_colors <- c(
    "Interaction"          = "#E67E22",
    "Signature"            = "#2ECC71",
    "Biological processes" = "#5DADE2",
    "Immune phenotype"     = "gold",
    "Tumor phenotype"      = "#8E44AD",
    "Prognosis"            = "#C0392B"
  )
  edge_cols <- c(
    "convergent"  = "#27AE60",
    "divergent"   = "#C0392B",
    "association" = "grey60",
    "regulation"  = "#E67E22"
  )
  etype_lty <- c("Regulation" = "dashed", "Association" = "solid")
  etype_w   <- c("Regulation" = 1.2,      "Association" = 1.2)
  
  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE)
  
  set.seed(123)
  p <- ggraph(g, layout = "circle") +
    geom_edge_arc(
      aes(edge_colour = conc_plot, edge_linetype = type, edge_width = type),
      arrow = grid::arrow(length = unit(2.2, "mm"), type = "closed"),
      strength = 0.8,
      edge_alpha = 0.9
    ) +
    geom_node_point(aes(fill = layer), size = 7, shape = 21, color = "white") +
    geom_node_text(aes(label = name), vjust = -1.5, size = 3, fontface = "bold") +
    scale_fill_manual(values = layer_colors, name = "Dimension") +
    scale_edge_colour_manual(values = edge_cols, name = "Interpretation") +
    scale_edge_width_manual(values = etype_w, guide = "none") +
    scale_edge_linetype_manual(values = etype_lty, name = "Relation type") +
    ggtitle("Regulatory Circuitry — Signature ↔ Interaction") +
    theme_void() +
    coord_equal(clip = "off") +
    theme(
      legend.position = "left",
      legend.box.margin = margin(r = 30, l = 30, t = 20, b = 20),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1.5),
      plot.title.position = "panel",
      plot.margin = margin(t = 20, r = 40, b = 20, l = 55)
    )
  
  list(plot = p, nodes = nodes, edges = edges)
}

# =========================== Summary builder ==================================
.build_summary_html <- function(row, tcga_types = NULL) {
  # Defensive access
  pick <- function(nm) if (nm %in% names(row)) row[[nm]][1] else NA
  
  cancer_abbr <- pick("CTAB")
  cancer_full <- cancer_abbr
  if (!is.null(tcga_types) && all(c("Cancer_abbreviation","Cancer_names") %in% names(tcga_types))) {
    hit <- tcga_types$Cancer_names[match(cancer_abbr, tcga_types$Cancer_abbreviation)]
    if (!is.na(hit) && nzchar(hit)) cancer_full <- hit
  }
  
  # Core fields
  nom   <- pick("Nomenclature")
  sig   <- pick("Signatures")
  met   <- pick("Metabolism")
  path  <- pick("Pathways")
  mcd   <- pick("Metabolic_cell_death")
  mcd_txt <- if (!is.na(mcd) && nzchar(mcd) && tolower(mcd) != "unrelated")
    glue("<p>This signature is associated with <b>{mcd}</b>.</p>") else ""
  
  # Interactions
  inter <- pick("Meaningful_interaction")
  inter_sentence <- if (!is.na(inter) && nzchar(inter) && tolower(inter) != "no meaningful interaction")
    glue("<p>Clinically meaningful interaction(s): <b>{inter}</b>.</p>") else ""
  
  # RHO / Phenotypes (both sides)
  rho_sig <- suppressWarnings(as.numeric(pick("rho_signature")))
  rho_int <- suppressWarnings(as.numeric(pick("rho_interaction")))
  phen_sig <- pick("Phenotypic_layer_signature")
  phen_int <- pick("Phenotypic_layer_interaction")
  omic_sig <- pick("Omic_layer_signature")
  omic_int <- pick("Omic_layer_interaction")
  
  fmt_rho <- function(r) if (!is.na(r)) sprintf("%.3f", r) else "NA"
  
  corr_sig <- if (!is.na(phen_sig) && nzchar(phen_sig))
    glue("<p>Signature side: correlation between <b>{omic_sig}</b> and <b>{phen_sig}</b> (ρ = {fmt_rho(rho_sig)}).</p>")
  else ""
  
  corr_int <- if (!is.na(phen_int) && nzchar(phen_int))
    glue("<p>Interaction side: correlation between <b>{omic_int}</b> and <b>{phen_int}</b> (ρ = {fmt_rho(rho_int)}).</p>")
  else ""
  
  # Concordances (phenotype / immune / cox / survival)
  phen_cd <- pick("Phenotypic_concordance") %>% as.character() %>% {ifelse(. %in% c("convergent","divergent"), ., NA)}
  phen_line <- if (!is.na(phen_cd)) glue("<p>Phenotypic concordance: <b>{phen_cd}</b>.</p>") else ""
  
  imm_sig <- pick("Immune_classification_signature")
  imm_int <- pick("Immune_classification_interaction")
  imm_cd  <- pick("Immune_concordance") %>% as.character()
  imm_lines <- if (!is.na(imm_cd) && imm_cd %in% c("convergent","divergent")) {
    glue("<p>Immune: signature=<b>{imm_sig}</b>, interaction=<b>{imm_int}</b>, concordance=<b>{imm_cd}</b>.</p>")
  } else ""
  
  # Cox
  cox_block <- {
    items <- list()
    add <- function(lbl, sig_t, int_t, cd) {
      if (!is.na(sig_t) && tolower(sig_t) == "ns" &&
          !is.na(int_t) && tolower(int_t) == "ns") return(NULL)
      cd_txt <- if (!is.na(cd) && cd %in% c("convergent","divergent")) glue(" (concordance: {cd})") else ""
      items <<- append(items, glue("{lbl}: signature={sig_t}, interaction={int_t}{cd_txt}"))
      NULL
    }
    add("OS",  pick("Cox_OS_type_signature"),  pick("Cox_OS_type_interaction"),  pick("Cox_concordance_OS"))
    add("DSS", pick("Cox_DSS_type_signature"), pick("Cox_DSS_type_interaction"), pick("Cox_concordance_DSS"))
    add("DFI", pick("Cox_DFI_type_signature"), pick("Cox_DFI_type_interaction"), pick("Cox_concordance_DFI"))
    add("PFI", pick("Cox_PFI_type_signature"), pick("Cox_PFI_type_interaction"), pick("Cox_concordance_PFI"))
    if (length(items)) glue("<p>Cox regression: <b>{paste(items, collapse = '; ')}</b>.</p>") else ""
  }
  
  # Survival (fallback/also show)
  surv_block <- {
    items <- list()
    add <- function(lbl, sig_t, int_t, cd) {
      if (is.na(sig_t) && is.na(int_t)) return(NULL)
      cd_txt <- if (!is.na(cd) && cd %in% c("convergent","divergent")) glue(" (concordance: {cd})") else ""
      items <<- append(items, glue("{lbl}: signature={sig_t}, interaction={int_t}{cd_txt}"))
      NULL
    }
    add("OS",  pick("OS_worst_prognosis_group_signature"),  pick("OS_worst_prognosis_group_interaction"),  pick("Survival_concordance_OS"))
    add("DSS", pick("DSS_worst_prognosis_group_signature"), pick("DSS_worst_prognosis_group_interaction"), pick("Survival_concordance_DSS"))
    add("DFI", pick("DFI_worst_prognosis_group_signature"), pick("DFI_worst_prognosis_group_interaction"), pick("Survival_concordance_DFI"))
    add("PFI", pick("PFI_worst_prognosis_group_signature"), pick("PFI_worst_prognosis_group_interaction"), pick("Survival_concordance_PFI"))
    if (length(items)) glue("<p>Survival analysis: <b>{paste(items, collapse = '; ')}</b>.</p>") else ""
  }
  
  intro <- glue(
    "<p>The signature <b>{nom}</b> (components: <b>{sig}</b>) was identified in <b>{cancer_full}</b>, ",
    "involved in <b>{met}</b> and <b>{path}</b>.</p>"
  )
  
  HTML(paste0(
    intro,
    inter_sentence,
    mcd_txt,
    corr_sig,
    corr_int,
    phen_line,
    imm_lines,
    cox_block,
    surv_block,
    if ("Final_concordance_summary" %in% names(row) && nzchar(as.character(row$Final_concordance_summary)))
      glue("<p>Overall concordance: <b>{row$Final_concordance_summary}</b>.</p>") else "",
    "<p><i>Use the plot for topology; export the row for full statistics.</i></p>"
  ))
}

# ============================= SERVER =========================================
mod_meaningful_interaction_server <- function(id, dataset, tcga_types = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Populate choices from dataset
    observe({
      req(dataset())
      df <- dataset()
      validate(need("Nomenclature" %in% names(df), "Dataset must contain 'Nomenclature'."))
      ch <- sort(unique(na.omit(trimws(as.character(df$Nomenclature)))))
      updateSelectizeInput(session, "nomen_input", choices = ch, server = TRUE)
    })
    
    # Pick the row for the chosen Nomenclature (use first if duplicates)
    chosen_row <- eventReactive(input$run, {
      df <- dataset(); req(df, input$nomen_input)
      hit <- df %>% filter(.data$Nomenclature == input$nomen_input)
      validate(need(nrow(hit) > 0, "Nomenclature not found."))
      hit[1, , drop = FALSE]
    }, ignoreInit = TRUE)
    
    # Summary
    output$summaryText <- renderUI({
      row <- chosen_row(); req(row)
      .build_summary_html(row, tcga_types = tcga_types)
    })
    
    # Plot
    net_obj <- reactive({
      row <- chosen_row(); req(row)
      .build_network_from_row(row)
    })
    
    output$netplot <- renderPlot({
      obj <- net_obj(); req(obj$plot)
      obj$plot
    }, res = 120)
    
    # Row table
    output$row_table <- renderDT({
      row <- chosen_row(); req(row)
      datatable(row, rownames = FALSE, options = list(scrollX = TRUE, pageLength = 5))
    })
    
    # Downloads
    output$dl_row <- downloadHandler(
      filename = function() paste0("meaningful_interaction_row_", gsub("[^A-Za-z0-9]+","_", input$nomen_input %||% "signature"), ".csv"),
      content = function(file) readr::write_csv(chosen_row(), file)
    )
    
    output$dl_plot_png <- downloadHandler(
      filename = function() paste0("meaningful_interaction_plot_", gsub("[^A-Za-z0-9]+","_", input$nomen_input %||% "signature"), ".png"),
      content = function(file) {
        obj <- net_obj(); req(obj$plot)
        ggsave(file, obj$plot, width = 14, height = 8, dpi = 300, bg = "white")
      }
    )
    
    output$dl_plot_pdf <- downloadHandler(
      filename = function() paste0("meaningful_interaction_plot_", gsub("[^A-Za-z0-9]+","_", input$nomen_input %||% "signature"), ".pdf"),
      content = function(file) {
        obj <- net_obj(); req(obj$plot)
        ggsave(file, obj$plot, width = 14, height = 8, device = cairo_pdf, bg = "white")
      }
    )
    
    output$dl_summary <- downloadHandler(
      filename = function() paste0("meaningful_interaction_summary_", gsub("[^A-Za-z0-9]+","_", input$nomen_input %||% "signature"), ".txt"),
      content = function(file) {
        html <- as.character(isolate(.build_summary_html(chosen_row(), tcga_types = tcga_types)))
        txt  <- gsub("<[^>]+>", "", html)
        writeLines(txt, file, useBytes = TRUE)
      }
    )
  })
}
