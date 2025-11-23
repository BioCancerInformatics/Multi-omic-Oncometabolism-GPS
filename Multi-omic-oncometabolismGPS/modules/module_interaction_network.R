mod_interaction_network_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        h4("Select signature (Nomenclature)"),
        
        # Escolha de fonte de dados
        selectInput(
          ns("mode"),
          "Data source:",
          choices = c(
            "Meaningful interactions" = "meaningful",
            "Regulatory circuitry"    = "regulatory"
          ),
          selected = "meaningful"
        ),
        
        # Nomenclature (vai ser populado a partir do dataset ativo)
        selectizeInput(
          ns("nomen_input"),
          "Nomenclature:",
          choices = NULL,
          multiple = FALSE,
          options = list(placeholder = "Type to search…")
        ),
        
        actionButton(ns("run"), "Search", class = "btn-primary"),
        hr(),
        
        downloadButton(ns("downloadRow"),     "Download selected row (.csv)"),
        br(), br(),
        downloadButton(ns("downloadAll"),     "Download all rows (.csv)"),
        br(), br(),
        downloadButton(ns("downloadSummary"), "Download summary (.txt)"),
        br(), br(),
        downloadButton(ns("downloadPlot"),    "Download network plot (.pdf)")
      ),
      
      mainPanel(
        h3("Integrative summary for the selected interaction"),
        uiOutput(ns("summaryText")),
        hr(),
        h3("Regulatory circuitry network"),
        plotOutput(ns("netplot"), height = "550px"),
        hr(),
        h3("All interactions for this Nomenclature"),
        p("Click a row to update the textual summary above."),
        DT::DTOutput(ns("row_table"))
      )
    )
  )
}



`%||%` <- function(a, b) if (!is.null(a)) a else b

coalesce_col <- function(row, candidates) {
  for (cn in candidates) {
    if (cn %in% names(row)) {
      v <- row[[cn]]
      if (!all(is.na(v))) {
        if (is.character(v)) {
          if (any(stringr::str_squish(v) != "")) return(v)
        } else {
          return(v)
        }
      }
    }
  }
  NA
}

split_multi <- function(x, delim = "\\s*[+;/|]\\s*") {
  if (is.null(x) || all(is.na(x))) return(character(0))
  x %>%
    as.character() %>%
    stringr::str_replace_all("[(){}\\[\\]]", "") %>%
    stringr::str_replace_all("[\u200B-\u200D\uFEFF]", "") %>%
    stringr::str_split(delim) %>%
    purrr::pluck(1) %>%
    purrr::map_chr(~ stringr::str_squish(.x)) %>%
    purrr::discard(~ .x == "" | is.na(.x)) %>%
    unique()
}

pick_immune_node <- function(label) {
  lab <- tolower(trimws(as.character(label)))
  if (lab %in% "hot")      return("Hot immune profile")
  if (lab %in% "cold")     return("Cold immune profile")
  if (lab %in% "variable") return("Variable immune profile")
  NA_character_
}

phenotype_from_rho <- function(phenotype_label, rho) {
  if (is.null(rho) || length(rho) == 0L) {
    return(list(label = NA_character_, dir = NA_character_))
  }
  if (is.null(phenotype_label) || length(phenotype_label) == 0L) {
    return(list(label = NA_character_, dir = NA_character_))
  }
  rho_num <- suppressWarnings(as.numeric(rho[1]))
  if (is.na(rho_num)) {
    return(list(label = NA_character_, dir = NA_character_))
  }
  plab <- tolower(trimws(as.character(phenotype_label[1])))
  if (plab == "" || is.na(plab)) {
    return(list(label = NA_character_, dir = NA_character_))
  }
  dir <- ifelse(rho_num > 0, "pos",
                ifelse(rho_num < 0, "neg", NA_character_))
  if (stringr::str_detect(plab, "stemness"))
    return(list(
      label = if (rho_num > 0) "Lower stemness" else "Higher stemness",
      dir   = dir
    ))
  if (stringr::str_detect(plab, "tumor mutational burden"))
    return(list(
      label = if (rho_num > 0) "Higher TMB" else "Lower TMB",
      dir   = dir
    ))
  if (stringr::str_detect(plab, "microsatellite instability"))
    return(list(
      label = if (rho_num > 0) "Higher MSI" else "Lower MSI",
      dir   = dir
    ))
  list(label = NA_character_, dir = dir)
}

norm_cd <- function(v) {
  v <- trimws(tolower(as.character(v)))
  v <- ifelse(v %in% c("convergent","divergent","ns"), v, NA)
  factor(v, levels = c("convergent","divergent","ns"))
}

mk_edge <- function(from, to, type, concordance = NA_character_) {
  from <- stringr::str_squish(as.character(from))
  to   <- stringr::str_squish(as.character(to))
  if (is.na(from) || is.na(to) || from == "" || to == "") return(NULL)
  if (!is.na(concordance) && identical(concordance, "ns")) return(NULL)
  tibble::tibble(from = from, to = to, type = type, concordance = concordance)
}

map_side_type_to_node <- function(x) {
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("protective","low")) return("Favorable prognosis")
  if (x %in% c("risky","high"))     return("Worse prognosis")
  NA_character_
}

pair_to_conc <- function(sig_type, int_type) {
  if (any(is.na(c(sig_type, int_type)))) return(NA_character_)
  a <- tolower(sig_type); b <- tolower(int_type)
  if ((a %in% c("protective","low")) && (b %in% c("protective","low"))) return("convergent")
  if ((a %in% c("risky","high"))      && (b %in% c("risky","high")))     return("convergent")
  if ((a %in% c("protective","low")) && (b %in% c("risky","high")))      return("divergent")
  if ((a %in% c("risky","high"))      && (b %in% c("protective","low"))) return("divergent")
  NA_character_
}

clean_side <- function(x){
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("ns", "", NA)) return(NA_character_)
  x
}










##############






add_prog_edges_for_endpoint <- function(edges_list, sig_type, int_type,
                                        ASSOC, NODE_SIG, regulators, row){
  sig_type <- clean_side(sig_type)
  int_type <- clean_side(int_type)
  
  sig_node <- map_side_type_to_node(sig_type)
  int_node <- map_side_type_to_node(int_type)
  conc_ep  <- pair_to_conc(sig_type, int_type)
  
  if (!is.na(sig_node) && !is.na(conc_ep))
    edges_list <- append(edges_list, list(mk_edge(NODE_SIG, sig_node, ASSOC, conc_ep)))
  
  if (!is.na(int_node) && !is.na(conc_ep) && length(regulators))
    edges_list <- append(
      edges_list,
      lapply(regulators, function(rg) mk_edge(rg, int_node, ASSOC, conc_ep))
    )
  
  edges_list
}

build_edges_from_row <- function(x_row, NODE_SIG){
  # Immune
  IMM_sig_node <- pick_immune_node(
    coalesce_col(x_row, c("Immune_classification_signature","Immune_classification_sig"))
  )
  IMM_ent_node <- pick_immune_node(
    coalesce_col(x_row, c("Immune_classification_interaction","Immune_classification_int"))
  )
  
  # Pathways + MCD
  paths_sig <- split_multi(
    coalesce_col(x_row, c("Pathways_signature","Pathways")),
    delim = "\\s*[+;/|]\\s*"
  )
  mcd_sig   <- split_multi(
    coalesce_col(x_row, c("Metabolic_cell_death_signature","Metabolic_cell_death")),
    delim = "\\s*[+;/|]\\s*"
  )
  mcd_sig   <- mcd_sig[!tolower(mcd_sig) %in% "unrelated"]
  all_sig_paths <- unique(c(paths_sig, mcd_sig))
  
  # Fenótipo (RHO) – aceita ambos esquemas
  PH_sig <- phenotype_from_rho(
    coalesce_col(x_row, c("Phenotypic_layer_signature","Phenotypic_layer_sig")),
    coalesce_col(x_row, c("rho_signature","Correlation_rho_sig"))
  )
  PH_ent <- phenotype_from_rho(
    coalesce_col(x_row, c("Phenotypic_layer_interaction","Phenotypic_layer_int")),
    coalesce_col(x_row, c("rho_interaction","Correlation_rho_int"))
  )
  phen_conc <- if (!is.na(PH_sig$dir) && !is.na(PH_ent$dir)) {
    if (PH_sig$dir == PH_ent$dir) "convergent" else "divergent"
  } else NA_character_
  
  # Reguladores
  regulators <- split_multi(coalesce_col(x_row, c("Meaningful_interaction", "Interaction", "Nomenclature_int")))
  if (!length(regulators)) regulators <- "(no regulator)"
  
  ASSOC <- "Signature and Interaction association"
  edges_row <- list()
  
  # (A) Reguladores -> assinatura
  edges_row <- append(
    edges_row,
    lapply(regulators, function(rg) mk_edge(rg, NODE_SIG, "Regulation", NA_character_))
  )
  
  # (B) Assinatura/Reguladores -> pathways / MCD
  if (length(all_sig_paths)) {
    edges_row <- append(
      edges_row,
      lapply(all_sig_paths, function(pth) mk_edge(NODE_SIG, pth, ASSOC, "convergent"))
    )
    if (length(regulators)) {
      edges_row <- append(
        edges_row,
        purrr::flatten(
          lapply(all_sig_paths, function(pth)
            lapply(regulators, function(rg) mk_edge(rg, pth, ASSOC, "convergent"))
          )
        )
      )
    }
  }
  
  # (D) Fenótipo
  if (!is.na(PH_sig$label))
    edges_row <- append(edges_row, list(mk_edge(NODE_SIG, PH_sig$label, ASSOC, phen_conc)))
  if (!is.na(PH_ent$label))
    edges_row <- append(edges_row, lapply(regulators, function(rg) mk_edge(rg, PH_ent$label, ASSOC, phen_conc)))
  
  # (E) Imune
  if (!is.na(IMM_sig_node))
    edges_row <- append(
      edges_row,
      list(mk_edge(
        NODE_SIG, IMM_sig_node, ASSOC,
        as.character(norm_cd(x_row$Immune_concordance))
      ))
    )
  if (!is.na(IMM_ent_node))
    edges_row <- append(
      edges_row,
      lapply(
        regulators,
        function(rg) mk_edge(
          rg, IMM_ent_node, ASSOC,
          as.character(norm_cd(x_row$Immune_concordance))
        )
      )
    )
  
  # (F) Cox lado-específico – aceita variantes *_signature / *_sig, *_interaction / *_int
  edges_row <- add_prog_edges_for_endpoint(
    edges_row,
    sig_type = coalesce_col(x_row, c("Cox_OS_type_signature","Cox_OS_type_sig")),
    int_type = coalesce_col(x_row, c("Cox_OS_type_interaction","Cox_OS_type_int")),
    ASSOC    = ASSOC,
    NODE_SIG = NODE_SIG,
    regulators = regulators,
    row = x_row
  )
  edges_row <- add_prog_edges_for_endpoint(
    edges_row,
    sig_type = coalesce_col(x_row, c("Cox_DSS_type_signature","Cox_DSS_type_sig")),
    int_type = coalesce_col(x_row, c("Cox_DSS_type_interaction","Cox_DSS_type_int")),
    ASSOC    = ASSOC,
    NODE_SIG = NODE_SIG,
    regulators = regulators,
    row = x_row
  )
  edges_row <- add_prog_edges_for_endpoint(
    edges_row,
    sig_type = coalesce_col(x_row, c("Cox_PFI_type_signature","Cox_PFI_type_sig")),
    int_type = coalesce_col(x_row, c("Cox_PFI_type_interaction","Cox_PFI_type_int")),
    ASSOC    = ASSOC,
    NODE_SIG = NODE_SIG,
    regulators = regulators,
    row = x_row
  )
  edges_row <- add_prog_edges_for_endpoint(
    edges_row,
    sig_type = coalesce_col(x_row, c("Cox_DFI_type_signature","Cox_DFI_type_sig")),
    int_type = coalesce_col(x_row, c("Cox_DFI_type_interaction","Cox_DFI_type_int")),
    ASSOC    = ASSOC,
    NODE_SIG = NODE_SIG,
    regulators = regulators,
    row = x_row
  )
  
  dplyr::bind_rows(purrr::compact(edges_row))
}

build_network_from_matches <- function(matches){
  stopifnot(nrow(matches) >= 1)
  stopifnot("Signatures" %in% names(matches))
  
  NODE_SIG <- unique(matches$Signatures); stopifnot(length(NODE_SIG) == 1L)
  NODE_SIG <- NODE_SIG[[1]]
  
  edges <- purrr::map_dfr(
    seq_len(nrow(matches)),
    function(i) build_edges_from_row(matches[i, , drop = FALSE], NODE_SIG = NODE_SIG)
  ) %>%
    dplyr::filter(!is.na(to), !is.na(from)) %>%
    dplyr::mutate(
      type = stringr::str_squish(as.character(type)),
      type = ifelse(type == "Regulation", "Regulation", "Signature and Interaction association"),
      type = factor(type, levels = c("Regulation","Signature and Interaction association")),
      concordance = factor(concordance, levels = c("convergent","divergent"))
    ) %>%
    dplyr::distinct(from, to, type, concordance, .keep_all = TRUE)
  
  layer_colors <- c(
    "Interaction"          = "#D55E00",
    "Signature"            = "#009E73",
    "Biological processes" = "#0072B2",
    "Immune phenotype"     = "#D0AD00",
    "Tumor phenotype"      = "#CC79A7",
    "Prognosis"            = "#000000"
  )
  
  edge_colours <- c(
    "sig_assoc"   = "#000000",
    "convergent"  = "#009E73",
    "divergent"   = "#C0392B",
    "association" = "#A9A9A9",
    "regulation"  = "#D55E00"
  )
  
  edges_for_plot <- edges %>%
    dplyr::mutate(
      conc_plot = dplyr::case_when(
        as.character(type) == "Regulation" ~ "regulation",
        from == NODE_SIG                   ~ "sig_assoc",
        is.na(concordance)                 ~ "association",
        TRUE                               ~ as.character(concordance)
      )
    )
  
  regs_all <- unique(unlist(purrr::map(
    seq_len(nrow(matches)),
    function(i) split_multi(matches$Meaningful_interaction[i] %||% matches$Interaction[i])
  )))
  
  layer_of <- function(name){
    if (isTRUE(name %in% regs_all)) return("Interaction")
    if (name %in% c(NODE_SIG)) return("Signature")
    if (grepl("stemness|tmb|msi", tolower(name))) return("Tumor phenotype")
    if (grepl("immune profile|microenvironment", tolower(name))) return("Immune phenotype")
    if (name %in% c("Favorable prognosis","Worse prognosis")) return("Prognosis")
    "Biological processes"
  }
  
  nodes <- tibble::tibble(name = unique(c(edges$from, edges$to))) %>%
    dplyr::mutate(layer = vapply(name, layer_of, character(1))) %>%
    dplyr::mutate(
      pt_size = dplyr::case_when(
        layer %in% c("Interaction", "Signature") ~ 12,
        TRUE                                     ~ 8
      )
    )
  
  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges_for_plot, directed = TRUE)
  
  etype_lty <- c(
    "Regulation"                            = "dashed",
    "Signature and Interaction association" = "solid"
  )
  etype_w   <- c(
    "Regulation"                            = 1.0,
    "Signature and Interaction association" = 1.0
  )
  
  set.seed(123)
  p <- ggraph::ggraph(g, layout = "circle") +
    ggraph::geom_edge_arc(
      ggplot2::aes(edge_colour = conc_plot, edge_linetype = type, edge_width = type),
      arrow      = grid::arrow(length = grid::unit(2.0, "mm"), type = "closed"),
      strength   = 0.3,
      edge_alpha = 0.5,
      start_cap  = ggraph::circle(6, "mm"),
      end_cap    = ggraph::circle(6, "mm")
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill = layer, size = pt_size),
      shape = 21,
      color = "white"
    ) +
    ggplot2::scale_size_identity() +
    ggraph::geom_node_text(
      data = function(d) d %>% dplyr::filter(layer %in% c("Signature","Interaction")),
      ggplot2::aes(label = name, color = layer),
      nudge_y  = 0.15,
      size = 3.0,
      fontface = "bold"
    ) +
    ggraph::geom_node_text(
      data = function(d) d %>%
        dplyr::filter(layer %in% c("Biological processes","Tumor phenotype","Immune phenotype","Prognosis")),
      ggplot2::aes(label = name, color = layer),
      nudge_y  = -0.15,
      size = 3.0,
      fontface = "bold"
    ) +
    ggplot2::scale_color_manual(values = layer_colors, guide = "none") +
    ggplot2::scale_fill_manual(
      values = layer_colors,
      name = "Dimension",
      guide = ggplot2::guide_legend(override.aes = list(size = 7))
    ) +
    ggraph::scale_edge_colour_manual(
      values = edge_colours,
      name   = "Interpretation",
      breaks = c("sig_assoc","convergent","divergent"),
      labels = c("Signature association","Convergent Sig-Int association","Divergent Sig-Int association")
    ) +
    ggraph::scale_edge_width_manual(values = etype_w, guide = "none") +
    ggraph::scale_edge_linetype_manual(values = etype_lty, name = "Relation type") +
    ggplot2::theme_void() +
    ggplot2::coord_equal(expand = FALSE, clip = "off") +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 10, face = "bold"),
      legend.text  = ggplot2::element_text(size = 8),
      legend.position = "left",
      legend.justification = c(-0.8, 0.0),
      legend.box.margin = ggplot2::margin(r=10, l=10, t=10, b=10),
      plot.margin = ggplot2::margin(t=20, r=20, b=20, l=20)
    )
  
  list(plot = p, nodes = nodes, edges = edges_for_plot, NODE_SIG = NODE_SIG)
}







#####################







build_summary_html <- function(row, tcga_types = NULL){
  get1 <- function(r, nm) if (nm %in% names(r)) r[[nm]][1] else NA
  
  nomen <- coalesce_col(row, c("Nomenclature","Nomenclature_sig"))
  sig_lab  <- get1(row, "Signatures")
  ctab_val <- get1(row, "CTAB")
  metabolism <- get1(row, "Metabolism")
  pathways   <- get1(row, "Pathways")
  mcd        <- get1(row, "Metabolic_cell_death")
  inter      <- coalesce_col(row, c("Meaningful_interaction","Interaction"))
  
  cancer_full <- ctab_val
  if (!is.null(tcga_types) &&
      all(c("Cancer_abbreviation","Cancer_names") %in% names(tcga_types))) {
    hit <- tcga_types$Cancer_names[match(ctab_val, tcga_types$Cancer_abbreviation)]
    if (!is.na(hit) && nzchar(hit)) cancer_full <- hit
  }
  
  omic_sig   <- coalesce_col(row, c("Omic_layer_signature","Omic_layer_sig"))
  phen_sig   <- coalesce_col(row, c("Phenotypic_layer_signature","Phenotypic_layer_sig"))
  rho_sig    <- suppressWarnings(as.numeric(
    coalesce_col(row, c("rho_signature","Correlation_rho_sig"))
  ))
  omic_int   <- coalesce_col(row, c("Omic_layer_interaction","Omic_layer_int"))
  phen_int   <- coalesce_col(row, c("Phenotypic_layer_interaction","Phenotypic_layer_int"))
  rho_int    <- suppressWarnings(as.numeric(
    coalesce_col(row, c("rho_interaction","Correlation_rho_int"))
  ))
  phen_conc  <- get1(row, "Phenotypic_concordance")
  
  fmt_rho <- function(x) if (!is.na(x)) sprintf("%.3f", x) else "NA"
  
  cox_conc_OS  <- get1(row, "Cox_concordance_OS")
  cox_conc_DSS <- get1(row, "Cox_concordance_DSS")
  cox_conc_DFI <- get1(row, "Cox_concordance_DFI")
  cox_conc_PFI <- get1(row, "Cox_concordance_PFI")
  cox_agg      <- get1(row, "Cox_concordance_aggregated")
  
  surv_conc_OS  <- get1(row, "Survival_concordance_OS")
  surv_conc_DSS <- get1(row, "Survival_concordance_DSS")
  surv_conc_DFI <- get1(row, "Survival_concordance_DFI")
  surv_conc_PFI <- get1(row, "Survival_concordance_PFI")
  surv_agg      <- get1(row, "Survival_concordance_aggregated")
  
  immune_sig  <- coalesce_col(row, c("Immune_classification_signature","Immune_classification_sig"))
  immune_int  <- coalesce_col(row, c("Immune_classification_interaction","Immune_classification_int"))
  immune_conc <- get1(row, "Immune_concordance")
  
  final_conc <- get1(row, "Final_concordance_summary")
  
  intro <- paste0(
    "<p>The nomenclature <b>", nomen, "</b> corresponds to the signature <b>",
    sig_lab, "</b> identified in <b>", cancer_full, "</b>. ",
    "This signature is mapped to <b>", metabolism, "</b> and the pathway <b>",
    pathways, "</b>."
  )
  if (!is.na(mcd) && tolower(mcd) != "unrelated" && nzchar(mcd)) {
    intro <- paste0(
      intro, " It is also associated with the metabolic cell death category <b>",
      mcd, "</b>."
    )
  }
  intro <- paste0(intro, "</p>")
  
  inter_txt <- ""
  if (!is.na(inter) && nzchar(inter) && inter != "No meaningful interaction") {
    inter_txt <- paste0(
      "<p>This signature establishes clinically meaningful interactions with <b>",
      inter, "</b>.</p>"
    )
  }
  
  phen_txt <- ""
  if (!is.na(omic_sig) && !is.na(phen_sig) && !is.na(rho_sig)) {
    phen_txt <- paste0(
      "<p>At the signature level, <b>", omic_sig, "</b> is correlated with <b>",
      phen_sig, "</b> (ρ = ", fmt_rho(rho_sig), ")."
    )
    if (!is.na(omic_int) && !is.na(phen_int) && !is.na(rho_int)) {
      phen_txt <- paste0(
        phen_txt, " The interaction layer (", omic_int, ") is associated with <b>",
        phen_int, "</b> (ρ = ", fmt_rho(rho_int), ")."
      )
    }
    if (!is.na(phen_conc) && phen_conc %in% c("convergent","divergent")) {
      phen_txt <- paste0(
        phen_txt, " Overall, phenotype-level behaviour is <b>",
        phen_conc, "</b> between signature and interaction."
      )
    }
    phen_txt <- paste0(phen_txt, "</p>")
  }
  
  collapse_non_ns <- function(...) {
    x <- c(...)
    x <- x[!is.na(x) & x != "" & x != "NS"]
    if (!length(x)) return(NULL)
    paste(x, collapse = ", ")
  }
  
  cox_vec <- collapse_non_ns(
    if (!is.na(cox_conc_OS))  paste0("OS: ",  cox_conc_OS)  else NULL,
    if (!is.na(cox_conc_DSS)) paste0("DSS: ", cox_conc_DSS) else NULL,
    if (!is.na(cox_conc_DFI)) paste0("DFI: ", cox_conc_DFI) else NULL,
    if (!is.na(cox_conc_PFI)) paste0("PFI: ", cox_conc_PFI) else NULL
  )
  if (!is.null(cox_vec)) {
    cox_txt <- paste0(
      "<p>Cox-based prognosis concordance across clinical endpoints: <b>",
      cox_vec, "</b>."
    )
    if (!is.na(cox_agg) && cox_agg != "" && cox_agg != "NS") {
      cox_txt <- paste0(
        cox_txt, " Aggregated concordance across endpoints is <b>",
        cox_agg, "</b>."
      )
    }
    cox_txt <- paste0(cox_txt, "</p>")
  } else {
    cox_txt <- "<p>Cox-based concordance was not classified across endpoints.</p>"
  }
  
  surv_vec <- collapse_non_ns(
    if (!is.na(surv_conc_OS))  paste0("OS: ",  surv_conc_OS)  else NULL,
    if (!is.na(surv_conc_DSS)) paste0("DSS: ", surv_conc_DSS) else NULL,
    if (!is.na(surv_conc_DFI)) paste0("DFI: ", surv_conc_DFI) else NULL,
    if (!is.na(surv_conc_PFI)) paste0("PFI: ", surv_conc_PFI) else NULL
  )
  if (!is.null(surv_vec)) {
    surv_txt <- paste0(
      "<p>Survival-group concordance across endpoints is described as <b>",
      surv_vec, "</b>."
    )
    if (!is.na(surv_agg) && surv_agg != "" && surv_agg != "NS") {
      surv_txt <- paste0(
        surv_txt, " The aggregated survival concordance is <b>",
        surv_agg, "</b>."
      )
    }
    surv_txt <- paste0(surv_txt, "</p>")
  } else {
    surv_txt <- "<p>Survival-group concordance metrics were not classified.</p>"
  }
  
  immune_txt <- ""
  if ((!is.na(immune_sig) && immune_sig != "" && immune_sig != "NS") ||
      (!is.na(immune_int) && immune_int != "" && immune_int != "NS")) {
    immune_txt <- "<p>At the immune level,"
    if (!is.na(immune_sig) && immune_sig != "" && immune_sig != "NS") {
      immune_txt <- paste0(
        immune_txt, " the signature is classified as <b>", immune_sig, "</b>"
      )
    }
    if (!is.na(immune_int) && immune_int != "" && immune_int != "NS") {
      if (!is.na(immune_sig) && immune_sig != "" && immune_sig != "NS") {
        immune_txt <- paste0(immune_txt, ", while the interaction is classified as <b>",
                             immune_int, "</b>")
      } else {
        immune_txt <- paste0(
          immune_txt, " the interaction is classified as <b>", immune_int, "</b>"
        )
      }
    }
    if (!is.na(immune_conc) && immune_conc != "" && immune_conc != "NS") {
      immune_txt <- paste0(
        immune_txt, ". The immune concordance between both layers is <b>",
        immune_conc, "</b>"
      )
    }
    immune_txt <- paste0(immune_txt, ".</p>")
  }
  
  final_txt <- ""
  if (!is.na(final_conc) && final_conc != "" && final_conc != "NS") {
    final_txt <- paste0(
      "<p>Integrating metabolic, phenotypic, immune and prognostic dimensions, ",
      "the overall regulatory circuitry is summarized as <b>", final_conc,
      "</b>.</p>"
    )
  }
  
  HTML(paste0(
    intro,
    inter_txt,
    phen_txt,
    cox_txt,
    surv_txt,
    immune_txt,
    final_txt,
    "<p><i>Use the regulatory network below to explore how the signature and its interaction partners",
    " connect metabolic pathways, phenotypes, immune context and prognosis.</i></p>"
  ))
}


######################





mod_interaction_network_server <- function(id,
                                           meaningful_dataset,   # reactiveVal / reactive (Dataset S3)
                                           regulatory_dataset,   # reactiveVal / reactive (Dataset S4)
                                           tcga_types = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # -------- 1) Escolher data frame ativo conforme modo --------
    active_df <- reactive({
      if (input$mode == "meaningful") {
        req(meaningful_dataset())
        meaningful_dataset()
      } else {
        req(regulatory_dataset())
        regulatory_dataset()
      }
    })
    
    # Nome da coluna de Nomenclature em cada dataset
    get_nomen_col <- function(df) {
      if ("Nomenclature" %in% names(df))     return("Nomenclature")
      if ("Nomenclature_sig" %in% names(df)) return("Nomenclature_sig")
      stop("Dataset must contain a 'Nomenclature' or 'Nomenclature_sig' column.")
    }
    
    # -------- 2) Atualizar choices de Nomenclature quando muda modo ou dataset --------
    observe({
      df <- active_df(); req(df)
      ncol <- get_nomen_col(df)
      ch <- sort(unique(na.omit(as.character(df[[ncol]]))))
      updateSelectizeInput(
        session, "nomen_input",
        choices  = ch,
        selected = NULL,
        server   = TRUE
      )
    })
    
    # -------- 3) Filtrar todas as linhas da Nomenclature escolhida --------
    matches_all <- eventReactive(input$run, {
      df <- active_df(); req(df, input$nomen_input)
      ncol <- get_nomen_col(df)
      hit <- df %>% dplyr::filter(.data[[ncol]] == input$nomen_input)
      validate(need(nrow(hit) > 0, "Nomenclature not found in dataset."))
      
      if ("Signatures" %in% names(hit)) {
        stopifnot(length(unique(hit$Signatures)) == 1L)
      }
      hit
    }, ignoreInit = TRUE)
    
    # -------- 4) Tabela de interações --------
    output$row_table <- DT::renderDT({
      rows <- matches_all(); req(rows)
      DT::datatable(
        rows,
        rownames  = FALSE,
        selection = "single",
        options   = list(scrollX = TRUE, pageLength = 5)
      )
    })
    
    # -------- 5) Linha selecionada (para resumo e download) --------
    selected_row <- reactive({
      rows <- matches_all(); req(rows)
      sel  <- input$row_table_rows_selected
      if (!is.null(sel) && length(sel) == 1) {
        rows[sel, , drop = FALSE]
      } else {
        rows[1, , drop = FALSE]
      }
    })
    
    # -------- 6) Resumo HTML --------
    output$summaryText <- renderUI({
      row <- selected_row(); req(row)
      build_summary_html(row, tcga_types = tcga_types)
    })
    
    # -------- 7) Grafo (todas as interações daquela Nomenclature) --------
    net_obj <- reactive({
      rows <- matches_all(); req(rows)
      build_network_from_matches(rows)
    })
    
    output$netplot <- renderPlot({
      obj <- net_obj(); req(obj$plot)
      obj$plot
    }, res = 120)
    
    # -------- 8) Download do plot --------
    output$downloadPlot <- downloadHandler(
      filename = function() {
        nm <- input$nomen_input %||% "signature"
        mode_tag <- if (input$mode == "meaningful") "S3" else "S4"
        paste0("interaction_network_",
               gsub("[^A-Za-z0-9]+","_", nm), "_", mode_tag, ".pdf")
      },
      content = function(file) {
        obj <- net_obj(); req(obj$plot)
        ggplot2::ggsave(
          filename = file,
          plot     = obj$plot,
          device   = "pdf",
          width    = 10,
          height   = 8
        )
      }
    )
    
    # -------- 9) Download da linha selecionada --------
    output$downloadRow <- downloadHandler(
      filename = function() {
        nm <- input$nomen_input %||% "signature"
        mode_tag <- if (input$mode == "meaningful") "S3" else "S4"
        paste0("interaction_row_",
               gsub("[^A-Za-z0-9]+","_", nm), "_", mode_tag, ".csv")
      },
      content = function(file) {
        df <- selected_row(); req(df)
        readr::write_csv(df, file)
      }
    )
    
    # -------- 10) Download de todas as linhas daquela Nomenclature --------
    output$downloadAll <- downloadHandler(
      filename = function() {
        nm <- input$nomen_input %||% "signature"
        mode_tag <- if (input$mode == "meaningful") "S3" else "S4"
        paste0("interaction_all_",
               gsub("[^A-Za-z0-9]+","_", nm), "_", mode_tag, ".csv")
      },
      content = function(file) {
        df <- matches_all(); req(df)
        readr::write_csv(df, file)
      }
    )
    
    # -------- 11) Download do resumo em texto --------
    output$downloadSummary <- downloadHandler(
      filename = function() {
        nm <- input$nomen_input %||% "signature"
        mode_tag <- if (input$mode == "meaningful") "S3" else "S4"
        paste0("interaction_summary_",
               gsub("[^A-Za-z0-9]+","_", nm), "_", mode_tag, ".txt")
      },
      content = function(file) {
        row  <- selected_row(); req(row)
        html <- as.character(build_summary_html(row, tcga_types = tcga_types))
        txt  <- gsub("<[^>]+>", "", html)
        writeLines(txt, file, useBytes = TRUE)
      }
    )
  })
}
