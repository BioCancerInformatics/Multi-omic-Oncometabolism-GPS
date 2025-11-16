# =========================================================
# Module: Custom Signature Builder
# File: modules/module_custom_signatures.R
# =========================================================

mod_custom_signatures_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Custom Signature Builder"),
    fluidRow(
      column(
        width = 3,
        wellPanel(
          tags$h5("1) Base columns (fixed)"),
          div(
            style = "background:#f7f7f7; border:1px solid #ddd; border-radius:6px; padding:8px; margin-bottom:8px;",
            tags$b("Mandatory grouping:"),
            tags$ul(
              tags$li("Cancer_types"),
              tags$li("Metabolism")
            )
          ),
          tags$h5("2) Optional grouping columns"),
          uiOutput(ns("optional_cols_ui")),
          tags$hr(),
          numericInput(
            ns("rho_threshold"), 
            label = HTML("Absolute |Spearman ρ| threshold"),
            value = 0.10, min = 0, step = 0.01
          ),
          checkboxInput(
            ns("split_by_sign"),
            label = "Separate positive vs. negative correlations",
            value = TRUE
          ),
          actionButton(ns("build"), "Build Signatures", class = "btn-primary")
        )
      ),
      column(
        width = 9,
        wellPanel(
          fluidRow(
            column(
              6,
              tags$h5("Summary"),
              verbatimTextOutput(ns("summary_text"))
            ),
            column(
              6,
              div(
                style = "text-align:right; margin-top:24px;",
                downloadButton(ns("download_csv"), "Download CSV")
              )
            )
          ),
          DT::DTOutput(ns("signatures_table")) %>% shinycssloaders::withSpinner(type = 4)
        )
      )
    )
  )
}


mod_custom_signatures_server <- function(id, dataset_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # ---------------------- Helpers (adapted from your code) ----------------------
    format_gene_signature <- function(Target) {
      # Keeps gene symbols intact, wraps whole block in parentheses and joins with " + "
      formatted_Target <- vapply(Target, function(gene) {
        if (grepl("-", gene)) gene else gene
      }, FUN.VALUE = character(1))
      signature <- paste(formatted_Target, collapse = " + ")
      paste0("(", signature, ")")
    }
    
    find_common_interactions <- function(Interactions) {
      interaction_list <- strsplit(Interactions, " / ", fixed = TRUE)
      
      if (length(interaction_list) == 1) return(Interactions[1])
      
      # identical sets?
      if (all(vapply(interaction_list, function(x) identical(x, interaction_list[[1]]), logical(1)))) {
        return(Interactions[1])
      }
      
      common <- Reduce(intersect, interaction_list)
      if (length(common) > 0) paste(common, collapse = " / ") else "No common interactions"
    }
    
    # ---------------------- Available columns UI ----------------------
    base_cols <- c("Cancer_types", "Metabolism")
    
    observe({
      df <- dataset_reactive()
      req(df, nrow(df) > 0)
      
      # Candidate optional columns = all columns minus base_cols and the known "measure" columns
      known_measures <- c(
        "Target", "Interactions", "Correlation_rho", "Correlation_type",
        "Members", "Multiomics_Signature", "Common_interaction"
      )
      
      candidates <- setdiff(colnames(df), unique(c(base_cols, known_measures)))
      candidates <- sort(candidates)
      
      # Pre-select a few sensible defaults if present
      suggested <- intersect(
        c("Molecular_class", "Pathways", "Metabolic_cell_death", "Tumor_vs_normal",
          "Omic_layer", "Phenotypic_layer", "Cox_OS_type", "Cox_DSS_type",
          "Cox_DFI_type", "Cox_PFI_type", "OS_worst_prognosis_group",
          "DSS_worst_prognosis_group", "DFI_worst_prognosis_group",
          "PFI_worst_prognosis_group", "Immune_classification"),
        candidates
      )
      
      output$optional_cols_ui <- renderUI({
        selectizeInput(
          ns("optional_cols"),
          label = "Select optional grouping columns (order matters):",
          choices = candidates,
          selected = suggested,
          multiple = TRUE,
          options = list(plugins = list("remove_button", "drag_drop"), placeholder = "Choose columns…")
        )
      })
    })
    
    # ---------------------- Build signatures ----------------------
    make_signatures <- eventReactive(input$build, {
      df <- dataset_reactive()
      req(df, nrow(df) > 0)
      req(all(base_cols %in% colnames(df)))
      validate(need("Target" %in% colnames(df), "'Target' column is required."))
      validate(need("Interactions" %in% colnames(df), "'Interactions' column is required."))
      validate(need("Correlation_rho" %in% colnames(df), "'Correlation_rho' column is required."))
      
      rho_thr <- as.numeric(input$rho_threshold)
      validate(need(!is.na(rho_thr) && rho_thr >= 0, "Provide a valid non-negative ρ threshold."))
      
      # Filter by |rho|
      df_f <- dplyr::filter(df, abs(.data$Correlation_rho) >= rho_thr)
      
      # Build grouping columns
      opt_cols <- input$optional_cols %||% character(0)
      grouping_cols <- unique(c(base_cols, opt_cols))
      
      # If splitting by sign, add Correlation_type and process separately
      split_sign <- isTRUE(input$split_by_sign)
      
      aggregate_once <- function(d) {
        # Defensive: keep only columns we actually use to avoid accidental type surprises
        keep_cols <- unique(c(grouping_cols, "Target", "Interactions"))
        keep_cols <- keep_cols[keep_cols %in% colnames(d)]
        d <- d[, keep_cols, drop = FALSE]
        
        d %>%
          group_by(across(all_of(grouping_cols))) %>%
          summarise(
            Target       = format_gene_signature(.data$Target),
            Interactions = find_common_interactions(.data$Interactions),
            .groups = "drop"
          )
      }
      
      if (split_sign) {
        positive <- df_f %>% filter(.data$Correlation_rho > 0)
        negative <- df_f %>% filter(.data$Correlation_rho < 0)
        
        pos_tbl <- if (nrow(positive)) aggregate_once(positive) %>% mutate(Correlation_type = "positive") else NULL
        neg_tbl <- if (nrow(negative)) aggregate_once(negative) %>% mutate(Correlation_type = "negative") else NULL
        
        aggregated <- bind_rows(pos_tbl, neg_tbl)
      } else {
        aggregated <- aggregate_once(df_f)
      }
      
      validate(need(!is.null(aggregated) && nrow(aggregated) > 0, "No signatures formed with current settings."))
      
      # Members: number of items inside the parentheses separated by " + "
      members_count <- function(sig) {
        core <- gsub("[()]", "", sig)
        if (nzchar(core)) length(strsplit(core, " \\+ ", fixed = FALSE)[[1]]) else 0
      }
      
      aggregated <- aggregated %>%
        mutate(Members = vapply(.data$Target, members_count, integer(1), USE.NAMES = FALSE))
      
      # Arrange by cancer then members desc
      if ("Cancer_types" %in% colnames(aggregated)) {
        aggregated <- aggregated %>% arrange(.data$Cancer_types, desc(.data$Members))
      } else {
        aggregated <- aggregated %>% arrange(desc(.data$Members))
      }
      
      # Rename + column order
      aggregated <- aggregated %>%
        rename(
          Multiomics_Signature = .data$Target,
          Common_interaction   = .data$Interactions
        )
      
      # Desired order: key metadata first
      front_cols <- c(
        "Members", "Multiomics_Signature", "Common_interaction",
        base_cols, opt_cols
      )
      if (split_sign) front_cols <- c(front_cols, "Correlation_type")
      
      # Only keep existing columns and reorder
      front_cols <- front_cols[front_cols %in% colnames(aggregated)]
      other_cols <- setdiff(colnames(aggregated), front_cols)
      aggregated <- aggregated[, c(front_cols, other_cols), drop = FALSE]
      
      aggregated
    }, ignoreInit = TRUE)
    
    # ---------------------- Outputs ----------------------
    output$summary_text <- renderText({
      req(make_signatures())
      x <- make_signatures()
      paste0(
        "Rows: ", nrow(x),
        " | Columns: ", ncol(x),
        " | Grouping: ",
        paste(c("Cancer_types", "Metabolism", input$optional_cols %||% character(0)), collapse = " + "),
        if (isTRUE(input$split_by_sign)) " | Split: by correlation sign" else ""
      )
    })
    
    output$signatures_table <- DT::renderDT({
      req(make_signatures())
      DT::datatable(
        make_signatures(),
        rownames = FALSE,
        filter = "top",
        extensions = c("Buttons"),
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = "Bfrtip",
          buttons = c("copy", "csv", "excel")
        )
      )
    })
    
    output$download_csv <- downloadHandler(
      filename = function() {
        paste0("custom_signatures_", format(Sys.time(), "%Y%m%d-%H%M%S"), ".csv")
      },
      content = function(file) {
        req(make_signatures())
        readr::write_csv(make_signatures(), file)
      }
    )
  })
}
