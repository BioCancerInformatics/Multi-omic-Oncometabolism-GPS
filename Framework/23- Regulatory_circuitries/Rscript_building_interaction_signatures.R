library(dplyr)
library(stringr)
library(rio)

df001 <- import("/Oncometabolism_GPS/22- Nomenclature/Regulatory_circuitries_with_nomenclature.rds")
df002 <- import("/Oncometabolism_GPS/22- Nomenclature/All_signature_results.rds")

# --- Improved function to format gene signatures ---
format_gene_signature <- function(Target) {
  # Apply backticks only to genes with '-'
  formatted_Target <- sapply(Target, function(gene) {
    if (grepl("-", gene)) paste0("`", gene, "`") else gene
  })
  
  # If there's only one gene, return it as is (no parentheses)
  if (length(formatted_Target) == 1) {
    return(formatted_Target)
  } else {
    # Combine and wrap in parentheses for multi-gene signatures
    signature <- paste(formatted_Target, collapse = " + ")
    return(paste0("(", signature, ")"))
  }
}

# --- Aggregation function ---
aggregate_Target <- function(df001) {
  aggregate_and_format <- function(df) {
    df %>%
      group_by(
        Nomenclature, CTAB, Metabolism_signature, Pathways_signature, 
        Metabolic_cell_death_signature, Molecular_class_interaction, 
        Omic_layer_interaction, Phenotypic_layer_interaction, Phenotypic_concordance, Phenotype_block,
        Cox_OS_type_interaction, Cox_concordance_OS, Cox_DSS_type_interaction, Cox_concordance_DSS, 
        Cox_concordance_PFI, Cox_DFI_type_interaction, Cox_concordance_DFI, Cox_PFI_type_interaction,
        Cox_metrics_aggregated, Cox_concordance_aggregated, Cox_block, 
        OS_worst_prognosis_group_interaction, Survival_concordance_OS, DSS_worst_prognosis_group_interaction, 
        Survival_concordance_DSS, DFI_worst_prognosis_group_interaction, Survival_concordance_DFI, 
        PFI_worst_prognosis_group_interaction, Survival_concordance_PFI, Survival_metrics_aggregated,
        Survival_concordance_aggregated, Survival_block,
        Immune_classification_interaction, Immune_concordance, Immune_block, Final_concordance_summary 
      ) %>%
      summarise(
        Meaningful_interaction = format_gene_signature(Meaningful_interaction),
        .groups = "drop"
      )
  }
  
  positive <- df001 %>% filter(rho_interaction > 0)
  negative <- df001 %>% filter(rho_interaction < 0)
  
  positive_aggregated <- aggregate_and_format(positive) %>%
    mutate(Correlation_type = "positive")
  
  negative_aggregated <- aggregate_and_format(negative) %>%
    mutate(Correlation_type = "negative")
  
  aggregated_data <- bind_rows(positive_aggregated, negative_aggregated)
  return(aggregated_data)
}

# --- Apply the function to create new signatures from df001 ---
df003_new_signatures <- aggregate_Target(df001)

# --- Step 2: Compare with df002 to annotate matches and retrieve Nomenclature ---
df004_matched_signatures <- df003_new_signatures %>%
  left_join(
    df002 %>%
      select(Nomenclature, CTAB, Metabolism, Pathways, Metabolic_cell_death,
             Omic_layer, Phenotypic_layer, Signatures),
    by = c(
      "CTAB" = "CTAB",
      "Metabolism_signature" = "Metabolism",
      "Pathways_signature" = "Pathways",
      "Metabolic_cell_death_signature" = "Metabolic_cell_death",
      "Omic_layer_interaction" = "Omic_layer",
      "Phenotypic_layer_interaction" = "Phenotypic_layer",
      "Meaningful_interaction" = "Signatures"
    )
  )

df004_matched_signatures <- df004_matched_signatures %>%
  filter(!is.na(Nomenclature.y))

# --- Step 3 ---
df005_matched_signatures <- df004_matched_signatures %>%
  left_join(
    df002 %>%
      select(Nomenclature, Signatures, Omic_layer, Phenotypic_layer, Molecular_class,
             Tumor_vs_normal, Tumor_vs_normal_p.adj, Correlation_rho, Correlation_p.adj,
             Combined_outcome_HRC, Cox_OS_type, Cox_OS_p.value, Cox_DSS_type, Cox_DSS_p.value,
             Cox_DFI_type, Cox_DFI_p.value, Cox_PFI_type, Cox_PFI_p.value, Combined_outcome_SMC,
             OS_worst_prognosis_group, OS_p.value, DSS_worst_prognosis_group, DSS_p.value,
             DFI_worst_prognosis_group, DFI_p.value, PFI_worst_prognosis_group, PFI_p.value,
             Microenvironment_classification, Immune_classification
             ),
    by = c(
      "Nomenclature.x" = "Nomenclature"
    )
  )

# Rename columns and reorder them
df006_renamed <- df005_matched_signatures %>%
  rename(
    Nomenclature_sig = Nomenclature.x,
    Signatures = Signatures,
    
    Nomenclature_int = Nomenclature.y,
    Interaction = Meaningful_interaction,
    
    CTAB = CTAB,
    Metabolism = Metabolism_signature,
    Pathways = Pathways_signature,
    Metabolic_cell_death = Metabolic_cell_death_signature,
    
    Molecular_class_sig = Molecular_class,   
    Molecular_class_int = Molecular_class_interaction,
    
    Omic_layer_sig = Omic_layer,
    Phenotypic_layer_sig = Phenotypic_layer,
    
    Omic_layer_int = Omic_layer_interaction,
    Phenotypic_layer_int = Phenotypic_layer_interaction,
    
    Correlation_rho_sig = Correlation_rho,
    Correlation_p.adj_sig = Correlation_p.adj,
    
    Phenotypic_concordance = Phenotypic_concordance,
    Phenotype_block = Phenotype_block, 
    
    Combined_outcome_HRC_sig = Combined_outcome_HRC,
    
    Cox_OS_type_sig = Cox_OS_type,
    Cox_OS_p.value_sig = Cox_OS_p.value,
    
    Cox_OS_type_int = Cox_OS_type_interaction,
    Cox_concordance_OS = Cox_concordance_OS,
    
    Cox_DSS_type_sig = Cox_DSS_type, 
    Cox_DSS_p.value_sig = Cox_DSS_p.value,
    
    Cox_DSS_type_int = Cox_DSS_type_interaction,
    Cox_concordance_DSS = Cox_concordance_DSS,
    
    Cox_DFI_type_sig = Cox_DFI_type, 
    Cox_DFI_p.value_sig = Cox_DFI_p.value,
    
    Cox_DFI_type_int = Cox_DFI_type_interaction,
    Cox_concordance_DFI = Cox_concordance_DFI,
    
    Cox_PFI_type_sig = Cox_PFI_type, 
    Cox_PFI_p.value_sig = Cox_PFI_p.value,
    
    Cox_PFI_type_int = Cox_PFI_type_interaction,
    Cox_concordance_PFI = Cox_concordance_PFI,
    
    Cox_metrics_aggregated = Cox_metrics_aggregated,
    Cox_concordance_aggregated = Cox_concordance_aggregated,
    Cox_block = Cox_block,
    
    Combined_outcome_SMC_sig = Combined_outcome_SMC,
    
    OS_worst_prognosis_group_sig = OS_worst_prognosis_group,
    OS_p.value_sig = OS_p.value,
    
    OS_worst_prognosis_group_interaction = OS_worst_prognosis_group_interaction,
    Survival_concordance_OS = Survival_concordance_OS,
    
    DSS_worst_prognosis_group_sig = DSS_worst_prognosis_group,
    DSS_p.value_sig = DSS_p.value,
    
    DSS_worst_prognosis_group_interaction = DSS_worst_prognosis_group_interaction,
    Survival_concordance_DSS = Survival_concordance_DSS,
    
    DFI_worst_prognosis_group_sig = DFI_worst_prognosis_group,
    DFI_p.value_sig = DFI_p.value,
    
    DFI_worst_prognosis_group_interaction = DFI_worst_prognosis_group_interaction,
    Survival_concordance_DFI = Survival_concordance_DFI,
    
    PFI_worst_prognosis_group_sig = PFI_worst_prognosis_group,
    PFI_p.value_sig = PFI_p.value,
    
    PFI_worst_prognosis_group_interaction = PFI_worst_prognosis_group_interaction,
    Survival_concordance_PFI = Survival_concordance_PFI,
    
    Survival_concordance_aggregated = Survival_concordance_aggregated,
    Survival_block = Survival_block,
    
    Microenvironment_classification_sig = Microenvironment_classification,
    
    Immune_classification_sig = Immune_classification,
    
    Immune_classification_int = Immune_classification_interaction,
    
    Immune_concordance = Immune_concordance,
    Immune_block = Immune_block,
    Final_concordance_summary = Final_concordance_summary
    
  )
         
# --- Step 4 ---
df007_matched_signatures <- df006_renamed %>%
  left_join(
    df002 %>%
      select(Nomenclature, Tumor_vs_normal, Tumor_vs_normal_p.adj, Correlation_rho, Correlation_p.adj,
             Combined_outcome_HRC, Cox_OS_p.value, Cox_DSS_p.value,
             Cox_DFI_p.value, Cox_PFI_p.value, Combined_outcome_SMC,
             OS_p.value, DSS_p.value,
             DFI_p.value, PFI_p.value,
             Microenvironment_classification
      ),
    by = c(
      "Nomenclature_int" = "Nomenclature"
    )
  )


# Rename columns and reorder them
df008_renamed <- df007_matched_signatures %>%
  rename(
    Nomenclature_sig = Nomenclature_sig,
    Signatures = Signatures,
    Molecular_class_sig = Molecular_class_sig, 
    
    Interaction = Interaction,
    Nomenclature_int = Nomenclature_int,
    Molecular_class_int = Molecular_class_int,
    
    CTAB = CTAB,
    Metabolism = Metabolism,
    Pathways = Pathways,
    Metabolic_cell_death = Metabolic_cell_death,
    
    Omic_layer_sig = Omic_layer_sig,
    Phenotypic_layer_sig = Phenotypic_layer_sig,
    
    Omic_layer_int = Omic_layer_int,
    Phenotypic_layer_int = Phenotypic_layer_int,
    
    Correlation_rho_sig = Correlation_rho_sig,
    Correlation_p.adj_sig = Correlation_p.adj_sig,
    
    Correlation_rho_int = Correlation_rho,
    Correlation_p.adj_int = Correlation_p.adj,
    
    Phenotypic_concordance = Phenotypic_concordance,
    Phenotype_block = Phenotype_block, 
    
    Combined_outcome_HRC_sig = Combined_outcome_HRC_sig,
    Combined_outcome_HRC_int = Combined_outcome_HRC,
    
    Cox_OS_type_sig = Cox_OS_type_sig,
    Cox_OS_p.value_sig = Cox_OS_p.value_sig,
    
    Cox_OS_type_int = Cox_OS_type_int,
    Cox_OS_p.value_int = Cox_OS_p.value,
    Cox_concordance_OS = Cox_concordance_OS,
    
    Cox_DSS_type_sig = Cox_DSS_type_sig, 
    Cox_DSS_p.value_sig = Cox_DSS_p.value_sig,
    
    Cox_DSS_type_int = Cox_DSS_type_int,
    Cox_DSS_p.value_int = Cox_DSS_p.value,
    Cox_concordance_DSS = Cox_concordance_DSS,
    
    Cox_DFI_type_sig = Cox_DFI_type_sig, 
    Cox_DFI_p.value_sig = Cox_DFI_p.value_sig,
    
    Cox_DFI_type_int = Cox_DFI_type_int,
    Cox_DFI_p.value_int = Cox_DFI_p.value,
    Cox_concordance_DFI = Cox_concordance_DFI,
    
    Cox_PFI_type_sig = Cox_PFI_type_sig, 
    Cox_PFI_p.value_sig = Cox_PFI_p.value_sig,
    
    Cox_PFI_type_int = Cox_PFI_type_int,
    Cox_PFI_p.value_int = Cox_PFI_p.value,
    Cox_concordance_PFI = Cox_concordance_PFI,
    
    Cox_metrics_aggregated = Cox_metrics_aggregated,
    Cox_concordance_aggregated = Cox_concordance_aggregated,
    Cox_block = Cox_block,
    
    Combined_outcome_SMC = Combined_outcome_SMC,
    
    OS_worst_prognosis_group_sig = OS_worst_prognosis_group_sig,
    OS_p.value_sig = OS_p.value_sig,
    
    OS_worst_prognosis_group_int = OS_worst_prognosis_group_interaction,
    OS_p.value_int = OS_p.value,
    Survival_concordance_OS = Survival_concordance_OS,
    
    DSS_worst_prognosis_group_sig = DSS_worst_prognosis_group_sig,
    DSS_p.value_sig = DSS_p.value_sig,
    
    DSS_worst_prognosis_group_int = DSS_worst_prognosis_group_interaction,
    DSS_p.value_int = DSS_p.value,
    Survival_concordance_DSS = Survival_concordance_DSS,
    
    DFI_worst_prognosis_group_sig = DFI_worst_prognosis_group_sig,
    DFI_p.value_sig = DFI_p.value_sig,
    
    DFI_worst_prognosis_group_int = DFI_worst_prognosis_group_interaction,
    DFI_p.value_int = DFI_p.value,
    Survival_concordance_DFI = Survival_concordance_DFI,
    
    PFI_worst_prognosis_group_sig = PFI_worst_prognosis_group_sig,
    PFI_p.value_sig = PFI_p.value_sig,
    
    PFI_worst_prognosis_group_int = PFI_worst_prognosis_group_interaction,
    PFI_p.value_int = PFI_p.value,
    Survival_concordance_PFI = Survival_concordance_PFI,
    
    Survival_concordance_aggregated = Survival_concordance_aggregated,
    Survival_block = Survival_block,
    
    Microenvironment_classification_sig = Microenvironment_classification_sig,
    
    Microenvironment_classification_int = Microenvironment_classification,
    
    Immune_classification_sig = Immune_classification_sig,
    
    Immune_classification_int = Immune_classification_int,
    
    Immune_concordance = Immune_concordance,
    Immune_block = Immune_block,
    Final_concordance_summary = Final_concordance_summary
    
  ) %>%
  select(
    Nomenclature_sig, Signatures, Molecular_class_sig,
    Nomenclature_int, Interaction, Molecular_class_int,
    CTAB, Metabolism, Pathways, Metabolic_cell_death,
    Omic_layer_sig, Phenotypic_layer_sig,
    Omic_layer_int, Phenotypic_layer_int,
    Correlation_rho_sig, Correlation_p.adj_sig,
    Correlation_rho_int, Correlation_p.adj_int,
    Phenotypic_concordance,
    Combined_outcome_HRC_sig, Combined_outcome_HRC_int,
    Cox_OS_type_sig, Cox_OS_p.value_sig,
    Cox_OS_type_int, Cox_OS_p.value_int, Cox_concordance_OS,
    Cox_DSS_type_sig, Cox_DSS_p.value_sig,
    Cox_DSS_type_int, Cox_DSS_p.value_int, Cox_concordance_DSS,
    Cox_DFI_type_sig, Cox_DFI_p.value_sig,
    Cox_DFI_type_int, Cox_DFI_p.value_int, Cox_concordance_DFI,
    Cox_PFI_type_sig, Cox_PFI_p.value_sig,
    Cox_PFI_type_int, Cox_PFI_p.value_int, Cox_concordance_PFI,
    Cox_metrics_aggregated, Cox_concordance_aggregated,
    Combined_outcome_SMC,
    OS_worst_prognosis_group_sig, OS_p.value_sig,
    OS_worst_prognosis_group_int, OS_p.value_int, Survival_concordance_OS,
    DSS_worst_prognosis_group_sig, DSS_p.value_sig,
    DSS_worst_prognosis_group_int, DSS_p.value_int, Survival_concordance_DSS,
    DFI_worst_prognosis_group_sig, DFI_p.value_sig,
    DFI_worst_prognosis_group_int, DFI_p.value_int, Survival_concordance_DFI,
    PFI_worst_prognosis_group_sig, PFI_p.value_sig,
    PFI_worst_prognosis_group_int, PFI_p.value_int, Survival_concordance_PFI,
    Survival_concordance_aggregated,
    Microenvironment_classification_sig, Microenvironment_classification_int,
    Immune_classification_sig, Immune_classification_int,
    Immune_concordance, Final_concordance_summary
  )


##### Remove duplicated



# Exploratory Analysis ----------------------------------------------------

# --- Step 1: Identify if Signatures and Interactions are unique or multiple ---
df009_classified <- df008_renamed %>%
  mutate(
    Signature_count = str_count(Signatures, "\\+") + 1,
    Interaction_count = str_count(Interaction, "\\+") + 1,
    
    Signature_type = ifelse(Signature_count > 1, "Multiple", "Unique"),
    Interaction_type = ifelse(Interaction_count > 1, "Multiple", "Unique")
  )

# --- Step 2: Categorize each pair into the 4 groups ---
df009_classified <- df009_classified %>%
  mutate(
    Category = case_when(
      Signature_type == "Unique" & Interaction_type == "Unique" ~ "Unique in Both",
      Signature_type == "Multiple" & Interaction_type == "Unique" ~ "Multiple Signatures, Unique Interaction",
      Signature_type == "Unique" & Interaction_type == "Multiple" ~ "Unique Signature, Multiple Interactions",
      Signature_type == "Multiple" & Interaction_type == "Multiple" ~ "Multiple in Both",
      TRUE ~ "Unclassified"
    )
  )


df009_classified$Circuitries_id <- paste(
  sub("\\..*", "", df009_classified$Nomenclature_sig),
  sub("\\..*", "", df009_classified$Nomenclature_int),
  sep = " / "
)

df009_classified <- df009_classified[, c("Circuitries_id", setdiff(names(df009_classified), "Circuitries_id"))]


export(df009_classified, "Dataset_S5.tsv")

# --- Step 3: Create separate DataFrames for each category ---
df010_unique_in_both <- df009_classified %>% filter(Category == "Unique in Both")
df011_multiple_sig_unique_int <- df009_classified %>% filter(Category == "Multiple Signatures, Unique Interaction")
df012_unique_sig_multiple_int <- df009_classified %>% filter(Category == "Unique Signature, Multiple Interactions")
df013_multiple_in_both <- df009_classified %>% filter(Category == "Multiple in Both")

# --- Step 4: Summary table (optional, for quick overview) ---
df014_category_summary <- df009_classified %>%
  count(Category, sort = TRUE)
