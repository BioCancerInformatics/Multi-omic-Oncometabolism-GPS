## ---- Packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(rio)       # import()
  library(dplyr)
  library(stringr)
  library(tibble)
})

setwd("E:/Oncometabolism_GPS/21- Interactions_and_signatures_cross-analysis/")

## =============================================================================
## Multi-omic signature × entity validation with per-endpoint concordance
## Outputs:
##   - Per-metric concordance columns:
##       Cox_concordance_OS, Cox_concordance_DSS, Cox_concordance_PFI, Cox_concordance_DFI
##       Survival_concordance_OS, Survival_concordance_DSS, Survival_concordance_PFI, Survival_concordance_DFI
##   - Optional aggregate summaries (kept for compatibility):
##       Cox_results, Clinical_cox_concordance, Survival_results, Clinical_survival_concordance
## =============================================================================

## ---- Load data ---------------------------------------------------------------
df001 <- import("/Oncometabolism_GPS/20- Checking_signature_results/df_gene_protein_mirna_lncrna_signatures.rds")
df002 <- import("/Oncometabolism_GPS/11- Checking_results_and_annotation/3- Checking_results_and_annotation/df063_gene_protein_mirna_lncrna_results.rds")

## ---- Clinical variables ------------------------------------------------------
cox_cols   <- c("Cox_OS_type", "Cox_DSS_type", "Cox_PFI_type", "Cox_DFI_type")
surv_cols  <- c("OS_worst_prognosis_group", "DSS_worst_prognosis_group", "DFI_worst_prognosis_group", "PFI_worst_prognosis_group")
clinical_cols <- c(cox_cols, surv_cols)

## ---- Helpers -----------------------------------------------------------------
is_valid_val <- function(x) !is.na(x) & x != "NS" & x != "No data" & x != ""

safe_as_numeric <- function(x) suppressWarnings(as.numeric(x))

is_valid_num <- function(x) {
  nx <- safe_as_numeric(x)
  is.finite(nx)
}

sign_nonzero <- function(x) {
  # returns -1 for <0, +1 for >0, NA for 0 or invalid
  nx <- safe_as_numeric(x)
  if (!is.finite(nx) || nx == 0) return(NA_real_)
  if (nx > 0)  return(1)
  if (nx < 0)  return(-1)
  NA_real_
}

# Per-endpoint concordance (single metric, categorical equality)
concordance_per_metric <- function(row_val, match_val) {
  if (!is_valid_val(row_val) || !is_valid_val(match_val)) {
    return("NS")
  } else if (row_val == match_val) {
    return("Convergent")
  } else {
    return("Divergent")
  }
}

# Aggregate concordance across endpoints given named vector in {"Convergent","Divergent","NS"}
aggregate_concordance <- function(vec_named) {
  vec <- vec_named[vec_named != "NS"]
  if (length(vec) == 0) {
    list(metrics = NA_character_, category = "NS")
  } else {
    eq   <- sum(vec == "Convergent")
    diff <- sum(vec == "Divergent")
    catg <- dplyr::case_when(
      eq > diff  ~ "Convergent",
      diff > eq  ~ "Divergent",
      TRUE       ~ "Mixed"
    )
    metrics <- names(vec_named)[vec_named != "NS"] %>% paste(collapse = " / ")
    list(metrics = metrics, category = catg)
  }
}

# ---- NEW: Phenotypic_layer concordance based on RHO sign ----------------------
# Rules:
# 1) If Phenotypic_layer differs -> "NS"
# 2) If same Phenotypic_layer:
#    - If either RHO is invalid or zero -> "NS"
#    - If signs opposite -> "Divergent"
#    - If signs equal -> "Convergent"
phenotypic_concordance_rho <- function(phen_sig, phen_ent, rho_sig, rho_ent) {
  if (!is_valid_val(phen_sig) || !is_valid_val(phen_ent)) return("NS")
  if (!identical(as.character(phen_sig), as.character(phen_ent))) return("NS")
  s_sig <- sign_nonzero(rho_sig)
  s_ent <- sign_nonzero(rho_ent)
  if (is.na(s_sig) || is.na(s_ent)) return("NS")
  if (s_sig != s_ent) return("Divergent")
  "Convergent"
}

## ---- Pre-filter: keep rows with at least one valid clinical metric -----------
df001_filtered <- df001 %>% filter(if_any(all_of(clinical_cols), is_valid_val))
df002_filtered <- df002 %>% filter(if_any(all_of(clinical_cols), is_valid_val))

## ---- Main loop ---------------------------------------------------------------
total_rows <- nrow(df001_filtered)
pb <- utils::txtProgressBar(min = 0, max = total_rows, style = 3)
result_list <- list()

for (i in seq_len(total_rows)) {
  row <- df001_filtered[i, ]
  
  # Split common interactions (one signature may map to multiple entities)
  interactions <- str_split(row$Common_interaction, " / ", simplify = TRUE) %>% str_trim()
  
  for (entity in interactions) {
    
    # Restrict candidate entity by cancer type
    match_row <- df002_filtered %>%
      filter(Target == entity, Cancer_types == row$Cancer_types)
    
    if (nrow(match_row) > 0) {
      match_row <- match_row[1, , drop = FALSE]  # first occurrence only
      
      # ---- Per-endpoint concordance: COX ----
      cox_conc_OS  <- concordance_per_metric(row$Cox_OS_type,  match_row$Cox_OS_type)
      cox_conc_DSS <- concordance_per_metric(row$Cox_DSS_type, match_row$Cox_DSS_type)
      cox_conc_PFI <- concordance_per_metric(row$Cox_PFI_type, match_row$Cox_PFI_type)
      cox_conc_DFI <- concordance_per_metric(row$Cox_DFI_type, match_row$Cox_DFI_type)
      
      cox_vec <- c(OS = cox_conc_OS, DSS = cox_conc_DSS, PFI = cox_conc_PFI, DFI = cox_conc_DFI)
      cox_agg <- aggregate_concordance(cox_vec)
      
      # ---- Per-endpoint concordance: SURVIVAL worst-prognosis group ----
      surv_conc_OS  <- concordance_per_metric(row$OS_worst_prognosis_group,  match_row$OS_worst_prognosis_group)
      surv_conc_DSS <- concordance_per_metric(row$DSS_worst_prognosis_group, match_row$DSS_worst_prognosis_group)
      surv_conc_PFI <- concordance_per_metric(row$PFI_worst_prognosis_group, match_row$PFI_worst_prognosis_group)
      surv_conc_DFI <- concordance_per_metric(row$DFI_worst_prognosis_group, match_row$DFI_worst_prognosis_group)
      
      surv_vec <- c(OS = surv_conc_OS, DSS = surv_conc_DSS, PFI = surv_conc_PFI, DFI = surv_conc_DFI)
      surv_agg <- aggregate_concordance(surv_vec)
      
      # ---- Immune classification concordance (unchanged) ----------------------
      immune_sig <- row$Immune_classification
      immune_ent <- match_row$Immune_classification
      immune_concordance <-
        if (!is_valid_val(immune_sig) || !is_valid_val(immune_ent)) {
          "NS"
        } else if (immune_sig == immune_ent) {
          "Convergent"
        } else {
          "Divergent"
        }
      
      # # ---- Microenvironment classification concordance (unchanged) ----------------------
      # Microenv_sig <- row$Microenvironment_classification
      # Microenv_ent <- match_row$Microenvironment_classification
      # Microenv_concordance <-
      #   if (!is_valid_val(Microenv_sig) || !is_valid_val(Microenv_ent)) {
      #     "NS"
      #   } else if (Microenv_sig == Microenv_ent) {
      #     "Convergent"
      #   } else {
      #     "Divergent"
      #   }
      
      # ---- Phenotypic layer concordance (NEW logic using RHO) -----------------
      phen_sig <- row$Phenotypic_layer
      phen_ent <- match_row$Phenotypic_layer
      
      # Expect RHO columns present on both tables; adapt names here if needed.
      rho_sig  <- row$Correlation_rho
      rho_ent  <- match_row$Correlation_rho
      
      phenotypic_conc <- phenotypic_concordance_rho(phen_sig, phen_ent, rho_sig, rho_ent)
      
      # ---- Output row ---------------------------------------------------------
      output_row <- tibble(
        # Signature-side metadata
        Signatures = row$Multiomics_Signature,
        Molecular_class_signature = row$Molecular_class,
        Metabolism_signature = row$Metabolism,
        Pathways_signature = row$Pathways,
        Metabolic_cell_death_signature = row$Metabolic_cell_death,
        CTAB = row$Cancer_types,
        Omic_layer_signature = row$Omic_layer,
        Phenotypic_layer_signature = phen_sig,
        Correlation_rho_signature = row$Correlation_rho,
        Common_interaction = row$Common_interaction,
        # Microenvironment_classification_signature = row$Microenvironment_classification,
        Immune_classification_signature = immune_sig,
        
        # Signature clinical metrics (reference)
        Cox_OS_type_signature  = row$Cox_OS_type,
        Cox_DSS_type_signature = row$Cox_DSS_type,
        Cox_PFI_type_signature = row$Cox_PFI_type,
        Cox_DFI_type_signature = row$Cox_DFI_type,
        OS_worst_prognosis_group_signature  = row$OS_worst_prognosis_group,
        DSS_worst_prognosis_group_signature = row$DSS_worst_prognosis_group,
        DFI_worst_prognosis_group_signature = row$DFI_worst_prognosis_group,
        PFI_worst_prognosis_group_signature = row$PFI_worst_prognosis_group,
        
        # Entity-side metadata
        Meaningful_interaction = entity,
        Molecular_class_interaction = match_row$Molecular_class,
        Metabolism_interaction = match_row$Metabolism,
        Pathways_interaction = match_row$Pathways,
        Metabolic_cell_death_interaction = match_row$Metabolic_cell_death,
        Omic_layer_interaction = match_row$Omic_layer,
        Phenotypic_layer_interaction = phen_ent,
        Correlation_rho_interaction = row$Correlation_rho,
        # Microenvironment_classification_entity = match_row$Microenvironment_classification,
        Immune_classification_interaction = immune_ent,
        
        # NEW: Phenotypic layer concordance based on RHO
        rho_signature = safe_as_numeric(rho_sig),
        rho_interaction    = safe_as_numeric(rho_ent),
        rho_sign_signature = sign_nonzero(rho_sig),
        rho_sign_interaction    = sign_nonzero(rho_ent),
        Phenotypic_concordance = phenotypic_conc,
        
        # Entity clinical metrics (reference)
        Cox_OS_type_interaction  = match_row$Cox_OS_type,
        Cox_DSS_type_interaction = match_row$Cox_DSS_type,
        Cox_PFI_type_interaction = match_row$Cox_PFI_type,
        Cox_DFI_type_interaction = match_row$Cox_DFI_type,
        OS_worst_prognosis_group_interaction  = match_row$OS_worst_prognosis_group,
        DSS_worst_prognosis_group_interaction = match_row$DSS_worst_prognosis_group,
        DFI_worst_prognosis_group_interaction = match_row$DFI_worst_prognosis_group,
        PFI_worst_prognosis_group_interaction = match_row$PFI_worst_prognosis_group,
        
        # Per-endpoint concordance (kept)
        Cox_concordance_OS  = cox_conc_OS,
        Cox_concordance_DSS = cox_conc_DSS,
        Cox_concordance_PFI = cox_conc_PFI,
        Cox_concordance_DFI = cox_conc_DFI,
        Survival_concordance_OS  = surv_conc_OS,
        Survival_concordance_DSS = surv_conc_DSS,
        Survival_concordance_PFI = surv_conc_PFI,
        Survival_concordance_DFI = surv_conc_DFI,
        
        # Optional aggregate summaries
        Cox_metrics_aggregated = cox_agg$metrics,
        Cox_concordance_aggregated = cox_agg$category,
        Survival_metrics_aggregated = surv_agg$metrics,
        Survival_concordance_aggregated = surv_agg$category,
        
        # # Microenvironment concordance
        # Microenv_concordance = Microenv_concordance,
        
        # Immune concordance (unchanged)
        Immune_concordance = immune_concordance
      )
      
      result_list[[length(result_list) + 1]] <- output_row
    } # if match_row
  }   # for entity
  
  utils::setTxtProgressBar(pb, i)
}
close(pb)

## ---- Final outputs -----------------------------------------------------------
df003_output <- dplyr::bind_rows(result_list)
df004_output_filtered <- df003_output %>% dplyr::distinct()

## 1) Columns to remove (hard redundancies)
drop_cols <- c("Correlation_rho_signature", "Correlation_rho_interaction")

df005_organized <- df004_output_filtered %>% select(-any_of(drop_cols))

## 2) Analysis-ready ordered selection
analysis_cols <- c(
  # Keys
  "Signatures", "CTAB", "Meaningful_interaction",
  # Signature metadata
  "Molecular_class_signature", "Omic_layer_signature",
  "Metabolism_signature", "Pathways_signature", "Metabolic_cell_death_signature",
  # Entity metadata
  "Molecular_class_interaction", "Omic_layer_interaction",
  "Metabolism_interaction", "Pathways_interaction", "Metabolic_cell_death_interaction",
  # Phenotype (RHO-based)
  "Phenotypic_layer_signature", "Phenotypic_layer_interaction",
  "rho_signature", "rho_interaction", "Phenotypic_concordance",
  # Microenvironment & Immune
  "Microenvironment_classification_signature", "Microenvironment_classification_interaction", 
  # "Microenv_concordance",
  "Immune_classification_signature", "Immune_classification_interaction", "Immune_concordance",
  # Cox per-endpoint
  "Cox_OS_type_signature", "Cox_OS_type_interaction", "Cox_concordance_OS", 
  "Cox_DSS_type_signature", "Cox_DSS_type_interaction", "Cox_concordance_DSS", 
  "Cox_PFI_type_signature", "Cox_PFI_type_interaction", "Cox_concordance_DFI",
  "Cox_DFI_type_signature", "Cox_DFI_type_interaction", "Cox_concordance_PFI", 
  # Survival per-endpoint
  "OS_worst_prognosis_group_signature", "OS_worst_prognosis_group_interaction", "Survival_concordance_OS", 
  "DSS_worst_prognosis_group_signature", "DSS_worst_prognosis_group_interaction", "Survival_concordance_DSS", 
  "DFI_worst_prognosis_group_signature", "DFI_worst_prognosis_group_interaction", "Survival_concordance_DFI",
  "PFI_worst_prognosis_group_signature", "PFI_worst_prognosis_group_interaction", "Survival_concordance_PFI", 
  # Optional aggregates (comment these four lines out if you don’t want them)
  "Cox_metrics_aggregated", "Cox_concordance_aggregated",
  "Survival_metrics_aggregated", "Survival_concordance_aggregated"
)

df005_organized <- df005_organized %>%
  # keep only columns that exist (robust to slight naming changes)
  select(any_of(analysis_cols))

# ## 3) Audit-ready = analysis + full raw traces
# audit_extras <- c(
#   # Interaction context
#   "Common_interaction",
#   # RHO signs (QC)
#   "rho_sign_signature", "rho_sign_interaction",
#   # Raw signature clinical labels
#   "Cox_OS_type_signature", "Cox_DSS_type_signature",
#   "Cox_PFI_type_signature", "Cox_DFI_type_signature",
#   "OS_worst_prognosis_group_signature", "DSS_worst_prognosis_group_signature",
#   "DFI_worst_prognosis_group_signature", "PFI_worst_prognosis_group_signature",
#   # Raw entity clinical labels
#   "Cox_OS_type_interaction", "Cox_DSS_type_interaction",
#   "Cox_PFI_type_interaction", "Cox_DFI_type_interaction",
#   "OS_worst_prognosis_group_interaction", "DSS_worst_prognosis_group_interaction",
#   "DFI_worst_prognosis_group_interaction", "PFI_worst_prognosis_group_interaction"
# )

# 
# df_audit <- df005_full %>%
#   select(any_of(c(analysis_cols, audit_extras)))
## 4) (Optional) Consistent renaming for aggregates
# If you prefer shorter names:
# df_analysis <- df_analysis %>%
#   rename(
#     Cox_results = Cox_metrics_aggregated,
#     Clinical_cox_concordance = Cox_concordance_aggregated,
#     Survival_results = Survival_metrics_aggregated,
#     Clinical_survival_concordance = Survival_concordance_aggregated
#   )
# df_audit <- df_audit %>%
#   rename(
#     Cox_results = Cox_metrics_aggregated,
#     Clinical_cox_concordance = Cox_concordance_aggregated,
#     Survival_results = Survival_metrics_aggregated,
#     Clinical_survival_concordance = Survival_concordance_aggregated
#   )


## --- helpers ------------------------------------------------------------------
# Collapse a vector of endpoint-level concordances to a single block label.
# Input values should be in {"Convergent","Divergent","Mixed","NS", NA}.
collapse_block <- function(x) {
  x <- x[!is.na(x) & x != "NS"]
  if (length(x) == 0) return("NS")
  # If all equal (all Convergent OR all Divergent OR all Mixed)
  if (length(unique(x)) == 1) return(unique(x))
  # Otherwise there is heterogeneity across endpoints
  "Mixed"
}

# Build a compact, human-readable summary string from named block labels,
# skipping blocks that are "NS".
compose_summary <- function(named_blocks) {
  nb <- named_blocks[named_blocks != "NS"]
  if (length(nb) == 0) return("NS")
  # If all Convergent:
  if (all(nb == "Convergent")) return("Only Convergent")
  # If all Divergent:
  if (all(nb == "Divergent")) return("Only Divergent")
  # Otherwise, compose detailed summary in a stable order:
  order_vec <- c("Immune", "Microenvironment", "Phenotype", "Cox", "Survival")
  nb <- nb[intersect(order_vec, names(nb))]
  paste0(names(nb), ": ", nb, collapse = "; ")
}

## --- computation --------------------------------------------------------------
df006_final <- df005_organized %>%
  rowwise() %>%
  mutate(
    # Aggregate Cox across endpoints
    Cox_block = collapse_block(c_across(any_of(c(
      "Cox_concordance_OS", "Cox_concordance_DSS",
      "Cox_concordance_PFI","Cox_concordance_DFI"
    )))),
    # Aggregate Survival across endpoints
    Survival_block = collapse_block(c_across(any_of(c(
      "Survival_concordance_OS","Survival_concordance_DSS",
      "Survival_concordance_PFI","Survival_concordance_DFI"
    )))),
    # Single-field blocks
    Immune_block = coalesce(Immune_concordance, "NS"),
    # Microenvironment_block = coalesce(Microenv_concordance, "NS"),
    Phenotype_block = coalesce(Phenotypic_concordance, "NS"),
    # Compose final per-row summary
    Final_concordance_summary = compose_summary(c(
      Immune = Immune_block,
      # Microenvironment = Microenvironment_block,
      Phenotype = Phenotype_block,
      Cox = Cox_block,
      Survival = Survival_block
    ))
  ) %>%
  ungroup()

colnames(df006_final)

# Salvar resultado
saveRDS(df006_final, "Regulatory_circuitries_results.rds")

###### Combined_outcome Survival column 

# Define semantic categories
inconclusive_terms <- c("NS", "No data")

# All known biological terms
biological_terms <- c(
  "Low", "High", "WT", "MT", "Normal", "Deleted", "Duplicated",
  "Duplicated and Deleted", "Duplicated and Normal",
  "Deleted and Duplicated", "Normal and Duplicated",
  "Deleted and Normal", "Normal and Deleted"
)

# Function to classify a single row
classify_combined_outcome <- function(values) {
  values <- na.omit(values)
  values <- unique(trimws(values))
  values <- values[values != ""]  # Remove empty strings
  
  # Partition
  has_inconclusive <- any(values %in% inconclusive_terms)
  has_biological <- any(values %in% biological_terms)
  
  all_inconclusive <- all(values %in% inconclusive_terms)
  all_same <- length(unique(values[values %in% biological_terms])) == 1 && !has_inconclusive
  
  if (length(values) == 0 || all_inconclusive) {
    return("Meaningless")
  } else if (all_same) {
    return(values[values %in% biological_terms][1])  # Return the biological term
  } else if (has_biological && !has_inconclusive) {
    return("Meaningful Mixed")
  } else if (has_biological && has_inconclusive) {
    return(paste("Meaningful", unique(values[values %in% biological_terms][1])))
  } else {
    return("Undefined")
  }
}

# Apply to df006
df007_outcome_SMC <- df006_final %>%
  rowwise() %>%
  mutate(Combined_outcome_SMC = classify_combined_outcome(c_across(c(
    OS_worst_prognosis_group_interaction,
    DSS_worst_prognosis_group_interaction,
    DFI_worst_prognosis_group_interaction,
    PFI_worst_prognosis_group_interaction
  )))) %>%
  ungroup()


###### Combined_outcome Cox column 

# Define semantic categories
inconclusive_terms <- c("NS", "No data")

# All known biological terms
biological_terms <- c(
  "Protective", "Risky"
)

# Function to classify a single row
classify_combined_outcome <- function(values) {
  values <- na.omit(values)
  values <- unique(trimws(values))
  values <- values[values != ""]  # Remove empty strings
  
  # Partition
  has_inconclusive <- any(values %in% inconclusive_terms)
  has_biological <- any(values %in% biological_terms)
  
  all_inconclusive <- all(values %in% inconclusive_terms)
  all_same <- length(unique(values[values %in% biological_terms])) == 1 && !has_inconclusive
  
  if (length(values) == 0 || all_inconclusive) {
    return("Meaningless")
  } else if (all_same) {
    return(values[values %in% biological_terms][1])  # Return the biological term
  } else if (has_biological && !has_inconclusive) {
    return("Meaningful Mixed")
  } else if (has_biological && has_inconclusive) {
    return(paste("Meaningful", unique(values[values %in% biological_terms][1])))
  } else {
    return("Undefined")
  }
}

# Apply to df009
df008_outcome_HRC <- df007_outcome_SMC %>%
  rowwise() %>%
  mutate(Combined_outcome_HRC = classify_combined_outcome(c_across(c(
    Cox_OS_type_interaction,
    Cox_DSS_type_interaction,
    Cox_DFI_type_interaction,
    Cox_PFI_type_interaction
  )))) %>%
  ungroup()


# Salvar resultado
saveRDS(df008_outcome_HRC, "Regulatory_circuitries_with_HRC_SMC.rds")

# Step 1: Summarize validated interactions from df006_final
df009_interaction_summary <- df006_final %>%
  group_by(Signatures, Omic_layer_signature, Metabolism_signature, Pathways_signature, CTAB) %>%
  summarise(Meaningful_interaction = str_c(unique(Meaningful_interaction), collapse = " / "), .groups = "drop")

# Step 2: Rename to match if needed
df009_interaction_summary <- df009_interaction_summary %>%
  rename(Omic_layer = Omic_layer_signature,
         Metabolism = Metabolism_signature,
         Pathways = Pathways_signature)

df001 <- df001 %>%
  rename(Signatures = Multiomics_Signature,
         CTAB = Cancer_types)

# Step 3: Left join and fill NAs with "No meaningful interaction"
df010_final_result <- df001 %>%
  left_join(df009_interaction_summary, by = c("Signatures", "Omic_layer", "Metabolism", "Pathways", "CTAB")) %>%
  mutate(Meaningful_interaction = replace_na(Meaningful_interaction, "No meaningful interaction"))

# Reorder the columns
New_column_order_01 <- c("Members", "Signatures", "Molecular_class", "Common_interaction", "Meaningful_interaction", "Metabolism", "Pathways", "Metabolic_cell_death", "CTAB", 
                         "Tumor_vs_normal",  "Tumor_vs_normal_p.adj", "Omic_layer", "Phenotypic_layer", "Correlation_rho", "Correlation_p.adj", "Cox_OS_type", "Cox_OS_p.value", 
                         "Cox_DSS_type", "Cox_DSS_p.value", "Cox_DFI_type", "Cox_DFI_p.value", "Cox_PFI_type", "Cox_PFI_p.value", "OS_worst_prognosis_group", "OS_p.value", 
                         "DSS_worst_prognosis_group", "DSS_p.value", "DFI_worst_prognosis_group", "DFI_p.value", "PFI_worst_prognosis_group", "PFI_p.value",
                         "Microenvironment_classification", "Microenvironment_score_details", "Immune_classification", "Immune_score_details")

df010_final_result <- df010_final_result[New_column_order_01]

# Salvar resultado
saveRDS(df010_final_result, "df_gene_protein_mirna_lncrna_signatures.rds")


# Transcript signatures #######################################################

## ---- Load data ---------------------------------------------------------------
df011 <- import("/Oncometabolism_GPS/20- Checking_signature_results/df_transcript_signatures.rds")
df012 <- import("/Oncometabolism_GPS/11- Checking_results_and_annotation/3- Checking_results_and_annotation/df063_gene_protein_mirna_lncrna_results.rds")

## ---- Clinical variables ------------------------------------------------------
cox_cols   <- c("Cox_OS_type", "Cox_DSS_type", "Cox_PFI_type", "Cox_DFI_type")
surv_cols  <- c("OS_worst_prognosis_group", "DSS_worst_prognosis_group", "DFI_worst_prognosis_group", "PFI_worst_prognosis_group")
clinical_cols <- c(cox_cols, surv_cols)

## ---- Helpers -----------------------------------------------------------------
is_valid_val <- function(x) !is.na(x) & x != "NS" & x != "No data" & x != ""

safe_as_numeric <- function(x) suppressWarnings(as.numeric(x))

is_valid_num <- function(x) {
  nx <- safe_as_numeric(x)
  is.finite(nx)
}

sign_nonzero <- function(x) {
  # returns -1 for <0, +1 for >0, NA for 0 or invalid
  nx <- safe_as_numeric(x)
  if (!is.finite(nx) || nx == 0) return(NA_real_)
  if (nx > 0)  return(1)
  if (nx < 0)  return(-1)
  NA_real_
}

# Per-endpoint concordance (single metric, categorical equality)
concordance_per_metric <- function(row_val, match_val) {
  if (!is_valid_val(row_val) || !is_valid_val(match_val)) {
    return("NS")
  } else if (row_val == match_val) {
    return("Convergent")
  } else {
    return("Divergent")
  }
}

# Aggregate concordance across endpoints given named vector in {"Convergent","Divergent","NS"}
aggregate_concordance <- function(vec_named) {
  vec <- vec_named[vec_named != "NS"]
  if (length(vec) == 0) {
    list(metrics = NA_character_, category = "NS")
  } else {
    eq   <- sum(vec == "Convergent")
    diff <- sum(vec == "Divergent")
    catg <- dplyr::case_when(
      eq > diff  ~ "Convergent",
      diff > eq  ~ "Divergent",
      TRUE       ~ "Mixed"
    )
    metrics <- names(vec_named)[vec_named != "NS"] %>% paste(collapse = " / ")
    list(metrics = metrics, category = catg)
  }
}

# ---- NEW: Phenotypic_layer concordance based on RHO sign ----------------------
# Rules:
# 1) If Phenotypic_layer differs -> "NS"
# 2) If same Phenotypic_layer:
#    - If either RHO is invalid or zero -> "NS"
#    - If signs opposite -> "Divergent"
#    - If signs equal -> "Convergent"
phenotypic_concordance_rho <- function(phen_sig, phen_ent, rho_sig, rho_ent) {
  if (!is_valid_val(phen_sig) || !is_valid_val(phen_ent)) return("NS")
  if (!identical(as.character(phen_sig), as.character(phen_ent))) return("NS")
  s_sig <- sign_nonzero(rho_sig)
  s_ent <- sign_nonzero(rho_ent)
  if (is.na(s_sig) || is.na(s_ent)) return("NS")
  if (s_sig != s_ent) return("Divergent")
  "Convergent"
}

## ---- Pre-filter: keep rows with at least one valid clinical metric -----------
df011_filtered <- df011 %>% filter(if_any(all_of(clinical_cols), is_valid_val))
df012_filtered <- df012 %>% filter(if_any(all_of(clinical_cols), is_valid_val))

## ---- Main loop ---------------------------------------------------------------
total_rows <- nrow(df011_filtered)
pb <- utils::txtProgressBar(min = 0, max = total_rows, style = 3)
result_list <- list()

for (i in seq_len(total_rows)) {
  row <- df011_filtered[i, ]
  
  # Split common interactions (one signature may map to multiple entities)
  interactions <- str_split(row$Common_interaction, " / ", simplify = TRUE) %>% str_trim()
  
  for (entity in interactions) {
    
    # Restrict candidate entity by cancer type
    match_row <- df012_filtered %>%
      filter(Target == entity, Cancer_types == row$Cancer_types)
    
    if (nrow(match_row) > 0) {
      match_row <- match_row[1, , drop = FALSE]  # first occurrence only
      
      # ---- Per-endpoint concordance: COX ----
      cox_conc_OS  <- concordance_per_metric(row$Cox_OS_type,  match_row$Cox_OS_type)
      cox_conc_DSS <- concordance_per_metric(row$Cox_DSS_type, match_row$Cox_DSS_type)
      cox_conc_PFI <- concordance_per_metric(row$Cox_PFI_type, match_row$Cox_PFI_type)
      cox_conc_DFI <- concordance_per_metric(row$Cox_DFI_type, match_row$Cox_DFI_type)
      
      cox_vec <- c(OS = cox_conc_OS, DSS = cox_conc_DSS, PFI = cox_conc_PFI, DFI = cox_conc_DFI)
      cox_agg <- aggregate_concordance(cox_vec)
      
      # ---- Per-endpoint concordance: SURVIVAL worst-prognosis group ----
      surv_conc_OS  <- concordance_per_metric(row$OS_worst_prognosis_group,  match_row$OS_worst_prognosis_group)
      surv_conc_DSS <- concordance_per_metric(row$DSS_worst_prognosis_group, match_row$DSS_worst_prognosis_group)
      surv_conc_PFI <- concordance_per_metric(row$PFI_worst_prognosis_group, match_row$PFI_worst_prognosis_group)
      surv_conc_DFI <- concordance_per_metric(row$DFI_worst_prognosis_group, match_row$DFI_worst_prognosis_group)
      
      surv_vec <- c(OS = surv_conc_OS, DSS = surv_conc_DSS, PFI = surv_conc_PFI, DFI = surv_conc_DFI)
      surv_agg <- aggregate_concordance(surv_vec)
      
      # ---- Immune classification concordance (unchanged) ----------------------
      immune_sig <- row$Immune_classification
      immune_ent <- match_row$Immune_classification
      immune_concordance <-
        if (!is_valid_val(immune_sig) || !is_valid_val(immune_ent)) {
          "NS"
        } else if (immune_sig == immune_ent) {
          "Convergent"
        } else {
          "Divergent"
        }
      
      # # ---- Microenvironment classification concordance (unchanged) ----------------------
      # Microenv_sig <- row$Microenvironment_classification
      # Microenv_ent <- match_row$Microenvironment_classification
      # Microenv_concordance <-
      #   if (!is_valid_val(Microenv_sig) || !is_valid_val(Microenv_ent)) {
      #     "NS"
      #   } else if (Microenv_sig == Microenv_ent) {
      #     "Convergent"
      #   } else {
      #     "Divergent"
      #   }
      
      # ---- Phenotypic layer concordance (NEW logic using RHO) -----------------
      phen_sig <- row$Phenotypic_layer
      phen_ent <- match_row$Phenotypic_layer
      
      # Expect RHO columns present on both tables; adapt names here if needed.
      rho_sig  <- row$Correlation_rho
      rho_ent  <- match_row$Correlation_rho
      
      phenotypic_conc <- phenotypic_concordance_rho(phen_sig, phen_ent, rho_sig, rho_ent)
      
      # ---- Output row ---------------------------------------------------------
      output_row <- tibble(
        # Signature-side metadata
        Signatures = row$Multiomics_Signature,
        Molecular_class_signature = row$Molecular_class,
        Metabolism_signature = row$Metabolism,
        Pathways_signature = row$Pathways,
        Metabolic_cell_death_signature = row$Metabolic_cell_death,
        CTAB = row$Cancer_types,
        Omic_layer_signature = row$Omic_layer,
        Phenotypic_layer_signature = phen_sig,
        Correlation_rho_signature = row$Correlation_rho,
        Common_interaction = row$Common_interaction,
        # Microenvironment_classification_signature = row$Microenvironment_classification,
        Immune_classification_signature = immune_sig,
        
        # Signature clinical metrics (reference)
        Cox_OS_type_signature  = row$Cox_OS_type,
        Cox_DSS_type_signature = row$Cox_DSS_type,
        Cox_PFI_type_signature = row$Cox_PFI_type,
        Cox_DFI_type_signature = row$Cox_DFI_type,
        OS_worst_prognosis_group_signature  = row$OS_worst_prognosis_group,
        DSS_worst_prognosis_group_signature = row$DSS_worst_prognosis_group,
        DFI_worst_prognosis_group_signature = row$DFI_worst_prognosis_group,
        PFI_worst_prognosis_group_signature = row$PFI_worst_prognosis_group,
        
        # Entity-side metadata
        Meaningful_interaction = entity,
        Molecular_class_interaction = match_row$Molecular_class,
        Metabolism_interaction = match_row$Metabolism,
        Pathways_interaction = match_row$Pathways,
        Metabolic_cell_death_interaction = match_row$Metabolic_cell_death,
        Omic_layer_interaction = match_row$Omic_layer,
        Phenotypic_layer_interaction = phen_ent,
        Correlation_rho_interaction = row$Correlation_rho,
        # Microenvironment_classification_entity = match_row$Microenvironment_classification,
        Immune_classification_interaction = immune_ent,
        
        # NEW: Phenotypic layer concordance based on RHO
        rho_signature = safe_as_numeric(rho_sig),
        rho_interaction    = safe_as_numeric(rho_ent),
        rho_sign_signature = sign_nonzero(rho_sig),
        rho_sign_interaction    = sign_nonzero(rho_ent),
        Phenotypic_concordance = phenotypic_conc,
        
        # Entity clinical metrics (reference)
        Cox_OS_type_interaction  = match_row$Cox_OS_type,
        Cox_DSS_type_interaction = match_row$Cox_DSS_type,
        Cox_PFI_type_interaction = match_row$Cox_PFI_type,
        Cox_DFI_type_interaction = match_row$Cox_DFI_type,
        OS_worst_prognosis_group_interaction  = match_row$OS_worst_prognosis_group,
        DSS_worst_prognosis_group_interaction = match_row$DSS_worst_prognosis_group,
        DFI_worst_prognosis_group_interaction = match_row$DFI_worst_prognosis_group,
        PFI_worst_prognosis_group_interaction = match_row$PFI_worst_prognosis_group,
        
        # Per-endpoint concordance (kept)
        Cox_concordance_OS  = cox_conc_OS,
        Cox_concordance_DSS = cox_conc_DSS,
        Cox_concordance_PFI = cox_conc_PFI,
        Cox_concordance_DFI = cox_conc_DFI,
        Survival_concordance_OS  = surv_conc_OS,
        Survival_concordance_DSS = surv_conc_DSS,
        Survival_concordance_PFI = surv_conc_PFI,
        Survival_concordance_DFI = surv_conc_DFI,
        
        # Optional aggregate summaries
        Cox_metrics_aggregated = cox_agg$metrics,
        Cox_concordance_aggregated = cox_agg$category,
        Survival_metrics_aggregated = surv_agg$metrics,
        Survival_concordance_aggregated = surv_agg$category,
        
        # # Microenvironment concordance
        # Microenv_concordance = Microenv_concordance,
        
        # Immune concordance (unchanged)
        Immune_concordance = immune_concordance
      )
      
      result_list[[length(result_list) + 1]] <- output_row
    } # if match_row
  }   # for entity
  
  utils::setTxtProgressBar(pb, i)
}
close(pb)


## ---- Final outputs -----------------------------------------------------------
df013_output <- dplyr::bind_rows(result_list)
df014_output_filtered <- df013_output %>% dplyr::distinct()

## 1) Columns to remove (hard redundancies)
drop_cols <- c("Correlation_rho_signature", "Correlation_rho_interaction")

df015_organized <- df014_output_filtered %>% select(-any_of(drop_cols))

## 2) Analysis-ready ordered selection
analysis_cols <- c(
  # Keys
  "Signatures", "CTAB", "Meaningful_interaction",
  # Signature metadata
  "Molecular_class_signature", "Omic_layer_signature",
  "Metabolism_signature", "Pathways_signature", "Metabolic_cell_death_signature",
  # Entity metadata
  "Molecular_class_interaction", "Omic_layer_interaction",
  "Metabolism_interaction", "Pathways_interaction", "Metabolic_cell_death_interaction",
  # Phenotype (RHO-based)
  "Phenotypic_layer_signature", "Phenotypic_layer_interaction",
  "rho_signature", "rho_interaction", "Phenotypic_concordance",
  # Microenvironment & Immune
  "Microenvironment_classification_signature", "Microenvironment_classification_interaction", 
  # "Microenv_concordance",
  "Immune_classification_signature", "Immune_classification_interaction", "Immune_concordance",
  # Cox per-endpoint
  "Cox_OS_type_signature", "Cox_OS_type_interaction", "Cox_concordance_OS", 
  "Cox_DSS_type_signature", "Cox_DSS_type_interaction", "Cox_concordance_DSS", 
  "Cox_PFI_type_signature", "Cox_PFI_type_interaction", "Cox_concordance_DFI",
  "Cox_DFI_type_signature", "Cox_DFI_type_interaction", "Cox_concordance_PFI", 
  # Survival per-endpoint
  "OS_worst_prognosis_group_signature", "OS_worst_prognosis_group_interaction", "Survival_concordance_OS", 
  "DSS_worst_prognosis_group_signature", "DSS_worst_prognosis_group_interaction", "Survival_concordance_DSS", 
  "DFI_worst_prognosis_group_signature", "DFI_worst_prognosis_group_interaction", "Survival_concordance_DFI",
  "PFI_worst_prognosis_group_signature", "PFI_worst_prognosis_group_interaction", "Survival_concordance_PFI", 
  # Optional aggregates (comment these four lines out if you don’t want them)
  "Cox_metrics_aggregated", "Cox_concordance_aggregated",
  "Survival_metrics_aggregated", "Survival_concordance_aggregated"
)

df015_organized <- df015_organized %>%
  # keep only columns that exist (robust to slight naming changes)
  select(any_of(analysis_cols))

# ## 3) Audit-ready = analysis + full raw traces
# audit_extras <- c(
#   # Interaction context
#   "Common_interaction",
#   # RHO signs (QC)
#   "rho_sign_signature", "rho_sign_interaction",
#   # Raw signature clinical labels
#   "Cox_OS_type_signature", "Cox_DSS_type_signature",
#   "Cox_PFI_type_signature", "Cox_DFI_type_signature",
#   "OS_worst_prognosis_group_signature", "DSS_worst_prognosis_group_signature",
#   "DFI_worst_prognosis_group_signature", "PFI_worst_prognosis_group_signature",
#   # Raw entity clinical labels
#   "Cox_OS_type_interaction", "Cox_DSS_type_interaction",
#   "Cox_PFI_type_interaction", "Cox_DFI_type_interaction",
#   "OS_worst_prognosis_group_interaction", "DSS_worst_prognosis_group_interaction",
#   "DFI_worst_prognosis_group_interaction", "PFI_worst_prognosis_group_interaction"
# )

# 
# df_audit <- df005_full %>%
#   select(any_of(c(analysis_cols, audit_extras)))
## 4) (Optional) Consistent renaming for aggregates
# If you prefer shorter names:
# df_analysis <- df_analysis %>%
#   rename(
#     Cox_results = Cox_metrics_aggregated,
#     Clinical_cox_concordance = Cox_concordance_aggregated,
#     Survival_results = Survival_metrics_aggregated,
#     Clinical_survival_concordance = Survival_concordance_aggregated
#   )
# df_audit <- df_audit %>%
#   rename(
#     Cox_results = Cox_metrics_aggregated,
#     Clinical_cox_concordance = Cox_concordance_aggregated,
#     Survival_results = Survival_metrics_aggregated,
#     Clinical_survival_concordance = Survival_concordance_aggregated
#   )


## --- helpers ------------------------------------------------------------------
# Collapse a vector of endpoint-level concordances to a single block label.
# Input values should be in {"Convergent","Divergent","Mixed","NS", NA}.
collapse_block <- function(x) {
  x <- x[!is.na(x) & x != "NS"]
  if (length(x) == 0) return("NS")
  # If all equal (all Convergent OR all Divergent OR all Mixed)
  if (length(unique(x)) == 1) return(unique(x))
  # Otherwise there is heterogeneity across endpoints
  "Mixed"
}

# Build a compact, human-readable summary string from named block labels,
# skipping blocks that are "NS".
compose_summary <- function(named_blocks) {
  nb <- named_blocks[named_blocks != "NS"]
  if (length(nb) == 0) return("NS")
  # If all Convergent:
  if (all(nb == "Convergent")) return("Only Convergent")
  # If all Divergent:
  if (all(nb == "Divergent")) return("Only Divergent")
  # Otherwise, compose detailed summary in a stable order:
  order_vec <- c("Immune", "Microenvironment", "Phenotype", "Cox", "Survival")
  nb <- nb[intersect(order_vec, names(nb))]
  paste0(names(nb), ": ", nb, collapse = "; ")
}

## --- computation --------------------------------------------------------------
df016_final <- df015_organized %>%
  rowwise() %>%
  mutate(
    # Aggregate Cox across endpoints
    Cox_block = collapse_block(c_across(any_of(c(
      "Cox_concordance_OS", "Cox_concordance_DSS",
      "Cox_concordance_PFI","Cox_concordance_DFI"
    )))),
    # Aggregate Survival across endpoints
    Survival_block = collapse_block(c_across(any_of(c(
      "Survival_concordance_OS","Survival_concordance_DSS",
      "Survival_concordance_PFI","Survival_concordance_DFI"
    )))),
    # Single-field blocks
    Immune_block = coalesce(Immune_concordance, "NS"),
    # Microenvironment_block = coalesce(Microenv_concordance, "NS"),
    Phenotype_block = coalesce(Phenotypic_concordance, "NS"),
    # Compose final per-row summary
    Final_concordance_summary = compose_summary(c(
      Immune = Immune_block,
      # Microenvironment = Microenvironment_block,
      Phenotype = Phenotype_block,
      Cox = Cox_block,
      Survival = Survival_block
    ))
  ) %>%
  ungroup()

colnames(df016_final)

# Salvar resultado
saveRDS(df016_final, "Regulatory_circuitries_transcript_results.rds")

###### Combined_outcome Survival column 

# Define semantic categories
inconclusive_terms <- c("NS", "No data")

# All known biological terms
biological_terms <- c(
  "Low", "High", "WT", "MT", "Normal", "Deleted", "Duplicated",
  "Duplicated and Deleted", "Duplicated and Normal",
  "Deleted and Duplicated", "Normal and Duplicated",
  "Deleted and Normal", "Normal and Deleted"
)

# Function to classify a single row
classify_combined_outcome <- function(values) {
  values <- na.omit(values)
  values <- unique(trimws(values))
  values <- values[values != ""]  # Remove empty strings
  
  # Partition
  has_inconclusive <- any(values %in% inconclusive_terms)
  has_biological <- any(values %in% biological_terms)
  
  all_inconclusive <- all(values %in% inconclusive_terms)
  all_same <- length(unique(values[values %in% biological_terms])) == 1 && !has_inconclusive
  
  if (length(values) == 0 || all_inconclusive) {
    return("Meaningless")
  } else if (all_same) {
    return(values[values %in% biological_terms][1])  # Return the biological term
  } else if (has_biological && !has_inconclusive) {
    return("Meaningful Mixed")
  } else if (has_biological && has_inconclusive) {
    return(paste("Meaningful", unique(values[values %in% biological_terms][1])))
  } else {
    return("Undefined")
  }
}

# Apply to df006
df017_outcome_SMC <- df016_final %>%
  rowwise() %>%
  mutate(Combined_outcome_SMC = classify_combined_outcome(c_across(c(
    OS_worst_prognosis_group_interaction,
    DSS_worst_prognosis_group_interaction,
    DFI_worst_prognosis_group_interaction,
    PFI_worst_prognosis_group_interaction
  )))) %>%
  ungroup()


###### Combined_outcome Cox column 

# Define semantic categories
inconclusive_terms <- c("NS", "No data")

# All known biological terms
biological_terms <- c(
  "Protective", "Risky"
)

# Function to classify a single row
classify_combined_outcome <- function(values) {
  values <- na.omit(values)
  values <- unique(trimws(values))
  values <- values[values != ""]  # Remove empty strings
  
  # Partition
  has_inconclusive <- any(values %in% inconclusive_terms)
  has_biological <- any(values %in% biological_terms)
  
  all_inconclusive <- all(values %in% inconclusive_terms)
  all_same <- length(unique(values[values %in% biological_terms])) == 1 && !has_inconclusive
  
  if (length(values) == 0 || all_inconclusive) {
    return("Meaningless")
  } else if (all_same) {
    return(values[values %in% biological_terms][1])  # Return the biological term
  } else if (has_biological && !has_inconclusive) {
    return("Meaningful Mixed")
  } else if (has_biological && has_inconclusive) {
    return(paste("Meaningful", unique(values[values %in% biological_terms][1])))
  } else {
    return("Undefined")
  }
}

# Apply to df009
df018_outcome_HRC <- df017_outcome_SMC %>%
  rowwise() %>%
  mutate(Combined_outcome_HRC = classify_combined_outcome(c_across(c(
    Cox_OS_type_interaction,
    Cox_DSS_type_interaction,
    Cox_DFI_type_interaction,
    Cox_PFI_type_interaction
  )))) %>%
  ungroup()


# Salvar resultado
saveRDS(df018_outcome_HRC, "Regulatory_circuitries_transcript_with_HRC_SMC.rds")

# Step 1: Summarize validated interactions from df004_final_output_filtered
df019_interaction_summary <- df016_final %>%
  group_by(Signatures, Omic_layer_signature, Metabolism_signature, Pathways_signature, CTAB) %>%
  summarise(Meaningful_interaction = str_c(unique(Meaningful_interaction), collapse = " / "), .groups = "drop")

# Step 2: Rename to match if needed
df019_interaction_summary <- df019_interaction_summary %>%
  rename(Omic_layer = Omic_layer_signature,
         Metabolism = Metabolism_signature,
         Pathways = Pathways_signature)


# Step 3: Left join and fill NAs with "No meaningful interaction"
df020_final_result <- df011 %>%
  left_join(df019_interaction_summary, by = c("Signatures", "Omic_layer", "Metabolism", "Pathways", "CTAB")) %>%
  mutate(Meaningful_interaction = replace_na(Meaningful_interaction, "No meaningful interaction"))

# Reorder the columns
New_column_order_02 <- c("Members", "Signatures", "TranscriptGene", "Molecular_class", "Common_interaction", "Meaningful_interaction", "Metabolism", "Pathways", "Metabolic_cell_death", "CTAB", 
                         "Tumor_vs_normal",  "Tumor_vs_normal_p.adj", "Omic_layer", "Phenotypic_layer", "Correlation_rho", "Correlation_p.adj", "Cox_OS_type", "Cox_OS_p.value", 
                         "Cox_DSS_type", "Cox_DSS_p.value", "Cox_DFI_type", "Cox_DFI_p.value", "Cox_PFI_type", "Cox_PFI_p.value", "OS_worst_prognosis_group", "OS_p.value", 
                         "DSS_worst_prognosis_group", "DSS_p.value", "DFI_worst_prognosis_group", "DFI_p.value", "PFI_worst_prognosis_group", "PFI_p.value",
                         "Microenvironment_classification", "Microenvironment_score_details", "Immune_classification", "Immune_score_details")

df020_final_result <- df020_final_result[New_column_order_02]

# Salvar resultado
saveRDS(df020_final_result, "df_transcript_signatures.rds")


# # 📂 Carregar dados
# df007 <- import("/Cancer_metabolism_project_02/20- Checking_signature_results/df_transcript_signatures.rds")
# df008 <- import("/Cancer_metabolism_project_02/11- Checking_results_and_annotation/3- Checking_results_and_annotation/df063_gene_protein_mirna_lncrna_results.rds")
# 
# # 🧬 Métricas clínicas
# cox_cols <- c("Cox_OS_type", "Cox_DSS_type", "Cox_PFI_type", "Cox_DFI_type")
# surv_cols <- c("OS_worst_prognosis_group", "DSS_worst_prognosis_group", "DFI_worst_prognosis_group", "PFI_worst_prognosis_group")
# clinical_cols <- c(cox_cols, surv_cols)
# 
# # 🧼 Filtragem e limpeza
# df007_filtered <- df007 %>%
#   filter(if_any(all_of(clinical_cols), ~ . != "NS" & . != "No data"))
# 
# df008_filtered <- df008 %>%
#   filter(if_any(all_of(clinical_cols), ~ . != "NS" & . != "No data"))
# 
# # 🔄 Inicialização
# total_rows <- nrow(df007_filtered)
# pb <- txtProgressBar(min = 0, max = total_rows, style = 3)
# result_list <- list()
# 
# # 🔁 Loop principal
# for (i in seq_len(total_rows)) {
#   row <- df007_filtered[i, ]
#   interactions <- str_split(row$Common_interaction, " / ", simplify = TRUE) %>% str_trim()
#   
#   for (entity in interactions) {
#     match_row <- df008_filtered %>%
#       filter(Target == entity, Cancer_types == row$Cancer_types)
#     
#     if (nrow(match_row) > 0) {
#       match_row <- match_row[1, ]  # usar apenas a primeira correspondência
#       
#       # Verificar métricas válidas compartilhadas
#       valid_shared_metrics <- clinical_cols[sapply(clinical_cols, function(col) {
#         row_val <- row[[col]]
#         match_val <- match_row[[col]]
#         !(row_val %in% c("NS", "No data") | match_val %in% c("NS", "No data")) &&
#           row_val != "" && match_val != ""
#       })]
#       
#       if (length(valid_shared_metrics) > 0) {
#         
#         # Função para classificar concordância
#         classify_concordance <- function(cols) {
#           shared <- intersect(valid_shared_metrics, cols)
#           if (length(shared) == 0) return(list(metrics = NA_character_, category = "NS"))
#           vals <- sapply(shared, function(col) {
#             if (row[[col]] == match_row[[col]]) "equal" else "diff"
#           })
#           metrics_str <- shared %>%
#             str_extract("(OS|DSS|PFI|DFI)") %>%
#             unique() %>%
#             paste(collapse = " / ")
#           eq <- sum(vals == "equal")
#           diff <- sum(vals == "diff")
#           concordance <- case_when(
#             eq > diff ~ "Convergent",
#             diff > eq ~ "Divergent",
#             eq == diff ~ "Mixed"
#           )
#           list(metrics = metrics_str, category = concordance)
#         }
#         
#         cox_result <- classify_concordance(cox_cols)
#         surv_result <- classify_concordance(surv_cols)
#         
#         # Concordância imune
#         immune_sig <- row$Immune_classification
#         immune_ent <- match_row$Immune_classification
#         immune_concordance <- if (immune_sig %in% c("NS", "No data") | immune_ent %in% c("NS", "No data")) {
#           "NS"
#         } else if (immune_sig == immune_ent) {
#           "Convergent"
#         } else {
#           "Divergent"
#         }
#         
#         # Concordância fenotípica
#         phen_sig <- row$Phenotypic_layer
#         phen_ent <- match_row$Phenotypic_layer
#         phenotypic_concordance <- if (phen_sig %in% c("NS", "No data") | phen_ent %in% c("NS", "No data")) {
#           "NS"
#         } else if (phen_sig == phen_ent) {
#           "Convergent"
#         } else {
#           "Divergent"
#         }
#         
#         # Criar linha de saída
#         output_row <- tibble(
#           Multiomics_Signature = row$Multiomics_Signature,
#           Molecular_class_Signature = row$Molecular_class,
#           Metabolism_Signature = row$Metabolism,
#           Pathways_Signature = row$Pathways,
#           Metabolic_cell_death_Signature = row$Metabolic_cell_death,
#           Omic_layer_signature = row$Omic_layer,
#           Phenotypic_layer_signature = phen_sig,
#           Cancer_types = row$Cancer_types,
#           Common_interaction = row$Common_interaction,
#           Cox_OS_type_signature = row$Cox_OS_type,
#           Cox_DSS_type_signature = row$Cox_DSS_type,
#           Cox_PFI_type_signature = row$Cox_PFI_type,
#           Cox_DFI_type_signature = row$Cox_DFI_type,
#           OS_worst_prognosis_group_signature = row$OS_worst_prognosis_group,
#           DSS_worst_prognosis_group_signature = row$DSS_worst_prognosis_group,
#           DFI_worst_prognosis_group_signature = row$DFI_worst_prognosis_group,
#           PFI_worst_prognosis_group_signature = row$PFI_worst_prognosis_group,
#           Microenvironment_classification_signature = row$Microenvironment_classification,
#           Immune_classification_signature = immune_sig,
#           Interaction_validated = entity,
#           Molecular_class_entity = match_row$Molecular_class,
#           Metabolism_entity = match_row$Metabolism,
#           Pathways_entity = match_row$Pathways,
#           Metabolic_cell_death_entity = match_row$Metabolic_cell_death,
#           Phenotypic_layer_entity = phen_ent,
#           Omic_layer_entity = match_row$Omic_layer,
#           Cox_OS_type_entity = match_row$Cox_OS_type,
#           Cox_DSS_type_entity = match_row$Cox_DSS_type,
#           Cox_PFI_type_entity = match_row$Cox_PFI_type,
#           Cox_DFI_type_entity = match_row$Cox_DFI_type,
#           OS_worst_prognosis_group_entity = match_row$OS_worst_prognosis_group,
#           DSS_worst_prognosis_group_entity = match_row$DSS_worst_prognosis_group,
#           DFI_worst_prognosis_group_entity = match_row$DFI_worst_prognosis_group,
#           PFI_worst_prognosis_group_entity = match_row$PFI_worst_prognosis_group,
#           Microenvironment_classification_entity = match_row$Microenvironment_classification,
#           Immune_classification_entity = immune_ent,
#           Phenotypic_concordance = phenotypic_concordance,
#           Cox_results = cox_result$metrics,
#           Clinical_cox_concordance = cox_result$category,
#           Survival_results = surv_result$metrics,
#           Clinical_survival_concordance = surv_result$category,
#           Immune_concordance = immune_concordance
#         )
#         
#         result_list[[length(result_list) + 1]] <- output_row
#       }
#     }
#   }
#   
#   setTxtProgressBar(pb, i)
# }
# close(pb)
# 
# # 💾 Output final
# df009_final_output <- bind_rows(result_list)
# df010_final_output_filtered <- df009_final_output %>% distinct()
# 
# 
# ###### Combined_outcome Survival column 
# 
# # Define semantic categories
# inconclusive_terms <- c("NS", "No data")
# 
# # All known biological terms
# biological_terms <- c(
#   "Low", "High", "WT", "MT", "Normal", "Deleted", "Duplicated",
#   "Duplicated and Deleted", "Duplicated and Normal",
#   "Deleted and Duplicated", "Normal and Duplicated",
#   "Deleted and Normal", "Normal and Deleted"
# )
# 
# # Function to classify a single row
# classify_combined_outcome <- function(values) {
#   values <- na.omit(values)
#   values <- unique(trimws(values))
#   values <- values[values != ""]  # Remove empty strings
#   
#   # Partition
#   has_inconclusive <- any(values %in% inconclusive_terms)
#   has_biological <- any(values %in% biological_terms)
#   
#   all_inconclusive <- all(values %in% inconclusive_terms)
#   all_same <- length(unique(values[values %in% biological_terms])) == 1 && !has_inconclusive
#   
#   if (length(values) == 0 || all_inconclusive) {
#     return("Meaningless")
#   } else if (all_same) {
#     return(values[values %in% biological_terms][1])  # Return the biological term
#   } else if (has_biological && !has_inconclusive) {
#     return("Meaningful Mixed")
#   } else if (has_biological && has_inconclusive) {
#     return(paste("Meaningful", unique(values[values %in% biological_terms][1])))
#   } else {
#     return("Undefined")
#   }
# }
# 
# # Apply to df004
# df011_outcome_SMC <- df010_final_output_filtered %>%
#   rowwise() %>%
#   mutate(Combined_outcome_SMC = classify_combined_outcome(c_across(c(
#     OS_worst_prognosis_group_entity,
#     DSS_worst_prognosis_group_entity,
#     DFI_worst_prognosis_group_entity,
#     PFI_worst_prognosis_group_entity
#   )))) %>%
#   ungroup()
# 
# 
# 
# ###### Combined_outcome Cox column 
# 
# # Define semantic categories
# inconclusive_terms <- c("NS", "No data")
# 
# # All known biological terms
# biological_terms <- c(
#   "Protective", "Risky"
# )
# 
# # Function to classify a single row
# classify_combined_outcome <- function(values) {
#   values <- na.omit(values)
#   values <- unique(trimws(values))
#   values <- values[values != ""]  # Remove empty strings
#   
#   # Partition
#   has_inconclusive <- any(values %in% inconclusive_terms)
#   has_biological <- any(values %in% biological_terms)
#   
#   all_inconclusive <- all(values %in% inconclusive_terms)
#   all_same <- length(unique(values[values %in% biological_terms])) == 1 && !has_inconclusive
#   
#   if (length(values) == 0 || all_inconclusive) {
#     return("Meaningless")
#   } else if (all_same) {
#     return(values[values %in% biological_terms][1])  # Return the biological term
#   } else if (has_biological && !has_inconclusive) {
#     return("Meaningful Mixed")
#   } else if (has_biological && has_inconclusive) {
#     return(paste("Meaningful", unique(values[values %in% biological_terms][1])))
#   } else {
#     return("Undefined")
#   }
# }
# 
# # Apply to df009
# df012_outcome_HRC <- df011_outcome_SMC %>%
#   rowwise() %>%
#   mutate(Combined_outcome_HRC = classify_combined_outcome(c_across(c(
#     Cox_OS_type_entity,
#     Cox_DSS_type_entity,
#     Cox_DFI_type_entity,
#     Cox_PFI_type_entity
#   )))) %>%
#   ungroup()
# 
# 
# # Salvar resultado
# export(df012_outcome_HRC, "Interaction_results_transcript.tsv")
# 
# # Step 1: Summarize validated interactions from df_003_final_output_filtered
# df011_interaction_summary <- df010_final_output_filtered %>%
#   group_by(Multiomics_Signature, Omic_layer_signature, Metabolism_Signature, Pathways_Signature, Cancer_types) %>%
#   summarise(Meaningful_interaction = str_c(unique(Interaction_validated), collapse = " / "), .groups = "drop")
# 
# # Step 2: Rename to match if needed
# df011_interaction_summary <- df011_interaction_summary %>%
#   rename(Omic_layer = Omic_layer_signature,
#          Metabolism = Metabolism_Signature,
#          Pathways = Pathways_Signature)
# 
# # Step 3: Left join and fill NAs with "No meaningful interaction"
# df012_final_result <- df007 %>%
#   left_join(df011_interaction_summary, by = c("Multiomics_Signature", "Omic_layer", "Metabolism", "Pathways", "Cancer_types")) %>%
#   mutate(Meaningful_interaction = replace_na(Meaningful_interaction, "No meaningful interaction"))
# 
# # Reorder the columns
# New_column_order_02 <- c("Members", "Multiomics_Signature", "TranscriptGene", "Molecular_class", "Common_interaction", "Meaningful_interaction", "Metabolism", "Pathways", "Metabolic_cell_death", "Cancer_types", 
#                          "Tumor_vs_normal",  "Tumor_vs_normal_p.adj", "Omic_layer", "Phenotypic_layer", "Correlation_rho", "Correlation_p.adj", "Cox_OS_type", "Cox_OS_p.value", 
#                          "Cox_DSS_type", "Cox_DSS_p.value", "Cox_DFI_type", "Cox_DFI_p.value", "Cox_PFI_type", "Cox_PFI_p.value", "OS_worst_prognosis_group", "OS_p.value", 
#                          "DSS_worst_prognosis_group", "DSS_p.value", "DFI_worst_prognosis_group", "DFI_p.value", "PFI_worst_prognosis_group", "PFI_p.value",
#                          "Microenvironment_classification", "Microenvironment_score_details", "Immune_classification", "Immune_score_details")
# 
# df012_final_result <- df012_final_result[New_column_order_02]
# 
# # Salvar resultado
# saveRDS(df012_final_result, "df_transcript_signatures.rds")



# 
# ###### Network of Validated interactions
# 
# library(dplyr)
# library(tidygraph)
# library(ggraph)
# library(igraph)
# 
# # 📌 Filter data for PAAD only
# df_paad <- df004_final_output_filtered %>%
#   filter(Cancer_types == "PAAD")
# 
# # 🔝 Identify top 3 signatures with the most validated interactions
# top_signatures <- df_paad %>%
#   count(Multiomics_Signature, name = "n_interactions") %>%
#   arrange(desc(n_interactions)) %>%
#   slice_head(n = 3) %>%
#   pull(Multiomics_Signature)
# 
# # 🧼 Filter df_paad to keep only interactions from top 3 signatures
# df_paad_top <- df_paad %>%
#   filter(Multiomics_Signature %in% top_signatures)
# 
# # ✅ Prepare edge list
# edges <- df_paad_top %>%
#   select(from = Multiomics_Signature,
#          to = Interaction_validated,
#          Omic_layer_signature,
#          Cancer_types) %>%
#   distinct()
# 
# # ✅ Create node list
# nodes <- unique(c(edges$from, edges$to)) %>%
#   as.data.frame() %>%
#   rename(name = ".")
# 
# # ✅ Add node type (Signature or Target)
# nodes <- nodes %>%
#   mutate(type = ifelse(name %in% df_paad_top$Multiomics_Signature, 
#                        "Signature", 
#                        "Target"))
# 
# # ✅ Create graph object
# graph <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
# 
# # ✅ Plot network
# ggraph(graph, layout = "fr") +
#   geom_edge_link(aes(color = Omic_layer_signature), alpha = 0.6) +
#   geom_node_point(aes(color = type), size = 4) +
#   geom_node_text(aes(label = name), repel = TRUE, size = 3) +
#   scale_color_manual(values = c("Signature" = "darkblue", "Target" = "orange")) +
#   theme_void() +
#   labs(title = "Top 3 Multiomics Signatures in PAAD",
#        subtitle = "Network of Validated Interactions")


# =========================== PACKAGES ========================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(tibble)
  library(purrr)
  library(igraph); library(tidygraph); library(ggraph); library(ggplot2)
  library(grid)
})

# =========================== INPUT ===========================================
df <- df006_final  # your real data.frame

# Select one regulatory circuitry by its Nomenclature label:
target_nomenclature <- "PUT_YOUR_NOMENCLATURE_HERE"  # <-- change this string

# Safety checks
stopifnot("Nomenclature" %in% names(df))

matches <- df %>% dplyr::filter(Nomenclature == target_nomenclature)

if (nrow(matches) == 0) {
  stop(paste0("No row found for Nomenclature = '", target_nomenclature, "'"))
}
if (nrow(matches) > 1) {
  message("⚠️ Multiple rows match this Nomenclature — using the first.")
}

x <- matches[1, , drop = FALSE]  # this replaces the old df[row_id, ]

# df <- df006_final  # your real data.frame
# row_id <- 1L
# stopifnot(row_id >= 1L, row_id <= nrow(df))
# x <- df[row_id, , drop = FALSE]

# # Case label (use your current columns)
# case_label <- x$Signatures

# =========================== HELPERS =========================================
norm_cd <- function(v) {
  v <- trimws(tolower(as.character(v)))
  v <- ifelse(v %in% c("convergent","divergent","ns"), v, NA_character_)
  factor(v, levels = c("convergent","divergent","ns"))
}

mk_edge <- function(from, to, type, concordance = NA_character_) {
  if (is.na(from) || is.na(to) || from == "" || to == "") return(NULL)
  # ignore explicit NS
  if (!is.na(concordance) && identical(concordance, "ns")) return(NULL)
  tibble(from = from, to = to, type = type, concordance = concordance)
}

split_multi <- function(x) {
  if (is.null(x) || all(is.na(x))) return(character(0))
  x %>%
    as.character() %>%
    str_replace_all("[(){}\\[\\]]", "") %>%
    str_split("\\s*[+,;/]\\s*") %>%
    pluck(1) %>%
    discard(~ .x == "" | is.na(.x)) %>%
    unique()
}

pick_immune_node <- function(label) {
  lab <- tolower(trimws(as.character(label)))
  if (lab %in% c("hot"))      return("Hot immune profile")
  if (lab %in% c("cold")) return("Cold immune profile")
  if (lab %in% c("variable"))            return("Variable immune profile")
  return(NA_character_)
}

# ---------- Phenotype-from-RHO (returns node label + sign direction)
phenotype_from_rho <- function(phenotype_label, rho) {
  if (is.na(rho)) return(list(label = NA_character_, dir = NA_character_))
  rho_num <- suppressWarnings(as.numeric(rho))
  if (is.na(rho_num)) return(list(label = NA_character_, dir = NA_character_))
  plab <- tolower(trimws(as.character(phenotype_label)))
  if (plab == "" || is.na(plab)) return(list(label = NA_character_, dir = NA_character_))
  
  dir <- ifelse(rho_num > 0, "pos", ifelse(rho_num < 0, "neg", NA_character_))
  
  # Stemness
  if (str_detect(plab, "stemness")) {
    lab <- if (rho_num > 0) "Lower stemness" else if (rho_num < 0) "Higher stemness" else NA_character_
    return(list(label = lab, dir = dir))
  }
  
  # TMB
  if (str_detect(plab, "tumor mutational burden")) {
    lab <- if (rho_num > 0) "Higher TMB" else if (rho_num < 0) "Lower TMB" else NA_character_
    return(list(label = lab, dir = dir))
  }
  
  # MSI
  if (str_detect(plab, "microsatellite instability")) {
    lab <- if (rho_num > 0) "Higher MSI" else if (rho_num < 0) "Lower MSI" else NA_character_
    return(list(label = lab, dir = dir))
  }
  
  list(label = NA_character_, dir = NA_character_)
}

# ---------- Prognosis from keywords (“Risk”/“Protective”)
prog_node_from_keywords <- function(text) {
  if (is.null(text) || all(is.na(text))) return(NA_character_)
  tx <- tolower(paste(na.omit(as.character(text)), collapse = " | "))
  if (tx == "") return(NA_character_)
  if (str_detect(tx, "\\brisk\\b"))       return("Worse prognosis")
  if (str_detect(tx, "\\bprotective\\b")) return("Favorable prognosis")
  return(NA_character_)
}

# Collect candidate text for each side
collect_signature_text <- function(row) {
  cols <- grep("signature", names(row), ignore.case = TRUE, value = TRUE)
  cols <- union(cols, c("Signatures"))
  unlist(row[ , intersect(cols, names(row)), drop = TRUE])
}
collect_regulator_text <- function(row) {
  cols <- grep("interaction", names(row), ignore.case = TRUE, value = TRUE)
  cols <- union(cols, c("Meaningful_interaction"))
  unlist(row[ , intersect(cols, names(row)), drop = TRUE])
}

# ---------- Robust string compare for pathways
.norm_str <- function(z) {
  z <- as.character(z)
  z <- stringr::str_squish(tolower(z))
  ifelse(is.na(z) | z == "" | z == "na", NA_character_, z)
}
same_or_diff <- function(a, b) {
  aa <- .norm_str(a); bb <- .norm_str(b)
  if (is.na(aa) || is.na(bb)) return(NA_character_)
  if (identical(aa, bb)) "convergent" else "divergent"
}

# ======================= CANONICAL LABELS / DYNAMIC MAPPING ==================
NODE_SIG <- x$Nomenclature

# Immune nodes (signature/entity)
IMM_sig_node <- pick_immune_node(x$Immune_classification_signature)
IMM_ent_node <- pick_immune_node(x$Immune_classification_interaction)

# Pathway labels (signature/entity) – fallback to generic if blank (for node label)
met_path_sig <- if (!is.na(x$Pathways_signature) && nzchar(x$Pathways_signature)) as.character(x$Pathways_signature) else "Metabolic pathway"
met_path_ent <- if (!is.na(x$Pathways_interaction) && nzchar(x$Pathways_interaction)) as.character(x$Pathways_interaction) else met_path_sig

# Concordance for pathways based on RAW equality
raw_path_sig <- x$Pathways_signature
raw_path_ent <- x$Pathways_interaction
path_conc    <- same_or_diff(raw_path_sig, raw_path_ent)  # "convergent"/"divergent"/NA

# Phenotype from RHO (signature & interaction), using the signature phenotype label
PH_sig <- phenotype_from_rho(x$Phenotypic_layer_signature, x$rho_signature)
PH_ent <- phenotype_from_rho(x$Phenotypic_layer_signature, x$rho_interaction)

# Concordance for phenotype from RHO sign
phen_conc <- NA_character_
if (!is.na(PH_sig$dir) && !is.na(PH_ent$dir)) {
  phen_conc <- if (PH_sig$dir == PH_ent$dir) "convergent" else "divergent"
}

# Prognosis from keywords per side
prog_sig_node <- prog_node_from_keywords(collect_signature_text(x))
prog_ent_node <- prog_node_from_keywords(collect_regulator_text(x))
prog_conc <- NA_character_
if (!is.na(prog_sig_node) && !is.na(prog_ent_node)) {
  prog_conc <- if (prog_sig_node == prog_ent_node) "convergent" else "divergent"
}

# Concordances from Cox (explicit)
Cox_OS_cd  <- norm_cd(x$Cox_concordance_OS)
Cox_DSS_cd <- norm_cd(x$Cox_concordance_DSS)
Cox_PFI_cd <- norm_cd(x$Cox_concordance_PFI)
Cox_DFI_cd <- norm_cd(x$Cox_concordance_DFI)

# ======================= REGULATOR NODES (REAL NAMES) ========================
regulators <- split_multi(x$Meaningful_interaction)
if (!length(regulators)) regulators <- "(no regulator)"

# ============================= EDGES =========================================
edges <- list()

# Unified association type label
ASSOC <- "Signature and Interaction association"

# (A) Regulation: each regulator -> signature
edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, NODE_SIG, "Regulation", NA_character_)))

# (B) Signature -> pathway (association) with pathway interpretation
if (!is.na(met_path_sig) && nzchar(met_path_sig)) {
  edges <- append(edges, list(mk_edge(NODE_SIG, met_path_sig, ASSOC, path_conc)))
}

# (C) Regulator -> pathway (association) with same interpretation
if (!is.na(met_path_ent) && nzchar(met_path_ent)) {
  edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, met_path_ent, ASSOC, path_conc)))
}

# (D) Phenotype from RHO (signature & interaction) with concordance by sign equality
if (!is.na(PH_sig$label))
  edges <- append(edges, list(mk_edge(NODE_SIG, PH_sig$label, ASSOC, phen_conc)))
if (!is.na(PH_ent$label))
  edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, PH_ent$label, ASSOC, phen_conc)))

# (E) Immune classification + concordance edges (ignore NS)
if (!is.na(IMM_sig_node))
  edges <- append(edges, list(mk_edge(NODE_SIG, IMM_sig_node, ASSOC, as.character(norm_cd(x$Immune_concordance)))))
if (!is.na(IMM_ent_node))
  edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, IMM_ent_node, ASSOC, as.character(norm_cd(x$Immune_concordance)))))

# (F) Prognosis from keywords (“Risk”/“Protective”) – concordance by equality
if (!is.na(prog_sig_node))
  edges <- append(edges, list(mk_edge(NODE_SIG, prog_sig_node, ASSOC, prog_conc)))
if (!is.na(prog_ent_node))
  edges <- append(edges, lapply(regulators, function(rg) mk_edge(rg, prog_ent_node, ASSOC, prog_conc)))

# (G) Cox-based prognosis (explicit convergent/divergent; NS ignored)
prog_node_from_cd <- function(cd) {
  if (is.na(cd)) return(NA_character_)
  if (cd == "convergent") return("Favorable prognosis")
  if (cd == "divergent")  return("Worse prognosis")
  return(NA_character_)
}
edges <- append(edges, list(
  mk_edge(NODE_SIG, prog_node_from_cd(Cox_OS_cd),  ASSOC, as.character(Cox_OS_cd)),
  mk_edge(NODE_SIG, prog_node_from_cd(Cox_DSS_cd), ASSOC, as.character(Cox_DSS_cd)),
  mk_edge(NODE_SIG, prog_node_from_cd(Cox_PFI_cd), ASSOC, as.character(Cox_PFI_cd)),
  mk_edge(NODE_SIG, prog_node_from_cd(Cox_DFI_cd), ASSOC, as.character(Cox_DFI_cd))
))

edges <- bind_rows(compact(edges)) %>%
  filter(!is.na(to), !is.na(from)) %>%
  mutate(
    type        = factor(type, levels = c("Regulation", ASSOC)),
    concordance = factor(concordance, levels = c("convergent","divergent"))
  )

# Use a separate plotting colour field: regulation gets a fixed label, associations keep concordance
edges_for_plot <- edges %>%
  mutate(conc_plot = ifelse(as.character(type) == "Regulation",
                            "regulation",
                            as.character(concordance)))


# ============================= NODES =========================================
layer_of <- function(name) {
  if (name %in% regulators)                                return("Interaction")
  if (name %in% c(NODE_SIG))                               return("Signature")
  if (grepl("stemness|tmb|msi", tolower(name)))            return("Tumor phenotype")
  if (grepl("immune profile", tolower(name)))              return("Immune phenotype")
  if (grepl("microenvironment", tolower(name)))            return("Immune phenotype")
  if (name %in% c("Favorable prognosis","Worse prognosis"))return("Prognosis")
  return("Metabolism") # pathways/metabolism default
}

nodes <- tibble(name = unique(c(edges$from, edges$to))) %>%
  mutate(layer = vapply(name, layer_of, character(1)))

# ============================ AESTHETICS =====================================
layer_colors <- c(
  "Interaction"      = "#E67E22",
  "Signature"        = "#2ECC71",
  "Metabolism"       = "#5DADE2",
  "Immune phenotype" = "gold",
  "Tumor phenotype"  = "#8E44AD",
  "Prognosis"        = "#C0392B"
)

concordance_cols <- c(
  "convergent" = "#27AE60",  # supports/reinforces
  "divergent"  = "#C0392B"   # opposes/modulates
)

etype_lty <- c(
  "Regulation" = "dashed",
  "Signature and Interaction association" = "solid"
)

etype_w <- c(
  "Regulation" = 1.2,
  "Signature and Interaction association" = 1.2
)


# ================================ GRAPH ======================================
g <- tidygraph::tbl_graph(nodes = nodes, edges = edges_for_plot, directed = TRUE)

# Color map includes 'regulation' (hidden from legend via breaks)
edge_colours <- c(
  "convergent" = "#27AE60",
  "divergent"  = "#C0392B",
  "regulation" = "#E67E22"   # fixed colour for regulation edges
)

# Two relation types: dashed for Regulation; solid for the unified association
etype_lty <- c(
  "Regulation"                          = "dashed",
  "Signature and Interaction association" = "solid"
)
etype_w <- c(
  "Regulation"                          = 1.2,
  "Signature and Interaction association" = 1.2
)

set.seed(123)
p <- ggraph::ggraph(g, layout = "circle") +
  # Single edge layer; color mapped to conc_plot
  ggraph::geom_edge_arc(
    ggplot2::aes(edge_colour = conc_plot, edge_linetype = type, edge_width = type),
    arrow      = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
    strength   = 0.8,
    edge_alpha = 0.9
  ) +
  # Nodes
  ggraph::geom_node_point(ggplot2::aes(fill = layer), size = 7, shape = 21, color = "white") +
  ggraph::geom_node_text(ggplot2::aes(label = name), vjust = -1.5, size = 3, fontface = "bold") +
  
  # Scales
  ggplot2::scale_fill_manual(values = layer_colors, name = "Dimension") +
  ggraph::scale_edge_colour_manual(
    values       = edge_colours,
    name         = "Interpretation",
    breaks       = c("convergent","divergent"),  # hide 'regulation' from legend
    na.translate = FALSE,                        # drop NA from legend
    na.value     = "grey70"                      # draw association edges lacking interpretation
  ) +
  ggraph::scale_edge_width_manual(values = etype_w, guide = "none") +
  ggraph::scale_edge_linetype_manual(values = etype_lty, name = "Relation type") +
  
  # Theme
  ggplot2::ggtitle(paste0("Regulatory Circuitry — Signature ↔ Regulator")) +
  ggplot2::theme_void() +
  ggplot2::coord_equal(clip = "off") +
  ggplot2::theme(
    legend.position = "left",
    legend.box.margin = ggplot2::margin(r = 30, l = 30, t = 20, b = 20),
    plot.title = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5, vjust = 1.5),
    plot.title.position = "panel",
    plot.margin = ggplot2::margin(t = 20, r = 40, b = 20, l = 55)
  )

print(p)


# ============================== EXPORT (PDF only) ============================
outfile_pdf <- sprintf("Regulatory_Circuitry_%s.pdf", gsub("[^A-Za-z0-9]+","_", case_label))

ggplot2::ggsave(
  filename = outfile_pdf,
  plot     = p,
  width    = 16,
  height   = 9,
  units    = "in",
  device   = cairo_pdf,   # keeps text + vectors editable
  bg       = "white"
)

message("Saved PDF:", outfile_pdf)
