##### RSCRIPT TO SURVIVAL ANALYSIS
##### AUTHOR: Higor Almeida Cordeiro Nogueira
##### UPDATED: 20/09/2024

# Load required libraries
library("rio")
library("UCSCXenaShiny")
library("dplyr")
library("survival")

# Set directory
setwd("F:/Higor/Cancer_metabolism_analysis_06/8- Survival_Analysis/lncRNA/CNV/")

# Import table
table <- import("F:/Higor/Cancer_metabolism_analysis_06/7- Cox_Analysis/lncRNA/CNV/Cox_results.tsv")

# Initialize new columns for each survival metric
survival_metrics <- c("OS", "DSS", "DFI", "PFI")
for (metric in survival_metrics) {
  table[[paste0(metric, "_log_rank_chisq")]] <- NA
  table[[paste0(metric, "_p_val")]] <- NA
  table[[paste0(metric, "_worst_prognosis_group")]] <- NA
}

# Function to compare survival groups
compare_groups <- function(data, group1, group2, metric) {
  data_filtered <- data %>% filter(group %in% c(group1, group2))
  if (nrow(data_filtered) == 0) return(NA)
  
  surv_diff_pair <- tryCatch({
    survdiff(Surv(data_filtered[[paste0(metric, ".time")]], data_filtered[[metric]]) ~ group, data = data_filtered)
  }, error = function(e) {
    warning(paste("Failed to calculate survival differences for", group1, "vs", group2, ":", e$message))
    return(NULL)
  })
  
  if (!is.null(surv_diff_pair) && is.list(surv_diff_pair) && !is.null(surv_diff_pair$chisq)) {
    return(1 - pchisq(surv_diff_pair$chisq, 1))
  } else {
    return(NA)
  }
}

# Function to perform survival analysis
perform_survival_analysis <- function(genes, cancer_type, metric) {
  cat("Running survival analysis for gene:", genes, "and cancer_type:", cancer_type, "on metric:", metric, "\n")
  
  # Load survival data with error handling
  data_surv <- tryCatch({
    tcga_surv_get(
      genes,   
      TCGA_cohort = cancer_type,
      profile = "cnv",
      TCGA_cli_data = dplyr::inner_join(
        load_data("tcga_clinical"), 
        load_data("tcga_surv"),  
        by = "sample"
      )
    )
  }, error = function(e) {
    warning(paste("Failed to load survival data for gene:", genes, "and cancer_type:", cancer_type, ":", e$message))
    return(NULL)
  })
  
  if (is.null(data_surv)) {
    warning(paste("No survival data for gene:", genes, "and cancer_type:", cancer_type))
    return(NULL)
  }
  
  # Create survival plot
  plot_surv <- tryCatch({
    tcga_surv_plot(
      data_surv,
      time = paste0(metric, ".time"),
      status = metric,
      profile = "cnv",
      palette = "aaas"
    )
  }, error = function(e) {
    warning(paste("Failed to create survival plot for", genes, "in", cancer_type, ":", e$message))
    return(NULL)
  })
  
  # Process data for survival analysis
  data_surv <- data_surv %>%
    mutate(group = case_when(
      .data$value == 0 ~ "Normal",
      .data$value > 0 ~ "Duplicated",
      .data$value < 0 ~ "Deleted",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(group))
  
  # Calculate log-rank test
  surv_diff <- tryCatch({
    survdiff(Surv(data_surv[[paste0(metric, ".time")]], data_surv[[metric]]) ~ group, data = data_surv)
  }, error = function(e) {
    warning(paste("Failed to calculate survival differences for", metric, "on gene:", genes, "in cancer_type:", cancer_type, ":", e$message))
    return(NULL)
  })
  
  if (is.null(surv_diff) || !is.list(surv_diff) || is.null(surv_diff$chisq)) {
    warning(paste("survdiff failed or returned invalid data for", metric, "on gene", genes, "in cancer_type", cancer_type))
    return(list(chisq = NA, p_val = NA, worst_prognosis_group = "NS"))
  }
  
  # Calculate p-value from log-rank test
  p_val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  chisq_log_rank <- surv_diff$chisq
  
  if (p_val >= 0.05) {
    return(list(chisq = chisq_log_rank, p_val = p_val, worst_prognosis_group = "NS"))
  }
  
  # Summarize survival times across groups, handling missing thresholds properly
  group_surv_multiple <- tryCatch({
    plot_surv[["plot"]][["data"]] %>%
      group_by(group) %>%
      summarize(
        time_40 = if (any(surv <= 0.4, na.rm = TRUE)) min(time[surv <= 0.4], na.rm = TRUE) else NA,
        time_50 = if (any(surv <= 0.5, na.rm = TRUE)) min(time[surv <= 0.5], na.rm = TRUE) else NA,
        time_60 = if (any(surv <= 0.6, na.rm = TRUE)) min(time[surv <= 0.6], na.rm = TRUE) else NA,
        time_70 = if (any(surv <= 0.7, na.rm = TRUE)) min(time[surv <= 0.7], na.rm = TRUE) else NA,
        time_80 = if (any(surv <= 0.8, na.rm = TRUE)) min(time[surv <= 0.8], na.rm = TRUE) else NA
      )
  }, error = function(e) {
    warning(paste("Failed to summarize survival data for", metric, ":", e$message))
    return(NULL)
  })
  
  # Check if group_surv_multiple is NULL or has valid survival times for at least one group
  if (is.null(group_surv_multiple) || all(is.na(group_surv_multiple %>% select(starts_with("time"))))) {
    return(list(chisq = chisq_log_rank, p_val = p_val, worst_prognosis_group = "NS"))
  }
  
  # Get the specific groups present in the data
  present_groups <- unique(group_surv_multiple$group)
  
  # Logic to handle cases with two or three groups
  if (length(present_groups) == 2) {
    # Two-group case: Determine which specific groups are present and assign names explicitly
    if (all(c("Normal", "Deleted") %in% present_groups)) {
      # Compare Normal and Deleted
      normal_wins <- sum(
        (group_surv_multiple %>% filter(group == "Normal") %>% select(starts_with("time"))) <
          (group_surv_multiple %>% filter(group == "Deleted") %>% select(starts_with("time"))),
        na.rm = TRUE
      )
      
      deleted_wins <- sum(
        (group_surv_multiple %>% filter(group == "Deleted") %>% select(starts_with("time"))) <
          (group_surv_multiple %>% filter(group == "Normal") %>% select(starts_with("time"))),
        na.rm = TRUE
      )
      
      worst_prognosis_group <- if (normal_wins > deleted_wins) "Normal" else "Deleted"
      
    } else if (all(c("Normal", "Duplicated") %in% present_groups)) {
      # Compare Normal and Duplicated
      normal_wins <- sum(
        (group_surv_multiple %>% filter(group == "Normal") %>% select(starts_with("time"))) <
          (group_surv_multiple %>% filter(group == "Duplicated") %>% select(starts_with("time"))),
        na.rm = TRUE
      )
      
      duplicated_wins <- sum(
        (group_surv_multiple %>% filter(group == "Duplicated") %>% select(starts_with("time"))) <
          (group_surv_multiple %>% filter(group == "Normal") %>% select(starts_with("time"))),
        na.rm = TRUE
      )
      
      worst_prognosis_group <- if (normal_wins > duplicated_wins) "Normal" else "Duplicated"
      
    } else if (all(c("Deleted", "Duplicated") %in% present_groups)) {
      # Compare Deleted and Duplicated
      deleted_wins <- sum(
        (group_surv_multiple %>% filter(group == "Deleted") %>% select(starts_with("time"))) <
          (group_surv_multiple %>% filter(group == "Duplicated") %>% select(starts_with("time"))),
        na.rm = TRUE
      )
      
      duplicated_wins <- sum(
        (group_surv_multiple %>% filter(group == "Duplicated") %>% select(starts_with("time"))) <
          (group_surv_multiple %>% filter(group == "Deleted") %>% select(starts_with("time"))),
        na.rm = TRUE
      )
      
      worst_prognosis_group <- if (deleted_wins > duplicated_wins) "Deleted" else "Duplicated"
    }
    
  } else if (length(present_groups) == 3) {
    # Three-group case: Calculate win counts as before
    
    normal_wins <- sum(
      (group_surv_multiple %>% filter(group == "Normal") %>% select(starts_with("time"))) <
        (group_surv_multiple %>% filter(group != "Normal") %>% select(starts_with("time")) %>% apply(1, min, na.rm = TRUE)),
      na.rm = TRUE
    )
    
    deleted_wins <- sum(
      (group_surv_multiple %>% filter(group == "Deleted") %>% select(starts_with("time"))) <
        (group_surv_multiple %>% filter(group != "Deleted") %>% select(starts_with("time")) %>% apply(1, min, na.rm = TRUE)),
      na.rm = TRUE
    )
    
    duplicated_wins <- sum(
      (group_surv_multiple %>% filter(group == "Duplicated") %>% select(starts_with("time"))) <
        (group_surv_multiple %>% filter(group != "Duplicated") %>% select(starts_with("time")) %>% apply(1, min, na.rm = TRUE)),
      na.rm = TRUE
    )
    
    # Debugging output to verify win counts
    cat("Win counts - Normal:", normal_wins, "Deleted:", deleted_wins, "Duplicated:", duplicated_wins, "\n")
    
    # Identify the two groups with the most occurrences of the lowest survival times
    win_counts <- c(Normal = normal_wins, Deleted = deleted_wins, Duplicated = duplicated_wins)
    top_groups <- names(sort(win_counts, decreasing = TRUE))[1:2]
    
    # Pairwise comparison between the top two groups
    pairwise_p_val <- compare_groups(data_surv, top_groups[1], top_groups[2], metric)
    
    # Determine worst prognosis group
    if (pairwise_p_val < 0.05) {
      worst_prognosis_group <- if (win_counts[top_groups[1]] > win_counts[top_groups[2]]) top_groups[1] else top_groups[2]
    } else {
      worst_prognosis_group <- paste(top_groups, collapse = " and ")
    }
  }
  
  # Final return
  return(list(chisq = chisq_log_rank, p_val = p_val, worst_prognosis_group = worst_prognosis_group))
}



# Loop through genes and cancer types for analysis
for (i in seq_len(nrow(table))) {
  gene <- table$genes[i]
  cancer_type <- table$cancer_types[i]
  
  for (metric in survival_metrics) {
    result <- perform_survival_analysis(gene, cancer_type, metric)
    
    if (!is.null(result)) {
      table[[paste0(metric, "_log_rank_chisq")]][i] <- result$chisq
      table[[paste0(metric, "_p_val")]][i] <- result$p_val
      table[[paste0(metric, "_worst_prognosis_group")]][i] <- result$worst_prognosis_group
    }
  }
}

# Salvar a tabela atualizada
export(table, "Updated_Cox_results_with_survival_analysis.tsv")

cat("Analysis complete. Results saved.\n")



