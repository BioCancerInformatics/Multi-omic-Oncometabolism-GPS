##### RSCRIPT TO SURVIVAL ANALYSIS
##### AUTHOR: Higor Almeida Cordeiro Nogueira
##### UPDATED: 20/09/2024

# Load required libraries
library("rio")
library("UCSCXenaShiny")
library("dplyr")
library("survival")

# Set directory
setwd("F:/Higor/Cancer_metabolism_analysis_06/8- Survival_Analysis/lncRNA/Expression/")

# Import table
table <- import("F:/Higor/Cancer_metabolism_analysis_06/7- Cox_Analysis/lncRNA/Expression/Cox_results.tsv")

# Initialize new columns for each survival metric
survival_metrics <- c("OS", "DSS", "DFI", "PFI")
for (metric in survival_metrics) {
  table[[paste0(metric, "_log_rank_chisq")]] <- NA
  table[[paste0(metric, "_p_val")]] <- NA
  table[[paste0(metric, "_worst_prognosis_group")]] <- NA
}

# Function to perform survival analysis for a given metric
perform_survival_analysis <- function(genes, cancer_type, metric) {
  cat("Running survival analysis for gene:", genes, "and cancer_type:", cancer_type, "on metric:", metric, "\n")
  
  # Load survival data with error handling
  data_surv <- tryCatch({
    tcga_surv_get( 
      genes,   
      TCGA_cohort = cancer_type,
      profile = "mRNA",
      TCGA_cli_data = dplyr::inner_join(
        load_data("tcga_clinical"), 
        load_data("tcga_surv"),  
        by = "sample"
      ))
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
      time = paste0(metric, ".time"),  # Use the appropriate survival time metric
      status = metric,  # Use the appropriate survival status metric
      cutoff_mode = "Custom", 
      cutpoint = c(50, 50), 
      profile = "mRNA",
      palette = "aaas"
    )
  }, error = function(e) {
    warning(paste("Failed to create survival plot for", genes, "in", cancer_type, ":", e$message))
    return(NULL)
  })
  
  # Process data for survival analysis
  data_surv <- data_surv %>%
    arrange(.data$value) %>%
    mutate(per_rank = percent_rank(.data$value) * 100) %>%
    mutate(group = case_when(
      .data$per_rank > 50 ~ "High",
      .data$per_rank <= 50 ~ "Low",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(group))
  
  # Calculate log-rank test
  surv_diff <- tryCatch({
    survdiff(Surv(data_surv[[paste0(metric, ".time")]], data_surv[[metric]]) ~ group, data = data_surv)
  }, error = function(e) {
    warning(paste("Failed to calculate survival differences for", metric, ":", e$message))
    return(NULL)
  })
  
  if (is.null(surv_diff)) return(NULL)
  
  # Extract p-value and chi-squared statistic
  p_val <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)
  chisq_log_rank <- surv_diff$chisq
  
  # Determine worst prognosis group considering multiple survival points
  group_surv_multiple <- tryCatch({
    plot_surv[["plot"]][["data"]] %>%
      group_by(group) %>%
      summarize(
        time_40 = min(time[surv <= 0.4], na.rm = TRUE),
        time_50 = min(time[surv <= 0.5], na.rm = TRUE),
        time_60 = min(time[surv <= 0.6], na.rm = TRUE),
        time_70 = min(time[surv <= 0.7], na.rm = TRUE),
        time_80 = min(time[surv <= 0.8], na.rm = TRUE),  # 0.8 survival point
        time_90 = min(time[surv <= 0.9], na.rm = TRUE)   # Added 0.9 survival point
      )
  }, error = function(e) {
    warning(paste("Failed to summarize survival data for", metric, ":", e$message))
    return(NULL)
  })
  
  worst_prognosis_group <- if (!is.null(group_surv_multiple)) {
    if (p_val >= 0.05) {
      "NS"  # Not significant
    } else {
      # Count how many times "Low" reaches the minimum time at each survival point
      low_wins <- sum(
        !is.na(group_surv_multiple$time_40) && group_surv_multiple$time_40[group_surv_multiple$group == "Low"] < group_surv_multiple$time_40[group_surv_multiple$group == "High"],
        !is.na(group_surv_multiple$time_50) && group_surv_multiple$time_50[group_surv_multiple$group == "Low"] < group_surv_multiple$time_50[group_surv_multiple$group == "High"],
        !is.na(group_surv_multiple$time_60) && group_surv_multiple$time_60[group_surv_multiple$group == "Low"] < group_surv_multiple$time_60[group_surv_multiple$group == "High"],
        !is.na(group_surv_multiple$time_70) && group_surv_multiple$time_70[group_surv_multiple$group == "Low"] < group_surv_multiple$time_70[group_surv_multiple$group == "High"],
        !is.na(group_surv_multiple$time_80) && group_surv_multiple$time_80[group_surv_multiple$group == "Low"] < group_surv_multiple$time_80[group_surv_multiple$group == "High"],
        !is.na(group_surv_multiple$time_90) && group_surv_multiple$time_90[group_surv_multiple$group == "Low"] < group_surv_multiple$time_90[group_surv_multiple$group == "High"]  # Added 0.9 comparison
      )
      
      high_wins <- sum(
        !is.na(group_surv_multiple$time_40) && group_surv_multiple$time_40[group_surv_multiple$group == "High"] < group_surv_multiple$time_40[group_surv_multiple$group == "Low"],
        !is.na(group_surv_multiple$time_50) && group_surv_multiple$time_50[group_surv_multiple$group == "High"] < group_surv_multiple$time_50[group_surv_multiple$group == "Low"],
        !is.na(group_surv_multiple$time_60) && group_surv_multiple$time_60[group_surv_multiple$group == "High"] < group_surv_multiple$time_60[group_surv_multiple$group == "Low"],
        !is.na(group_surv_multiple$time_70) && group_surv_multiple$time_70[group_surv_multiple$group == "High"] < group_surv_multiple$time_70[group_surv_multiple$group == "Low"],
        !is.na(group_surv_multiple$time_80) && group_surv_multiple$time_80[group_surv_multiple$group == "High"] < group_surv_multiple$time_80[group_surv_multiple$group == "Low"],
        !is.na(group_surv_multiple$time_90) && group_surv_multiple$time_90[group_surv_multiple$group == "High"] < group_surv_multiple$time_90[group_surv_multiple$group == "Low"]  # Added 0.9 comparison
      )
      
      if (low_wins > high_wins) {
        "Low"
      } else if (high_wins > low_wins) {
        "High"
      } else {
        # Tie-breaking criterion at 0.5
        ifelse(group_surv_multiple$time_50[group_surv_multiple$group == "Low"] < group_surv_multiple$time_50[group_surv_multiple$group == "High"],
               "Low", "High")
      }
    }
  } else {
    NA
  }
  
  return(list(chisq = chisq_log_rank, p_val = p_val, worst_prognosis_group = worst_prognosis_group))
}

# Loop to analyze each gene and cancer_type for all survival metrics
for (i in 1:nrow(table)) {
  cancer_type <- table$cancer_type[i]
  genes <- table$genes[i]
  
  cat("Processing gene:", genes, "for cancer_type:", cancer_type, "\n")
  
  for (metric in survival_metrics) {
    cat("Analyzing metric:", metric, "for gene:", genes, "and cancer_type:", cancer_type, "\n")
    
    result <- perform_survival_analysis(genes, cancer_type, metric)
    
    if (!is.null(result)) {
      table[[paste0(metric, "_log_rank_chisq")]][i] <- result$chisq
      table[[paste0(metric, "_p_val")]][i] <- result$p_val
      table[[paste0(metric, "_worst_prognosis_group")]][i] <- result$worst_prognosis_group
    }
    
    cat("Completed metric:", metric, "for gene:", genes, "and cancer_type:", cancer_type, "\n")
  }
  
  cat("Finished processing gene:", genes, "for cancer_type:", cancer_type, "\n")
}

# Save the updated table with results for all metrics
export(table, "multi_metric_survival_analysis_results.tsv")

cat("Analysis complete. Results saved.\n")