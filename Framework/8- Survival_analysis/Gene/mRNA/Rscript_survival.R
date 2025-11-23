##### RSCRIPT TO SURVIVAL ANALYSIS
##### AUTHOR: Higor Almeida Cordeiro Nogueira
##### UPDATED: 20/09/2024

# Load required libraries
library("rio")
library("UCSCXenaShiny")
library("dplyr")
library("survival")

# Set directory
setwd("F:/Higor/Metabolic_gene_analysis_05/7- Survival_Analysis/mRNA/")

# Import table
table <- import("F:/Higor/Metabolic_gene_analysis_05/6- Cox_Analysis/mRNA/Cox_results.tsv")

# # Filter the table to include only rows where genotypic is "mRNA"
# table <- table %>% filter(genotypic == "mRNA")
# 
# genes <- c("TP53", "AOC3", "AXL", "AOC3")
# tumor <- c("PAAD", "BRCA", "PAAD", "OV")
# table <- data.frame(genes, tumor)

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
  
  # Determine worst prognosis group
  group_surv_50 <- tryCatch({
    plot_surv[["plot"]][["data"]] %>%
      group_by(group) %>%
      summarize(time_50 = min(time[surv <= 0.7], na.rm = TRUE))
  }, error = function(e) {
    warning(paste("Failed to summarize survival data for", metric, ":", e$message))
    return(NULL)
  })
  
  worst_prognosis_group <- if (!is.null(group_surv_50)) {
    ifelse(p_val >= 0.05, "NS", 
           ifelse(group_surv_50$time_50[group_surv_50$group == "Low"] < group_surv_50$time_50[group_surv_50$group == "High"],
                  "Low", "High"))
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

