### Rscript

### Libraries
library(rio)
library(dplyr)
library(tidyr)
library(stringr)

### Setting working directory
setwd("E:/Oncometabolism_GPS/21- Meaningful_interactions/")

### Data frame
df001_gene_protein_mirna_lncrna_signatures <- import("df_gene_protein_mirna_lncrna_signatures.rds")
df002_transcript_signatures <- import("df_transcript_signatures.rds")
df002_transcript_signatures <- df002_transcript_signatures %>% 
  select(-TranscriptGene)

### All signatures
df003_all_signatures <- rbind(df001_gene_protein_mirna_lncrna_signatures, df002_transcript_signatures)


# -----------------------------------------------
# ⚙️ Definir mapeamento
# -----------------------------------------------

map_rcd <- c(
  "Unrelated"                      = 0,
  "Alkaliptosis"                   = 1,
  "Apoptosis"                      = 2,
  "Autophagy_dependent_cell_death" = 3,
  "Cuproptosis"                    = 4,
  "Entotic_cell_death"             = 5,
  "Ferroptosis"                    = 6,
  "Lysosome_dependent_cell_death"  = 7,
  "Necroptosis"                    = 8,
  "Oxeiptosis"                     = 9,
  "Parthanatos"                    = 10,
  "Pyroptosis"                     = 11
)


map_metabolism <- c(
  "Amino acid metabolism"                = 1,
  "Carbohydrate metabolism"              = 2,
  "Energy metabolism"                    = 3,
  "Lipid metabolism"                     = 4,
  "Metabolism of cofactors and vitamins" = 5,
  "Metabolism of other amino acids"      = 6,
  "Nucleotide metabolism"                = 7
)

### Mapping Metabolic Cell Death 
df004_mapping <- df003_all_signatures %>%
  mutate(
    # Codificação para Metabolic_cell_death
    MCD = map_rcd[Metabolic_cell_death]
    )
    
### Mapping INTERACTION
df004_mapping <- df004_mapping %>%
  mutate(
    # Codificação para Interaction
    RCC = if_else(
      Meaningful_interaction == "No meaningful interaction", 0, 1
    )
  )

# -----------------------------------------------
# 🛠️ Passos para numerar corretamente os Pathways
# -----------------------------------------------

# 1. Atribuir código do metabolismo
df004_mapping <- df004_mapping %>%
  mutate(
    MET = map_metabolism[Metabolism]
  )

# 2. Criar tabela auxiliar com Pathways únicos numerados dentro de cada Metabolism
pathway_codes <- df004_mapping %>%
  distinct(Metabolism, Pathways) %>%  # pegar só combinações únicas
  arrange(Metabolism, Pathways) %>%
  group_by(Metabolism) %>%
  mutate(
    PATH = row_number()        # numerar Pathways únicos
  ) %>%
  ungroup()

# 3. Fazer o join de volta no dataframe original
df004_mapping <- df004_mapping %>%
  left_join(pathway_codes, by = c("Metabolism", "Pathways"))

# Mapeamento textual para o código TNC
map_tnc <- c(
  "No data"        = 0,
  "Unchanged"      = 1,
  "Underexpression" = 2,
  "Overexpression"  = 3
)

# -----------------------------------------------
# 🛠️ Aplicar os mapeamentos para SCS e TNC
# -----------------------------------------------

# Atualizar o dataframe df001
df004_mapping <- df004_mapping %>%
  mutate(
    # SCS: sinal da correlação
    SCS = ifelse(Correlation_rho < 0, "N", "P"),
    
    # TNC: tumor vs. tecido normal
    TNC = map_tnc[Tumor_vs_normal],
    
    # # LSC: 1 if LassoCox_selection is present (not NA / "NA" / empty); 0 otherwise
    # LSC = ifelse(
    #   is.na(LassoCox_selection) |
    #     toupper(trimws(as.character(LassoCox_selection))) %in% c("", "NA", "N/A"),
    #   0L, 1L
    # )
  )


# Mapeamento da Omic_layer para OFC (multi-omic feature code)
map_omic_layer <- c(
  "Protein expression" = 1,
  "Mutation"           = 2,
  "CNV"                = 3,
  "miRNA expression"   = 4,
  "Transcript expression" = 5,
  "Gene expression"    = 6,
  "Methylation"        = 7
)

# Mapeamento da Phenotypic_layer para PFC (phenotypic feature code)
map_phenotypic_layer <- c(
  "Tumor mutational burden"    = 1,
  "Microsatellite instability" = 2,
  "Stemness"                   = 3
)

# -----------------------------------------------
# 🛠️ Aplicar os mapeamentos
# -----------------------------------------------

# Atualizar o dataframe
df004_mapping <- df004_mapping %>%
  mutate(
    OFC = map_omic_layer[Omic_layer],
    PFC = map_phenotypic_layer[Phenotypic_layer]
  )


# Função vetorizada para HRC
gerar_letra_hrc <- function(cox_type) {
  case_when(
    is.na(cox_type) ~ "A",
    cox_type %in% c("NS", "No data") ~ "A",
    cox_type == "Risky" ~ "B",
    cox_type == "Protective" ~ "C",
    TRUE ~ "A" # fallback para casos inesperados
  )
}

# Função vetorizada para SMC
gerar_letra_smc <- function(worst_group, ofc) {
  case_when(
    is.na(worst_group) ~ "A",
    worst_group %in% c("NS", "No data") ~ "A",
    
    # Para Protein, miRNA, Transcript, mRNA, Methylation
    ofc %in% c(1, 4, 5, 6, 7) & worst_group == "High" ~ "B",
    ofc %in% c(1, 4, 5, 6, 7) & worst_group == "Low" ~ "C",
    
    # Para Mutation
    ofc == 2 & worst_group == "MT" ~ "B",
    ofc == 2 & worst_group == "WT" ~ "C",
    
    # Para CNV
    ofc == 3 & worst_group == "Deleted" ~ "B",
    ofc == 3 & worst_group %in% c("Duplicated", "Normal") ~ "C",
    ofc == 3 & worst_group %in% c(
      "Deleted and Duplicated", "Duplicated and Deleted",
      "Normal and Deleted", "Deleted and Normal",
      "Normal and Duplicated", "Duplicated and Normal"
    ) ~ "D",
    
    TRUE ~ "A" # fallback
  )
}

# -----------------------------------------------
# ⚙️ Aplicar geração de HRC e SMC
# -----------------------------------------------

df004_mapping <- df004_mapping %>%
  mutate(
    # Geração do HRC concatenado
    HRC = paste0(
      "1", gerar_letra_hrc(Cox_DSS_type),
      "2", gerar_letra_hrc(Cox_DFI_type),
      "3", gerar_letra_hrc(Cox_PFI_type),
      "4", gerar_letra_hrc(Cox_OS_type)
    ),
    
    # Geração do SMC concatenado
    SMC = paste0(
      "1", gerar_letra_smc(DSS_worst_prognosis_group, OFC),
      "2", gerar_letra_smc(DFI_worst_prognosis_group, OFC),
      "3", gerar_letra_smc(PFI_worst_prognosis_group, OFC),
      "4", gerar_letra_smc(OS_worst_prognosis_group, OFC)
    )
  )

# Mapeamento da coluna do Microambiente Tumoral (TMC)
map_tmc <- c(
  "anti-tumoral" = 1,
  "dual" = 2,        # função dupla pró e anti-tumoral
  "pro-tumoral" = 3,
  "NS" = 4
)

# Mapeamento da coluna de Infiltração Linfocitária (TIC)
map_tic <- c(
  "Hot" = 1,         # alta infiltração
  "Variable" = 2,    # infiltração variável
  "Cold" = 3,        # baixa infiltração
  "NS" = 4
)

# -----------------------------------------------
# 🛠️ Aplicar os mapeamentos
# -----------------------------------------------

# Atualizar o dataframe
df004_mapping <- df004_mapping %>%
  mutate(
    TMC = map_tmc[Microenvironment_classification],
    TIC = map_tic[Immune_classification]
  )

# #########################
# 
# #### No effect value signatures:
# filtered_rows_meaningless <- df004_mapping %>%
#   filter(HRC == "1A2A3A4A", SMC == "1A2A3A4A")
# 
# #### Select all rows that are meaningful by "HRC_series" and "SMC_series"
# # Select all rows where the value in "HRC_series" an  ""SMC_series) is different that "1A2A3A4A" 
# filtered_rows_meaningful <- df004_mapping %>%
#   filter(HRC != "1A2A3A4A", SMC != "1A2A3A4A")
# 
# setwd("E:/Oncometabolism_GPS/22- Nomenclature/")
# 
# export(filtered_rows_meaningful, "Meaningful_filtered_signatures.tsv") 

##########
########## Making the ranking templates for HRC and SMC components
##########
# Creating the mapping template (df12) as a named vector for easy lookup

HRC_rank_template <- c('A' = 0, 'B' = 1, 'C' = 1, '1' = 2, '2' = 2, '3' = 2, '4' = 1)

# Create a data frame
HRC_map <- data.frame(
  HRC = names(HRC_rank_template),
  HRC_Rank = as.vector(HRC_rank_template)
)

# Define a function to rank each series
HRC_rank_series <- function(series_value) {
  # Split the series into individual parts, e.g., "1A2A3D4A" -> c("1A", "2A", "3D", "4A")
  split_series <- strsplit(series_value, "(?<=[A-D])", perl = TRUE)[[1]]
  
  total_rank <- 0
  
  # Loop through each part of the series (e.g., "1A", "2A", etc.)
  for (part in split_series) {
    # Extract the number and letter
    num <- substr(part, 1, 1)  # First character is the number
    letter <- substr(part, 2, 2)  # Second character is the letter
    
    # Check if letter is "A"; if so, skip adding the rank for this number-letter pair
    if (letter == "A") {
      next  # Skip to the next part
    }
    
    # Add the ranks from the template
    total_rank <- total_rank + HRC_rank_template[num] + HRC_rank_template[letter]
  }
  
  return(total_rank)
}

# Apply the ranking function to the "series" column in df12
df005_ranking <- df004_mapping %>%
  mutate(HRC_Rank = sapply(HRC, HRC_rank_series))

# Create the mapping template (df13) as a named vector for easy lookup
SMC_rank_template <- c('A' = 0, 'B' = 1, 'C' = 1,  'D' = 1, '1' = 2, '2' = 2, '3' = 2, '4' = 1)

# Create a data frame
SMC_map <- data.frame(
  SMC = names(SMC_rank_template),
  SMC_Rank = as.vector(SMC_rank_template)
)

# Define a function to rank each series
SMC_rank_series <- function(series_value) {
  # Split the series into individual parts, e.g., "1A2A3D4A" -> c("1A", "2A", "3D", "4A")
  split_series <- strsplit(series_value, "(?<=[A-D])", perl = TRUE)[[1]]
  
  total_rank <- 0
  
  # Loop through each part of the series (e.g., "1A", "2A", etc.)
  for (part in split_series) {
    # Extract the number and letter
    num <- substr(part, 1, 1)  # First character is the number
    letter <- substr(part, 2, 2)  # Second character is the letter
    
    # Check if letter is "A"; if so, skip adding the rank for this number-letter pair
    if (letter == "A") {
      next  # Skip to the next part
    }
    
    # Add the ranks from the template
    total_rank <- total_rank + SMC_rank_template[num] + SMC_rank_template[letter]
  }
  
  return(total_rank)
}

# Apply the ranking function to the "series" column
df005_ranking <- df005_ranking %>%
  mutate(SMC_Rank = sapply(SMC, SMC_rank_series))

#######
### Signature components rank maps!
#######
# Genomic Feature Ranking
OFC_map <- data.frame(
  OFC = c("1", "2", "3", "4", "5", "6", "7"),
  OFC_Rank = c(7, 6, 5, 4, 3, 2, 1) # Protein Expression = 7,
  # Mutations = 6, Copy Number Variation = 5, miRNA Expression = 4, 
  # Transcript Expression = 3, mRNA Expression = 2, CpG Methylation =1
)

# Phenotypic Feature Ranking
PFC_map <- data.frame(
  PFC = c("1", "2", "3"),
  PFC_Rank = c(3, 2, 1) # TMB=3, MSI=2, TSM=1
)

# Metabolic cell death Ranking
MCD_map <- data.frame(
  MCD = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"),
  MCD_Rank = c("0", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1", "1")
)

# Meaningful interaction Ranking (Regulatory circuitry)
RCC_map <- data.frame(
  RCC = c("0", "1"),
  RCC_Rank = c("0", "2")
)

# Spearman Correlation Sign (SCS) Feature Ranking
# Ranking SCS features based on their correlation signs.
SCS_map <- data.frame(
  SCS = c("P", "N"),
  SCS_Rank = c(2, 1) # P=Positive=2 , N=Negative=1
)

# TCGA vs. GTEx expression Code (TNC) Feature Ranking
# Ranking TNC features based on their expression codes.
# NOTE: MODIFIED revised on 01/12/2024
TNC_map <- data.frame(
  TNC = c("0", "1", "2", "3"),
  TNC_Rank = c(0, 1, 2, 2) # No_data; unchanged=1, underexpressed=3, overexpressed=3 NOTE: MODIFIED revised on 01/12/2024
)

# TMC and TIC Feature Ranking
# A combined ranking for TMC and TIC features.
TMC_TIC_map <- data.frame(
  TMC = c("1", "2", "3", "4"),
  TIC = c("1", "2", "3", "4"),
  TMC_Rank = c(7, 4, 1, 0), 
  TIC_Rank = c(7, 4, 1, 0) # "Hot", "Variable", "Cold", "Anti-tumoral, "Dual", "Pro-tumoral", "NS"
)


# Rename columns
df005_ranking <- df005_ranking %>%
  rename(
    CTAB = CTAB,
    HRC_Template = HRC,
    SMC_Template = SMC,
    HRC = HRC_Rank,
    SMC = SMC_Rank
  )

# Create the USI (Unique Series Identifier) DISTRIBUTION variable based on CTAB and Signature
df006_USI <- df005_ranking %>%
  group_by(CTAB) %>%  # Group by CTAB
  mutate(USI = row_number()) %>%  # Create the series number within each CTAB group
  ungroup() %>%  # Ungroup to finish the transformation
  relocate(USI, .after = CTAB)  # Move GSI column to be after the CTAB column

# Distribution of signatures per cancer type
# Summarize the DISTRIBUTION number of Signatures per CTAB and add total row
summary_signatures_GSI <- df006_USI %>%
  group_by(CTAB) %>%
  summarise(num_signatures = n()) %>%  # Count the number of rows (Signatures) for each CTAB
  ungroup() %>%
  # Add a row for the total number of CTABs and Signatures
  bind_rows(
    summarise(., CTAB = "Total", num_signatures = sum(num_signatures))
  )

# Reorder columns
df007_reorder <- df006_USI %>%
  select(
    Members,
    Signatures,
    Molecular_class,
    Common_interaction,
    Meaningful_interaction,
    Metabolism,
    Pathways,
    Metabolic_cell_death,
    CTAB,
    USI,
    MET,
    PATH,
    MCD,
    RCC,
    OFC,
    PFC,
    SCS,
    TNC,
    HRC,
    SMC,
    TMC,
    TIC,
    # LSC,
    HRC_Template,
    SMC_Template,
    Tumor_vs_normal,
    Tumor_vs_normal_p.adj,
    Omic_layer,
    Phenotypic_layer,
    Correlation_rho,
    Correlation_p.adj,
    Cox_OS_type,
    Cox_OS_p.value,
    Cox_DSS_type,
    Cox_DSS_p.value,
    Cox_DFI_type,
    Cox_DFI_p.value,
    Cox_PFI_type,
    Cox_PFI_p.value,
    OS_worst_prognosis_group,
    OS_p.value,
    DSS_worst_prognosis_group,
    DSS_p.value,
    DFI_worst_prognosis_group,
    DFI_p.value,
    PFI_worst_prognosis_group,
    PFI_p.value,
    Microenvironment_classification,
    Microenvironment_score_details,
    Immune_classification,
    Immune_score_details,
    # LassoCox_selection,
    # pi_hat,
    
  )

# Create the "Nomenclature" column and place it after the "Signatures" column
df008_nomenclature <- df007_reorder %>%
  mutate(
    Nomenclature = apply(df007_reorder[, 9:22], 1, function(row) {
      row <- trimws(row)  # Remove espaços em branco à esquerda/direita de todos os elementos
      paste(c(paste(row[1], row[2], sep = "-"), row[3:14]), collapse = ".")
    })
  ) %>%
  relocate(Nomenclature, .after = Signatures)


# Define a function to calculate the rank for a string like "1B2C3A4C"
calculate_rank <- function(code_string, letter_scores, digit_scores) {
  # Split into pairs like "1B", "2C", ...
  parts <- unlist(regmatches(code_string, gregexpr(".{2}", code_string)))
  
  total_score <- 0
  for (part in parts) {
    digit <- substr(part, 1, 1)
    letter <- substr(part, 2, 2)
    
    if (letter != "A") {
      total_score <- total_score + digit_scores[digit] + letter_scores[letter]
    }
  }
  
  return(total_score)
}


### HRC ###

# Letter and digit scores
letters_HRC <- c("A", "B", "C")
letter_scores_HRC <- c("A" = 0, "B" = 1, "C" = 1)
digit_scores <- c("1" = 2, "2" = 2, "3" = 2, "4" = 1)

# Generate combinations
combs_HRC <- expand.grid(L1 = letters_HRC, 
                         L2 = letters_HRC, 
                         L3 = letters_HRC, 
                         L4 = letters_HRC,
                         stringsAsFactors = FALSE)

# Create codes like "1A2B3C4A"
combs_HRC$Array_series_class <- paste0("1", combs_HRC$L1, 
                                       "2", combs_HRC$L2, 
                                       "3", combs_HRC$L3, 
                                       "4", combs_HRC$L4)

# Apply the function to each code
combs_HRC$HRC_Rank <- sapply(combs_HRC$Array_series_class, calculate_rank, 
                             letter_scores = letter_scores_HRC, 
                             digit_scores = digit_scores)

# Final table
df_HRC_map <- combs_HRC[, c("Array_series_class", "HRC_Rank")]


### SMC ###

# Letter and digit scores
letters_SMC <- c("A", "B", "C", "D")
letter_scores_SMC <- c("A" = 0, "B" = 1, "C" = 1, "D" = 1)
digit_scores <- c("1" = 2, "2" = 2, "3" = 2, "4" = 1)

# Generate combinations
combs_SMC <- expand.grid(L1 = letters_SMC, 
                         L2 = letters_SMC, 
                         L3 = letters_SMC, 
                         L4 = letters_SMC,
                         stringsAsFactors = FALSE)

# Create codes like "1A2B3C4D"
combs_SMC$Array_series_class <- paste0("1", combs_SMC$L1, 
                                       "2", combs_SMC$L2, 
                                       "3", combs_SMC$L3, 
                                       "4", combs_SMC$L4)

# Apply the function
combs_SMC$SMC_Rank <- sapply(combs_SMC$Array_series_class, calculate_rank, 
                             letter_scores = letter_scores_SMC, 
                             digit_scores = digit_scores)

# Final table
df_SMC_map <- combs_SMC[, c("Array_series_class", "SMC_Rank")]

# Function to map ranks and add new ranking columns at the right end of the dataframe
assign_ranks <- function(df, col_name, map_df, map_col, rank_col) {
  # Match the values in the target column with the mapping column and extract the rank
  rank_values <- map_df[[rank_col]][match(df[[col_name]], map_df[[map_col]])]
  # Create the new column name
  rank_col_name <- paste0(col_name, "_ranking")
  # Add the new ranking column at the right end of the dataframe
  df[[rank_col_name]] <- rank_values
  return(df)
}

# Apply the function for each column with its respective map

# GFC ranking
df009_ranking <- assign_ranks(df008_nomenclature, "OFC", OFC_map, "OFC", "OFC_Rank")

# PFC ranking
df009_ranking <- assign_ranks(df009_ranking, "PFC", PFC_map, "PFC", "PFC_Rank")

# MCD ranking
df009_ranking <- assign_ranks(df009_ranking, "MCD", MCD_map, "MCD", "MCD_Rank")

# RC ranking
df009_ranking <- assign_ranks(df009_ranking, "RCC", RCC_map, "RCC", "RCC_Rank")

# SCS ranking
df009_ranking <- assign_ranks(df009_ranking, "SCS", SCS_map, "SCS", "SCS_Rank")

# TNC ranking
df009_ranking <- assign_ranks(df009_ranking, "TNC", TNC_map, "TNC", "TNC_Rank")

# HRC ranking
df009_ranking <- assign_ranks(df009_ranking, "HRC_Template", df_HRC_map, "Array_series_class", "HRC_Rank")

# SMC ranking
df009_ranking <- assign_ranks(df009_ranking, "SMC_Template", df_SMC_map, "Array_series_class", "SMC_Rank")

# TMC ranking
df009_ranking <- assign_ranks(df009_ranking, "TMC", TMC_TIC_map, "TMC", "TMC_Rank")

# TIC ranking
df009_ranking <- assign_ranks(df009_ranking, "TIC", TMC_TIC_map, "TIC", "TIC_Rank")

#######
# Ensure the Final_Rank column is numeric
df009_ranking$Final_Rank <- rowSums(
  sapply(df009_ranking[, c("OFC_ranking", "PFC_ranking", "MCD_ranking", "RCC_ranking", 
                           "SCS_ranking", "TNC_ranking", "TMC_ranking", "TIC_ranking", 
                           "HRC_Template_ranking", "SMC_Template_ranking")], as.numeric),
  na.rm = TRUE
)

###### Combined_outcome SMC column 

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

# Apply to df009
df010_outcome <- df009_ranking %>%
  rowwise() %>%
  mutate(Combined_outcome_SMC = classify_combined_outcome(c_across(c(
    OS_worst_prognosis_group,
    DSS_worst_prognosis_group,
    DFI_worst_prognosis_group,
    PFI_worst_prognosis_group
  )))) %>%
  ungroup()


###### Combined_outcome_HRC column 

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
df010_outcome_HRC <- df010_outcome %>%
  rowwise() %>%
  mutate(Combined_outcome_HRC = classify_combined_outcome(c_across(c(
    Cox_OS_type,
    Cox_DSS_type,
    Cox_DFI_type,
    Cox_PFI_type
  )))) %>%
  ungroup()

# Reorder columns
df011_reorder <- df010_outcome_HRC %>%
  select(
    Members,
    Signatures,
    Nomenclature,
    Molecular_class,
    Common_interaction,
    Meaningful_interaction,
    Metabolism,
    Pathways,
    Metabolic_cell_death,
    CTAB,
    Tumor_vs_normal,
    Tumor_vs_normal_p.adj,
    Omic_layer,
    Phenotypic_layer,
    Correlation_rho,
    Correlation_p.adj,
    Combined_outcome_HRC,
    Cox_OS_type,
    Cox_OS_p.value,
    Cox_DSS_type,
    Cox_DSS_p.value,
    Cox_DFI_type,
    Cox_DFI_p.value,
    Cox_PFI_type,
    Cox_PFI_p.value,
    Combined_outcome_SMC,
    OS_worst_prognosis_group,
    OS_p.value,
    DSS_worst_prognosis_group,
    DSS_p.value,
    DFI_worst_prognosis_group,
    DFI_p.value,
    PFI_worst_prognosis_group,
    PFI_p.value,
    Microenvironment_classification,
    Microenvironment_score_details,
    Immune_classification,
    Immune_score_details,
    # LassoCox_selection,
    # pi_hat,
    USI,
    MET,
    PATH,
    MCD,
    MCD_ranking,
    RCC,
    RCC_ranking,
    OFC,
    OFC_ranking,
    PFC,
    PFC_ranking,
    SCS,
    SCS_ranking,
    TNC,
    TNC_ranking,
    HRC_Template,
    HRC,
    SMC_Template,
    SMC,
    TMC,
    TMC_ranking,
    TIC,
    TIC_ranking,
    # LSC,
    # LSC_ranking,
    Final_Rank
  )


######
# # Select rows where HRC_Rank and SMC_Rank are 11, and TMC_Rank and TIC_Rank are 7
# top_ranked_signatures <- df011_reorder[df011_reorder$SMC == 11 &
#                                 df011_reorder$TMC_ranking == 7 &
#                                 df011_reorder$TIC_ranking == 7, ]
# 
# top_ranked_signatures_v2 <- df011_reorder[
#   df011_reorder$SMC == 11 & 
#     df011_reorder$TMC_ranking == 7 & 
#     df011_reorder$TIC_ranking == 7 &
#     df011_reorder$RC_ranking == 2, ]
# 
# top_ranked_signatures_v3 <- df011_reorder[
#   df011_reorder$SMC == 11 & 
#     df011_reorder$TMC_ranking == 7 & 
#     df011_reorder$TIC_ranking == 7 &
#     df011_reorder$RC_ranking == 2 &
#     df011_reorder$MCD_ranking == 1, ]
# 
# # Select rows where HRC_Rank and SMC_Rank are 11, and TMC_Rank and TIC_Rank are 7
# worst_ranked_signatures <- df011_reorder[df011_reorder$SMC == 11 & 
#                                            df011_reorder$TMC_ranking == 1 & 
#                                            df011_reorder$TIC_ranking == 1, ]
# 
# # Select rows where HRC_Rank and SMC_Rank are 11, and TMC_Rank and TIC_Rank are 7
# worst_ranked_signatures_v2 <- df011_reorder[df011_reorder$SMC == 11 & 
#                                               df011_reorder$TMC_ranking == 1 & 
#                                               df011_reorder$TIC_ranking == 1 & 
#                                      !(df011_reorder$TNC_ranking %in% c(0, 1)),
# ]
# 
# # # Select rows where LSC is 1, and TMC_Rank and TIC_Rank are 7
# # LSC_ranked_signatures <- df011_reorder[df011_reorder$LSC == 1, ]
# 
# saveRDS(top_ranked_signatures, "top_ranked_signatures.rds")
# 
# saveRDS(top_ranked_signatures_v2, "top_ranked_signatures_v2.rds")
# 
# saveRDS(top_ranked_signatures_v3, "top_ranked_signatures_v3.rds")
# 
# saveRDS(worst_ranked_signatures, "worst_ranked_signatures.rds")
# 
# # Preparing top_ranked_signatures_v3 to Shiny app
# df012_shiny_filtered <- top_ranked_signatures_v3 %>% 
#   select(-USI, -MET, -PATH, -MCD, -RC, -OFC, -PFC, -SCS, -TNC, -HRC, -SMC, -TMC, -TIC, -HRC_Template, 
#          -SMC_Template, -OFC_ranking, -PFC_ranking, -MCD_ranking, -RC_ranking, -SCS_ranking, -TNC_ranking, 
#          -TMC_ranking, -TIC_ranking, -Immune_score_details, -Microenvironment_score_details)
# 
# saveRDS(df012_shiny_filtered, "top_ranked_signatures_shiny.rds")
# 
# # # Preparing LSC_ranked_signatures to Shiny app
# # df012_LSC_shiny_filtered <- LSC_ranked_signatures %>% 
# #   select(-USI, -MET, -PATH, -MCD, -INT, -OFC, -PFC, -SCS, -TNC, -HRC, -SMC, -TMC, -TIC, -LSC, -HRC_Template, 
# #          -SMC_Template, -OFC_ranking, -PFC_ranking, -MCD_ranking, -INT_ranking, -SCS_ranking, -TNC_ranking, 
# #          -TMC_ranking, -TIC_ranking, -LSC_ranking, -Immune_score_details, -Microenvironment_score_details)
# 
# # saveRDS(df012_LSC_shiny_filtered, "LSC_ranked_signatures_shiny.rds")

# Preparing data set to Shiny app
df013_shiny_filtered <- df011_reorder %>% 
  select(-USI, -MET, -PATH, -MCD, -RCC, -OFC, -PFC, -SCS, -TNC, -HRC, -SMC, -TMC, -TIC, -HRC_Template, 
         -SMC_Template, -OFC_ranking, -PFC_ranking, -MCD_ranking, -RCC_ranking, -SCS_ranking, -TNC_ranking, 
         -TMC_ranking, -TIC_ranking, -Immune_score_details, -Microenvironment_score_details)

saveRDS(df013_shiny_filtered, "All_signature_results.rds")

# #### Select all rows that are meaningful by "HRC_series" and "SMC_series"
# # Select all rows where the value in "HRC_series" an  ""SMC_series) is different that "1A2A3A4A" 
# filtered_rows_meaningful <- df011_reorder %>%
#   filter(HRC_Template != "1A2A3A4A", SMC_Template != "1A2A3A4A")
# 
# # Preparing data set to Shiny app
# filtered_rows_meaningful <- filtered_rows_meaningful %>% 
#   select(-USI, -MET, -PATH, -MCD, -RC, -OFC, -PFC, -SCS, -TNC, -HRC, -SMC, -TMC, -TIC, -HRC_Template, 
#          -SMC_Template, -OFC_ranking, -PFC_ranking, -MCD_ranking, -RC_ranking, -SCS_ranking, -TNC_ranking, 
#          -TMC_ranking, -TIC_ranking, -Immune_score_details, -Microenvironment_score_details)
# 
# saveRDS(filtered_rows_meaningful, "Meaningful_signatures.rds")

# Create a new folder named "New_Folder"
dir.create("Omic_layer_specific_signatures")

setwd("E:/Oncometabolism_GPS/22- Nomenclature/Omic_layer_specific_signatures/")

# Creating separate data frames - omic-specific
df014_Gene_expression_results <- df013_shiny_filtered %>% filter(Omic_layer == "Gene expression")
saveRDS(df014_Gene_expression_results, "Gene_expression_results.rds")

df015_Methylation_results <- df013_shiny_filtered %>% filter(Omic_layer == "Methylation")
saveRDS(df015_Methylation_results, "Methylation_results.rds")

df016_miRNA_expression_results <- df013_shiny_filtered %>% filter(Omic_layer == "miRNA expression")
saveRDS(df016_miRNA_expression_results, "miRNA_expression_results.rds")

df017_Protein_expression_results <- df013_shiny_filtered %>% filter(Omic_layer == "Protein expression")
saveRDS(df017_Protein_expression_results, "Protein_expression_results.rds")

df018_Mutation_results <- df013_shiny_filtered %>% filter(Omic_layer == "Mutation")
saveRDS(df018_Mutation_results, "Mutation_results.rds")

df019_CNV_results <- df013_shiny_filtered %>% filter(Omic_layer == "CNV")
saveRDS(df019_CNV_results, "CNV_results.rds")

df020_Transcript_results <- df013_shiny_filtered %>% filter(Omic_layer == "Transcript expression")

df021_transcript_signatures <- import("/Oncometabolism_GPS/21- Meaningful_interactions/df_transcript_signatures.rds")

df021_transcript_signatures <- df021_transcript_signatures %>% 
  select(Signatures, TranscriptGene)

df021_transcript_signatures <- df021_transcript_signatures %>% 
  distinct()

df022_merged <- df020_Transcript_results %>%
  left_join(
    df021_transcript_signatures %>% select(Signatures, TranscriptGene),
    by = c("Signatures" = "Signatures")
  )

# Move the last column to the 3rd position
df022_merged <- df022_merged %>%
  relocate(last_col(), .before = 3)

saveRDS(df022_merged, "Transcript_expression_results.rds")

######## IT IS NECESSARY TO REVIEW

# 📂 Step 1: Load interaction results from RDS and TSV sources
df023_Interactions <- import("/Oncometabolism_GPS/21- Meaningful_interactions/Regulatory_circuitries_results.rds")
df024_Interactions_transcript <- import("/Oncometabolism_GPS/21- Meaningful_interactions/Regulatory_circuitries_transcript_results.rds")

# 🔗 Step 2: Combine both interaction datasets into a unified data frame
df025_All_interactions <- rbind(df023_Interactions, df024_Interactions_transcript)

# # 🔎 Step 3: Select relevant columns from the filtered signature annotation table
# df026_selected <- df011_reorder %>%
#   select(Signatures, CTAB, Omic_layer, Phenotypic_layer, Metabolism, Pathways, Metabolic_cell_death,
#          Nomenclature, Combined_outcome_HRC, HRC_Template, HRC, Combined_outcome_SMC, SMC_Template, SMC, 
#          TMC_ranking, TIC_ranking)
# 
# # 🏷️ Step 4: Rename columns in df026 to match key columns in df025 for joining
# df026_selected <- df026_selected %>%
#   rename(
#     Multiomics_Signature = Signatures,
#     Metabolism_Signature = Metabolism,
#     Pathways_Signature = Pathways,
#     Phenotypic_layer_signature = Phenotypic_layer,
#     Metabolic_cell_death_Signature = Metabolic_cell_death,
#     Cancer_types = CTAB,
#     Omic_layer_signature = Omic_layer,
#     Combined_outcome_HRC_Signature = Combined_outcome_HRC,
#     Combined_outcome_SMC_Signature = Combined_outcome_SMC
#   )
# 
# # 🔗 Step 5: Merge interaction results with signature annotations by multi-key join
# df027_merged <- df025_All_interactions %>%
#   left_join(df026_selected, 
#             by = c("Multiomics_Signature", "Metabolism_Signature", "Pathways_Signature",
#                    "Metabolic_cell_death_Signature", "Phenotypic_layer_signature", "Cancer_types", "Omic_layer_signature"))
# 
# # 🏷️ Step 6: Rename merged columns to distinguish between interaction and signature-derived data
# df028_rename <- df027_merged %>%
#   rename(
#     Signatures = Multiomics_Signature,
#     Metabolism = Metabolism_Signature,
#     Pathways = Pathways_Signature,
#     Metabolic_cell_death = Metabolic_cell_death_Signature,
#     CTAB = Cancer_types,
#     Common_interaction = Common_interaction,
#     Meaningful_Interaction = Interaction_validated,
#     Molecular_class_interaction = Molecular_class_entity,
#     Metabolism_interaction = Metabolism_entity, 
#     Pathways_interaction = Pathways_entity,
#     Metabolic_cell_death_interaction = Metabolic_cell_death_entity,
#     Omic_layer_interaction = Omic_layer_entity,
#     Cox_OS_type_interaction = Cox_OS_type_entity,
#     Cox_DSS_type_interaction = Cox_DSS_type_entity,
#     Cox_DFI_type_interaction = Cox_DFI_type_entity,
#     Cox_PFI_type_interaction = Cox_PFI_type_entity,
#     OS_worst_prognosis_group_interaction = OS_worst_prognosis_group_entity,
#     DSS_worst_prognosis_group_interaction = DSS_worst_prognosis_group_entity,
#     DFI_worst_prognosis_group_interaction = DFI_worst_prognosis_group_entity,
#     PFI_worst_prognosis_group_interaction = PFI_worst_prognosis_group_entity,
#     Microenvironment_classification_interaction = Microenvironment_classification_entity,
#     Immune_classification_interaction = Immune_classification_entity,
#     Combined_outcome_SMC_interaction = Combined_outcome_SMC, 
#     Combined_outcome_HRC_interaction = Combined_outcome_HRC,
#     HRC_template_signature = HRC_Template,
#     SMC_template_signature = SMC_Template,
#     HRC_signature = HRC,
#     SMC_signature = SMC,
#     TMC_ranking_signature = TMC_ranking,
#     TIC_ranking_signature = TIC_ranking,
#     Cox_concordance = Clinical_cox_concordance,
#     Survival_concordance = Clinical_survival_concordance,
#     Cox_metrics_concordance = Cox_results,
#     Survival_metrics_concordance = Survival_results,
#     Phenotypic_association_concordance = Phenotypic_concordance
#   )
# 
# # 📐 Step 7: Reorganize columns for logical and publication-ready output
# df029_ordered <- df028_rename %>%
#   select(
#     Nomenclature, Signatures, Molecular_class_Signature, Common_interaction,
#     Metabolism, Pathways, Metabolic_cell_death,
#     Omic_layer_signature, Phenotypic_layer_signature, CTAB,
#     Combined_outcome_HRC_Signature, Cox_OS_type_signature, Cox_DSS_type_signature,
#     Cox_DFI_type_signature, Cox_PFI_type_signature, Combined_outcome_SMC_Signature,
#     OS_worst_prognosis_group_signature, DSS_worst_prognosis_group_signature, 
#     DFI_worst_prognosis_group_signature, PFI_worst_prognosis_group_signature,
#     Microenvironment_classification_signature, Immune_classification_signature,
#     Meaningful_Interaction, Molecular_class_interaction, Omic_layer_interaction, Phenotypic_association_concordance,
#     Combined_outcome_HRC_interaction, Cox_OS_type_interaction, Cox_DSS_type_interaction,
#     Cox_DFI_type_interaction, Cox_PFI_type_interaction, Cox_concordance, Cox_metrics_concordance, Combined_outcome_SMC_interaction,
#     OS_worst_prognosis_group_interaction, DSS_worst_prognosis_group_interaction,
#     DFI_worst_prognosis_group_interaction, PFI_worst_prognosis_group_interaction, Survival_concordance, Survival_metrics_concordance,
#     Microenvironment_classification_interaction, Immune_classification_interaction, Immune_concordance
#   )
# 
# # 💾 Step 8: Save final structured and ordered data frame to RDS file
# setwd("E:/Cancer_metabolism_project_02/22- Nomenclature/")
# saveRDS(df029_ordered, "Meaningful_interaction.rds")
# 
# 

# 🔗 Step 2: Merge Nomenclature column 
df026_regulatory_circuitry <- merge(
  df025_All_interactions,
  df008_nomenclature[, c("Signatures", "CTAB", "Metabolism", "Pathways", 
            "Metabolic_cell_death", "Omic_layer", 
            "Phenotypic_layer", "Nomenclature")],
  by.x = c("Signatures", "CTAB", "Metabolism_signature", "Pathways_signature",
           "Metabolic_cell_death_signature", "Omic_layer_signature", 
           "Phenotypic_layer_signature"),
  by.y = c("Signatures", "CTAB", "Metabolism", "Pathways", 
           "Metabolic_cell_death", "Omic_layer", 
           "Phenotypic_layer"),
  all.x = TRUE
)

df027_regulatory_circuitry_filtered <- df026_regulatory_circuitry %>% 
  select(-Metabolism_interaction, -Pathways_interaction, -Metabolic_cell_death_interaction) %>% 
  distinct()

# 💾 Step 8: Save final structured and ordered data frame to RDS file
setwd("E:/Oncometabolism_GPS/22- Nomenclature/")
saveRDS(df027_regulatory_circuitry_filtered, "Regulatory_circuitries_with_nomenclature.rds")

################################################################################


# ## Named vectors already defined:
# ## map_metabolism, map_rcd
# 
# # A. Tabela de códigos para Metabolism
# df_codes_metabolism <- enframe(map_metabolism, name = "Metabolism", value = "Metabolism_code") %>%
#   arrange(Metabolism_code)
# 
# # B. Tabela de códigos para Pathways (numerados dentro de cada Metabolism)
# df_codes_metab_pathways <- df004_mapping %>%
#   distinct(Metabolism, Pathways) %>%
#   inner_join(df_codes_metabolism, by = "Metabolism") %>%
#   group_by(Metabolism, Metabolism_code) %>%
#   arrange(Pathways, .by_group = TRUE) %>%
#   mutate(Pathway_code = row_number()) %>%
#   ungroup()
# 
# 
# # C. Tabela de códigos para "Metabolic cell death" (RCD)
# df_codes_rcd <- enframe(map_rcd, name = "Metabolic_cell_death", value = "RCD_code") %>%
#   arrange(RCD_code)
# 
# 
# df_codes_tnc <- enframe(map_tnc, name = "TNC", value = "TNC_code") %>%
#   arrange(TNC_code)
# 
# df_codes_omic_layer <- enframe(map_omic_layer, name = "Omic_layer", value = "OFC") %>%
#   arrange(OFC)
# 
# df_codes_phenotypic_layer <- enframe(map_phenotypic_layer, name = "Phenotypic_layer", value = "PFC") %>%
#   arrange(PFC)
# 
# 
# df_codes_tmc <- enframe(map_tmc, name = "Microenvironment_classification", value = "TMC_code") %>%
#   arrange(TMC_code)
# 
# df_codes_tic <- enframe(map_tic, name = "Immune_infiltration", value = "TIC_code") %>%
#   arrange(TIC_code)
# 
# 
# 
# # Endpoints index used in HRC/SMC strings
# df_endpoint_key <- tibble::tribble(
#   ~Pos, ~Endpoint,
#   1L,  "DSS",
#   2L,  "DFI",
#   3L,  "PFI",
#   4L,  "OS"
# )
# 
# # If you already have map_omic_layer -> OFC (1..7), categorize OFC for SMC rules:
# # Expression layers = 1,4,5,6,7 ; Mutation = 2 ; CNV = 3
# df_ofc_category <- tibble::tribble(
#   ~OFC, ~OFC_category,
#   1L,  "Expression",
#   4L,  "Expression",
#   5L,  "Expression",
#   6L,  "Expression",
#   7L,  "Expression",
#   2L,  "Mutation",
#   3L,  "CNV"
# )
# 
# # Letter → meaning for HRC (Cox-based risk code)
# df_hrc_letter_key <- tibble::tribble(
#   ~Letter, ~Meaning,                                           ~Triggers,
#   "A",     "Neutral/indeterminate",                            "NA, NS, No data, or unexpected value",
#   "B",     "Risky (hazard ↑ with higher feature)",             "Cox_type == 'Risky'",
#   "C",     "Protective (hazard ↓ with higher feature)",        "Cox_type == 'Protective'"
# )
# 
# # A) General meaning for SMC letters
# df_smc_letter_key <- tibble::tribble(
#   ~Letter, ~Meaning,
#   "A", "Neutral/indeterminate (NA, NS, No data, or fallback)",
#   "B", "Worse prognosis group for this OFC category",
#   "C", "Better prognosis group for this OFC category",
#   "D", "Mixed CNV pattern (worse-than-ambiguous set for CNV only)"
# )
# 
# # B) Operational rules used by your gerar_letra_smc()
# # We encode the rule space explicitly
# df_smc_rules <- bind_rows(
#   # Expression layers (OFC: 1,4,5,6,7)
#   tibble(
#     OFC_category = "Expression",
#     Worst_group  = c("High","Low"),
#     Letter       = c("B","C"),
#     Interpretation = c("High expression = worse (Risky)","Low expression = better (Protective)")
#   ),
#   # Mutation (OFC: 2)
#   tibble(
#     OFC_category = "Mutation",
#     Worst_group  = c("MT","WT"),
#     Letter       = c("B","C"),
#     Interpretation = c("Mutant = worse (Risky)","Wild-type = better (Protective)")
#   ),
#   # CNV (OFC: 3)
#   tibble(
#     OFC_category = "CNV",
#     Worst_group  = c("Deleted","Duplicated","Normal",
#                      "Deleted and Duplicated","Duplicated and Deleted",
#                      "Normal and Deleted","Deleted and Normal",
#                      "Normal and Duplicated","Duplicated and Normal"),
#     Letter       = c("B","C","C",  # Deleted->B; Duplicated/Normal->C
#                      rep("D",6)),  # mixed patterns -> D
#     Interpretation = c(
#       "Deletion = worse (Risky)",
#       "Duplication = better (Protective)",
#       "Normal = better (Protective)",
#       rep("Mixed CNV states", 6)
#     )
#   )
# )
# 
# 
