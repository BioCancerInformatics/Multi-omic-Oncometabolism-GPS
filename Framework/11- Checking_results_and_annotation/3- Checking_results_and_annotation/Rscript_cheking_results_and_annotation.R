##### Checking results and annotation

library(dplyr)
library(tidyr)
library(rio)

setwd("E:/Higor/Cancer_metabolism_project/10- Tumor_immune_and_microenvironment_classification/")

##### CODING GENE RESULTS

# List of gene results to be imported
filepaths_gene <- c(
  "Gene/CNV/Immune_classification.tsv",
  "Gene/Methylation/Immune_classification.tsv",
  "Gene/mRNA/Immune_classification.tsv",
  "Gene/Mutation/Immune_classification.tsv"
)

# Import and combine all files 
df001_gene_results <- filepaths_gene %>%
  lapply(rio::import) %>% # Importa cada arquivo da lista
  bind_rows()             # Combina os arquivos em um único data.frame


# Remove specific columns
df002_gene_results_filtered <- df001_gene_results %>%
  select(-log10_correlation, -Expression_p.sgnif, -log10_p.adj, -immune_cells, -cor, -padj, -Infiltrate_profile)

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
df002_gene_results_filtered <- df002_gene_results_filtered %>%
  distinct()

# Replace empty strings with "no data"
df002_gene_results_filtered <- df002_gene_results_filtered %>%
  mutate(across(
    c(Expression, Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI,   
      OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group, 
      Classification, Microenvironment_score_details, Immune_classification, Immune_score_details),
    ~ ifelse(. == "", "No data", .)
  ))


# Create the column "Molecular class" with "Coding gene" in all lines
df002_gene_results_filtered <- df002_gene_results_filtered %>%
  mutate(Molecular_class = "Coding gene")

# Remover as linhas onde a coluna 'genes' é igual a "MTFMT"
df002_gene_results_filtered <- df002_gene_results_filtered[df002_gene_results_filtered$genes != "MTFMT", ]

# Import pathways information
df003_metabolic_pathways <- import("/Higor/Cancer_metabolism_project/1- Searching_for_target_coding_genes/Selected_Metabolic_Pathways_KEGG.tsv")

# Replace occurrences 
df003_metabolic_pathways <- df003_metabolic_pathways %>%
  mutate(Genes = recode(Genes, 
                        "PEDS1-UBE2V1" = "TMEM189-UBE2V1",
                        "BPNT2" = "IMPAD1",
                        "PLAAT3" = "PLA2G16",
                        "TAFAZZIN" = "TAZ"))


# Merge entre 'df002_gene_results_filtered' e 'df003_metabolic_pathways' com base no gene
df004_merge_pathways <- merge(
  df002_gene_results_filtered, df003_metabolic_pathways, 
  by.x = "genes", by.y = "Genes", 
  all.x = TRUE, # Inclui todos os genes de 'df002_gene_results_filtered'
  all.y = FALSE # Exclui genes de 'df003_metabolic_pathways' que não estão em 'df002_gene_results_filtered'
)

# Import miRNA target
df005_mirna_target <- import("E:/Cancer_metabolism_project/2- Searching_for_target_miRNA/Target_mirna_pubmed_id.tsv")

# Remove duplication of genes that interact with miRNA and concatenate the miRNA associated with the common genes
df005_mirna_target <- df005_mirna_target %>%
  group_by(Target_genes) %>%
  summarise(Mature_mirna_id = paste(unique(na.omit(Mature_mirna_id)), collapse = " / ")) %>%
  ungroup()

# Add the column "Interactions" to "df004_merge_pathways"
df006_merge_interactions <- df004_merge_pathways %>%
  left_join(df005_mirna_target, by = c("genes" = "Target_genes")) %>%
  rename(Interactions = Mature_mirna_id)

# Import lncRNA target
df007_lncRNA_target <- import("/Higor/Cancer_metabolism_project/3- Searching_for_target_lncRNA/lncRNA-related_genes.tsv")

# vetor de genes a remover
genes_remover <- c("ABTB1", "AGO2", "ANKRD40CL", "ARFRP1", "ARK2C", "BLACAT1", "BMP1", 
                   "C1QTNF1", "C20orf204", "C2orf92", "CALB2", "CCL7", "CIR1", "CXCL1", 
                   "ERVH48-1", "GLRX3", "GREP1", "IL7R", "MAP3K20", "MICAL2", "NBDY", 
                   "NBPF4", "NDRG1", "NLRP3", "PEG10", "POU3F3", "PPT2-EGFL8", "PRAC2", 
                   "PSEN1", "PTGS2", "RC3H2", "ROBO2", "SCAMP1", "SENP3-EIF4A1", "SFTA3", 
                   "SLC25A15", "SLCO4A1", "SMARCC2", "SMIM30", "SMIM31", "SPAAR", "TATDN1", 
                   "TCTN2", "TINCR", "TMEM235", "TMEM238L", "TUG1", "UFC1", "YWHAE", 
                   "ZNF800", "ZNF883")

df007_lncRNA_target <- df007_lncRNA_target %>%
  filter(!`LncRNA name` %in% genes_remover)

# Remove duplication of genes that interact with lncRNA and concatenate the lncRNA associated with the common genes
df007_lncRNA_target <- df007_lncRNA_target %>%
  group_by(`Interaction target`) %>%
  summarise(`LncRNA name` = paste(unique(na.omit(`LncRNA name`)), collapse = " / ")) %>%
  ungroup()

# Add lncRNA information to the "interactions" column
df008_merge_all_interactions <- df006_merge_interactions %>%
  left_join(df007_lncRNA_target, by = c("genes" = "Interaction target")) %>%
  mutate(Interactions = ifelse(is.na(`LncRNA name`), 
                               Interactions, 
                               ifelse(is.na(Interactions), 
                                      `LncRNA name`, 
                                      paste(Interactions, `LncRNA name`, sep = " / ")))) %>%
  select(-`LncRNA name`)  # Remover a coluna auxiliar após a adição

# Import metabolic cell death
df009_metabolic_cell_death <- import("/Higor/Cancer_metabolism_project/11- Checking_results_and_annotation/1- Searching_for_metabolic_cell_death_gene/Metabolic_cell_death_target.tsv")

df009_metabolic_cell_death <- df009_metabolic_cell_death %>%
  distinct(deathtype, gene, .keep_all = TRUE)

# Merge entre 'df008_merge_all_interactions' e 'df009_metabolic_cell_death' com base no gene
df010_merge_metabolic_cell_death <- merge(
  df008_merge_all_interactions, df009_metabolic_cell_death, 
  by.x = "genes", by.y = "gene", 
  all.x = TRUE, 
  all.y = FALSE 
)

# Import driver genes
df011_driver_gene <- import("/Higor/Cancer_metabolism_project/11- Checking_results_and_annotation/2- Searching_for_driver_genes/Driver_gene_target.tsv") 

df012_merge_driver_gene <- df010_merge_metabolic_cell_death %>%
  mutate(`Driver gene` = ifelse(genes %in% df011_driver_gene$`Gene Symbol`, "Yes", "Not"))

# Reorder the columns
New_column_order_01 <- c("genes", "Molecular_class", "Interactions", "Metabolism", "Pathways", "deathtype", "Driver gene", "cancer_types", "Expression",  "Expression_p.adj",  
                         "Genotypic_var", "Phenotypic_var", "correlation", "Correlation_p.adj", "Type_Cox_OS", "p.value_Cox_OS", "Type_Cox_DSS", 
                         "p.value_Cox_DSS", "Type_Cox_DFI", "p.value_Cox_DFI", "Type_Cox_PFI", "p.value_Cox_PFI", "OS_log_rank_chisq",              
                         "OS_p_val", "OS_worst_prognosis_group", "DSS_log_rank_chisq", "DSS_p_val", "DSS_worst_prognosis_group",      
                         "DFI_log_rank_chisq", "DFI_p_val", "DFI_worst_prognosis_group", "PFI_log_rank_chisq", "PFI_p_val", "PFI_worst_prognosis_group", 
                         "Classification", "Microenvironment_score_details", "Immune_classification",          
                         "Immune_score_details")

df013_verified_gene_results <- df012_merge_driver_gene[New_column_order_01]

# Change column names
df013_verified_gene_results <- df013_verified_gene_results %>%
  rename(
    Target = genes,
    Molecular_class = Molecular_class,
    Interactions = Interactions,
    Metabolism = Metabolism,
    Pathways = Pathways,
    Metabolic_cell_death = deathtype, 
    Driver_gene = 'Driver gene',
    Cancer_types = cancer_types,
    Tumor_vs_normal = Expression,
    Tumor_vs_normal_p.adj = Expression_p.adj,
    Omic_layer = Genotypic_var,
    Phenotypic_layer = Phenotypic_var,
    Correlation_rho = correlation,
    Correlation_p.adj = Correlation_p.adj,
    Cox_OS_type = Type_Cox_OS,
    Cox_OS_p.value = p.value_Cox_OS,
    Cox_DSS_type = Type_Cox_DSS,
    Cox_DSS_p.value = p.value_Cox_DSS,
    Cox_DFI_type = Type_Cox_DFI,
    Cox_DFI_p.value = p.value_Cox_DFI,
    Cox_PFI_type = Type_Cox_PFI,
    Cox_PFI_p.value = p.value_Cox_PFI,
    OS_log_rank_chisq = OS_log_rank_chisq,
    OS_p.Value = OS_p_val,
    OS_worst_prognosis_group = OS_worst_prognosis_group,
    DSS_log_rank_chisq = DSS_log_rank_chisq,
    DSS_p.value = DSS_p_val,
    DSS_worst_prognosis_group = DSS_worst_prognosis_group,
    DFI_log_rank_chisq = DFI_log_rank_chisq,
    DFI_p.value = DFI_p_val,
    DFI_worst_prognosis_group = DFI_worst_prognosis_group,
    PFI_log_rank_chisq = PFI_log_rank_chisq,
    PFI_p.value = PFI_p_val,
    PFI_worst_prognosis_group = PFI_worst_prognosis_group,
    Microenvironment_classification = Classification,
    Microenvironment_score_details = Microenvironment_score_details,
    Immune_classification = Immune_classification,
    Immune_score_details = Immune_score_details
  )


# Change names
df013_verified_gene_results <- df013_verified_gene_results %>%
  mutate(Phenotypic_layer = case_when(
    Phenotypic_layer == "MSI" ~ "Microsatellite instability",
    Phenotypic_layer == "TMB" ~ "Tumor mutational burden",
    TRUE ~ Phenotypic_layer # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df013_verified_gene_results <- df013_verified_gene_results %>%
  mutate(Immune_classification = case_when(
    Immune_classification == "Frio" ~ "Cold",
    Immune_classification == "Quente" ~ "Hot",
    Immune_classification == "Variável" ~ "Variable",
    TRUE ~ Immune_classification # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df013_verified_gene_results <- df013_verified_gene_results %>%
  mutate(Metabolic_cell_death = case_when(
    is.na(Metabolic_cell_death) ~ "Unrelated",  # Replace NA with "Unrelated"
    TRUE ~ Metabolic_cell_death  # Keep existing values unchanged
  ))


##### PROTEIN RESULTS

# Lista de arquivos a serem importados
filepaths_protein <- c(
  "Protein/Immune_classification.tsv"
)

# Importar e combinar todos os arquivos
df014_protein_results <- filepaths_protein %>%
  lapply(rio::import) %>% # Importa cada arquivo da lista
  bind_rows()             # Combina os arquivos em um único data.frame

# Remove specific columns
df015_protein_results_filtered <- df014_protein_results %>%
  select(-log10_correlation, -Expression_p.sgnif, -log10_p.adj, -immune_cells, -cor, -padj, -Infiltrate_profile)

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
df015_protein_results_filtered <- df015_protein_results_filtered %>%
  distinct()

# Replace empty strings with "no data"
df015_protein_results_filtered <- df015_protein_results_filtered %>%
  mutate(across(
    c(Expression, Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI,   
      OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group,
      Classification, Microenvironment_score_details, Immune_classification, Immune_score_details),
    ~ ifelse(. == "", "No data", .)
  ))

# Create the column "Molecular class" with "Coding gene" in all lines
df015_protein_results_filtered <- df015_protein_results_filtered %>%
  mutate(Molecular_class = "Protein")

# Definir um vetor com os pares de substituição
Change_name <- c("COMPLEXIISUBUNIT30" = "SDHB", "THYMIDILATESYNTHASE" = "TYMS")

# Substituir as observações no dataframe 
df015_protein_results_filtered$genes <- recode(df015_protein_results_filtered$genes, !!!Change_name)

# Merge entre 'All_results_gene' e 'Pathways' com base no gene
df016_merge_pathways <- merge(
  df015_protein_results_filtered, df003_metabolic_pathways, 
  by.x = "genes", by.y = "Genes", 
  all.x = TRUE, # Inclui todos os genes de 'All_results_gene'
  all.y = FALSE # Exclui genes de 'Pathways' que não estão em 'All_results_gene'
)

# Add the column "Interactions" to "Aggregate_gene_pathways"
df017_merge_interactions <- df016_merge_pathways %>%
  left_join(df005_mirna_target, by = c("genes" = "Target_genes")) %>%
  rename(Interactions = Mature_mirna_id)

# Add lncRNA information to the "interactions" column
df018_all_merge_interactions <- df017_merge_interactions %>%
  left_join(df007_lncRNA_target, by = c("genes" = "Interaction target")) %>%
  mutate(Interactions = ifelse(is.na(`LncRNA name`), 
                               Interactions, 
                               ifelse(is.na(Interactions), 
                                      `LncRNA name`, 
                                      paste(Interactions, `LncRNA name`, sep = " / ")))) %>%
  select(-`LncRNA name`)  # Remover a coluna auxiliar após a adição

# Merge entre 'df018_all_merge_interactions' e 'df009_metabolic_cell_death' com base no gene
df019_merge_metabolic_cell_death <- merge(
  df018_all_merge_interactions, df009_metabolic_cell_death, 
  by.x = "genes", by.y = "gene", 
  all.x = TRUE, 
  all.y = FALSE 
)

df020_merge_driver_gene <- df019_merge_metabolic_cell_death %>%
  mutate(`Driver gene` = ifelse(genes %in% df011_driver_gene$`Gene Symbol`, "Yes", "Not"))

# Definir um vetor com os pares de substituição
Change_name <- c("SDHB" = "COMPLEXIISUBUNIT30", "TYMS" = "THYMIDILATESYNTHASE")

# Substituir as observações no dataframe
df020_merge_driver_gene$genes <- recode(df020_merge_driver_gene$genes, !!!Change_name)

df021_verified_protein_results <- df020_merge_driver_gene[New_column_order_01]

# Change column names
df021_verified_protein_results <- df021_verified_protein_results %>%
  rename(
    Target = genes,
    Molecular_class = Molecular_class,
    Interactions = Interactions,
    Metabolism = Metabolism,
    Pathways = Pathways,
    Metabolic_cell_death = deathtype, 
    Driver_gene = 'Driver gene',
    Cancer_types = cancer_types,
    Tumor_vs_normal = Expression,
    Tumor_vs_normal_p.adj = Expression_p.adj,
    Omic_layer = Genotypic_var,
    Phenotypic_layer = Phenotypic_var,
    Correlation_rho = correlation,
    Correlation_p.adj = Correlation_p.adj,
    Cox_OS_type = Type_Cox_OS,
    Cox_OS_p.value = p.value_Cox_OS,
    Cox_DSS_type = Type_Cox_DSS,
    Cox_DSS_p.value = p.value_Cox_DSS,
    Cox_DFI_type = Type_Cox_DFI,
    Cox_DFI_p.value = p.value_Cox_DFI,
    Cox_PFI_type = Type_Cox_PFI,
    Cox_PFI_p.value = p.value_Cox_PFI,
    OS_log_rank_chisq = OS_log_rank_chisq,
    OS_p.Value = OS_p_val,
    OS_worst_prognosis_group = OS_worst_prognosis_group,
    DSS_log_rank_chisq = DSS_log_rank_chisq,
    DSS_p.value = DSS_p_val,
    DSS_worst_prognosis_group = DSS_worst_prognosis_group,
    DFI_log_rank_chisq = DFI_log_rank_chisq,
    DFI_p.value = DFI_p_val,
    DFI_worst_prognosis_group = DFI_worst_prognosis_group,
    PFI_log_rank_chisq = PFI_log_rank_chisq,
    PFI_p.value = PFI_p_val,
    PFI_worst_prognosis_group = PFI_worst_prognosis_group,
    Microenvironment_classification = Classification,
    Microenvironment_score_details = Microenvironment_score_details,
    Immune_classification = Immune_classification,
    Immune_score_details = Immune_score_details
  )


# Change names
df021_verified_protein_results <- df021_verified_protein_results %>%
  mutate(Phenotypic_layer = case_when(
    Phenotypic_layer == "MSI" ~ "Microsatellite instability",
    Phenotypic_layer == "TMB" ~ "Tumor mutational burden",
    TRUE ~ Phenotypic_layer # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df021_verified_protein_results <- df021_verified_protein_results %>%
  mutate(Immune_classification = case_when(
    Immune_classification == "Frio" ~ "Cold",
    Immune_classification == "Quente" ~ "Hot",
    Immune_classification == "Variável" ~ "Variable",
    TRUE ~ Immune_classification # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df021_verified_protein_results <- df021_verified_protein_results %>%
  mutate(Metabolic_cell_death = case_when(
    is.na(Metabolic_cell_death) ~ "Unrelated",  # Replace NA with "Unrelated"
    TRUE ~ Metabolic_cell_death  # Keep existing values unchanged
  ))


##### CODING GENE TRANSCRIPT RESULTS

# Lista de arquivos a serem importados
filepaths_transcript_gene <- c(
  "Transcript/Coding gene/Immune_classification.tsv"
)

# Importar e combinar todos os arquivos
df022_gene_transcript_results <- filepaths_transcript_gene %>%
  lapply(rio::import) %>% # Importa cada arquivo da lista
  bind_rows()             # Combina os arquivos em um único data.frame

# Criar a coluna "Molecular Class" com o valor "coding gene" em todas as linhas
df023_gene_transcript_results_filtered <- df022_gene_transcript_results %>%
  mutate(Molecular_class = "Coding Transcript Isoform")

# Remover as colunas específicas
df023_gene_transcript_results_filtered <- df023_gene_transcript_results_filtered %>%
  select(-log10_correlation, -Expression_p.sgnif, -log10_p.adj, -immune_cells, -cor, -padj, -Infiltrate_profile)

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
df023_gene_transcript_results_filtered <- df023_gene_transcript_results_filtered %>%
  distinct()

# Substituir strings vazias por "no data" 
df023_gene_transcript_results_filtered <- df023_gene_transcript_results_filtered %>%
  mutate(across(
    c(Expression, Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI,   
      OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group,
      Classification, Microenvironment_score_details, Immune_classification, Immune_score_details),
    ~ ifelse(. == "", "No data", .)
  ))

# Import transcript isoforms
df024_gene_transcript_target <- import("/Higor/Cancer_metabolism_project/4- Searching_for_target_transcript_isoforms/Coding gene/gene_info.csv")

# Adicionar a coluna "Interactions" ao 'gene_transcript_results'
df025_gene_transcript_annotation <- df023_gene_transcript_results_filtered %>%
  left_join(
    df024_gene_transcript_target %>% select(Transcript_ID, Display_Name),  # Seleciona apenas as colunas necessárias
    by = c("genes" = "Transcript_ID")  # Faz a correspondência correta entre 'genes' e 'Transcript_ID'
  ) %>%
  rename(TranscriptGene = Display_Name)  # Renomeia a coluna para 'Interactions'

# Remover as linhas onde Interactions é "MTFMT"
df025_gene_transcript_annotation <- df025_gene_transcript_annotation[df025_gene_transcript_annotation$TranscriptGene != "MTFMT", ]

# Replace occurrences of "PEDS1-UBE2V1" with "TMEM189-UBE2V1"
df025_gene_transcript_annotation <- df025_gene_transcript_annotation %>%
  mutate(TranscriptGene = ifelse(TranscriptGene == "PEDS1-UBE2V1", "TMEM189-UBE2V1", TranscriptGene))

# Merge entre 'All_results_gene' e 'Pathways_lnc_gene' com base no gene
df026_merge_pathways <- merge(
  df025_gene_transcript_annotation, df003_metabolic_pathways, 
  by.x = "TranscriptGene", by.y = "Genes", 
  all.x = TRUE, # Inclui todos os genes de 'All_results_gene'
  all.y = FALSE # Exclui genes de 'Pathways' que não estão em 'All_results_gene'
)

# Add the column "Interactions" to "df026_merge_pathways"
df027_merge_interactions <- df026_merge_pathways %>%
  left_join(df005_mirna_target, by = c("TranscriptGene" = "Target_genes")) %>%
  rename(Interactions = Mature_mirna_id)

# Add lncRNA information to the "interactions" column
df028_merge_all_interactions <- df027_merge_interactions %>%
  left_join(df007_lncRNA_target, by = c("TranscriptGene" = "Interaction target")) %>%
  mutate(Interactions = ifelse(is.na(`LncRNA name`), 
                               Interactions, 
                               ifelse(is.na(Interactions), 
                                      `LncRNA name`, 
                                      paste(Interactions, `LncRNA name`, sep = " / ")))) %>%
  select(-`LncRNA name`)  # Remover a coluna auxiliar após a adição

# Merge entre 'df028_merge_all_interactions' e 'df009_metabolic_cell_death' com base no gene
df029_merge_metabolic_cell_death <- merge(
  df028_merge_all_interactions, df009_metabolic_cell_death, 
  by.x = "TranscriptGene", by.y = "gene", 
  all.x = TRUE, 
  all.y = FALSE 
)

df030_merge_driver_gene <- df029_merge_metabolic_cell_death %>%
  mutate(`Driver gene` = ifelse(TranscriptGene %in% df011_driver_gene$`Gene Symbol`, "Yes", "Not"))

# Reorder the columns
New_column_order_02 <- c("genes", "TranscriptGene", "Molecular_class", "Interactions", "Metabolism", "Pathways", "deathtype", "Driver gene", "cancer_types", "Expression",  "Expression_p.adj",  
                         "Genotypic_var", "Phenotypic_var", "correlation", "Correlation_p.adj", "Type_Cox_OS", "p.value_Cox_OS", "Type_Cox_DSS", 
                         "p.value_Cox_DSS", "Type_Cox_DFI", "p.value_Cox_DFI", "Type_Cox_PFI", "p.value_Cox_PFI", "OS_log_rank_chisq",              
                         "OS_p_val", "OS_worst_prognosis_group", "DSS_log_rank_chisq", "DSS_p_val", "DSS_worst_prognosis_group",      
                         "DFI_log_rank_chisq", "DFI_p_val", "DFI_worst_prognosis_group", "PFI_log_rank_chisq", "PFI_p_val", "PFI_worst_prognosis_group", 
                         "Classification", "Microenvironment_score_details", "Immune_classification",          
                         "Immune_score_details")

df031_verified_gene_transcript_results <- df030_merge_driver_gene[New_column_order_02]

# Change column names
df031_verified_gene_transcript_results <- df031_verified_gene_transcript_results %>%
  rename(
    Target = genes,
    TranscriptGene = TranscriptGene,
    Molecular_class = Molecular_class,
    Interactions = Interactions,
    Metabolism = Metabolism,
    Pathways = Pathways,
    Metabolic_cell_death = deathtype, 
    Driver_gene = 'Driver gene',
    Cancer_types = cancer_types,
    Tumor_vs_normal = Expression,
    Tumor_vs_normal_p.adj = Expression_p.adj,
    Omic_layer = Genotypic_var,
    Phenotypic_layer = Phenotypic_var,
    Correlation_rho = correlation,
    Correlation_p.adj = Correlation_p.adj,
    Cox_OS_type = Type_Cox_OS,
    Cox_OS_p.value = p.value_Cox_OS,
    Cox_DSS_type = Type_Cox_DSS,
    Cox_DSS_p.value = p.value_Cox_DSS,
    Cox_DFI_type = Type_Cox_DFI,
    Cox_DFI_p.value = p.value_Cox_DFI,
    Cox_PFI_type = Type_Cox_PFI,
    Cox_PFI_p.value = p.value_Cox_PFI,
    OS_log_rank_chisq = OS_log_rank_chisq,
    OS_p.Value = OS_p_val,
    OS_worst_prognosis_group = OS_worst_prognosis_group,
    DSS_log_rank_chisq = DSS_log_rank_chisq,
    DSS_p.value = DSS_p_val,
    DSS_worst_prognosis_group = DSS_worst_prognosis_group,
    DFI_log_rank_chisq = DFI_log_rank_chisq,
    DFI_p.value = DFI_p_val,
    DFI_worst_prognosis_group = DFI_worst_prognosis_group,
    PFI_log_rank_chisq = PFI_log_rank_chisq,
    PFI_p.value = PFI_p_val,
    PFI_worst_prognosis_group = PFI_worst_prognosis_group,
    Microenvironment_classification = Classification,
    Microenvironment_score_details = Microenvironment_score_details,
    Immune_classification = Immune_classification,
    Immune_score_details = Immune_score_details
  )

# Change names
df031_verified_gene_transcript_results <- df031_verified_gene_transcript_results %>%
  mutate(Phenotypic_layer = case_when(
    Phenotypic_layer == "MSI" ~ "Microsatellite instability",
    Phenotypic_layer == "TMB" ~ "Tumor mutational burden",
    TRUE ~ Phenotypic_layer # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df031_verified_gene_transcript_results <- df031_verified_gene_transcript_results %>%
  mutate(Immune_classification = case_when(
    Immune_classification == "Frio" ~ "Cold",
    Immune_classification == "Quente" ~ "Hot",
    Immune_classification == "Variável" ~ "Variable",
    TRUE ~ Immune_classification # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df031_verified_gene_transcript_results <- df031_verified_gene_transcript_results %>%
  mutate(Metabolic_cell_death = case_when(
    is.na(Metabolic_cell_death) ~ "Unrelated",  # Replace NA with "Unrelated"
    TRUE ~ Metabolic_cell_death  # Keep existing values unchanged
  ))


##### MIRNA RESULTS

# List of gene results to be imported
filepaths_mirna <- c(
  "miRNA/Immune_classification.tsv"
)

# Import and combine all files
df032_mirna_results <- filepaths_mirna %>%
  lapply(rio::import) %>% # Importa cada arquivo da lista
  bind_rows()             # Combina os arquivos em um único data.frame

# Remove specific columns
df033_mirna_results_filtered <- df032_mirna_results %>%
  select(-log10_correlation, -Expression_p.sgnif, -log10_p.adj, -immune_cells, -cor, -padj, -Infiltrate_profile)

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
df033_mirna_results_filtered <- df033_mirna_results_filtered %>%
  distinct()

# Replace empty strings with "no data" 
df033_mirna_results_filtered <- df033_mirna_results_filtered %>%
  mutate(across(
    c(Expression, Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI,   
      OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group,
      Classification, Microenvironment_score_details, Immune_classification, Immune_score_details),
    ~ ifelse(. == "", "No data", .)
  ))

# Create the column "Molecular class" with "miRNA" in all lines
df033_mirna_results_filtered <- df033_mirna_results_filtered %>%
  mutate(Molecular_class = "miRNA")

# Import miRNA target
df034_mirna_target <- import("E:/Cancer_metabolism_project/2- Searching_for_target_miRNA/Target_mirna_pubmed_id.tsv")

# Substituir todas as ocorrências de "ATP5MGL" por "ATP5L2"
df034_mirna_target$Target_genes <- gsub("\\bATP5MGL\\b", "ATP5L2", df034_mirna_target$Target_genes)

# Remover as linhas onde Target_genes é "MTFMT" ou "pten"
df034_mirna_target <- df034_mirna_target[!df034_mirna_target$Target_genes %in% c("MTFMT", "pten"), ]

# Remove duplication of miRNA that interact with genes and concatenate the genes associated with the common miRNA
df034_mirna_target <- df034_mirna_target %>%
  group_by(Mature_mirna_id) %>%
  summarise(Target_genes = paste(unique(na.omit(Target_genes)), collapse = " / ")) %>%
  ungroup()

# Add the column "Interactions" 
df035_merge_interactions <- df033_mirna_results_filtered %>%
  left_join(df034_mirna_target, by = c("genes" = "Mature_mirna_id")) %>%
  rename(Interactions = Target_genes)

# Transformar interações concatenadas em linhas individuais
df035_merge_interactions <- df035_merge_interactions %>%
  separate_rows(Interactions, sep = " / ") # Divide as interações e cria novas linhas

# Merge entre 'gene_results' e 'Pathways' com base no gene
df036_merge_pathways <- merge(
  df035_merge_interactions, df003_metabolic_pathways, 
  by.x = "Interactions", by.y = "Genes", 
  all.x = TRUE, # Inclui todos os genes de 'gene_results'
  all.y = FALSE # Exclui genes de 'Pathways' que não estão em 'All_results_gene'
)


# Store the original Target column before modifications
df036_merge_pathways$genes_original <- df036_merge_pathways$genes

# Apply transformations to Target column
df036_merge_pathways$genes <- df036_merge_pathways$genes %>%
  sub("^hsa-", "", .) %>%       # Remove "hsa-" if present
  sub("miR", "MIR", .) %>%      # Replace "miR" with "MIR"
  gsub("-", "", .) %>%          # Remove all hyphens
  sub("(5p|3p)$", "", .)        # Remove only "-5p" or "-3p" at the end

# Merge entre 'df037_merge_metabolic_cell_death' e 'df009_metabolic_cell_death' com base no gene
df037_merge_metabolic_cell_death <- merge(
  df036_merge_pathways, df009_metabolic_cell_death, 
  by.x = "genes", by.y = "gene", 
  all.x = TRUE, 
  all.y = FALSE 
)

# Merge driver gene info 
df038_merge_driver_gene <- df037_merge_metabolic_cell_death %>%
  mutate(`Driver gene` = ifelse(genes %in% df011_driver_gene$`Gene Symbol`, "Yes", "Not"))

# return the original values (hsa)
df038_merge_driver_gene$genes <- df038_merge_driver_gene$genes_original

# Agrupar pelas colunas desejadas e concatenar 'interactions', mantendo as demais colunas
df039_merge_interactions <- df038_merge_driver_gene %>%
  group_by(genes, Molecular_class, Metabolism, Pathways, cancer_types, Genotypic_var, Phenotypic_var, deathtype) %>%
  summarise(
    Interactions = paste(unique(Interactions), collapse = " / "),
    across(everything(), ~ first(.)),  # Mantém os valores das outras colunas
    .groups = "drop"
  )

df040_verified_mirna_results <- df039_merge_interactions[New_column_order_01]

# Change column names
df040_verified_mirna_results <- df040_verified_mirna_results %>%
  rename(
    Target = genes,
    Molecular_class = Molecular_class,
    Interactions = Interactions,
    Metabolism = Metabolism,
    Pathways = Pathways,
    Metabolic_cell_death = deathtype, 
    Driver_gene = 'Driver gene',
    Cancer_types = cancer_types,
    Tumor_vs_normal = Expression,
    Tumor_vs_normal_p.adj = Expression_p.adj,
    Omic_layer = Genotypic_var,
    Phenotypic_layer = Phenotypic_var,
    Correlation_rho = correlation,
    Correlation_p.adj = Correlation_p.adj,
    Cox_OS_type = Type_Cox_OS,
    Cox_OS_p.value = p.value_Cox_OS,
    Cox_DSS_type = Type_Cox_DSS,
    Cox_DSS_p.value = p.value_Cox_DSS,
    Cox_DFI_type = Type_Cox_DFI,
    Cox_DFI_p.value = p.value_Cox_DFI,
    Cox_PFI_type = Type_Cox_PFI,
    Cox_PFI_p.value = p.value_Cox_PFI,
    OS_log_rank_chisq = OS_log_rank_chisq,
    OS_p.Value = OS_p_val,
    OS_worst_prognosis_group = OS_worst_prognosis_group,
    DSS_log_rank_chisq = DSS_log_rank_chisq,
    DSS_p.value = DSS_p_val,
    DSS_worst_prognosis_group = DSS_worst_prognosis_group,
    DFI_log_rank_chisq = DFI_log_rank_chisq,
    DFI_p.value = DFI_p_val,
    DFI_worst_prognosis_group = DFI_worst_prognosis_group,
    PFI_log_rank_chisq = PFI_log_rank_chisq,
    PFI_p.value = PFI_p_val,
    PFI_worst_prognosis_group = PFI_worst_prognosis_group,
    Microenvironment_classification = Classification,
    Microenvironment_score_details = Microenvironment_score_details,
    Immune_classification = Immune_classification,
    Immune_score_details = Immune_score_details
  )


# Change names
df040_verified_mirna_results <- df040_verified_mirna_results %>%
  mutate(Phenotypic_layer = case_when(
    Phenotypic_layer == "MSI" ~ "Microsatellite instability",
    Phenotypic_layer == "TMB" ~ "Tumor mutational burden",
    TRUE ~ Phenotypic_layer # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df040_verified_mirna_results <- df040_verified_mirna_results %>%
  mutate(Immune_classification = case_when(
    Immune_classification == "Frio" ~ "Cold",
    Immune_classification == "Quente" ~ "Hot",
    Immune_classification == "Variável" ~ "Variable",
    TRUE ~ Immune_classification # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df040_verified_mirna_results <- df040_verified_mirna_results %>%
  mutate(Metabolic_cell_death = case_when(
    is.na(Metabolic_cell_death) ~ "Unrelated",  # Replace NA with "Unrelated"
    TRUE ~ Metabolic_cell_death  # Keep existing values unchanged
  ))


##### LNCRNA RESULTS

# Lista de arquivos a serem importados
filepaths_lncrna <- c(
  "lncRNA/CNV/Immune_classification.tsv",
  "lncRNA/Expression/Immune_classification.tsv",
  "lncRNA/Methylation/Immune_classification.tsv",
  "lncRNA/Mutation/Immune_classification.tsv"
)

# Importar e combinar todos os arquivos
df041_lncrna_results <- filepaths_lncrna %>%
  lapply(rio::import) %>% # Importa cada arquivo da lista
  bind_rows()             # Combina os arquivos em um único data.frame

# Import lncRNA target
df042_lncRNA_target <- import("/Higor/Cancer_metabolism_project/3- Searching_for_target_lncRNA/lncRNA-related_genes.tsv")

# vetor de genes a remover
genes_remover <- c("ABTB1", "AGO2", "ANKRD40CL", "ARFRP1", "ARK2C", "BLACAT1", "BMP1", 
                   "C1QTNF1", "C20orf204", "C2orf92", "CALB2", "CCL7", "CIR1", "CXCL1", 
                   "ERVH48-1", "GLRX3", "GREP1", "IL7R", "MAP3K20", "MICAL2", "NBDY", 
                   "NBPF4", "NDRG1", "NLRP3", "PEG10", "POU3F3", "PPT2-EGFL8", "PRAC2", 
                   "PSEN1", "PTGS2", "RC3H2", "ROBO2", "SCAMP1", "SENP3-EIF4A1", "SFTA3", 
                   "SLC25A15", "SLCO4A1", "SMARCC2", "SMIM30", "SMIM31", "SPAAR", "TATDN1", 
                   "TCTN2", "TINCR", "TMEM235", "TMEM238L", "TUG1", "UFC1", "YWHAE", 
                   "ZNF800", "ZNF883")

df042_lncRNA_target <- df042_lncRNA_target %>%
  filter(!`LncRNA name` %in% genes_remover)

# Remove duplication of genes that interact with lncRNA and concatenate the lncRNA associated with the common genes
df042_lncRNA_target <- df042_lncRNA_target %>%
  group_by(`LncRNA name`) %>%
  summarise(`Interaction target` = paste(unique(na.omit(`Interaction target`)), collapse = " / ")) %>%
  ungroup()

# Filtrando xxx para manter apenas as linhas onde 'genes' está presente em 'lncrna' de yyy
df043_lncrna_results_filtered <- df041_lncrna_results %>%
  filter(genes %in% df042_lncRNA_target$`LncRNA name`)

# Remover as colunas específicas
df043_lncrna_results_filtered <- df043_lncrna_results_filtered %>%
  select(-log10_correlation, -Expression_p.sgnif, -log10_p.adj, -immune_cells, -cor, -padj, -Infiltrate_profile)

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
df043_lncrna_results_filtered <- df043_lncrna_results_filtered %>%
  distinct()

# Substituir strings vazias por "no data" 
df043_lncrna_results_filtered <- df043_lncrna_results_filtered %>%
  mutate(across(
    c(Expression, Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI,   
      OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group,
      Classification, Microenvironment_score_details, Immune_classification, Immune_score_details),
    ~ ifelse(. == "", "No data", .)
  ))

# Criar a coluna "Molecular Class" com o valor "coding gene" em todas as linhas
df043_lncrna_results_filtered <- df043_lncrna_results_filtered %>%
  mutate(Molecular_class = "lncRNA")

# Add lncRNA information to the "interactions" column
df044_merge_interactions <- df043_lncrna_results_filtered %>%
  mutate(Interactions = NA_character_) %>%  # Cria a coluna 'Interactions' com valores NA
  left_join(df042_lncRNA_target, by = c("genes" = "LncRNA name")) %>%  # Faz a junção
  mutate(Interactions = ifelse(is.na(Interactions), 
                               `Interaction target`, 
                               paste(Interactions, `Interaction target`, sep = " / "))) %>%
  select(-`Interaction target`)  # Remove a coluna 'Interaction target' após o uso

# Transformar interações concatenadas em linhas individuais
df044_merge_interactions <- df044_merge_interactions %>%
  separate_rows(Interactions, sep = " / ") # Divide as interações e cria novas linhas

# Merge entre 'gene_results' e 'Pathways' com base no gene
df045_merge_pathways <- merge(
  df044_merge_interactions, df003_metabolic_pathways, 
  by.x = "Interactions", by.y = "Genes", 
  all.x = TRUE, # Inclui todos os genes de 'gene_results'
  all.y = FALSE # Exclui genes de 'Pathways' que não estão em 'All_results_gene'
)

# Merge entre 'df037_merge_metabolic_cell_death' e 'df009_metabolic_cell_death' com base no gene
df046_merge_metabolic_cell_death <- merge(
  df045_merge_pathways, df009_metabolic_cell_death, 
  by.x = "genes", by.y = "gene", 
  all.x = TRUE, 
  all.y = FALSE 
)

# Merge driver gene info 
df047_merge_driver_gene <- df046_merge_metabolic_cell_death %>%
  mutate(`Driver gene` = ifelse(genes %in% df011_driver_gene$`Gene Symbol`, "Yes", "Not"))


# Agrupar pelas colunas desejadas e concatenar 'interactions', mantendo as demais colunas
df048_merge_interactions <- df047_merge_driver_gene %>%
  group_by(genes, Molecular_class, Metabolism, Pathways, cancer_types, Genotypic_var, Phenotypic_var, deathtype) %>%
  summarise(
    Interactions = paste(unique(Interactions), collapse = " / "),
    across(everything(), ~ first(.)),  # Mantém os valores das outras colunas
    .groups = "drop"
  )

df049_verified_lncrna_results <- df048_merge_interactions[New_column_order_01]

# Change column names
df049_verified_lncrna_results <- df049_verified_lncrna_results %>%
  rename(
    Target = genes,
    Molecular_class = Molecular_class,
    Interactions = Interactions,
    Metabolism = Metabolism,
    Pathways = Pathways,
    Metabolic_cell_death = deathtype, 
    Driver_gene = 'Driver gene',
    Cancer_types = cancer_types,
    Tumor_vs_normal = Expression,
    Tumor_vs_normal_p.adj = Expression_p.adj,
    Omic_layer = Genotypic_var,
    Phenotypic_layer = Phenotypic_var,
    Correlation_rho = correlation,
    Correlation_p.adj = Correlation_p.adj,
    Cox_OS_type = Type_Cox_OS,
    Cox_OS_p.value = p.value_Cox_OS,
    Cox_DSS_type = Type_Cox_DSS,
    Cox_DSS_p.value = p.value_Cox_DSS,
    Cox_DFI_type = Type_Cox_DFI,
    Cox_DFI_p.value = p.value_Cox_DFI,
    Cox_PFI_type = Type_Cox_PFI,
    Cox_PFI_p.value = p.value_Cox_PFI,
    OS_log_rank_chisq = OS_log_rank_chisq,
    OS_p.Value = OS_p_val,
    OS_worst_prognosis_group = OS_worst_prognosis_group,
    DSS_log_rank_chisq = DSS_log_rank_chisq,
    DSS_p.value = DSS_p_val,
    DSS_worst_prognosis_group = DSS_worst_prognosis_group,
    DFI_log_rank_chisq = DFI_log_rank_chisq,
    DFI_p.value = DFI_p_val,
    DFI_worst_prognosis_group = DFI_worst_prognosis_group,
    PFI_log_rank_chisq = PFI_log_rank_chisq,
    PFI_p.value = PFI_p_val,
    PFI_worst_prognosis_group = PFI_worst_prognosis_group,
    Microenvironment_classification = Classification,
    Microenvironment_score_details = Microenvironment_score_details,
    Immune_classification = Immune_classification,
    Immune_score_details = Immune_score_details
  )


# Change names
df049_verified_lncrna_results <- df049_verified_lncrna_results %>%
  mutate(Phenotypic_layer = case_when(
    Phenotypic_layer == "MSI" ~ "Microsatellite instability",
    Phenotypic_layer == "TMB" ~ "Tumor mutational burden",
    TRUE ~ Phenotypic_layer # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df049_verified_lncrna_results <- df049_verified_lncrna_results %>%
  mutate(Immune_classification = case_when(
    Immune_classification == "Frio" ~ "Cold",
    Immune_classification == "Quente" ~ "Hot",
    Immune_classification == "Variável" ~ "Variable",
    TRUE ~ Immune_classification # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df049_verified_lncrna_results <- df049_verified_lncrna_results %>%
  mutate(Metabolic_cell_death = case_when(
    is.na(Metabolic_cell_death) ~ "Unrelated",  # Replace NA with "Unrelated"
    TRUE ~ Metabolic_cell_death  # Keep existing values unchanged
  ))


##### LNCRNA TRANSCRIPT RESULTS

# Lista de arquivos a serem importados
filepaths_transcript_lnc <- c(
  "Transcript/lncRNA/Immune_classification.tsv"
)

# Importar e combinar todos os arquivos
df050_lncrna_transcript_results <- filepaths_transcript_lnc %>%
  lapply(rio::import) %>% # Importa cada arquivo da lista
  bind_rows()             # Combina os arquivos em um único data.frame

# Criar a coluna "Molecular Class" com o valor "coding gene" em todas as linhas
df051_lncrna_transcript_results_filtered <- df050_lncrna_transcript_results %>%
  mutate(Molecular_class = "Long non-coding Transcript Isoforms")

# Remover as colunas específicas
df051_lncrna_transcript_results_filtered <- df051_lncrna_transcript_results_filtered %>%
  select(-log10_correlation, -Expression_p.sgnif, -log10_p.adj, -immune_cells, -cor, -padj, -Infiltrate_profile)

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
df051_lncrna_transcript_results_filtered <- df051_lncrna_transcript_results_filtered %>%
  distinct()

# Substituir strings vazias por "no data" 
df051_lncrna_transcript_results_filtered <- df051_lncrna_transcript_results_filtered %>%
  mutate(across(
    c(Expression, Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI,   
      OS_worst_prognosis_group, DSS_worst_prognosis_group, DFI_worst_prognosis_group, PFI_worst_prognosis_group,
      Classification, Microenvironment_score_details, Immune_classification, Immune_score_details),
    ~ ifelse(. == "", "No data", .)
  ))

# Import target
df052_lncrna_transcript_target <- import("E:/Higor/Cancer_metabolism_project/4- Searching_for_target_transcript_isoforms/lncRNA/gene_info.csv")

df052_lncrna_transcript_target <- df052_lncrna_transcript_target %>%
  semi_join(df042_lncRNA_target, by = c("Display_Name" = "LncRNA name"))

# vetor de genes a remover
genes_remover <- c("ABTB1", "AGO2", "ANKRD40CL", "ARFRP1", "ARK2C", "BLACAT1", "BMP1", 
                   "C1QTNF1", "C20orf204", "C2orf92", "CALB2", "CCL7", "CIR1", "CXCL1", 
                   "ERVH48-1", "GLRX3", "GREP1", "IL7R", "MAP3K20", "MICAL2", "NBDY", 
                   "NBPF4", "NDRG1", "NLRP3", "PEG10", "POU3F3", "PPT2-EGFL8", "PRAC2", 
                   "PSEN1", "PTGS2", "RC3H2", "ROBO2", "SCAMP1", "SENP3-EIF4A1", "SFTA3", 
                   "SLC25A15", "SLCO4A1", "SMARCC2", "SMIM30", "SMIM31", "SPAAR", "TATDN1", 
                   "TCTN2", "TINCR", "TMEM235", "TMEM238L", "TUG1", "UFC1", "YWHAE", 
                   "ZNF800", "ZNF883")

df052_lncrna_transcript_target <- df052_lncrna_transcript_target %>%
  filter(!Display_Name %in% genes_remover)

# Remover duplicatas baseadas na coluna Transcript_ID
df052_lncrna_transcript_target <- df052_lncrna_transcript_target %>%
  distinct(Transcript_ID, .keep_all = TRUE)  # Remove duplicatas, mantendo as outras colunas


df051_lncrna_transcript_results_filtered <- df051_lncrna_transcript_results_filtered %>%
  semi_join(df052_lncrna_transcript_target, by = c("genes" = "Transcript_ID"))

# Adicionar a coluna "TranscriptGene" ao 'transcript_lnc_results'
df053_lncrna_transcript_annotation <- df051_lncrna_transcript_results_filtered %>%
  left_join(
    df052_lncrna_transcript_target %>% select(Transcript_ID, Display_Name),  # Seleciona apenas as colunas necessárias
    by = c("genes" = "Transcript_ID")  # Faz a correspondência correta entre 'genes' e 'Transcript_ID'
  ) %>%
  rename(TranscriptGene = Display_Name)  # Renomeia a coluna para 'TranscriptGene'

# Remover as linhas onde Interactions é "MTFMT"
df053_lncrna_transcript_annotation <- df053_lncrna_transcript_annotation[df053_lncrna_transcript_annotation$TranscriptGene != "MTFMT", ]

# Filtrando xxx para manter apenas as linhas onde 'genes' está presente em 'lncrna' de yyy
df054_lncrna_transcript_filtered <- df053_lncrna_transcript_annotation %>%
  filter(TranscriptGene %in% df042_lncRNA_target$`LncRNA name`)

# Merge entre 'All_results_gene' e 'Pathways' com base no gene
df055_merge_pathways <- merge(
  df007_lncRNA_target, df003_metabolic_pathways, 
  by.x = "Interaction target", by.y = "Genes", 
  all.x = TRUE, # Inclui todos os genes de 'All_results_gene'
  all.y = FALSE # Exclui genes de 'Pathways' que não estão em 'All_results_gene'
)

# Expandir as observações da coluna "Interaction target" em linhas separadas
df055_merge_pathways <- df055_merge_pathways %>%
  separate_rows(`LncRNA name`, sep = " / ")

# Change column names
df056_merge_interactions <- df055_merge_pathways %>%
  rename(
    Genes = `LncRNA name`,
    Interactions = `Interaction target`,
  )

# Selecionar apenas as colunas de interesse do data frame 'lncrna_mirna_pathways'
df056_merge_interactions <- df056_merge_interactions %>%
  select(Genes, Metabolism, Pathways)

# Merge entre 'All_results_gene' e 'Pathways_lnc_gene' com base no gene
df057_merge_pathways <- merge(
  df054_lncrna_transcript_filtered, df056_merge_interactions, 
  by.x = "TranscriptGene", by.y = "Genes", 
  all.x = TRUE, # Inclui todos os genes de 'All_results_gene'
  all.y = FALSE # Exclui genes de 'Pathways' que não estão em 'All_results_gene'
)

# Add the column "Interactions" to "df057_merge_pathways"
df058_merge_interactions <- df057_merge_pathways %>%
  left_join(df042_lncRNA_target, by = c("TranscriptGene" = "LncRNA name")) %>%
  rename(Interactions = 'Interaction target')

# Merge entre 'df028_merge_all_interactions' e 'df009_metabolic_cell_death' com base no gene
df059_merge_metabolic_cell_death <- merge(
  df058_merge_interactions, df009_metabolic_cell_death, 
  by.x = "TranscriptGene", by.y = "gene", 
  all.x = TRUE, 
  all.y = FALSE 
)

df060_merge_driver_gene <- df059_merge_metabolic_cell_death %>%
  mutate(`Driver gene` = ifelse(TranscriptGene %in% df011_driver_gene$`Gene Symbol`, "Yes", "Not"))

# Reorder the columns
df061_verified_lncrna_transcript_results <- df060_merge_driver_gene[New_column_order_02]

# Change column names
df061_verified_lncrna_transcript_results <- df061_verified_lncrna_transcript_results %>%
  rename(
    Target = genes,
    TranscriptGene = TranscriptGene,
    Molecular_class = Molecular_class,
    Interactions = Interactions,
    Metabolism = Metabolism,
    Pathways = Pathways,
    Metabolic_cell_death = deathtype, 
    Driver_gene = 'Driver gene',
    Cancer_types = cancer_types,
    Tumor_vs_normal = Expression,
    Tumor_vs_normal_p.adj = Expression_p.adj,
    Omic_layer = Genotypic_var,
    Phenotypic_layer = Phenotypic_var,
    Correlation_rho = correlation,
    Correlation_p.adj = Correlation_p.adj,
    Cox_OS_type = Type_Cox_OS,
    Cox_OS_p.value = p.value_Cox_OS,
    Cox_DSS_type = Type_Cox_DSS,
    Cox_DSS_p.value = p.value_Cox_DSS,
    Cox_DFI_type = Type_Cox_DFI,
    Cox_DFI_p.value = p.value_Cox_DFI,
    Cox_PFI_type = Type_Cox_PFI,
    Cox_PFI_p.value = p.value_Cox_PFI,
    OS_log_rank_chisq = OS_log_rank_chisq,
    OS_p.Value = OS_p_val,
    OS_worst_prognosis_group = OS_worst_prognosis_group,
    DSS_log_rank_chisq = DSS_log_rank_chisq,
    DSS_p.value = DSS_p_val,
    DSS_worst_prognosis_group = DSS_worst_prognosis_group,
    DFI_log_rank_chisq = DFI_log_rank_chisq,
    DFI_p.value = DFI_p_val,
    DFI_worst_prognosis_group = DFI_worst_prognosis_group,
    PFI_log_rank_chisq = PFI_log_rank_chisq,
    PFI_p.value = PFI_p_val,
    PFI_worst_prognosis_group = PFI_worst_prognosis_group,
    Microenvironment_classification = Classification,
    Microenvironment_score_details = Microenvironment_score_details,
    Immune_classification = Immune_classification,
    Immune_score_details = Immune_score_details
  )

# Change names
df061_verified_lncrna_transcript_results <- df061_verified_lncrna_transcript_results %>%
  mutate(Phenotypic_layer = case_when(
    Phenotypic_layer == "MSI" ~ "Microsatellite instability",
    Phenotypic_layer == "TMB" ~ "Tumor mutational burden",
    TRUE ~ Phenotypic_layer # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df061_verified_lncrna_transcript_results <- df061_verified_lncrna_transcript_results %>%
  mutate(Immune_classification = case_when(
    Immune_classification == "Frio" ~ "Cold",
    Immune_classification == "Quente" ~ "Hot",
    Immune_classification == "Variável" ~ "Variable",
    TRUE ~ Immune_classification # Mantém os valores originais caso não sejam alterados
  ))

# Change names
df061_verified_lncrna_transcript_results <- df061_verified_lncrna_transcript_results %>%
  mutate(Metabolic_cell_death = case_when(
    is.na(Metabolic_cell_death) ~ "Unrelated",  # Replace NA with "Unrelated"
    TRUE ~ Metabolic_cell_death  # Keep existing values unchanged
  ))


##### ALL gene_protein_mirna_lncrna_results

df062_all_gene_protein_mirna_lncrna_results <- rbind(df013_verified_gene_results, df021_verified_protein_results, df040_verified_mirna_results, df049_verified_lncrna_results)

df063_filtered_gene_protein_mirna_lncrna_results <- df062_all_gene_protein_mirna_lncrna_results %>%
  mutate(across(c(Tumor_vs_normal, Tumor_vs_normal_p.adj), as.character))

df063_filtered_gene_protein_mirna_lncrna_results[df063_filtered_gene_protein_mirna_lncrna_results$Omic_layer %in% c("Methylation", "CNV", "Mutation", "Protein expression"), 
                                                 c("Tumor_vs_normal", "Tumor_vs_normal_p.adj")] <- "No data"

df063_filtered_gene_protein_mirna_lncrna_results <- df063_filtered_gene_protein_mirna_lncrna_results %>%
  mutate(Interactions = case_when(
    is.na(Interactions) ~ "No interactions",  # Replace NA with "Unrelated"
    TRUE ~ Interactions  # Keep existing values unchanged
  ))

df063_filtered_gene_protein_mirna_lncrna_results$Pathways <- gsub(" - Homo sapiens \\(human\\)", "", df063_filtered_gene_protein_mirna_lncrna_results$Pathways)

df063_filtered_gene_protein_mirna_lncrna_results <- df063_filtered_gene_protein_mirna_lncrna_results %>%
  arrange(Cancer_types, Target)

df063_filtered_gene_protein_mirna_lncrna_results <- df063_filtered_gene_protein_mirna_lncrna_results %>% 
  distinct()

setwd("E:/Higor/Cancer_metabolism_project/11- Checking_results_and_annotation/3- Checking_results_and_annotation/")

export(df063_filtered_gene_protein_mirna_lncrna_results, "df063_gene_protein_mirna_lncrna_results.tsv")

saveRDS(df063_filtered_gene_protein_mirna_lncrna_results, "df063_gene_protein_mirna_lncrna_results.rds")

##### ALL transcript_results

df064_all_transcript_results <- rbind(df031_verified_gene_transcript_results,  df061_verified_lncrna_transcript_results)

df065_filtered_transcript_results <- df064_all_transcript_results %>%
  mutate(across(c(Tumor_vs_normal, Tumor_vs_normal_p.adj), as.character))

df065_filtered_transcript_results[df065_filtered_transcript_results$Omic_layer %in% c("Methylation", "CNV", "Mutation", "Protein expression"), 
                                  c("Tumor_vs_normal", "Tumor_vs_normal_p.adj")] <- "No data"

df065_filtered_transcript_results$Pathways <- gsub(" - Homo sapiens \\(human\\)", "", df065_filtered_transcript_results$Pathways)

df065_filtered_transcript_results <- df065_filtered_transcript_results %>%
  mutate(Interactions = case_when(
    is.na(Interactions) ~ "No interactions",  # Replace NA with "Unrelated"
    TRUE ~ Interactions  # Keep existing values unchanged
  ))

df065_filtered_transcript_results <- df065_filtered_transcript_results %>%
  arrange(Cancer_types, Target)

df065_filtered_transcript_results <- df065_filtered_transcript_results %>% 
  distinct()

setwd("E:/Higor/Cancer_metabolism_project/11- Checking_results_and_annotation/3- Checking_results_and_annotation/")

export(df065_filtered_transcript_results, "df065_transcript_results.tsv")

saveRDS(df065_filtered_transcript_results, "df065_transcript_results.rds")

# Creating separate data frames
df066_coding_gene_results <- df063_filtered_gene_protein_mirna_lncrna_results %>% filter(Molecular_class == "Coding gene")
saveRDS(df066_coding_gene_results, "coding_gene_results.rds")

df067_lncRNA_results <- df063_filtered_gene_protein_mirna_lncrna_results %>% filter(Molecular_class == "lncRNA")
saveRDS(df067_lncRNA_results, "lncRNA_results.rds")

df068_miRNA_results <- df063_filtered_gene_protein_mirna_lncrna_results %>% filter(Molecular_class == "miRNA")
saveRDS(df068_miRNA_results, "miRNA_results.rds")

df069_protein_results <- df063_filtered_gene_protein_mirna_lncrna_results %>% filter(Molecular_class == "Protein")
saveRDS(df069_protein_results, "protein_results.rds")

df070_transcript_coding_gene <- df065_filtered_transcript_results %>% filter(Molecular_class == "Coding Transcript Isoform")
saveRDS(df070_transcript_coding_gene, "transcript_coding_gene.rds")

df071_transcript_lncrna <- df065_filtered_transcript_results %>% filter(Molecular_class == "Long non-coding Transcript Isoforms")
saveRDS(df071_transcript_lncrna, "transcript_lncrna.rds")
