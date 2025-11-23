#### RSCRIPT TO ANALYSE RESULTS
#### HIGOR ALMEIDA, M.sC


# Carregar as bibliotecas necessárias
library(ggplot2)
library(dplyr)
library(rio)

setwd("F:/Higor/Cancer_metabolism_analysis_06/14- Results analysis/")

# Carregar a tabela
results <- import("/Higor/Cancer_metabolism_analysis_06/13- Drug Response/DGIdb 5.0 drug-gene interaction analysis/All_signature_pathways_interactions.rds")




# Criar data frame de resultados únicos
df072_unique_gene_protein_mirna_lncrna_results <- df063_filtered_gene_protein_mirna_lncrna_results %>%
  select(Target, Molecular_class, Cancer_types, Omic_layer, Phenotypic_layer)

search_your_target <- search_your_target %>%
  distinct()


df073_unique_gene_protein_mirna_lncrna <- df063_filtered_gene_protein_mirna_lncrna_results %>%
  select(Target, Molecular_class)

search_your_target <- search_your_target %>%
  distinct()

saveRDS(search_your_target, "search_your_target.rds")


# Filtering out rows where any of the selected columns contain "NS"
filtered_df <- All_results %>%
  filter(OS_Worst_Prognosis_Group != "NS" &
           DSS_Worst_Prognosis_Group != "NS" &
           DFI_Worst_Prognosis_Group != "NS" &
           PFI_Worst_Prognosis_Group != "NS")

# Criar data frame dos alvos disponíveis para busca no shiny
filtered_df <- filtered_df %>%
  select(Target, Molecular_Class)

filtered_df <- filtered_df %>%
  distinct()

# Filtering for desired classes
desired_classes <- c("Coding gene", "lncRNA", "miRNA")

filtered_df <- filtered_df %>% filter(Molecular_Class %in% desired_classes)

# Transform Target column
filtered_df$Target <- filtered_df$Target %>%
  sub("^hsa-", "", .) %>%       # Remove "hsa-" if present
  sub("miR", "MIR", .) %>%      # Replace "miR" with "MIR"
  gsub("-", "", .) %>%          # Remove all hyphens
  sub("(5p|3p)$", "", .)      # Remove only "-5p" or "-3p" at the end

export(filtered_df, "meaningful_targets_ survival.tsv")
























# Filtrando Linhas Onde Todas as Colunas de Interesse Possuem Valores Válidos
results_filtered <- results %>%
  filter(if_all(c(Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI, OS_Worst_Prognosis_Group, DSS_Worst_Prognosis_Group, 
                  DFI_Worst_Prognosis_Group, PFI_Worst_Prognosis_Group, 
                  Microenvironment_Classification, Immune_Classification, Annotation), 
                ~ . != "NS" & . != "No data" &. != "No Drug interaction information"))

results_filtered_2 <- results %>%
  filter(if_all(c(OS_Worst_Prognosis_Group, DSS_Worst_Prognosis_Group, 
                  DFI_Worst_Prognosis_Group, PFI_Worst_Prognosis_Group, 
                  Microenvironment_Classification, Immune_Classification, Annotation), 
                ~ . != "NS" & . != "No data" &. != "No Drug interaction information"))

results_filtered_3 <- results %>%
  filter(if_all(c(Type_Cox_OS, Type_Cox_DSS, Type_Cox_DFI, Type_Cox_PFI, OS_Worst_Prognosis_Group, DSS_Worst_Prognosis_Group, 
                  DFI_Worst_Prognosis_Group, PFI_Worst_Prognosis_Group), 
                ~ . != "NS" & . != "No data"),
         Microenvironment_Classification == "anti-tumoral",
         Immune_Classification == "Hot")

results_filtered_4 <- results %>%
  filter(if_all(c(OS_Worst_Prognosis_Group, DSS_Worst_Prognosis_Group, 
                  DFI_Worst_Prognosis_Group, PFI_Worst_Prognosis_Group, 
                  Microenvironment_Classification, Immune_Classification), 
                ~ . != "NS" & . != "No data"))

# Define the directories
directories <- c("/Higor/Cancer_metabolism_analysis_06/15- Cancer GPS Pathways Shiny/Exploring pan-cancer metabolism/", 
                 "/Higor/Cancer_metabolism_analysis_06/14- Results analysis/")

# Save the RDS file in multiple directories
lapply(directories, function(dir) {
  saveRDS(results_filtered_2, file.path(dir, "the_best_signature.rds"))
})

# Create a new folder named "New_Folder"
dir.create("Metabolism_signatures")

setwd("F:/Higor/Cancer_metabolism_analysis_06/14- Results analysis/Metabolism_signatures/")

# Filtrar assinaturas específicas para cada camada ômica e salvar
# ----------------------------------------------------------------
Lipid_signature <- immune_prognostic_results_filtered %>% filter(Metabolism == "Lipid metabolism")
saveRDS(Lipid_signature, "the_best_Lipid_signature.rds")

Nucleotide_signature <- immune_prognostic_results_filtered %>% filter(Metabolism == "Nucleotide metabolism")
saveRDS(Nucleotide_signature, "the_best_Nucleotide_signature.rds")

Energy_signature <- immune_prognostic_results_filtered %>% filter(Metabolism == "Energy metabolism")
saveRDS(Energy_signature, "the_best_Energy_signature.rds")

Amino_acid_signature <- immune_prognostic_results_filtered %>% filter(Metabolism == "Amino acid metabolism")
saveRDS(Amino_acid_signature, "the_best_Amino_acid_signature.rds")

Carbohydrate_mutation_signature <- immune_prognostic_results_filtered %>% filter(Metabolism == "Carbohydrate metabolism")
saveRDS(Carbohydrate_mutation_signature, "the_best_Carbohydrate_mutation_signature.rds")

cofactors_and_vitamins_signature <- immune_prognostic_results_filtered %>% filter(Metabolism == "Metabolism of cofactors and vitamins")
saveRDS(cofactors_and_vitamins_signature, "the_best_cofactors_and_vitamins_signature.rds")

other_amino_acids_signature <- immune_prognostic_results_filtered %>% filter(Metabolism == "Metabolism of other amino acids")
saveRDS(other_amino_acids_signature, "the_best_other_amino_acids_signature.rds")

setwd("F:/Higor/Cancer_metabolism_analysis_06/13- Results analysis/Signatures/")

# Lista das colunas a serem mantidas
colunas_para_manter <- c(
  "Members", "Multiomics_Signatures", "Molecular_Class", "Cancer_Types", 
  "Omic_Layer", "Phenotypic_Layer", "correlation_type", 
  "Type_Cox_OS", "Type_Cox_DSS", "Type_Cox_DFI", "Type_Cox_PFI", 
  "OS_Worst_Prognosis_Group", "DSS_Worst_Prognosis_Group", 
  "DFI_Worst_Prognosis_Group", "PFI_Worst_Prognosis_Group", 
  "Microenvironment_Classification", "Immune_Classification"
)

# Mantendo apenas as colunas selecionadas
results_filtered <- results[, colunas_para_manter]

# Remover linhas duplicadas completamente, mas respeitar variações nas colunas
results_filtered <- results_filtered %>%
  distinct()

# Filtrar assinaturas específicas para cada camada ômica e salvar
# ----------------------------------------------------------------
# Filtrando Linhas Onde Todas as Colunas de Interesse Possuem Valores Válidos
immune_prognostic_results_filtered <- results_filtered %>%
  filter(if_all(c(OS_Worst_Prognosis_Group, DSS_Worst_Prognosis_Group, 
                  DFI_Worst_Prognosis_Group, PFI_Worst_Prognosis_Group, 
                  Microenvironment_Classification, Immune_Classification), 
                ~ . != "NS" & . != "No data"))

# Create a new folder named "New_Folder"
dir.create("Omic_signatures")

setwd("F:/Higor/Cancer_metabolism_analysis_06/13- Results analysis/Signatures/Omic_signatures/")

# Filtrar assinaturas específicas para cada camada ômica e salvar
# ----------------------------------------------------------------
omic_gene_signature <- immune_prognostic_results_filtered %>% filter(Omic_Layer == "Gene expression")
saveRDS(omic_gene_signature, "the_best_coding_gene_signature.rds")

omic_methylation_signature <- immune_prognostic_results_filtered %>% filter(Omic_Layer == "Methylation")
saveRDS(omic_methylation_signature, "the_best_methylation_signature.rds")

omic_miRNA_signature <- immune_prognostic_results_filtered %>% filter(Omic_Layer == "miRNA expression")
saveRDS(omic_miRNA_signature, "the_best_miRNA_signature.rds")

omic_protein_signature <- immune_prognostic_results_filtered %>% filter(Omic_Layer == "Protein expression")
saveRDS(omic_protein_signature, "the_best_protein_signature.rds")

omic_mutation_signature <- immune_prognostic_results_filtered %>% filter(Omic_Layer == "Mutation")
saveRDS(omic_mutation_signature, "the_best_mutation_signature.rds")

omic_cnv_signature <- immune_prognostic_results_filtered %>% filter(Omic_Layer == "CNV")
saveRDS(omic_cnv_signature, "the_best_cnv_signature.rds")

omic_transcript_signature <- immune_prognostic_results_filtered %>% filter(Omic_Layer == "Transcript expression")
saveRDS(omic_transcript_signature, "the_best_transcript_signature.rds")



























# Análise da quantidade de observações únicas
Multiomics_Signatures <- results_filtered %>%
  group_by(Multiomics_Signatures) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

Multiomics_Signatures_cancer <- results_filtered %>%
  group_by(Cancer_Types, Multiomics_Signatures) %>%
  summarise(quantidade = n(), .groups = "drop") %>%
  arrange(Cancer_Types, desc(quantidade))

# Análise da quantidade de observações únicas
Cancer_Types <- results_filtered_2 %>%
  group_by(Cancer_Types) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

# Análise da quantidade de observações únicas
Molecular_Class <- results_filtered %>%
  group_by(Molecular_Class) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

# Análise da quantidade de observações únicas
Tumor_vs_Normal <- results_filtered %>%
  group_by(Tumor_vs_Normal) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

# Análise da quantidade de observações únicas
Omic_Layer <- results_filtered %>%
  group_by(Omic_Layer) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

Omic_Cancer <- results_filtered %>%
  group_by(Cancer_Types, Omic_Layer) %>%
  summarise(quantidade = n(), .groups = "drop") %>%
  arrange(Cancer_Types, desc(quantidade))

# Análise da quantidade de observações únicas
Phenotypic_Layer <- results_filtered %>%
  group_by(Phenotypic_Layer) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

Phenotypic_Cancer <- results_filtered %>%
  group_by(Cancer_Types, Phenotypic_Layer) %>%
  summarise(quantidade = n(), .groups = "drop") %>%
  arrange(Cancer_Types, desc(quantidade))

# Análise da quantidade de observações únicas
Microenvironment_Classification <- results_filtered %>%
  group_by(Microenvironment_Classification) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

Microenvironment_Cancer <- results_filtered %>%
  group_by(Cancer_Types, Microenvironment_Classification) %>%
  summarise(quantidade = n(), .groups = "drop") %>%
  arrange(Cancer_Types, desc(quantidade))

Microenvironment_Cancer_Target <- results_filtered %>%
  group_by(Cancer_Types, Microenvironment_Classification, Target) %>%
  summarise(quantidade = n(), .groups = "drop") %>%
  arrange(Cancer_Types, Microenvironment_Classification, desc(quantidade))

# Análise da quantidade de observações únicas
Immune_Classification <- results_filtered %>%
  group_by(Immune_Classification) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

Immune_Cancer <- results_filtered %>%
  group_by(Cancer_Types, Microenvironment_Classification) %>%
  summarise(quantidade = n(), .groups = "drop") %>%
  arrange(Cancer_Types, desc(quantidade))

Immune_Cancer_Target <- results_filtered %>%
  group_by(Cancer_Types, Microenvironment_Classification, Target) %>%
  summarise(quantidade = n(), .groups = "drop") %>%
  arrange(Cancer_Types, Microenvironment_Classification, desc(quantidade))

# Filtrando Linhas com Pelo Menos uma Observação Válida nas Colunas de Prognóstico
prognostic_results_filtered_df1 <- results_filtered %>%
  filter(if_any(c(OS_Worst_Prognosis_Group, DSS_Worst_Prognosis_Group, 
                  DFI_Worst_Prognosis_Group, PFI_Worst_Prognosis_Group), 
                ~ . != "NS" & . != "No data"))


# Filtrando Linhas Onde Todas as Colunas de Prognóstico Possuem Valores Válidos
prognostic_results_filtered_df2 <- results_filtered %>%
  filter(if_all(c(OS_Worst_Prognosis_Group, DSS_Worst_Prognosis_Group, 
                  DFI_Worst_Prognosis_Group, PFI_Worst_Prognosis_Group), 
                ~ . != "NS" & . != "No data"))

















# Análise da quantidade de observações únicas
Multiomics_Signatures <- immune_prognostic_results_filtered %>%
  group_by(Multiomics_Signatures) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))

Multiomics_Signatures_cancer <- immune_prognostic_results_filtered %>%
  group_by(Cancer_Types, Multiomics_Signatures) %>%
  summarise(quantidade = n(), .groups = "drop") %>%
  arrange(Cancer_Types, desc(quantidade))


# Análise da quantidade de observações únicas
Cancer_Types <- immune_prognostic_results_filtered %>%
  group_by(Cancer_Types) %>%
  summarise(quantidade = n()) %>%
  arrange(desc(quantidade))
