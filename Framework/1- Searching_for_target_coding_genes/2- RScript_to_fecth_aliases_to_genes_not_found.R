### RSTUDIO CODE TO FECTH ALIASES TO GENES NOT FOUND IN XENA BROWSER - PhD PROJECT
### Creator: HIGOR ALMEIDA CORDEIRO NOGUEIRA 
### Last Version - 10/03/2025

# Definir o diretório de trabalho
setwd("F:/Higor/Cancer_metabolism_analysis_07/1- Searching for coding genes/")

# Load required packages
library(httr)
library(jsonlite)
library(readr)
library(ggplot2)
library(KEGGREST)
library(rio)

# Read the file with the list of all metabolic pathway genes
gene_set <- import("Selected_Metabolic Pathways_Without_Duplicated_Genes.tsv")

# Create data frame with list of genes not found in Xena Browser
Genes_not_found <- import("Genes not found.txt")

### CODE SNIPPET TO GET ALIASES OF NOT FOUND GENES IN XENA BROWSER

# Define a function to search for KEGG gene symbols
search_kegg_symbol <- function(symbol) {
  result <- keggFind(database = "genes", query = symbol)
  if (length(result) > 0) {
    return(result)
  } else {
    return(NA)
  }
}

# Apply the function to each gene symbol in the dataframe
Genes_not_found$KEGG_Symbol <- sapply(Genes_not_found$Genes, search_kegg_symbol)

# export(Genes_not_found, "Aliases_gene_not_found.csv")

Aliases_genes_not_found <- import("Aliases_gene_not_found.csv")

Aliases_found <- import("Aliases_found.xlsx")

# Função para buscar alias em Aliases_found e filtrar os presentes
search_alias_filtered <- function(row, aliases_found) {
  # Extrair os aliases da linha e separar por ", "
  all_aliases <- unlist(strsplit(row["Aliases_found"], ", "))
  
  # Filtrar apenas os aliases que estão no Aliases_found$Aliases_found
  found_aliases <- all_aliases[all_aliases %in% aliases_found$Aliases_found]
  
  # Retornar os aliases encontrados ou NA
  if (length(found_aliases) > 0) {
    return(paste(found_aliases, collapse = ", "))
  } else {
    return(NA)
  }
}

# Aplicar a função para cada linha de Aliases_genes_not_found
Aliases_genes_not_found$Aliases_found <- apply(Aliases_genes_not_found, 1, search_alias_filtered, aliases_found = Aliases_found)

export(Aliases_genes_not_found, "Aliases_found_genes_not_found.tsv")

# Função para substituir genes no gene_set
replace_genes <- function(genes_column, replacements) {
  sapply(genes_column, function(gene) {
    replacement <- replacements$Aliases_found[replacements$Genes_not_found == gene]
    if (length(replacement) > 0 && !is.na(replacement)) {
      return(replacement)
    } else {
      return(gene)
    }
  })
}

# Substituir genes no gene_set
gene_set$Genes <- replace_genes(gene_set$Genes, Aliases_genes_not_found)

export(gene_set, "Uptaded_Target_genes.tsv")
