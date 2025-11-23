# Carregar pacotes necessários
library(biomaRt)
library(httr)
library(jsonlite)
library(rio)
setwd("F:/Higor/Cancer_metabolism_analysis_06/4- Searching for target transcript isoforms/lncRNA/")

# Função para obter informações do gene e adicionar ao dataframe
get_gene_info <- function(gene_symbol) {
  base_url <- "https://rest.ensembl.org/"
  endpoint <- paste0("lookup/symbol/human/", URLencode(gene_symbol, reserved = TRUE), "?expand=1")
  
  tryCatch({
    response <- GET(paste0(base_url, endpoint), add_headers(Accept = "application/json"))
    if (http_status(response)$category == "Success") {
      return(fromJSON(content(response, "text"), flatten = TRUE))
    } else {
      cat("Erro ao buscar informações do gene", gene_symbol, "\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("Erro de conexão, tentando novamente para o gene", gene_symbol, "\n")
    Sys.sleep(5)  # Pausa de 5 segundos antes de tentar novamente
    return(get_gene_info(gene_symbol))
  })
}

# Definir caminhos de arquivos
csv_file <- "gene_info.csv"
rds_file <- "progresso_analise.rds"

# Carregar o progresso anterior ou iniciar novo processo
if (file.exists(rds_file)) {
  progresso_analise <- readRDS(rds_file)
} else {
  genes <- import("/Higor/Cancer_metabolism_analysis_06/3- Searching for target lncRNA/Target_lncrna.tsv")
  progresso_analise <- list(Genes = genes$Target_lncrna, combined_gene_df = data.frame())
  cat("Arquivo de progresso não encontrado. Iniciando do zero.\n")
}

# Processamento dos genes
analyzed_genes <- unique(progresso_analise$combined_gene_df$Gene_ID)
genes_to_analyze <- setdiff(progresso_analise$Genes, analyzed_genes)

if (length(genes_to_analyze) > 0) {
  for (gene_symbol in genes_to_analyze) {
    cat("Analisando o gene:", gene_symbol, "\n")
    gene_info <- get_gene_info(gene_symbol)
    
    if (!is.null(gene_info)) {
      gene_df <- data.frame(
        Display_Name = gene_info$display_name,
        Transcript_ID = gene_info$Transcript$id,
        Gene_ID = gene_info$id,
        Biotype = gene_info$biotype,
        Start = gene_info$start,
        End = gene_info$end,
        Chromosome = gene_info$seq_region_name,
        Description = gene_info$description,
        Canonical_Transcript = gene_info$canonical_transcript,
        stringsAsFactors = FALSE
      )
      
      progresso_analise$combined_gene_df <- rbind(progresso_analise$combined_gene_df, gene_df)
      
      # Atualizar o arquivo CSV em lotes
      if (nrow(progresso_analise$combined_gene_df) %% 10 == 0) {
        write.table(progresso_analise$combined_gene_df, csv_file, append = TRUE, sep = ",", col.names = !file.exists(csv_file), row.names = FALSE)
        progresso_analise$combined_gene_df <- data.frame()  # Resetar o dataframe após salvar
      }
    }
    
    saveRDS(progresso_analise, rds_file)
  }
  # Certifique-se de gravar quaisquer dados restantes no DataFrame após o loop
  if (nrow(progresso_analise$combined_gene_df) > 0) {
    write.table(progresso_analise$combined_gene_df, csv_file, append = TRUE, sep = ",", col.names = !file.exists(csv_file), row.names = FALSE)
  }
} else {
  cat("Nenhum gene novo para analisar.\n")
}

# Exibir o dataframe final (se necessário)
print(progresso_analise$combined_gene_df)

info <- import("gene_info.csv")
