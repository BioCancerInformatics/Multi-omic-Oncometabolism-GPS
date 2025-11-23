### RSTUDIO CODE TO FECTH GENE AND PATHWAYS FROM KEGG - PhD PROJECT
### HIGOR ALMEIDA, M.sC 
### LAST VERSION - 11/06/2024

# Definir o diretório de trabalho
setwd("F:/Higor/Metabolic_gene_analysis_05/1- Searching_for_metabolic_genes")

### ETAPA 01: INSTALAR E CARREGAR AS BIBLIOTECAS NECESSÁRIAS

# Vetor com os nomes dos pacotes necessários ao longo do código
library_packages <- c("BiocManager", "KEGGREST", "tidyverse", "dplyr", "rio")   

# Carregar as bibliotecas necessárias
lapply(library_packages, function(pkg) {                                        
  if (!requireNamespace(pkg, quietly = TRUE)) {                                 
    warning(paste0("Pacote ", pkg, " não encontrado."))                         
  } else {
    library(pkg, character.only = TRUE)                                         
  }
})

### ETAPA 02: BUSCAR VIAS METABÓLICAS HUMANAS DO KEGG

# Obter todas as vias metabólicas humanas
all_pathways <- keggList("pathway", "hsa")                                      # Esse código usa a função keggList do pacote KEGGREST para obter uma lista de vias metabólicas humanas do KEGG (Kyoto Encyclopedia of Genes and Genomes), que é uma base de dados de informações sobre vias metabólicas, genomas e outras informações biológicas.O operador <- atribui o resultado da função keggList("pathway", "hsa") à variável all_pathways. 

# Converter para um data frame para facilitar a manipulação
ID_all_pathways_df <- tibble(KEGG_ID = names(all_pathways), Pathway_Name = as.character(all_pathways)) # names(all_pathways): Extrai os nomes (identificadores) das vias. as.character(all_pathways): Converte as descrições das vias para um vetor de caracteres.

# Função para processar cada KEGG ID
process_kegg_id <- function(kegg_id) {                                          # Define uma função chamada process_kegg_id que aceita um argumento kegg_id.
  pathway_info <- try(keggGet(kegg_id), silent = TRUE)                          # keggGet(kegg_id): Chama a função keggGet para obter informações da via metabólica associada ao kegg_id.try(..., silent = TRUE): Usa try para capturar erros silenciosamente. Se keggGet falhar, try retornará um objeto de erro em vez de interromper a execução.
  if (inherits(pathway_info, "try-error") || is.null(pathway_info[[1]]$GENE)) return(NULL) # inherits(pathway_info, "try-error"): Verifica se pathway_info é um objeto de erro (ou seja, se keggGet falhou). is.null(pathway_info[[1]]$GENE): Verifica se a lista de genes ($GENE) em pathway_info é nula. return(NULL): Se qualquer uma das condições acima for verdadeira, a função retorna NULL e termina a execução.
  
  gene_info <- pathway_info[[1]]$GENE
  if (is.null(gene_info)) return(NULL)
  
  gene_info_df <- tibble(gene_info = gene_info) %>%                             # tibble(gene_info = gene_info): Cria um tibble (um tipo de data frame) com uma coluna gene_info contendo as informações dos genes.
    separate(gene_info, into = c("GeneID", "Description"), sep = ";\\s*", remove = TRUE, convert = TRUE) %>%    # separate(...): Separa a coluna gene_info em duas novas colunas, GeneID e Description, usando o ponto-e-vírgula seguido por espaços (";\\s*") como separador. into = c("GeneID", "Description"): Define os nomes das novas colunas. remove = TRUE: Remove a coluna original gene_info. convert = TRUE: Converte as novas colunas para tipos de dados apropriados.
    mutate(GeneSymbol = word(Description, 1),                                   # mutate(...): Cria novas colunas ou modifica as existentes. GeneSymbol = word(Description, 1): Extrai a primeira palavra de Description como GeneSymbol.
           GeneName = str_remove_all(Description, "\\[.*?\\]"),                 # GeneName = str_remove_all(Description, "\\[.*?\\]"): Remove qualquer texto entre colchetes em Description para criar GeneName.
           KEGG_ID = as.character(kegg_id),                                     # KEGG_ID = as.character(kegg_id): Adiciona a coluna KEGG_ID com o valor do kegg_id.
           Pathway_Name = as.character(pathway_info[[1]]$NAME)) %>%             # Pathway_Name = as.character(pathway_info[[1]]$NAME): Adiciona a coluna Pathway_Name com o nome da via metabólica.
    select(-Description)                                                        # select(-Description): Remove a coluna Description do tibble. 
  
  return(gene_info_df)                                                          # Retorna o tibble gene_info_df contendo as informações formatadas dos genes.
}

# Extrair os KEGG IDs para processamento
kegg_ids_to_process <- ID_all_pathways_df$KEGG_ID                               # Este código extrai os IDs KEGG de todas as vias metabólicas presentes no data frame ID_all_pathways_df e os armazena na variável kegg_ids_to_process.

# Processar todos os KEGG IDs e combinar os dados
combined_data <- bind_rows(lapply(kegg_ids_to_process, process_kegg_id)) %>%    # kegg_ids_to_process: Um vetor contendo os IDs KEGG que você deseja processar. lapply: Aplica a função process_kegg_id a cada elemento de kegg_ids_to_process. O resultado é uma lista onde cada elemento é o data frame retornado pela função process_kegg_id.Aqui, bind_rows pega a lista de data frames gerada por lapply e os combina em um único data frame combined_data.
  distinct() %>%                                                                # O operador pipe (%>%) é usado para encadear operações em uma sequência clara e legível. distinct é uma função do pacote dplyr que remove linhas duplicadas de um data frame.
  drop_na()                                                                     # drop_na é uma função do pacote tidyr (parte do tidyverse) que remove todas as linhas contendo valores NA (nulos).

# Lista dos códigos KEGG dos pathways desejados
pathways_desejados <- c("hsa00010", "hsa00020", "hsa00030", "hsa00040", "hsa00051", "hsa00052", "hsa00053", "hsa00500", "hsa00520", "hsa00620", "hsa00630", "hsa00640", 
                        "hsa00650", "hsa00660", "hsa00562", "hsa00190", "hsa00195", "hsa00196", "hsa00710", "hsa00720", "hsa00680", "hsa00910", "hsa00920", "hsa00061", 
                        "hsa00062", "hsa00071", "hsa00073", "hsa00100", "hsa00120", "hsa00121", "hsa00140", "hsa00561", "hsa00564", "hsa00565", "hsa00600", "hsa00590", 
                        "hsa00591", "hsa00592", "hsa01040", "hsa00230", "hsa00240", "hsa00250", "hsa00260", "hsa00270", "hsa00280", "hsa00290", "hsa00300", "hsa00310", 
                        "hsa00220", "hsa00330", "hsa00340", "hsa00350", "hsa00360", "hsa00380", "hsa00400", "hsa00410", "hsa00430", "hsa00440", "hsa00450", "hsa00460", 
                        "hsa00470", "hsa00480", "hsa00730", "hsa00740", "hsa00750", "hsa00760", "hsa00770", "hsa00780", "hsa00785", "hsa00790", "hsa00670", "hsa00830", 
                        "hsa00860", "hsa00130")

# Filtrar apenas os pathways desejados
combined_data_selecionado <- combined_data %>% filter(KEGG_ID %in% pathways_desejados) # combined_data: O data frame inicial com todos os dados combinados. %>%: O operador pipe que passa combined_data para a função filter(). filter(KEGG_ID %in% pathways_desejados): Filtra as linhas de combined_data para incluir apenas aquelas onde KEGG_ID está presente no vetor pathways_desejados. %in%: É um operador que verifica se os elementos do lado esquerdo estão presentes no vetor do lado direito. 

# Definir a função para atribuir o metabolismo com base no código KEGG
atribuir_metabolismo <- function(kegg_id) {                                     # atribuir_metabolismo <- function(kegg_id) define a função atribuir_metabolismo que recebe um argumento kegg_id.
  metabolismo <- switch(kegg_id,                                                # switch é uma função condicional que avalia o valor de kegg_id e retorna um resultado com base em correspondências específicas.Kegg_id: O valor a ser avaliado pela função switch. Se kegg_id for hsa00010, switch retorna "Carbohydrate metabolism". 
                        hsa00010 = "Carbohydrate metabolism",
                        hsa00020 = "Carbohydrate metabolism",
                        hsa00030 = "Carbohydrate metabolism",
                        hsa00040 = "Carbohydrate metabolism",
                        hsa00051 = "Carbohydrate metabolism",
                        hsa00052 = "Carbohydrate metabolism",
                        hsa00053 = "Carbohydrate metabolism",
                        hsa00500 = "Carbohydrate metabolism",
                        hsa00520 = "Carbohydrate metabolism",
                        hsa00620 = "Carbohydrate metabolism",
                        hsa00630 = "Carbohydrate metabolism",
                        hsa00640 = "Carbohydrate metabolism",
                        hsa00650 = "Carbohydrate metabolism",
                        hsa00660 = "Carbohydrate metabolism",
                        hsa00562 = "Carbohydrate metabolism",
                        hsa00190 = "Energy metabolism",
                        hsa00195 = "Energy metabolism",
                        hsa00196 = "Energy metabolism",
                        hsa00710 = "Energy metabolism",
                        hsa00720 = "Energy metabolism",
                        hsa00680 = "Energy metabolism",
                        hsa00910 = "Energy metabolism",
                        hsa00920 = "Energy metabolism",
                        hsa00061 = "Lipid metabolism",
                        hsa00062 = "Lipid metabolism",
                        hsa00071 = "Lipid metabolism",
                        hsa00073 = "Lipid metabolism",
                        hsa00100 = "Lipid metabolism",
                        hsa00120 = "Lipid metabolism",
                        hsa00121 = "Lipid metabolism",
                        hsa00140 = "Lipid metabolism",
                        hsa00561 = "Lipid metabolism",
                        hsa00564 = "Lipid metabolism",
                        hsa00565 = "Lipid metabolism",
                        hsa00600 = "Lipid metabolism",
                        hsa00590 = "Lipid metabolism",
                        hsa00591 = "Lipid metabolism",
                        hsa00592 = "Lipid metabolism",
                        hsa01040 = "Lipid metabolism",
                        hsa00230 = "Nucleotide metabolism",
                        hsa00240 = "Nucleotide metabolism",
                        hsa00250 = "Amino acid metabolism",
                        hsa00260 = "Amino acid metabolism",
                        hsa00270 = "Amino acid metabolism",
                        hsa00280 = "Amino acid metabolism",
                        hsa00290 = "Amino acid metabolism",
                        hsa00300 = "Amino acid metabolism",
                        hsa00310 = "Amino acid metabolism",
                        hsa00220 = "Amino acid metabolism",
                        hsa00330 = "Amino acid metabolism",
                        hsa00340 = "Amino acid metabolism",
                        hsa00350 = "Amino acid metabolism",
                        hsa00360 = "Amino acid metabolism",
                        hsa00380 = "Amino acid metabolism",
                        hsa00400 = "Amino acid metabolism",
                        hsa00410 = "Metabolism of other amino acids",
                        hsa00430 = "Metabolism of other amino acids",
                        hsa00440 = "Metabolism of other amino acids",
                        hsa00450 = "Metabolism of other amino acids",
                        hsa00460 = "Metabolism of other amino acids",
                        hsa00470 = "Metabolism of other amino acids",
                        hsa00480 = "Metabolism of other amino acids",
                        hsa00730 = "Metabolism of cofactors and vitamins",
                        hsa00740 = "Metabolism of cofactors and vitamins",
                        hsa00750 = "Metabolism of cofactors and vitamins",
                        hsa00760 = "Metabolism of cofactors and vitamins",
                        hsa00770 = "Metabolism of cofactors and vitamins",
                        hsa00780 = "Metabolism of cofactors and vitamins",
                        hsa00785 = "Metabolism of cofactors and vitamins",
                        hsa00790 = "Metabolism of cofactors and vitamins",
                        hsa00670 = "Metabolism of cofactors and vitamins",
                        hsa00830 = "Metabolism of cofactors and vitamins",
                        hsa00860 = "Metabolism of cofactors and vitamins",
                        hsa00130 = "Metabolism of cofactors and vitamins",
                        NA)
  return(metabolismo)
}

# Adicionar a coluna "Metabolism" ao dataframe combinado
combined_data_selecionado <- combined_data_selecionado %>%                      # o operador pipe passa o data frame combined_data_selecionado como entrada para a função mutate().
  mutate(Metabolism = sapply(KEGG_ID, atribuir_metabolismo))                    # mutate é uma função do pacote dplyr usada para criar novas colunas ou modificar colunas existentes em um data frame. sapply: É uma função de aplicação do pacote base do R que aplica uma função a cada elemento de um vetor ou lista e simplifica o resultado. KEGG_ID: Refere-se à coluna KEGG_ID do data frame combined_data_selecionado. atribuir_metabolismo: É a função definida anteriormente que atribui uma categoria de metabolismo com base no código KEGG. sapply(KEGG_ID, atribuir_metabolismo): Aplica a função atribuir_metabolismo a cada elemento da coluna KEGG_ID e retorna um vetor com os resultados.

# Reordenar e renomear algumas colunas
combined_data_selecionado <- combined_data_selecionado %>%
  select(KEGG_ID, Metabolism, Pathway_Name, GeneID, GeneSymbol, GeneName) %>%
  rename(Pathways = Pathway_Name, Genes = GeneID, Genes_name = GeneSymbol, Isoforms = GeneName)

# Exportar os dados para um arquivo
export(combined_data_selecionado, "Selected_Metabolic_Pathways_KEGG.tsv")

# Remover duplicatas na coluna "Genes" e combinar as observações das outras colunas relacionadas
Target_genes <- combined_data_selecionado %>%
  group_by(Genes) %>%
  summarize(
    KEGG_ID = paste(unique(KEGG_ID), collapse = " / "),
    Metabolism = paste(unique(Metabolism), collapse = " / "),
    Pathways = paste(unique(Pathways), collapse = " / "),
    Genes_name = paste(unique(Genes_name), collapse = " / "),
    Isoforms = paste(unique(Isoforms), collapse = " / ")
  ) %>%
  ungroup()

# Exportar os dados para um arquivo
export(Target_genes, "Selected_Metabolic Pathways_Without_Duplicated_Genes.tsv", row.names = FALSE)
