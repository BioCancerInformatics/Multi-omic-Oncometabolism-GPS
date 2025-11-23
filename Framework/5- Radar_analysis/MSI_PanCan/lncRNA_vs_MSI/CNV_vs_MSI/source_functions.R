#### UCSCXenaShiny Pancan genotype-phenotype correlations
## UCSC Xena Issues 
## Authors: Juan Carlo S Silva, Emanuell Rodrigues de Souza,Leonardo Henrique da Silva
## Ana Beatriz Garcia, Enrique Medina-Acostai
## Updated: 05/05/2024
## version:v10

## source_functions snippet
## CNV vs MSI

## Function 1: Define a wrapper function
vis_gene_wrapper <- function(sig) {
  
  p = vis_gene_msi_cor(Gene = sig, data_type = "cnv", cor_method = "spearman")
  p[["data"]] <- p[["data"]] %>%
    mutate(padj = p.adjust(p.value, method = "holm"))

  p[["data"]] <- p[["data"]] %>%
    filter(padj < 0.05)
  return(p)
}

## Function 2: Define data handler function
gene_analysis = function(genes){
  
  print(paste("Gene:", genes))
  p = lapply(genes, vis_gene_wrapper)
  return(p)
  
}

## Function 3: Reshaping
pivot_wider_manager = function(df) {

  df %>%
    filter(padj< 0.05) %>%
    select(cancer, cor) %>%
    pivot_wider(names_from = cancer,
                values_from = cor)

}

# Function 4: Run gene MSI
get_msi_data <- function(genes) {
  p <- lapply(genes, function(gene) {
    tryCatch({
      result <- gene_analysis(gene)
      return(result)
    }, error = function(e) {
      warning(paste("Error retrieving data for gene", gene, ":", e$message))
      return(NULL)
    })
  })
  names(p) <- genes
  p <- unlist(p, recursive = FALSE)
  return(p)
}

## Function 5: Get gene signatures
format_gene_signatures = function(nodes, network_table, correlation_status){
  
  if(correlation_status == 'all'){
    network_table = network_table
  }
  if(correlation_status == 'positive'){
    network_table = network_table[network_table$correlation > 0,]
  }
  if(correlation_status == 'negative'){
    network_table = network_table[network_table$correlation < 0,]
  }
  
  genes_by_cancer = aggregate(genes ~ cancer_types,
                              data = network_table,
                              FUN = function(x) paste(unique(x),
                                                      collapse = " + "))
  status_gene_sigs = paste0(correlation_status, '_gene_signatures')
  
  colnames(genes_by_cancer)[2] = status_gene_sigs
  nodes = left_join(nodes, genes_by_cancer, by=join_by(nodes == cancer_types))
  
  # Format gene signatures
  nodes[status_gene_sigs][is.na(nodes[status_gene_sigs])] = ''
  nodes[status_gene_sigs][nodes[status_gene_sigs] != ''] = paste0(
    "(", nodes[status_gene_sigs][nodes[status_gene_sigs] != ''], ")")
  
  # Find gene name containing '-'
  nodes[,ncol(nodes)] = gsub("([[:alpha:]]+-[[:alpha:]]+)", 
                               "'\\1'",
                             nodes[,ncol(nodes)])
  
  return(nodes)
  
}

## Function 6: Set nodes 
get_gene_cancer_node_table = function(network_table){
  
  ##### Set up nodes  ----
  gene_nodes = as.data.frame(network_table$genes)
  gene_nodes['node_class'] = 'gene'
  gene_nodes['node_color'] = '#CCCCCC'
  gene_nodes['node_bordercolor'] = '#404040'
  
  cancer_nodes = as.data.frame(network_table$cancer_types)
  cancer_nodes['node_class'] = 'cancer'
  cancer_nodes['node_color'] = '#ffcc66'
  cancer_nodes['node_bordercolor'] = '#404040'
  
  names(gene_nodes)[1] = 'nodes'
  names(cancer_nodes)[1] = 'nodes'
  
  ###### Get nodes by row
  nodes = rbind(gene_nodes, cancer_nodes)
  freq_nodes = table(nodes$nodes)
  nodes = arrange(nodes, nodes)
  nodes = nodes[!duplicated(nodes$nodes),]
  nodes$node_size = log2((freq_nodes)+1)*10
  nodes$interactions = freq_nodes
  
  ###### Set node labels
  nodes['labels>2'] = nodes$nodes
  nodes['labels>2'][nodes$node_size < 15,] = ''
  
  ###### Get mean correlation values and direction
  mean_corr = aggregate(network_table$correlation,
                        by = list(network_table$genes),
                        FUN = mean)
  
  colnames(mean_corr) = c('genes', 'mean_correlation')
  nodes = left_join(nodes,
                    mean_corr,
                    by=join_by(nodes == genes))
  
  nodes['corr_direction'] = ifelse(nodes$mean_correlation > 0,
                                   "positive",
                                   "negative")
  
  return(nodes)
  
}

