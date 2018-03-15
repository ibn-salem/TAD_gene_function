filter_paralogs <- function(cispairs, bio_mart) {
  
  #retrieves paralogs for all human protein coding genes from ENSEMBL
  paralog_atb <- c("ensembl_gene_id", "hsapiens_paralog_ensembl_gene")
  paralogs <- getBM(filters = "biotype", values = "protein_coding", attributes = paralog_atb, mart = bio_mart)
  
  #builds a tibble with all human paralog gene pairs
  paralogs <- as.tibble(paralogs)
  paralogs <- paralogs %>% 
    filter(hsapiens_paralog_ensembl_gene != "") %>% 
    mutate("paralog" = "paralog")
  names(paralogs) <- c("g1_id", "g2_id", "paralog")
  
  #classifies gene pairs whether or not they are paralogs and removes those that are
  cispairs <- cispairs %>% 
    left_join(paralogs, by = c("g1_id", "g2_id")) 

'%>% 
    filter(is.na(paralog)) %>% 
    select(g1, g2, dist, g1_id, g2_id, go_sim_BP, go_sim_MF, go_sim_CC)'
  
  return(cispairs)
} 