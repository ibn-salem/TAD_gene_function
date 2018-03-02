compare_kegg_pathways <- function(cispairs, bio_mart){
  
  #retrieves KEGG pathway annotations for all human protein coding genes from ENSEMBL
  kegg_atb <- c("ensembl_gene_id", "kegg_enzyme")
  kegg_pws <- getBM(filters = "biotype", values = "protein_coding", attributes = kegg_atb, mart = bio_mart)
  kegg_pws <- as_tibble(kegg_pws) 
  kegg_pws <- kegg_pws %>% 
    filter(kegg_enzyme != "")
  
  #selects only the two gene ids from the cispair table  
  merge_kegg_cisp <- cispairs %>% 
    select("g1_id", "g2_id")
  
  #adds information about the KEGG annotations correlating with the g1 ids  
  names(kegg_pws) <- c("g1_id", "kegg_enzyme")
  merge_kegg_cisp <- left_join(merge_kegg_cisp, kegg_pws, by = c("g1_id"))
    
  #adds information about the KEGG annotations correlating with the g2 ids, removes all
  #rows that do no carry KEGG annotations in g1 and g2 columns, compares the two KEGG
  #annotations and filters for the rows where they are the same.
  names(kegg_pws) <- c("g2_id", "kegg_enzyme")
  merge_kegg_cisp <- merge_kegg_cisp %>% 
    left_join(kegg_pws, by = c("g2_id")) %>% 
    filter(!is.na(kegg_enzyme.x)) %>% 
    filter(!is.na(kegg_enzyme.y)) %>% 
    mutate(same_kegg_pathway = ifelse(kegg_enzyme.x == kegg_enzyme.y, "same", "not_same")) %>% 
    filter(same_kegg_pathway == "same") %>% 
    select("g1_id", "g2_id", "same_kegg_pathway")
  
  #adds a colum with the information whether or not a gene pair is connected via
  #the same KEGG annotations 
  cispairs <- cispairs %>% 
    left_join(merge_kegg_cisp, by = c("g1_id", "g2_id")) %>% 
    mutate(same_kegg_pathway = ifelse(is.na(same_kegg_pathway), "not same", "same"))
  
  return(cispairs)
}