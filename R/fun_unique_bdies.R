find_genes_surround_unique_bdies <- function(grl_bdies, gr_HumProtGenes, meta_BDY, window_BDY, hum_seqinfo){

  #Creates Tibble with all Boundaries, their frequency and additional metadata provided
  bdies <- as_tibble(as.data.frame(grl_bdies)) 
  bdies_frequency <- dplyr::count(bdies, start, seqnames) 
  bdies <- bdies %>% 
    left_join(bdies_frequency, by = c("start", "seqnames")) %>% 
    left_join(meta_BDY, by = "group_name") %>% 
    select(Origin, Type, group_name, Name, seqnames, start, end, n)

  
  #Creates GRangesobject out of all regions that sourround the boundaries in a given 
  #window. Thereby excludes the boundary coordinates and trims al Ranges that 
  #are out of bounds.
  left <- GRanges(seqnames = bdies$seqnames,
                          IRanges(start = bdies$start - window_BDY,
                                  end = bdies$start),
                          seqinfo = hum_seqinfo,
                          "group_name" = bdies$group_name,
                          "start_of_realbdy" = bdies$start)
  left <- trim(left)
  right <- GRanges(seqnames = bdies$seqnames,
                           IRanges(start = bdies$end,
                                   end = bdies$end + window_BDY),
                           seqinfo = hum_seqinfo,
                           "group_name" = bdies$group_name,
                           "start_of_realbdy" = bdies$start)
  right<- trim(right)
  all_surrounding_ranges <- c(left, right)
  
  #Finds all Genes overlapping the surrounding ranges and retreives their ENSEMBL Gene IDs.
  #Further retrieves frequency and metadata information from the "bdies" object by using
  #the original boundary and the original celltype of the queryHits
  hits_in_window <- as.tibble(as.data.frame(findOverlaps(all_surrounding_ranges, gr_HumProtGenes)))
  genes_in_window <- gr_HumProtGenes[hits_in_window$subjectHits] %>% 
    as.data.frame() %>% 
    as.tibble() %>% 
    select(seqnames, gene_id) %>% 
    mutate(start = all_surrounding_ranges$start_of_realbdy[hits_in_window$queryHits]) %>% 
    mutate(group_name = all_surrounding_ranges$group_name[hits_in_window$queryHits]) %>% 
    left_join(bdies, by = c("seqnames", "start", "group_name"))

  #Saves a list of all distinct Gene IDs in the window to the tidydata folder to use as 
  #background for the GO enrichment analysis
  all_geneID_in_window <- genes_in_window %>% 
    select(gene_id) %>% 
    distinct() %>% 
    as.matrix()
  myfile <- file.path("results/tidydata", paste0("AllGeneIDsInWindow_", window_BDY))
  write(all_geneID_in_window, myfile, sep = "<ENTER>")
  
  #Saves a list of genes in a window around a boundary that is uniqe to a celltype. Does
  #this for all celltypes to use it for the GO-Enrichment analysis
  celltypes <- genes_in_window %>% 
    select(group_name) %>% 
    distinct() %>% 
    as.data.frame()
  celltypes <- as.list(celltypes$group_name)
    
  for(celltype in celltypes) {
    myfile <- file.path("results/tidydata", paste0("GeneIDsInUniqueWindow_", celltype, "_", window_BDY))
    genes_in_window %>% 
      filter(group_name == celltype) %>% 
      filter(n == 1) %>% 
      select(gene_id) %>% 
      distinct() %>% 
      as.matrix() %>% 
      write(myfile, sep = "<ENTER>")
  }
  
  return(bdies)
}








