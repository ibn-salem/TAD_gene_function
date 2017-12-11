getAllRdmBoundaries <- function (grangeslist, size, seqinfo){

  grangeslist <- read_rds("Datasets/data_grangeslist_all_BDY.rds")
  
  size <- 40000
  seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
  

  
  
  rdm_bdy <- function(bdy_per_chr, bdysize, hum_seqinfo){
    
    bdy_per_chr <- bdy_counts %>% filter(cell_type == cell_type)
    
    bdy_per_chr <- bdy_per_chr %>%
      mutate(chrom_len = seqlengths(hum_seqinfo[seqnames])- bdysize)
    
    # sample the coordinates of random breakpoints

    starts <- map2(bdy_per_chr$chrom_len, bdy_per_chr$n, sample.int)
    seqnames <- rep(bdy_per_chr$seqnames, bdy_per_chr$n)
    
    rdm_bdy <- GRanges(seqnames,
                       IRanges(unlist(starts), width = size),
                       seqinfo = hum_seqinfo)
    return(rdm_bdy)
  }
  
  fun <- function(cell_type){
    bdy_counts_cell_type <- bdy_counts %>% filter(cell_type == cell_type)
    rdm_bdy <- rdm_bdy(bdy_counts_cell_type, size, seqinfo)
    mcols(rdm_bdy)$cell_type <- cell_type 
    return(rdm_bdy)
  }
  

  list_rdm_bdy <- c()
  
  for (i in 1:length(grangeslist)){
    #counts number of boundaries per celltype and chromosome
    bdy_counts <- grangeslist[i] %>% 
      as_tibble() %>% 
      rename_( "cell_type" = "group_name") %>% 
      group_by(cell_type, seqnames) %>% 
      summarise(n = n()) %>%
      mutate(seqnames = as.character(seqnames))
    
    
    #loop over celltypes and apply the function to randomly introduce boundaries
    #returns a list of GRangesobjects per celltype carrying the new boundaries
    cell_type <- unique(unlist(bdy_counts$cell_type))
    list_rdm_bdy <- c(list_rdm_bdy, map(cell_type, fun))
  }
  
  list_rdm_bdy <- GRangesList(list_rdm_bdy)
  names(list_rdm_bdy) <- names(grangeslist)
  
  return(list_rdm_bdy)

}
