# sample random boundaries for cell type
sample_rdm_bdies <- function(real_bdies_gr, bdysize, hum_seqinfo){
  
  #counts number of boundaries per chromosome
  bdy_per_chr <- real_bdies_gr %>% 
    as_tibble() %>% 
    group_by(seqnames) %>% 
    summarise(n = n()) %>%
    mutate(seqnames = as.character(seqnames),
           #shrink chromsize by one bdy length, to restrict sampling out of bounds 
           chrom_len = seqlengths(hum_seqinfo[seqnames])- bdysize)
  
  #sample the coordinates of random boundaries
  starts <- map2(bdy_per_chr$chrom_len, bdy_per_chr$n, sample.int)
  
  #provides chromosome names for random boundaries
  seqnames <- rep(bdy_per_chr$seqnames, bdy_per_chr$n)
  
  rdm_bdies_gr <- GRanges(seqnames,
                     IRanges(unlist(starts), width = bdysize),
                     seqinfo = hum_seqinfo)
  return(rdm_bdies_gr)
}