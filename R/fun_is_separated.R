is_separated <- function(cispairs, genes_granges, grangeslist, seqinfo){
  
  nams <- names(cispairs)
  nams <- c(nams, names(grangeslist))
  
  cisp_granges <- GRanges(seqnames = seqnames(genes_granges[cispairs$g1]),
                          IRanges(start = start(genes_granges[cispairs$g1]),
                                  end = end(genes_granges[cispairs$g2])),
                          seqinfo = seqinfo)
  
  
  for (i in 1:length(grangeslist)){
    cispairs <- cispairs %>% 
      cbind(overlapsAny(cisp_granges, grangeslist[i]))
  }
  
  names(cispairs) <- nams
  
  tidyDF <- cispairs %>% 
    as_tibble() %>%
    gather(key = "celltype", value = "separated", GM12878:hESC)
  
  return(tidyDF)
}