is_separated <- function(cispairs, genes_granges, grangeslist, seqinfo){
  
  #extracts the names of variables stored in cispairs and in the grangeslist 
  #object(celltypes in thsi case) and concatenates them
  nams <- names(cispairs)
  nams <- c(nams, names(grangeslist))
  
  #builds a new granges object with the tss' of g1 and g2 of each cispair as start and
  #end of a range
  cisp_granges <- GRanges(seqnames = seqnames(genes_granges[cispairs$g1]),
                          IRanges(start = start(genes_granges[cispairs$g1]),
                                  end = end(genes_granges[cispairs$g2])),
                          seqinfo = seqinfo)
  
  #test for all entries in the GRangesList wether or not there are any overlaps between
  #the cispair ranges and the ranges in Granges object and stores the logical output in
  #an additional column
  for (i in 1:length(grangeslist)){
    cispairs <- cispairs %>% 
      cbind(overlapsAny(cisp_granges, grangeslist[i]))
  }
  
  #renames the colums of the tibble created above
  names(cispairs) <- nams
  
  #transforms teh tibble to have all celltypes as a variable and therefore stored
  #in one column
  tidyDF <- cispairs %>% 
    as_tibble() %>%
    gather(key = "celltype", value = "separated", GM12878:hESC)
  
  return(tidyDF)
}