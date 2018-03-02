#builds a GRanges object containing boundaries by fixing the end coordinate of a given 
#GRanges Object as center and expanding the sequence to both sides to a given size

introduce_boundaries <- function(bedfilepath, size, seqinfo){

granges <- import(bedfilepath)      
granges <- as_tibble(as.data.frame(granges))
granges <- GRanges (seqnames = granges$seqnames,IRanges(start = granges$end,end = granges$end),
                    seqinfo = hum_seqinfo)
granges <- granges %>% 
  resize(fix = "center", width = size)
granges <- trim(granges)


grangeslist <- GRangesList(granges)
names(grangeslist) <- "hESC"

return(grangeslist)
}