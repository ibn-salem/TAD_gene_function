#builds a GRanges object containing boundaries by fixing the end coordinate of a given GRanges
#object as center and expanding the sequence to both sides to a given size

introduce_boundaries <- function(bedfile, size, seqinfo){

granges <- import(bedfile)      
granges <- as_tibble(granges)
granges <- GRanges (seqnames = granges$seqnames,IRanges(start = granges$end,end = granges$end),
                    seqinfo = hum_seqinfo)
granges <- granges %>% 
  resize(fix = "center", width = size) %>% 
  trim()

grangeslist <- GRangesList(granges)
names(grangeslist) <- "hESC"

return(grangeslist)
}