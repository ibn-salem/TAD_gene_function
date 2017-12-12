#builds a GRanges object containing boundaries by fixing the end coordinate of a given GRanges
#object as center and expanding the sequence to both sides to a given size

introduce_boundaries <- function(bedfilepath, size, seqinfo){

granges <- import(bedfilepath)      
granges <- as_tibble(granges)
granges <- GRanges (seqnames = granges$seqnames,IRanges(start = granges$end,end = granges$end),
                    seqinfo = hum_seqinfo)

#produces 6 out of bound ranges at chromsome ends of chr7,10,11,15,19,21. Therefore these
#boundaries are samller than 40k bp. Since they map at chromosome ends this is not
#considered to introduce a major problem when analysing genepairs separated by a boundary
#because only genes that are located on the same chromosome can form a genepair
granges <- granges %>% 
  resize(fix = "center", width = size) %>% 
  trim()

grangeslist <- GRangesList(granges)
names(grangeslist) <- "hESC"

return(grangeslist)
}