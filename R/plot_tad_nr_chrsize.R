plot_tad_nr_chrsize <- function (grobject){

#stores the sequence names
chr_names <- levels(seqnames(grobject))

#extracts the nr of tads of hEsc_Tad
nr_of_TAD <- c()
length_chromosome <- c()
for(entry in chr_names){
  s <- subset(grobject, seqnames == entry)
  
  #gets nr of Tads per chromosome
  nr_of_TAD <- c(nr_of_TAD, NROW(s))
  
  #gets length of chromosome
  length_chromosome <- c(length_chromosome, sum(width(s)))
}
df <- data.frame(chr_names, nr_of_TAD, length_chromosome)

level_order <- c( "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                  "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", 
                  "chr21", "chr22", "chrX")
x <- factor(df$chr_names, levels = level_order)

return(plot(x, df$nr_of_TAD))
}