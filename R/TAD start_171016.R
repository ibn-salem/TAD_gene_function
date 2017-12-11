#loading required libraries
library(GenomicRanges)
library(biomaRt)
library(rtracklayer)
library(genepair)
library(GOSemSim)

library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(vegan)
library(labdsv)
library(mclust)
library(devtools)
library(gplots)

library(tidyverse)
library(readxl)
library(magrittr)

#sourcing required functions
source("R Scripts/fun_Excel_to_Boundary.R")
source("R Scripts/jaccard_matrix_GRList.R")
source("R Scripts/plot_tad_nr_chrsize.R")
source("R Scripts/ensembl_to_granges.R")
source("R Scripts/fun_mart_to_granges.R")
source("R Scripts/fun_introduce_boundaries.R")
source("R Scripts/fun_getAllRdmBoundaries.R")

#contains information about the chromosome size etc.
hum_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)



#construction of a random boundary list and a real boudary list containing boundaries in all
#available celltypes

'  #introduces boundaries of 40000bp into dixon dataset
  hESC_TADs_File <- file.path("Datasets", "hESC.hg18.bed.hg38.bed")
  hESC_BDY <- introduce_boundaries(hESC_TADs_File, 40000, hum_seqinfo)

  #builds grangeslist of Schmitt boundaries
  TAD_file_Schmitt2016 <- file.path("Datasets", "mmc4.xlsx")
  schmitt_BDY <- xltoBDY(TAD_file_Schmitt2016, hum_seqinfo)
  all_BDY <- c(schmitt_BDY, hESC_BDY)

  #creates a GRangesList with randomly introduced boundaries for every celltype in the List
  all_rdm_BDY <- getAllRdmBoundaries(all_BDY, 40000, hum_seqinfo)
  names(all_rdm_BDY) <- names(all_BDY)'

  #write_rds(all_BDY, "Datasets/data_grangeslist_all_BDY.rds")
  #write_rds(all_rdm_BDY, "Datasets/data_grangeslist_all_rdm_BDY.rds")
  all_BDY <- read_rds("Datasets/data_grangeslist_all_BDY.rds")
  all_rdm_BDY <- read_rds("Datasets/data_grangeslist_all_rdm_BDY.rds")

  
  
#downloads all human genes from ENSEMBL and filters them for protein coding genes. 
#Keeps only the longest transcript of every GeneID and builds a GRangesobject out of
#all remaining genes
  
' ch38_query <- useMart(biomart ="ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl")
  ch38_atb <- c("ensembl_gene_id", "chromosome_name", "transcript_end", "transcript_start", 
              "transcription_start_site", "gene_biotype", "goslim_goa_accession")
  ch38 <- getBM(attributes = ch38_atb, mart = ch38_query)
  
  gr_ch38 <- mart_to_granges(ch38, hum_seqinfo)'
  
  #write_rds(gr_ch38, "Datasets/data_granges_AllFilteredGenes")
  gr_ch38 <- read_rds("Datasets/data_granges_AllFilteredGenes")

  
  
#retrieves all possible cis pairs of the genes stored in a genomic ranges object and
#adds two colums containing their gene ids
  
  cisp_ch38 <- getAllCisPairs(gr_ch38)
  cisp_ch38 <- cisp_ch38 %>% 
    as_tibble() %>% 
    mutate("g1_id" = gr_ch38$gene_id[cisp_ch38$g1]) %>% 
    mutate("g2_id" = gr_ch38$gene_id[cisp_ch38$g2])





#constructs a GOSemSim data class using the human genome and ENSEMBl Gene ids and an 
#ontology out of: ’BP’, ’MF’, ’CC’
hsGO <- godata('org.Hs.eg.db', keytype = "ENSEMBL", ont="BP")

#calculates the Graph-based similarity of GO Terms from all Gene Pairs
go_sim <- vector(mode = "character", length = length(cisp_ch38$g1))
i <- 1

while (i <= length(cisp_ch38$g1)){
  g <- geneSim(gr_ch38$gene_id[cisp_ch38$g1[i]],
             gr_ch38$gene_id[cisp_ch38$g2[i]],
             semData = hsGO,
             measure = "Resnik",
             combine = "max")
  if(typeof(g) == "logical"){
    go_sim[i] <-"NA"
    i <- i + 1
  } else{
    go_sim[i] <- g$geneSim
    i <- i + 1
  }
}

#write_rds(go_sim, "Datasets/gosimlist_bP.rds")
go_sim <- read_rds("Datasets/gosimlist_bP.rds")
cisp_ch38 <- cbind(cisp_ch38, "go_sim" = go_sim)

#converts the list "cisp_ch38" into a tibble
bdy_between_tss <- as_tibble(cisp_ch38)






#extracts the name from the dataframe above and from all celltypes 
#from the boundary list and concatenates them
nms <- names(bdy_between_tss)
nms <- c(nms, names(all_BDY))

#builds a new GRanges object with the seqnames from gene 1 of the genepair (assumption is that gene pairs
#that are found by "getAllCisPairs" are one the same chromosme only) and the coordinates of the tss of both 
#genes using g1 tss as start and g2 tss as end of a range.
cisp_gr <-GRanges(seqnames = seqnames(gr_ch38[cisp_ch38$g1]),
                  IRanges(start = start(gr_ch38[cisp_ch38$g1]),
                          end = start(gr_ch38[cisp_ch38$g2])),
                  seqinfo = hum_seqinfo)

#extends the dataframe from above in every loop with one colum for a celltype. 
#In this colum it is written, wether the genepair ist seperated by a boundary (TRUE)
#or not ( FALSE)
for(i in 1:length(all_BDY)){
  bdy_between_tss <- cbind(bdy_between_tss, overlapsAny(cisp_gr, all_BDY[i]))
}

#assigns names to all colums of the data.frame
names(bdy_between_tss) <- nms

#creates a data frame containing the gene pairs, their distance and whether or not they are
#separated in each cell type
tidyDF <- bdy_between_tss %>% 
  as_tibble() %>%
  gather(key = "celltype", value = "separated", GM12878:hESC)

write_rds(tidyDF, "Datasets/tidy_df_BP.rds")

# creates a dataframe with TAD ids per cell type and the number of genes that they carry 
TAD_id <- c()
n_genes <- c()
celltype <- c()
i <- 1

while (i <= length(all_TAD)){
  TAD_id <- c(TAD_id, 1:NROW(all_TAD[[i]]))
  celltype<- c(celltype, rep(names(all_TAD[i]), NROW(all_TAD[[i]])))
  ovl <- findOverlaps(all_TAD[[i]], gr_ch38)
  n_genes <- c(n_genes, countLnodeHits(ovl))
  i<-i+1
}

gpert <- tibble(TAD_id, n_genes, celltype)
write_rds(gpert, "Datasets/genes_per_TAD")


#plots the TAD sizes per celltype as boxplot
plot_tadsizes <- as.data.frame(width(all_TAD))
ggplot(plot_tadsizes, aes(x = group_name, y = value, fill = group_name)) + 
  labs(title="Distribution of TAD sizes",
       x="Cell- / Tissuetype",
       y="TAD size") +
  geom_boxplot() + 
  scale_y_log10()

#stores the number of TADs per celltype in a list and plots it as bargraph
nr_of_tad <- c()
i <- 1
while (i <= length(all_TAD)){
  nr_of_tad <- c(nr_of_tad, NROW(all_TAD[[i]]))
  i <- i + 1
}
celltypes <- names(all_TAD)
plot_nroftad <- as.data.frame(nr_of_tad)
plot_nroftad <- cbind(plot_nroftad, celltypes)
ggplot(plot_nroftad, aes(x = celltypes, y = nr_of_tad, fill = nr_of_tad)) + geom_col() 


#creating a heatmap of the jaccard coefficients of each entry pair 
#for the coordinates of a GRList object
jmat <- jaccard_matrix(schmitt_BDY)
my_palette <- colorRampPalette(c("grey", "red", "blue"))(n = 299)
heatmap.2(jmat, tracecol = NA, symm = TRUE, col = my_palette)

#converts the jaccard_matrix into a distance matrix and does pcomp
distmatrix <- vegdist(jmat, method ="bray", diag =TRUE)
pcomp <- pco(distmatrix, k = 2)
plot(pcomp$points, cex = 2, pch = 16 , col="blue")
jmat <- as.data.frame(jmat)
for(i in 1:length(names(jmat))){
  x <- pcomp$points[i,1]
  y <- pcomp$points[i,2]
  text(x,y, labels = names(jmat)[i], pos = 2)
}
