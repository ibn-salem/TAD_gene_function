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
source("R/fun_Excel_to_Boundary.R")
source("R/jaccard_matrix_GRList.R")
source("R/plot_tad_nr_chrsize.R")
source("R/fun_mart_to_granges.R")
source("R/fun_introduce_boundaries.R")
source("R/fun_is_separated.R")
source("R/fun_sample_rdm_bdies.R")
source("R/fun_get_go_sim.R")
source("R/fun_tidy_gosimList.R")

#contains information about the human chromosome size etc...
hum_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)


#construction of a random boundary list and a real boudary list containing boundaries in all
#available celltypes

  #introduces boundaries of 40000bp into dixon dataset
  hESC_TADs_File <- file.path("data", "hESC.hg18.bed.hg38.bed")
  hESC_BDY <- introduce_boundaries(hESC_TADs_File, 40000, hum_seqinfo)

  #builds grangeslist of Schmitt boundaries
  TAD_file_Schmitt2016 <- file.path("data", "mmc4.xlsx")
  schmitt_BDY <- xltoBDY(TAD_file_Schmitt2016, hum_seqinfo)
  
  #merges grangeslists of schmitt and dixon datasets
  all_BDY <- c(schmitt_BDY, hESC_BDY)

  #creates a GRangesList with randomly introduced boundaries for every celltype in the List
  all_rdm_BDY <- lapply(all_BDY, sample_rdm_bdies, 40000, hum_seqinfo)
  all_rdm_BDY <- GRangesList(all_rdm_BDY)
  
  #write_rds(all_BDY, "results/tidydata/all_real_bdies.rds")
  #write_rds(all_rdm_BDY, "results/tidydata/all_rdm_bdies.rds")
  
  all_BDY <- read_rds("results/tidydata/all_real_bdies.rds")
  all_rdm_BDY <- read_rds("results/tidydata/all_rdm_bdies.rds")

  
#downloads all human genes from ENSEMBL and filters them for protein coding genes. 
#Keeps only the longest transcript of every GeneID and builds a GRangesobject out of
#all remaining genes
  
  ch38_query <- useMart(biomart ="ENSEMBL_MART_ENSEMBL",
                  dataset = "hsapiens_gene_ensembl")
  ch38_atb <- c("ensembl_gene_id", "chromosome_name", "transcript_end", "transcript_start", 
              "transcription_start_site", "gene_biotype", "goslim_goa_accession")
  ch38 <- getBM(attributes = ch38_atb, mart = ch38_query)
  gr_ch38 <- mart_to_granges(ch38, hum_seqinfo)
  
  #write_rds(gr_ch38, "Datasets/data_granges_AllFilteredGenes.rds")
  gr_ch38 <- read_rds("results/tidydata/data_granges_AllFilteredGenes")
  
  
#retrieves all possible cis pairs of the genes stored in a genomic ranges object and
#adds two colums containing their gene ids
  
  cisp <- getAllCisPairs(gr_ch38)
  cisp <- cisp %>% 
    as_tibble() %>% 
    mutate("g1_id" = gr_ch38$gene_id[cisp$g1]) %>% 
    mutate("g2_id" = gr_ch38$gene_id[cisp$g2])


#construction of a data classes to compare semantic similarities among GO Terms
#via ensembl geneIDs and comparison of GO annotations for all annotated genepairs

  #constructs the data classes of all three GO annotations
  hsGO_BP <- godata('org.Hs.eg.db', keytype = "ENSEMBL", ont = "BP")
  hsGO_MF <- godata('org.Hs.eg.db', keytype = "ENSEMBL", ont = "MF")
  hsGO_CC <- godata('org.Hs.eg.db', keytype = "ENSEMBL", ont = "CC")
  
  #adds the gosimilarity annotations to the cispair tibble
  #WARNING: Takes forever(~ 12h) could not find a way to speed it up so far. 
  #In addition, the readout type is not uniform (list of lists with different lengths)
  cisp <- cisp %>% 
    mutate(go_sim_BP = get_go_sim(cisp, hsGO_BP)) %>% 
    mutate(go_sim_MF = get_go_sim(cisp, hsGO_MF)) %>% 
    mutate(go_sim_CC = get_go_sim(cisp, hsGO_CC))
  
  #write_rds(cisp, "results/tidydata/cisp_with_allgosim.rds")
  cisp <- read_rds("results/tidydata/cisp_with_allgosim.rds")
  
  #extracts only the geneSim value of each list of lists
  cisp <- cisp %>% 
    mutate(go_sim_BP = tidy_gosimList(go_sim_BP)) %>% 
    mutate(go_sim_MF = tidy_gosimList(go_sim_MF)) %>% 
    mutate(go_sim_CC = tidy_gosimList(go_sim_CC))
  
  
#checks for all celltypes, genepairs and random or real bdies wether two genes 
#are separated by a boundary or not and stores all in one tibble
  gene_sep_allbdy <- is_separated(cisp, gr_ch38, all_BDY, hum_seqinfo)
  gene_sep_allbdy <- gene_sep_allbdy %>% 
    mutate(real = TRUE)
  
  gene_sep_allrdmbdy <- is_separated(cisp, gr_ch38, all_rdm_BDY, hum_seqinfo)
  gene_sep_allrdmbdy <- gene_sep_allrdmbdy %>% 
    mutate(real = FALSE)
  
  mytibble <- rbind(gene_sep_allbdy, gene_sep_allrdmbdy)

  #write_rds(mytibble, "results/tidydata/data_mytibble.rds")
  mytibble <- read_rds("results/tidydata/data_mytibble.rds")
  

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
