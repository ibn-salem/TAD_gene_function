#loading required libraries
library(GenomicRanges)
library(biomaRt)
library(rtracklayer)
library(genepair)
library(GOSemSim)
library(RColorBrewer)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(vegan)
library(labdsv)
library(mclust)
library(devtools)
library(gplots)
library(tidyverse)
library(readxl)

#sourcing required functions
source("R/fun_Excel_to_Boundary.R")
source("R/fun_Excel_to_TAD.R")
source("R/jaccard_matrix_GRList.R")
source("R/plot_tad_nr_chrsize.R")
source("R/fun_mart_to_granges.R")
source("R/fun_introduce_boundaries.R")
source("R/fun_is_separated.R")
source("R/fun_sample_rdm_bdies.R")
source("R/fun_get_go_sim.R")
source("R/fun_tidy_gosimList.R")
source("R/fun_unique_bdies.R")
source("R/fun_filter_paralogs.R")
source("R/fun_compare_kegg_pathways.R")

#contains information about the human chromosome size etc.
hum_seqinfo <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)

#imports TAD coordinates and introduces boundaries of 40000bp into Dixon et al. IMR90 and 
#hESC dataset
imr90_Dixon_TAD_path <- file.path("data", "IMR90.hg18.bed.hg38.bed")
imr90_Dixon_BDY <- introduce_boundaries(imr90_Dixon_TAD_path, 40000, hum_seqinfo)
names(imr90_Dixon_BDY) <- "IMR90_Dixon"

hESC_TADs_File <- file.path("data", "hESC.hg18.bed.hg38.bed")
hESC_BDY <- introduce_boundaries(hESC_TADs_File, 40000, hum_seqinfo)
names(hESC_BDY) <- "hESC_Dixon"

#imports TAD coordinates and introduce boundaries of 40000bp into Rao et al. IMR90 an 
#GM12878 dataset
imr90_Rao_TAD_path <- file.path("data", "GSE63525_IMR90_Arrowhead_domainlist.bed.hg38.bed")
imr90_Rao_BDY <- introduce_boundaries(imr90_Rao_TAD_path, 40000, hum_seqinfo)
names(imr90_Rao_BDY) <- "IMR90_Rao"

gm12878_Rao_TAD_path <- file.path("data", "GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.bed.hg38.bed")
gm12878_Rao_BDY <- introduce_boundaries(gm12878_Rao_TAD_path, 40000, hum_seqinfo)
names(gm12878_Rao_BDY) <- "GM12878_Rao"

#Imports boundary coordinates of IMR90, hESC and GM12878 cells from Schmitt et al.
TAD_file_Schmitt2016 <- file.path("data", "mmc4.xlsx")
schmitt_BDY <- xltoBDY(TAD_file_Schmitt2016, hum_seqinfo)
  
#Merges all boundary annotations into one GRangesListobject 
all_BDY <- c(schmitt_BDY, imr90_Rao_BDY, gm12878_Rao_BDY, hESC_BDY, imr90_Dixon_BDY)
remove(schmitt_BDY, imr90_Rao_BDY, gm12878_Rao_BDY, hESC_BDY, imr90_Dixon_BDY)


#creates a GRangesList with randomly introduced boundaries for every celltype in the List
all_rdm_BDY <- lapply(all_BDY, sample_rdm_bdies, 40000, hum_seqinfo)
all_rdm_BDY <- GRangesList(rdm_BDY)


  
#downloads all human genes from ENSEMBL and filters them for protein coding genes. 
#Keeps only the longest transcript of every GeneID and builds a GRangesobject out of
#all remaining gene tss
ch38_Mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
gene_atb <- c("ensembl_gene_id",
              "chromosome_name",
              "transcript_end",
              "transcript_start",
              "transcription_start_site",
              "gene_biotype")
  
ch38 <- getBM(attributes = gene_atb, mart = ch38_Mart)
remove(gene_atb)
gr_ch38 <- mart_to_granges(ch38, hum_seqinfo)
  
#write_rds(gr_ch38, "results/tidydata/AllStableProteinGenes_Ch38.rds")
gr_ch38 <- read_rds("results/tidydata/AllStableProteinGenes_Ch38.rds")

#Builds all possible cis pairs out of the genes stored in a genomic ranges object and
#adds two colums containing their gene ids
cisp <- getAllCisPairs(gr_ch38)
cisp <- cisp %>% 
  as_tibble() %>% 
  mutate("g1_id" = gr_ch38$gene_id[cisp$g1]) %>%
  mutate("g2_id" = gr_ch38$gene_id[cisp$g2]) 


#constructs the data classes of all three GO annotations
hsGO_BP <- godata('org.Hs.eg.db', keytype = "ENSEMBL", ont = "BP")
hsGO_MF <- godata('org.Hs.eg.db', keytype = "ENSEMBL", ont = "MF")
hsGO_CC <- godata('org.Hs.eg.db', keytype = "ENSEMBL", ont = "CC")
  
#adds the gosimilarity annotations to the cispair tibble
#WARNING: Takes forever(~ 12h) could not find a way to speed it up so far.
start <- Sys.time()
cisp <- cisp %>% 
  mutate(go_sim_BP = get_go_sim(g1_id, g2_id, hsGO_BP)) %>% 
  mutate(go_sim_MF = get_go_sim(g1_id, g2_id, hsGO_MF)) %>% 
  mutate(go_sim_CC = get_go_sim(g1_id, g2_id, hsGO_CC))
end <- Sys.time()
print(end - start)
  
#write_rds(cisp, "results/tidydata/GenepairWithGoSimUnfiltered.rds")
cisp <- read_rds("results/tidydata/GenepairWithGoSimUnfiltered.rds")

#filters all gene pairs where all three similarities are 'NA', to reduce the size of the data tibble.
#This filtering would otherwise be done automaticly by ggplot during plotting process
variables <- names(cisp)
cisp <- cisp %>% 
  mutate("missing_values" = ifelse(is.na(go_sim_BP),
                                   ifelse(is.na(go_sim_MF),
                                          ifelse(is.na(go_sim_CC),
                                                 "missing",
                                                 "not_missing"
                                          ),
                                          "not_missing"
                                   ),
                                   "not_missing"
  )
  ) %>% 
  filter(missing_values == "not_missing") %>% 
  select(variables)

#downloads paralog and kegg information of all human genes from ensembl, filters the cisp tibble for
#all non paralog pairs and returns them, determines whether or not genepairs are connected via the 
#same KEGG pathway
cisp <- filter_paralogs(cisp, ch38_Mart) 
cisp <- compare_kegg_pathways(cisp, ch38_Mart)

#write_rds(cisp, "results/tidydata/CispWoParalogsAndWithKegg.rds")
cisp <- read_rds("results/tidydata/CispWoParalogsAndWithKegg.rds")


#checks for all celltypes, genepairs and random or real bdies wether two genes 
#are separated by a boundary or not and stores all in one tibble
gene_sep_allbdy <- is_separated(cisp, gr_ch38, all_BDY, hum_seqinfo)
gene_sep_allbdy <- gene_sep_allbdy %>% 
  mutate(real = TRUE)

gene_sep_allrdmbdy <- tibble()
for (number in 1:100) {

  rdm_BDY <- lapply(all_BDY, sample_rdm_bdies, 40000, hum_seqinfo)
  rdm_BDY <- GRangesList(rdm_BDY)
  
  gene_sep_rdmbdy <- is_separated(cisp, gr_ch38, rdm_BDY, hum_seqinfo)
  gene_sep_rdmbdy <- gene_sep_rdmbdy %>% 
    mutate(real = FALSE) %>% 
    mutate(replicate = number)
  
  gene_sep_allrdmbdy <- rbind(gene_sep_allrdmbdy, gene_sep_rdmbdy)
}


sample_rand_boundaries <- function(number){
  
  rdm_BDY <- lapply(all_BDY, sample_rdm_bdies, 40000, hum_seqinfo)
  rdm_BDY <- GRangesList(rdm_BDY)
  
  gene_sep_rdmbdy <- is_separated(cisp, gr_ch38, rdm_BDY, hum_seqinfo)
  separated <- gene_sep_rdmbdy$separated
  # gene_sep_rdmbdy <- gene_sep_rdmbdy %>% 
  #   mutate(real = FALSE) %>% 
  #   mutate(replicate = number)
  
  return(separated)
}

n_rand <- 10
randomList <- map(1:n_rand, sample_rand_boundaries)

large_separated <- unlist(randomList)

largeDF <- cisp[rep(1:nrow(cisp), length(large_separated) / nrow(cisp)), ]

largeDF <- largeDF %>%   
  mutate(
    separated = large_separated,
    real = FALSE,
    replicate = rep(1:n_rand, each = n()/n_rand )
    )

df <- df_mit_allen_scores_und_gruppen

medianDF <- df %>% 
  group(real, separated, replicate) %>% 
  summarize(
    n = n(),
    median = median(GO_sim_BP)
  ) %>% 
  group_by(real, separated) %>% 
  summarize(
    n = n(),
    mean_score = mean(meidan),
    sd_score = sd(meidan)
  )

p <- ggplot(df, aes(x = separted, y = mean_score)) + 
  facet_grid(.~real) +
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(ymin = mean_score - sd_score, ymax = mean_score + sd_score))

  
mytibble <- rbind(gene_sep_allbdy, gene_sep_allrdmbdy)
remove(gene_sep_allbdy, gene_sep_allrdmbdy)

#write_rds(mytibble, "results/tidydata/AllGpairsWoParal.rds")
mytibble <- read_rds("results/tidydata/AllGpairsWoParal.rds")

  
#creating GRangesList object of TADs from Dixon and Schmitt Datasets, count the Genes
#each TAD carries, the NR of TAD per Celltype
  
  #import schmitt tads from Excel-file
  schmitt_TAD <- xltoTAD(TAD_file_Schmitt2016, 10000000)
  
  #import dixon Tads from Bed-file
  dixon_TAD <- import(hESC_TADs_File)
  
  #build a GRangeslist out of schmitt and dixon GRanges
  dixon_TAD <- GRangesList(dixon_TAD)
  names(dixon_TAD) <- "hESC"
  all_TAD <- c(schmitt_TAD, dixon_TAD) 

  
  #counts genes that map into TADs for all celltypes
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

  
  #build a tibble with celltype, the tad ids and the nr of genes per tad id
  gpert <- tibble(celltype, TAD_id, n_genes)

  #stores NR of TAD per Celltypes and their sizes in the built tibble
  tad_per_ctype <- as.data.frame(all_TAD) %>%
    as_tibble() %>% 
    group_by(group_name) %>% 
    summarise(n = n())
  
  names(tad_per_ctype) <- c("celltype", "TAD_nr")
  all_TAD <- all_TAD %>% 
    as.data.frame() %>% 
    as.tibble()
  
  gpert <- gpert %>% 
    mutate(width = all_TAD$width) %>% 
    left_join(tad_per_ctype, by = c("celltype"))
  
  names(gpert)[4] <- "TAD_size"
  
  total_genes <- gr_ch38 %>% 
    as.data.frame() %>% 
    as.tibble() %>% 
    summarise(n = n()) %>% 
    as.integer()
  
  gpert <- gpert %>% 
    mutate(nr_of_all_Genes = total_genes) %>% 
    mutate(meanrange_of_pairs = mean(cisp$dist)) %>% 
    mutate(Nr_of_Pairs = nrow(cisp))
  
  write_rds(gpert, "results/tidydata/AllTadWithGenes.rds")


#loads metadata corresponding to the schmitt Boundaries. 
meta_schmitt_BDY<-read_xlsx("data/ctypes_metadata.xlsx") 

#Counts the Frequency of each boundary across all celltypes annotated in schmitt_BDY.
#Further retrieves all geneIDs of the genes that map into a defined window surround these
#boundaries an saves them in a file within the results/tidydata folder. Another
#file that simply contains all Genes around unique boundaries without caring about the
#celltype. Thes stored GeneID lists can be used to run a Gene Ontology enrichment (for
#example with the Webtool GOrilla)
all_bdies_freq <- find_genes_surround_unique_bdies(schmitt_BDY,
                                                   gr_ch38,
                                                   meta_schmitt_BDY,
                                                   20000,
                                                   hum_seqinfo)
write_rds(all_bdies_freq, "results/tidydata/AllBdiesWithFrequency.rds")

#creating a heatmap of the jaccard coefficients of each entry pair 
#for the coordinates of a GRList object
  jmat <- jaccard_matrix(schmitt_BDY)
  my_palette <- colorRampPalette(c("grey", "red"))(n = 100)
  png(file = "results/figures/HeatmapSchmittBdies.png", height = 380, width = 450)
  heatmap.2(jmat,
            tracecol = NA,
            symm = TRUE,
            col = my_palette,
            dendrogram = "row",
            lhei = c(1,8),
            key = FALSE)
  dev.off()

#converts the jaccard_matrix into a distance matrix and does pcomp
  distmatrix <- vegdist(jmat, method ="bray", diag =TRUE)
  pcomp <- pco(distmatrix, k = 2)
  
  png(file = "results/figures/ClusterSchmittBdies.png", width = 900, height = 720)
  plot(pcomp$points, cex = 2, pch = 16 , col="blue")
  jmat <- as.data.frame(jmat)
  for(i in 1:length(names(jmat))){
    x <- pcomp$points[i,1]
    y <- pcomp$points[i,2]
    text(x,y, labels = names(jmat)[i], pos = 1)
  }
  dev.off()
  