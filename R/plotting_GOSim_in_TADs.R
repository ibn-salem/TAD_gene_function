library(tidyverse)
library(ggsignif)
library(ggpubr)
library(RColorBrewer)
library(readxl)

source("R/fun_plot_binned_gosim.R")
source("R/fun_plot_gosim.R")
source("R/fun_plot_gosim_distance.R")

#load tibble with gene pairs and their go annotations
all_gosims <- read_rds("results/tidydata/AllGpairsWoParal.rds")
metadata <- read_xlsx("data/ctypes_metadata.xlsx")
names(metadata)[1] <- "celltype"

all_gosims <- all_gosims %>% 
  mutate(sep  = ifelse( separated == TRUE, "different TAD", "same TAD")) %>% 
  mutate(random = ifelse(real == TRUE, "Hi-C", "randomised")) %>% 
  select(celltype, "real" = random, "separated" = sep, go_sim_BP, go_sim_MF, go_sim_CC, dist, same_kegg_pathway) 


#adds a binning variable to genepairs according to their distance to each other
all_gosims <- all_gosims %>% 
  left_join(metadata, by = "celltype") %>% 
  filter(go_sim_BP != 'NA') %>% 
  mutate(bin = cut_interval(dist, n = 4)) %>% 
  mutate(bin = ifelse(bin == "[0,2.5e+05]",
                       "0 - 250kb", 
                       ifelse(bin == "(2.5e+05,5e+05]",
                              "250kb - 500kb",
                              ifelse(bin == "(5e+05,7.5e+05]",
                                     "500kb - 750kb",
                                     "750kb - 1000kb")
                              )
                       )
         ) 
#adds a grouping variable for the combinations of real and separated
all_gosims <- all_gosims %>% 
  mutate(xfactor = ifelse(real == "Hi-C", 
                           ifelse(separated == "same TAD",
                                  "Hi-C & same TAD",
                                  "Hi-C & different TAD" ),
                           ifelse(separated == "same TAD",
                                  "randomised & same TAD",
                                  "randomised & different TAD")
                           )
  )


all_gosims %>%
  group_by(separated, same_kegg_pathway) %>% 
  summarise(n = n())


gosim <- all_gosims %>%
  ggplot(aes(x = xfactor , y = go_sim_BP)) +
  facet_grid(.~Name+Dataset) +
  geom_boxplot(aes(fill = separated, color = real),
               outlier.shape = NA,
               width = 0.5,
               lwd = 1) +
  geom_signif(comparisons = list(c("Hi-C & same TAD", "Hi-C & different TAD")),
              y_position = 0.75,
              tip_length = 0.02,
              textsize = 4,
              vjust = -0.2) +
  geom_signif(comparisons = list(c("randomised & same TAD", "randomised & different TAD")),
              y_position = 0.75,
              tip_length = 0.02,
              textsize = 4,
              vjust = -0.2) +
  geom_signif(comparisons = list(c("Hi-C & same TAD", "randomised & same TAD")),
              y_position = 0.875,
              tip_length = 0.02,
              textsize = 4,
              vjust = -0.2) +
  scale_fill_manual(values = c("grey95", "red4")) +
  scale_color_manual(values = c("navyblue", "snow4")) +
  labs(y = "BP GO similarity scores") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size =15, color = "black"),
        axis.ticks.y = element_line(size = 1),
        strip.text = element_text(size = 15),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "bottom")

gosim

ggsave(filename = "GosimEscForPoster.pdf",
       path = "results/figures",
       device = c("pdf"),
       units = c("cm"),
       height = 10,
       width = 20)

gosim_dist <- all_gosims %>% 
  ggplot(aes(x = dist, y = go_sim_BP, color = separated)) +
  geom_smooth(lwd = 2) +    
  scale_color_manual(values = c("grey60", "red4")) +
  facet_grid(Name+Dataset~real) +
  labs(x = "Genomic Distance (bp)",
       y = "BP GO similarity scores") +
  theme_bw() +
  theme(axis.ticks = element_line(size = 1),
        axis.text.x = element_text(size = 20, angle = 45, hjust = 1, color = "black"),
        axis.title = element_text(size = 21, color = "black"),
        axis.text.y = element_text(size = 20, color = "black"),
        strip.text = element_text(size = 20, color = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, color = "black"),
        legend.position = "bottom")

gosim_dist

ggsave(filename = "GosimDistanceEscForPoster.pdf",
       path = "results/figures",
       device = c("pdf"),
       units = c("cm"),
       height = 12,
       width = 20)

binned_gosim <- all_gosims %>% 
  ggplot(aes(x = xfactor, y = go_sim_BP)) +
  geom_boxplot(aes(fill = separated, color = real),
               outlier.shape = NA,
               width = 0.5,
               lwd = 1) +
  geom_signif(comparisons = list(c("Hi-C & same TAD", "Hi-C & different TAD")),
              y_position = 0.75,
              tip_length = 0.02,
              textsize = 5,
              vjust = -0.2) +
  geom_signif(comparisons = list(c("randomised & same TAD", "randomised & different TAD")),
              y_position = 0.75,
              tip_length = 0.02,
              textsize = 5,
              vjust = -0.2) +
  geom_signif(comparisons = list(c("Hi-C & same TAD", "randomised & same TAD")),
              y_position = 0.875,
              tip_length = 0.02,
              textsize = 5,
              vjust = -0.2) +
  scale_fill_manual(values = c("grey95", "red4")) +
  scale_color_manual(values = c("navyblue", "snow4")) +
  facet_grid(bin ~ Name+Dataset) +
  labs(y = "BP GO similarity scores") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size =20, color = "black"),
        axis.ticks.y = element_line(size = 2),
        strip.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.position = "bottom") 

binned_gosim
ggsave(filename = "GosimBinnedEscForPoster.pdf",
       path = "results/figures",
       device = c("pdf"),
       units = c("cm"),
       height = 36,
       width = 22)




#bins and go sim scores plotted for all celltypes

part_cts <- c("hESC", "H1", "IMR90", "GM12878", "SX")
names(part_cts) <- "Example"
plot_gosim(part_cts, all_gosims, 600, 250)
dev.off()
plot_binned_gosim(part_cts, all_gosims, 620, 600)
dev.off()
plot_gosim_distance(part_cts, all_gosims, 400,680)
dev.off()


all_celltypes <- unique(all_gosims$celltype)

first_cts <- all_celltypes[1:8]
names(first_cts) <- "1-8"
plot_gosim(first_cts, all_gosims, 920, 250)
dev.off()
plot_binned_gosim(first_cts, all_gosims, 990, 600)
dev.off()
plot_gosim_distance(first_cts, all_gosims, 400, 980)
dev.off()

middle_cts <- all_celltypes[9:16]
names(middle_cts)<- "9-16"
plot_gosim(middle_cts, all_gosims, 920, 250)
dev.off()
plot_binned_gosim(middle_cts, all_gosims, 990, 600)
dev.off()
plot_gosim_distance(middle_cts, all_gosims, 400, 980)
dev.off()

last_cts <- as.list(all_celltypes[17:22])
names(last_cts) <-"17-22"
plot_gosim(last_cts, all_gosims, 920, 250)
dev.off()
plot_binned_gosim(last_cts, all_gosims, 740, 600)
dev.off()
plot_gosim_distance(last_cts, all_gosims, 400, 980)
dev.off()

  



  
#plots the number of Genes separated and not separated by a boundary into a facet grid
sept_plot <- all_gosims %>% 
  filter(real == "Hi-C") %>% 
  group_by(Name, celltype) %>% 
  count(separated)

ggplot(data = sept_plot, aes(x = Name, y = n, fill = separated)) +
  geom_bar(stat = 'identity', position = "dodge") +
  scale_fill_manual(values =c("red4", "grey60")) +
  labs(x = " ", y = " Nr. of Genepairs ") +
  theme_bw() + 
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 28),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(size = 28, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 28, color = "black"),
        axis.title = element_text(size = 30))

ggsave(filename = "NrOfGenepairsInCtypes.pdf",
       path = "results/figures",
       device = c("pdf"),
       units = c("cm") ,
       width = 35,
       height = 18)


#loads the dataframe containing annotations about the nr of genes in each tad
all_tad_annotations <- read_rds("results/tidydata/AllTadWithGenes.rds")

all_tad_annotations %>% 
  filter(celltype %in% c("hESC", "H1")) %>% 
  mutate(celltype = ifelse(celltype == "hESC","Embryonic Stem Cell (2012)","Embryonic Stem Cell (2016")) %>% 
  group_by(n_genes) %>% 
  ggplot(aes (x = n_genes, fill = celltype)) +
  geom_histogram(bins = 50) + 
  scale_fill_manual(values = c("red4", "orange2" )) +
  labs(y = "Frequency", x = "Nr of Genes in a TAD") +
  xlim(0,50) +
  facet_grid(.~celltype) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.x = element_text(size = 38),
        axis.title = element_text(size = 38),
        axis.text = element_text(size = 38, color = "black"))

ggsave(filename = "GenesPerTadDistribution.pdf",
       path = "results/figures",
       device = c("pdf"),
       units = c("cm") ,
       width = 44,
       height = 20)

x <- all_gosims %>% 
  select(celltype, Name) %>% 
  unique()

all_tad_annotations %>% 
  left_join(x, by = c("celltype")) %>% 
  group_by(Name) %>% 
  ggplot(aes(x = Name, fill = TAD_nr)) +
  geom_bar() +
  scale_fill_gradient(low = "grey60" , high = "firebrick4") +
  theme_bw() +
  labs(x = "", y = "", fill = "Nr. of TADs") + 
  theme(axis.text.x = element_text(size = 27, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size =27, color = "black"),
        axis.ticks = element_line(size = 2),
        legend.position=  "none")

ggsave(filename = "TADsPerCtype.pdf",
       path = "results/figures",
       device = c("pdf"),
       units = c("cm") ,
       width = 46,
       height = 23)
  

x <- all_tad_annotations %>% 
  group_by(celltype) %>% 
  mutate(nr_genes_in_TAD = sum(n_genes)) %>% 
  mutate(mean_TAD_size = mean(TAD_size))

x <- x %>% 
  select_("celltype", "TAD_nr","mean_TAD_size", "nr_of_all_Genes", "nr_genes_in_TAD", "Nr_of_Pairs", "meanrange_of_pairs")

xnames <- c("Celltype", "Nr_TAD","mean_size_TAD", "All_Genes", "Genes_in_TAD", "Nr_Pairs", "meanrange_Pairs")
names(x) <- xnames
x <- x %>% 
  mutate(mean_size_TAD = as.integer(mean_size_TAD)) %>% 
  mutate(meanrange_Pairs = as.integer(meanrange_Pairs))


x <- unique(x)

write_tsv(x, "results/tidydata/meta_for_report")

