library(tidyverse)

#loads tidy Dataframe with GO similarity scires of Molecular function
tidydf <- read_rds("results/tidydata/data_great_tibbl.rds")

ggplot(data = tidydf, aes(x = dist, y = go_sim)) + 
  geom_smooth(aes(color = separated, linetype = real)) +
  facet_wrap(~celltype, nrow = 3)


all_gosims <- read_rds("results/tidydata/data_mytibble.rds")

#plot of BP gosims
BP_gosims <- all_gosims %>% 
  filter(go_sim_BP != "NA")

ggplot(data = BP_gosims, aes(x = dist, y = go_sim_BP)) + 
  geom_smooth(aes(color = separated, linetype = real)) +
  facet_wrap(~celltype, nrow = 3)
  
#plot MF gosims
MF_gosims <- all_gosims %>% 
  filter(go_sim_MF != "NA")

ggplot(data = MF_gosims, aes(x = dist, y = go_sim_MF)) + 
  geom_smooth(aes(color = separated, linetype = real)) +
  facet_wrap(~celltype, nrow = 3)

#plot cc gosims
CC_gosims <- all_gosims %>% 
  filter(go_sim_CC != "NA")

ggplot(data = CC_gosims, aes(x = dist, y = go_sim_CC)) + 
  geom_smooth(aes(color = separated, linetype = real)) +
  facet_wrap(~celltype, nrow = 3)




#plots the number of Genes separated and not separated by a boundary into a facet grid
sept_plot <- tidydf %>% 
  group_by(celltype, real) %>% 
  count(separated)

g <- ggplot(data = sept_plot, aes(x = separated , y = n, fill = real)) +
  geom_bar(stat = 'identity', position = "dodge") +
  scale_fill_brewer(palette = "Pastel1") +
  geom_text(aes(label = n), vjust = 1, size = 1.9, position = position_dodge(width = 0.9)) +
  facet_wrap(~celltype, nrow = 3) +
  labs(x = " ", y = " Nr. of Genepairs ") +
  theme_minimal()
g 


#plots distribiution of similarity scores within a boxplot
ggplot(tidydf, aes(x = separated, y = go_sim, fill = real)) +
  geom_boxplot() +
  theme_bw() +
  facet_wrap(~celltype, nrow = 3)


#loads the dataframe containing annotations about the nr of genes in each tad
genes_per_tad <- read_rds("Datasets/genes_per_TAD")

#plots the frequency in which several numbers of genes per tad occur
gpt_plot <- genes_per_tad %>% 
  group_by(celltype, n_genes) %>% 
  summarize(n = n()) %>% 
  filter(n_genes <=100)

p <- ggplot(gpt_plot, aes(x = n_genes, y = n, fill = n)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~celltype, nrow = 5) +
  theme_bw() +
  labs(x = "Nr of Genes per TAD", y = "Frequency of Occurence")
p




