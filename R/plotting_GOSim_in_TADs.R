library(tidyverse)

all_gosims <- read_rds("results/tidydata/data_mytibble.rds")

  #plot of BP gosims

  BP_gosims <- all_gosims %>% 
    filter(go_sim_BP != "NA")
  
  #plots GOsim of all real boundaries into a geom_smooth plot to show correlation between
  #distance within genepairs and their similarity
  all_gosims %>% 
    filter(go_sim_BP != 'NA') %>%
    filter(real == TRUE) %>% 
    ggplot(aes(x = dist, y = go_sim_BP)) + 
      geom_smooth(aes(color = separated)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs( x = "Distance between Genepairs", y = "BP GO similarity") +
      facet_wrap(~celltype, nrow = 3)
    
  #plots BP GOsim distribution of real boundaries into boxplots
  all_gosims %>% 
    filter(go_sim_BP != 'NA') %>% 
    filter(real == TRUE) %>% 
    ggplot(aes(x = separated, y = go_sim_BP, fill = separated)) +
      geom_boxplot(outlier.shape = NA) +
      ylim(0, 0.75) +
      labs( x = "Separated by Boundary", y = "BP GO similarity") +
      facet_wrap(~celltype, nrow = 3) +
      theme_bw()

  all_gosims %>% 
    filter(go_sim_BP != 'NA') %>% 
    ggplot(aes(x = dist, y = go_sim_BP)) + 
    geom_smooth(aes(color = separated, linetype = real)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs( x = "Distance between Genepairs", y = "BP GO similarity") +
    facet_wrap(~celltype, nrow = 3)
  
    
  #plot MF gosims

  MF_gosims <- all_gosims %>% 
    filter(go_sim_MF != "NA")
  
  ggplot(data = MF_gosims, aes(x = dist, y = go_sim_MF)) + 
    geom_smooth(aes(color = separated, linetype = real)) +
    facet_wrap(~celltype, nrow = 3)
  
  ggplot(data = MF_gosims, aes(x = separated, y = go_sim_MF, fill = separated)) +
    geom_boxplot(outlier.shape = NA) +
    ylim(0, 0.75) +
    facet_wrap(~celltype, nrow = 3) +
    theme_bw()
  
  #plot cc gosims
  
  ggplot(all_gosims, aes(x = separated, y = go_sim_CC, fill = real)) +
    geom_boxplot(outlier.shape = NA) +
    facet_wrap(~celltype, nrow = 3)
  
  CC_gosims <- all_gosims %>% 
    filter(go_sim_CC != "NA")
  
  ggplot(data = CC_gosims, aes(x = dist, y = go_sim_CC)) + 
    geom_smooth(aes(color = separated, linetype = real)) +
    facet_wrap(~celltype, nrow = 3)

  ggplot(data = CC_gosims, aes(x = separated, y = go_sim_CC, fill = separated)) +
    geom_boxplot(outlier.shape = NA) +
    ylim(0, 0.75) +
    facet_wrap(~celltype, nrow = 3) +
    theme_bw()
  
#plots the number of Genes separated and not separated by a boundary into a facet grid
  sept_plot <- tidydf %>% 
    group_by(celltype, real) %>% 
    count(separated)

  ggplot(data = sept_plot, aes(x = separated , y = n, fill = real)) +
    geom_bar(stat = 'identity', position = "dodge") +
    scale_fill_brewer(palette = "Pastel1") +
    geom_text(aes(label = n), vjust = 1, size = 1.9, position = position_dodge(width = 0.9)) +
    facet_wrap(~celltype, nrow = 3) +
    labs(x = " ", y = " Nr. of Genepairs ") +
    theme_minimal()



#loads the dataframe containing annotations about the nr of genes in each tad
all_tad_annotations <- read_rds("results/tidydata/data_all_tad_annotations.rds")

  #plots the frequency in which several numbers of genes per tad occur
  genes_per_tad_plot <- all_tad_annotations %>% 
    group_by(celltype, n_genes) %>% 
    summarize(n = n()) %>% 
    filter(n_genes <=100)
  
  #density distribution of Genefrequency in each TAD and per Celltype
  pdf(file = "results/figures/Genes_per_TAD_densitydist.pdf", width = 12, height = 8)
  ggplot(genes_per_tad_plot, aes(x = n_genes, fill = "gray")) +
    scale_fill_manual(values = "gray") +
    geom_density() +
    facet_wrap(~celltype, nrow = 3) +
    theme_minimal()+
    theme(legend.position = "none") +
    labs(x = "Nr of Genes per TAD", y = "Density")
  dev.off()

  #histogram that shows distribuiton of Genefrequency in each TAD and per Celltype
  pdf(file = "results/figures/Genes_per_TAD_histogramdist.pdf", width = 12, height = 8)
  ggplot(genes_per_tad_plot, aes(x = n_genes, y = n, fill = n)) +
    geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(~celltype, nrow = 3) +
    theme_bw() +
    labs(x = "Nr of Genes per TAD", y = "Frequency of Occurence")
  dev.off()

  #boxplots to show distribution of TAD sizes per Celltype
  pdf(file = "results/figures/TADsize_distribution_per_celltype.pdf", width = 12, height = 8)
  ggplot(all_tad_annotations, aes(x = celltype, y = TAD_size)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    labs(x = "Celltype", y = "TAD sizes") +
    ylim(0,3500000)
  dev.off()
  
  
  pdf(file = "results/figures/NrOfTADs_per_Celltype.pdf", width = 12, height = 8)
  all_tad_annotations %>% 
    group_by(celltype) %>% 
    summarise(n = n()) %>% 
    ggplot( aes(x = celltype, y = n, fill = n)) +
      geom_bar(stat= "identity") +
      theme_bw()
  dev.off()
  


