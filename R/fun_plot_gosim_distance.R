plot_gosim_distance <- function(wanted_celltypes, annotation_tibble, width, height) {
  
  separated_cols <- c("chartreuse3", "steelblue1")
  
  gosim_dist <- annotation_tibble %>% 
    filter(celltype %in% wanted_celltypes)  %>% 
    ggplot(aes(x = dist, y = go_sim_BP, color= separated)) +
    geom_smooth() +    
    scale_color_manual(values = c("grey0", "royalblue2")) +
    facet_grid(Name~real)+
    labs(x = "Distance between genes of the same pair (kb)",
         y = "GO similarity of biological processes") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.title = element_text(size = 16),
          axis.text.y = element_text(size = 12),
          legend.position = "bottom",
          legend.text = element_text(size = 14),
          legend.title = element_blank())
  
  myfile <- file.path("results/figures/BP", paste0("gosim_distance_",names(wanted_celltypes)[1]))
  png(myfile, width = width, height = height)
  gosim_dist
}