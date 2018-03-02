plot_gosim <- function(wanted_celltypes, annotation_tibble, width, height) {
  
  separated_cols <- c("chartreuse3", "steelblue1")
  
  gosim <- annotation_tibble %>% 
    filter(celltype %in% wanted_celltypes)  %>% 
    ggplot(aes(x = xfactor, y = go_sim_BP)) +
    geom_boxplot(aes(fill = real, color = separated),
                 outlier.shape = NA,
                 width = 0.5,
                 lwd = 0.4) +
    geom_signif(comparisons = list(c("Hi-C & same TAD", "Hi-C & different TAD")),
                y_position = 0.75,
                tip_length = 0.02,
                textsize = 3,
                vjust = -0.2)+
    geom_signif(comparisons = list(c("randomised & same TAD", "randomised & different TAD")),
                y_position = 0.75,
                tip_length = 0.02,
                textsize = 3,
                vjust = -0.2)+
    geom_signif(comparisons = list(c("Hi-C & same TAD", "randomised & same TAD")),
                y_position = 0.875,
                textsize = 3,
                vjust = -0.2) +
    scale_fill_manual(values = separated_cols) +
    scale_color_manual(values = c("grey0", "grey62")) +
    facet_grid(.~Name) +
    labs(y = "GO similarity of biological processes") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size =10),
          legend.title = element_blank(),
          legend.text = element_text(size = 12),
          legend.position = "bottom")
  
  myfile <- file.path("results/figures/BP", paste0("gosim_",names(wanted_celltypes)[1]))
  png(myfile, width = width, height = height)
  gosim
}